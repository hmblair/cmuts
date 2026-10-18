/* phmm.c -- a banded pair HMM.
 *
 * The read is cut into segments. A read whose sequencer reports run-length errors has
 * one segment per run of one base, since the errors enter where a run ends. A read
 * without them has segments of up to SEGMENT_BASES bases, since nothing enters where
 * such a segment ends.
 *
 * The matrix has one outer row per segment, and row 0 before the first. A cell of an
 * outer row is a prefix length of the reference, and holds three states: the segment
 * is complete and the last event was a pairing, a deletion is open after it, or an
 * insertion is open after it.
 *
 * Between two outer rows lies the inner matrix of a segment, which pairs the bases of
 * the segment with the reference one at a time. Its rows are the bases of the segment,
 * and one more for a template base that an Ultima sequencer did not count. Each
 * row holds the prefix lengths within the band of its base, and nothing else. Its
 * states are the three of a pair HMM over bases, with the weights of that model. A
 * pairing pays the weight of its base after no gap, or the opening rate of the gap it
 * closes, times the emission of the comparison. A gap pays its extension rate per base
 * after the first.
 *
 * The reverse transcriptase copies the reference 3' to 5', so the read's 5'-most
 * pairing is where it reached the end state. The passes run 5' to 3', so the end
 * state is where every path starts: row 0 holds the termination rate of each base,
 * and the first pairing pays nothing more. Stopping is one arm of the decision at
 * each later pairing, beside the openings of the two gaps.
 *
 * A deletion that follows the last pairing of a segment leaves the inner matrix as the
 * open deletion state of the outer cell, and extends or closes in the next inner
 * matrix. A deletion inside a segment opens, extends, and closes within one inner
 * matrix. The two never describe the same path, since an open deletion extends only in
 * the inner matrix that closes it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "phmm.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "align.h"
#include "nuc.h"

/* The states of an outer cell. */
enum { CLOSED, DELETION, INSERTION, N_OUTER };

/* The states of an inner cell. */
enum { PAIRED, GAP_DELETION, GAP_INSERTION, N_INNER };

/* Emission of a comparison that carries no information. */
#define UNINFORMATIVE (1.0 / NUC_BASES)

/* A misread or modified base is equally likely to be any base the reference does not
 * hold. */
#define OTHER_BASES ((double)(NUC_BASES - 1))

/* Forward times backward sums to one on the first row within this tolerance. */
#define NORMALIZATION_TOLERANCE 1e-6

/* The most bases a segment of a read without run-length errors holds. The paths
 * through this many bases stay far above the smallest double. */
#define SEGMENT_BASES 32

/* An outer row is divided by its total only when the total leaves this range. The
 * posteriors of a segment divide by the scale of its row, which is one where the row
 * was not divided, so the result is the same either way. The range keeps the totals
 * far from the smallest and the largest double. */
#define SCALE_FLOOR   1e-50
#define SCALE_CEILING 1e50

/* The step over an inner row is written once and compiled twice, for a read base the
 * sequencer counted and for one it did not. The compiler makes the two copies only
 * when told to inline the step at both calls. */
#define ALWAYS_INLINE static inline __attribute__((always_inline))

void phmm_rates_set_transitions(phmm_rates *rates, const phmm_uniform_rates *uniform)
{
    rates->match_to_insertion     = 1.0 - uniform->extend_insertion;
    rates->match_to_deletion      = 1.0 - uniform->extend_deletion;
    rates->insertion_to_insertion = uniform->extend_insertion;
    rates->deletion_to_deletion   = uniform->extend_deletion;
}

/* ------------------------------------------------------------------------ */
/* The model                                                                 */
/* ------------------------------------------------------------------------ */

/* Returns a 0-based base clamped to a reference of len bases. */
static hts_pos_t clamp_base(size_t len, hts_pos_t b)
{
    if (b < 0) {
        b = 0;
    }
    if ((size_t)b >= len) {
        b = (hts_pos_t)len - 1;
    }

    return b;
}

/* Returns a per-base rate at a 0-based base, clamped to the reference. */
static double rate_at(const double *values, size_t len, hts_pos_t b)
{
    return values[clamp_base(len, b)];
}

/* Returns the weight a pairing of base b pays when no gap precedes it: what opening
 * either gap and stopping leave. */
static double after_no_gap_at(const phmm_rates *rates, size_t len, hts_pos_t b)
{
    return 1.0 - rate_at(rates->open_insertion, len, b)
               - rate_at(rates->open_deletion, len, b - 1)
               - rate_at(rates->termination, len, b);
}

static int grow(void **ptr, size_t *cap, size_t needed, size_t size)
{
    void *grown;

    if (needed <= *cap) {
        return 0;
    }

    grown = realloc(*ptr, needed * size);

    if (!grown) {
        return -1;
    }

    *ptr = grown;
    *cap = needed;

    return 0;
}

static int grow_model(phmm_model *model, size_t len)
{
    if (grow((void **)&model->runs, &model->cap_runs, len, sizeof *model->runs) < 0) {
        return -1;
    }

    return grow((void **)&model->bases, &model->cap_bases, len, sizeof *model->bases);
}

static phmm_base base_at(const phmm_rates *rates, const cm_fasta_record *ref,
                           hts_pos_t b)
{
    size_t len = ref->len;

    return (phmm_base){
        .after_no_gap    = after_no_gap_at(rates, len, b),
        .after_deletion  = rate_at(rates->open_deletion, len, b - 1)
                         * rates->match_to_deletion,
        .after_insertion = rate_at(rates->open_insertion, len, b),
        .modification    = rates->modification[b],
        .code            = (uint8_t)nuc_from_char(ref->seq[b]),
        .run_bit         = 0,
    };
}

static void fill_bases(phmm_model *model, const cm_fasta_record *ref)
{
    for (size_t b = 0; b < ref->len; b++) {
        model->bases[b] = base_at(model->rates, ref, (hts_pos_t)b);
    }

    for (size_t m = 0; m < model->n_runs; m++) {
        const run *r = &model->runs[m];

        /* The runs were split from the reference, and both buffers are sized to it.
         * NOLINTNEXTLINE(clang-analyzer-security.ArrayBound) */
        if (r->len == 1 && nuc_is_base(r->base)) {
            model->bases[r->start].run_bit = (uint8_t)(1U << r->base);
        }
    }
}

int phmm_model_prepare(phmm_model *model, const phmm_rates *rates,
                         const cm_fasta_record *ref, int band)
{
    if (grow_model(model, ref->len) < 0) {
        return -1;
    }

    model->rates  = rates;
    model->len    = ref->len;
    model->band   = band;
    model->n_runs = runs_of_text(ref->seq, ref->len, model->runs);

    fill_bases(model, ref);

    return 0;
}

void phmm_model_free(phmm_model *model)
{
    free(model->runs);
    free(model->bases);
    *model = (phmm_model){ 0 };
}

/* The comparison of one read base with each reference base, by the code of the
 * reference base, as lines in the modification rate m of the reference base. The
 * emission is intercept + slope * m. The part of the emission on which the template
 * differed is share * m. A read base agrees with its template base with chance
 * (1 - m)(1 - e) + m e / 3, and differs from it into one given base with chance
 * (1 - m) e / 3 + m (1 - e / 3) / 3, where e is the chance the base was misread. An
 * unnamed base on either side emits the uninformative chance and carries nothing. A
 * template base the sequencer did not count was never observed, so it emits one and
 * carries nothing. */
typedef struct {
    double intercept[NUC_COUNT];
    double slope[NUC_COUNT];
    double share[NUC_COUNT];
} comparisons;

/* Writes the comparison with a reference base of code theirs. */
static void set_comparison(comparisons *o, nuc theirs, double intercept, double slope,
                           double share)
{
    o->intercept[theirs] = intercept;
    o->slope[theirs]     = slope;
    o->share[theirs]     = share;
}

/* Writes the comparisons of a read base ours, counted or not, that was misread with
 * chance error. */
static void set_comparisons(comparisons *o, nuc ours, bool uncounted, double error)
{
    double other = error / OTHER_BASES;
    double shift = 1.0 - error - other;
    double blank = uncounted ? 1.0 : UNINFORMATIVE;

    for (nuc theirs = NUC_N; theirs < NUC_COUNT; theirs++) {
        set_comparison(o, theirs, blank, 0.0, 0.0);
    }

    if (uncounted || !nuc_is_base(ours)) {
        return;
    }

    for (nuc theirs = NUC_A; theirs < NUC_COUNT; theirs++) {
        set_comparison(o, theirs, other, shift / OTHER_BASES,
                       (1.0 - other) / OTHER_BASES);
    }

    set_comparison(o, ours, 1.0 - error, -shift, other);
}

/* The comparisons of a read base of each base at each quality score, and of such a
 * base in a read that stores no qualities, which is never misread. The last entry is
 * for a template base the sequencer did not count. */
struct phmm_emissions {
    comparisons at_quality[NUC_COUNT][PHRED_MAX + 1];
    comparisons without_quality[NUC_COUNT];
    comparisons uncounted;
};

phmm_emissions *phmm_emissions_create(const phred *quality)
{
    phmm_emissions *emissions = malloc(sizeof *emissions);

    if (!emissions) {
        return NULL;
    }

    for (nuc ours = NUC_N; ours < NUC_COUNT; ours++) {
        for (unsigned q = 0; q <= PHRED_MAX; q++) {
            set_comparisons(&emissions->at_quality[ours][q], ours, false,
                            phred_error(quality, (uint8_t)q));
        }

        set_comparisons(&emissions->without_quality[ours], ours, false, 0.0);
    }

    set_comparisons(&emissions->uncounted, NUC_N, true, 0.0);

    return emissions;
}

void phmm_emissions_destroy(phmm_emissions *emissions)
{
    free(emissions);
}

/* ------------------------------------------------------------------------ */
/* Scratch                                                                   */
/* ------------------------------------------------------------------------ */

/* Where one outer row lies and what its inner matrix needs. Row 0 stands before the
 * first segment, around where the CIGAR starts. */
typedef struct {
    int32_t start;                   /* first read base of the segment; unused on row 0 */
    int32_t len;                     /* bases of the segment; unused on row 0 */
    int32_t lo, hi;                  /* the prefixes the row's cells cover */
    int32_t last_inner_row;          /* the segment's bases, and the uncounted rows */
    const double *length;            /* chance of each template length, by offset */
    double  missed;                  /* chance a one-base run after it was not read */
    uint8_t missed_bases;            /* the bases that run can be, as bits by nuc */
    size_t  inner;                   /* offset of the inner matrix in scratch */
    size_t  first_span;              /* index of the span of inner row 0 */
} outer_row;

/* One inner row of a segment: the prefixes of its band, from from to to; where its
 * cells lie, as an offset in doubles from the start of the inner matrix of its
 * segment; the comparisons of its read base; and what the row allows. */
typedef struct {
    int32_t                 from;
    int32_t                 to;
    size_t                  cells;
    const comparisons      *comparisons;
    bool                    uncounted;  /* a template base the sequencer did not count */
    bool                    deletes;    /* a deletion can open inside after the row */
    bool                    inserts;    /* the read base of the row can be inserted */
} inner_span;

/* The band around one point of the CIGAR path, within the reference. It runs from the
 * band before the prefix reached once the base is taken, to the band after the prefix
 * reached once the deletion straight after the base is taken too. */
typedef struct {
    int32_t from;
    int32_t to;
} band_edges;

struct phmm_scratch {
    aln_place     *places;
    size_t         cap_places;
    band_edges    *bands;        /* one per place */
    size_t         cap_bands;
    outer_row    *rows;
    size_t         cap_rows;
    inner_span    *spans;        /* one per inner row of every inner matrix */
    size_t         cap_spans;
    double        *forward;      /* outer, the last row, widest * N_OUTER */
    double        *backward;     /* outer, two rows of widest * N_OUTER */
    size_t         cap_outer;
    double        *scale;        /* one per outer row */
    size_t         cap_scale;
    double        *inner;        /* the inner matrices of every segment; row 0 of each
                                    is the outer row before it */
    size_t         cap_inner;
    double        *inner_back;   /* two inner rows of the widest span, backward */
    size_t         cap_inner_back;
    phmm_position *window;
    size_t         cap_window;
};

phmm_scratch *phmm_scratch_create(void)
{
    return calloc(1, sizeof(phmm_scratch));
}

void phmm_scratch_destroy(phmm_scratch *scratch)
{
    if (!scratch) {
        return;
    }

    free(scratch->places);
    free(scratch->bands);
    free(scratch->rows);
    free(scratch->spans);
    free(scratch->forward);
    free(scratch->backward);
    free(scratch->scale);
    free(scratch->inner);
    free(scratch->inner_back);
    free(scratch->window);
    free(scratch);
}

/* The weights of the gaps, which hold along every reference. A deletion pays its
 * closing with the opening rate of the base it closes at, which that base carries. An
 * insertion pays its closing once, as the first inserted base, and the opening rate at
 * the base it closes at. Each inserted base pays the chance of being the base it is. */
typedef struct {
    double extend_deletion;
    double extend_insertion;  /* per inserted base after the first */
    double open_insertion;    /* the first inserted base */
} gap_weights;

/* Returns the weights of the gaps under rates. */
static gap_weights gap_weights_of(const phmm_rates *rates)
{
    return (gap_weights){
        .extend_deletion  = rates->deletion_to_deletion,
        .extend_insertion = rates->insertion_to_insertion * UNINFORMATIVE,
        .open_insertion   = rates->match_to_insertion * UNINFORMATIVE,
    };
}

/* ------------------------------------------------------------------------ */
/* One read                                                                  */
/* ------------------------------------------------------------------------ */

typedef struct {
    const phmm_model         *model;
    const phmm_rates         *rates;
    const phmm_emissions     *emissions;
    const cm_bam_record      *read;
    const run                *read_runs;  /* of the placed span */
    size_t                    n_read_runs;
    const phmm_length_errors *errors;     /* NULL for a read without them */
    phmm_scratch             *scratch;

    aln_span    span;
    gap_weights gaps;         /* the weights of the gaps under the rates */
    size_t      n_segments;
    int32_t     widest;       /* cells of the widest outer row */
    int32_t     widest_span;  /* cells of the widest inner row */
    size_t      n_spans;      /* inner rows of every inner matrix */
    int32_t     origin;       /* reference position of window value 0 */
    size_t      window_len;
} context;

static size_t n_rows(const context *ctx)
{
    return ctx->n_segments + 1;
}

static inline outer_row *row_of(const context *ctx, size_t k)
{
    return &ctx->scratch->rows[k];
}

/* Returns the forward cells of outer row k, at prefix lo. They are inner row 0 of the
 * next segment, which holds the same prefixes. The last outer row has a row of its
 * own. */
static inline double *outer_cells(const context *ctx, size_t k)
{
    const outer_row *next;

    if (k == ctx->n_segments) {
        return ctx->scratch->forward;
    }

    next = row_of(ctx, k + 1);

    return ctx->scratch->inner + next->inner
         + ctx->scratch->spans[next->first_span].cells;
}

/* Returns the cell at prefix x of a row of outer cells that starts at prefix lo. */
static inline double *outer_at(double *cells, int32_t lo, int32_t x)
{
    return cells + (size_t)(x - lo) * N_OUTER;
}

/* Returns the backward cells of outer row k at prefix x. The pass fills row k - 1 from
 * the inner matrix of row k, so it keeps two rows, and row k lies in the slot of its
 * parity. The cells of a row have the layout of inner cells, so that inner row 0 of a
 * segment can be written straight into them. */
static inline double *backward_at(const context *ctx, size_t k, int32_t x)
{
    const outer_row *row = row_of(ctx, k);

    return ctx->scratch->backward
         + ((k & 1) * (size_t)ctx->widest + (size_t)(x - row->lo)) * N_OUTER;
}

/* Returns the span of inner row j of row. */
static inline const inner_span *span_at(const context *ctx, const outer_row *row,
                                        int32_t j)
{
    return &ctx->scratch->spans[row->first_span + (size_t)j];
}

/* Returns the forward cells of the inner row of segment k that span describes. The
 * first cell is at prefix span->from. */
static inline double *inner_row_at(const context *ctx, size_t k, const inner_span *span)
{
    return ctx->scratch->inner + row_of(ctx, k)->inner + span->cells;
}

/* Returns the backward cells of inner row j. The pass keeps two rows, and row j lies
 * in the slot of its parity. */
static inline double *inner_back_row_at(const context *ctx, int32_t j)
{
    return ctx->scratch->inner_back
         + ((size_t)j & 1) * (size_t)ctx->widest_span * N_INNER;
}

/* Returns the cell at prefix x of a row whose cells start at span->from. */
static inline double *cell_at(double *cells, const inner_span *span, int32_t x)
{
    return cells + (size_t)(x - span->from) * N_INNER;
}

/* Returns the cell at prefix x of a row whose cells start at span->from, for reading. */
static inline const double *cell_in(const double *cells, const inner_span *span,
                                    int32_t x)
{
    return cells + (size_t)(x - span->from) * N_INNER;
}

/* ------------------------------------------------------------------------ */
/* One base against the reference                                            */
/* ------------------------------------------------------------------------ */

/* Returns the comparisons of read base q. Every base of an Ultima read takes one
 * score, since its quality string describes run lengths. */
static const comparisons *comparisons_of(const context *ctx, int32_t q)
{
    const phmm_emissions *emissions = ctx->emissions;
    nuc                   ours      = nuc_from_read(ctx->read->seq, q);

    if (ctx->errors) {
        return &emissions->at_quality[ours][(uint8_t)ctx->errors->phred];
    }
    if (!ctx->read->qual) {
        return &emissions->without_quality[ours];
    }

    return &emissions->at_quality[ours][ctx->read->qual[q]];
}

/* Returns the emission of the comparison of a read base with reference base base. */
static inline double emission_of(const comparisons *o, const phmm_base *base)
{
    return o->intercept[base->code] + o->slope[base->code] * base->modification;
}

/* Returns the part of that emission on which the template differed. */
static inline double modified_of(const comparisons *o, const phmm_base *base)
{
    return o->share[base->code] * base->modification;
}

/* ------------------------------------------------------------------------ */
/* Layout                                                                    */
/* ------------------------------------------------------------------------ */

/* Returns whether read base j of a segment can be inserted. The first base of the read
 * cannot be, nor can a template base the sequencer did not count. */
static bool may_insert(const outer_row *row, size_t k, int32_t j)
{
    return j < row->len && (j >= 1 || k >= 2);
}

/* Returns whether a deletion can open inside the segment after inner row j. */
static bool may_delete_inside(const outer_row *row, int32_t j)
{
    return j >= 1 && j < row->last_inner_row;
}

/* Clamps a prefix to the reference. */
static int32_t clamp_prefix(const context *ctx, int32_t x)
{
    if (x < 0) {
        return 0;
    }
    if ((size_t)x > ctx->model->len) {
        return (int32_t)ctx->model->len;
    }

    return x;
}

/* Writes the band around every point of the CIGAR path of the placed span. */
static void fill_bands(const context *ctx)
{
    size_t placed = (size_t)(ctx->span.end - ctx->span.begin);

    for (size_t prefix = 0; prefix <= placed; prefix++) {
        const aln_place *place = &ctx->scratch->places[prefix];
        band_edges      *band  = &ctx->scratch->bands[prefix];

        band->from = clamp_prefix(ctx, (int32_t)place->first - ctx->model->band);
        band->to   = clamp_prefix(ctx, (int32_t)place->last + ctx->model->band);
    }
}

/* Returns the band around where the CIGAR has reached after j bases of the segment
 * that starts at read base start. */
static const band_edges *band_at(const context *ctx, int32_t start, int32_t j)
{
    size_t prefix = (size_t)(start - ctx->span.begin) + (size_t)j;

    return &ctx->scratch->bands[prefix];
}

/* Writes the window of a row: the band around the prefix the CIGAR has reached after
 * the segment. */
static void lay_out_window(const context *ctx, outer_row *row)
{
    const band_edges *band = band_at(ctx, row->start, row->len);

    row->lo = band->from;
    row->hi = band->to;
}

/* Fills span with inner row j of the segment of row k: its band, its comparisons, and
 * what the row allows. The band is the band around the prefix the CIGAR has reached
 * after j bases of the segment. A row past the bases of the segment takes the band of
 * the last base. */
static void set_inner_span(inner_span *span, const context *ctx, size_t k, int32_t j)
{
    const outer_row  *row  = row_of(ctx, k);
    int32_t           real = j < row->len ? j : row->len;
    const band_edges *band = band_at(ctx, row->start, real);

    span->from        = band->from;
    span->to          = band->to;
    span->uncounted   = j >= row->len;
    span->deletes     = may_delete_inside(row, j);
    span->inserts     = may_insert(row, k, j);
    span->comparisons = span->uncounted ? &ctx->emissions->uncounted
                                        : comparisons_of(ctx, row->start + j);
}

/* Returns how many cells the row holds. */
static int32_t span_width(const inner_span *span)
{
    return span->to - span->from + 1;
}

/* Returns whether a read run continues into the clipped bases on one side, which
 * leaves its length unknown. */
static bool clipped_before(const context *ctx, const run *read_run)
{
    int32_t before = ctx->span.begin - 1;

    return before >= 0 && nuc_is_base(read_run->base)
        && nuc_from_read(ctx->read->seq, before) == read_run->base;
}

static bool clipped_after(const context *ctx, const run *read_run)
{
    int32_t after = ctx->span.end;

    return after < ctx->read->l_qseq && nuc_is_base(read_run->base)
        && nuc_from_read(ctx->read->seq, after) == read_run->base;
}

/* Returns how many bases past the read run the longest template length with any chance
 * reaches, which is how many inner rows the run needs past its bases. */
static int32_t uncounted_rows(const double *length)
{
    for (int offset = PHMM_MAX_LENGTH_ERROR; offset > 0; offset--) {
        if (length[PHMM_MAX_LENGTH_ERROR + offset] > 0.0) {
            return offset;
        }
    }

    return 0;
}

/* Writes the chance of each template length of the run that is segment k. A run whose
 * length the read does not give in full has the read length. */
static void lay_out_run_lengths(const context *ctx, outer_row *row, size_t k)
{
    const run *read_run = &ctx->read_runs[k - 1];
    bool       clipped  = (k == 1 && clipped_before(ctx, read_run))
                       || (k == ctx->n_segments && clipped_after(ctx, read_run));

    if (clipped || !nuc_is_base(read_run->base)) {
        return;
    }

    row->length         = ctx->errors->of_run[k - 1].length;
    row->last_inner_row = row->len + uncounted_rows(row->length);
}

/* The chances of the template lengths of a segment whose length the read gives. */
static const double READ_LENGTH[PHMM_LENGTH_OFFSETS] = {
    [PHMM_MAX_LENGTH_ERROR] = 1.0,
};

/* Points the row at the chance of each template length of segment k. A read without
 * the tags has the read length. */
static void lay_out_lengths(const context *ctx, outer_row *row, size_t k)
{
    row->length         = READ_LENGTH;
    row->last_inner_row = row->len;

    if (ctx->errors) {
        lay_out_run_lengths(ctx, row, k);
    }
}

/* Writes the chance that a one-base run is missing after the run that is segment k,
 * where that run is not the last of the read. */
static void lay_out_missed(const context *ctx, outer_row *row, size_t k)
{
    row->missed       = 0.0;
    row->missed_bases = 0;

    if (ctx->errors && k < ctx->n_segments) {
        row->missed       = ctx->errors->of_run[k - 1].missed;
        row->missed_bases = ctx->errors->of_run[k - 1].missed_bases;
    }
}

/* Returns how many segments the read is cut into. */
static size_t count_segments(const context *ctx)
{
    size_t placed = (size_t)(ctx->span.end - ctx->span.begin);

    if (ctx->errors) {
        return ctx->n_read_runs;
    }

    return (placed + SEGMENT_BASES - 1) / SEGMENT_BASES;
}

/* Writes the read bases of segment k. */
static void lay_out_segment(const context *ctx, outer_row *row, size_t k)
{
    int32_t end;

    if (ctx->errors) {
        row->start = ctx->read_runs[k - 1].start;
        row->len   = ctx->read_runs[k - 1].len;
        return;
    }

    row->start = ctx->span.begin + (int32_t)(k - 1) * SEGMENT_BASES;
    end        = row->start + SEGMENT_BASES;
    row->len   = (end < ctx->span.end ? end : ctx->span.end) - row->start;
}

/* Lays out every outer row. */
static void lay_out_rows(context *ctx)
{
    outer_row        *rows   = ctx->scratch->rows;
    const band_edges *start;
    int32_t           widest = 0;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        lay_out_segment(ctx, &rows[k], k);
        lay_out_window(ctx, &rows[k]);
        lay_out_lengths(ctx, &rows[k], k);
        lay_out_missed(ctx, &rows[k], k);
    }

    start = band_at(ctx, rows[1].start, 0);

    rows[0].lo             = start->from;
    rows[0].hi             = start->to;
    rows[0].last_inner_row = 0;

    ctx->n_spans = 0;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        rows[k].first_span = ctx->n_spans;
        ctx->n_spans      += (size_t)rows[k].last_inner_row + 1;
    }

    for (size_t k = 0; k < n_rows(ctx); k++) {
        int32_t width = rows[k].hi - rows[k].lo + 1;

        if (width > widest) {
            widest = width;
        }
    }

    ctx->widest = widest;
}

/* Writes the span of every inner row of every inner matrix, and where each lies.
 * Returns the doubles the inner matrices need in all. */
static size_t lay_out_spans(context *ctx)
{
    size_t  inner  = 0;
    int32_t widest = 0;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        outer_row *row   = row_of(ctx, k);
        size_t     cells = 0;

        row->inner = inner;

        for (int32_t j = 0; j <= row->last_inner_row; j++) {
            inner_span *span = &ctx->scratch->spans[row->first_span + (size_t)j];

            set_inner_span(span, ctx, k, j);
            span->cells = cells;
            cells      += (size_t)span_width(span) * N_INNER;

            if (span_width(span) > widest) {
                widest = span_width(span);
            }
        }

        inner += cells;
    }

    ctx->widest_span = widest;

    return inner;
}

static void lay_out_output(context *ctx)
{
    const outer_row *rows = ctx->scratch->rows;
    int32_t           lo   = rows[0].lo;
    int32_t           hi   = rows[0].hi;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        if (rows[k].lo < lo) {
            lo = rows[k].lo;
        }
        if (rows[k].hi > hi) {
            hi = rows[k].hi;
        }
    }

    ctx->origin     = lo - 1;
    ctx->window_len = (size_t)(hi - ctx->origin) + 1;
}

/* Returns the doubles one outer row holds. */
static size_t outer_row_doubles(const context *ctx)
{
    return (size_t)ctx->widest * N_OUTER;
}

/* Grows every buffer the passes need. The outer buffers are a pair, so they are grown
 * together. */
static int grow_buffers(context *ctx, size_t inner)
{
    phmm_scratch *s     = ctx->scratch;
    size_t          outer = 2 * outer_row_doubles(ctx);
    size_t          back  = 2 * (size_t)ctx->widest_span * N_INNER;

    if (outer > s->cap_outer) {
        double *forward  = realloc(s->forward, outer * sizeof *forward);
        double *backward = realloc(s->backward, outer * sizeof *backward);

        if (forward) {
            s->forward = forward;
        }
        if (backward) {
            s->backward = backward;
        }
        if (!forward || !backward) {
            return -1;
        }

        s->cap_outer = outer;
    }

    if (grow((void **)&s->scale, &s->cap_scale, n_rows(ctx), sizeof *s->scale) < 0
        || grow((void **)&s->inner, &s->cap_inner, inner, sizeof *s->inner) < 0
        || grow((void **)&s->inner_back, &s->cap_inner_back, back,
                sizeof *s->inner_back) < 0
        || grow((void **)&s->window, &s->cap_window, ctx->window_len,
                sizeof *s->window) < 0) {
        return -1;
    }

    return 0;
}

static phmm_status prepare(context *ctx)
{
    phmm_scratch      *s    = ctx->scratch;
    const cm_bam_record *read = ctx->read;
    size_t               inner;

    if (grow((void **)&s->places, &s->cap_places, (size_t)read->l_qseq + 1,
             sizeof *s->places) < 0
        || grow((void **)&s->bands, &s->cap_bands, (size_t)read->l_qseq + 1,
                sizeof *s->bands) < 0) {
        return PHMM_NO_MEMORY;
    }

    ctx->span = aln_places(read, s->places);

    if (ctx->span.begin >= ctx->span.end) {
        return PHMM_NO_PATH;
    }

    fill_bands(ctx);
    ctx->n_segments = count_segments(ctx);

    if (grow((void **)&s->rows, &s->cap_rows, n_rows(ctx), sizeof *s->rows) < 0) {
        return PHMM_NO_MEMORY;
    }

    lay_out_rows(ctx);
    lay_out_output(ctx);

    if (grow((void **)&s->spans, &s->cap_spans, ctx->n_spans, sizeof *s->spans) < 0) {
        return PHMM_NO_MEMORY;
    }

    inner = lay_out_spans(ctx);

    if (grow_buffers(ctx, inner) < 0) {
        return PHMM_NO_MEMORY;
    }

    memset(s->window, 0, ctx->window_len * sizeof *s->window);

    return PHMM_OK;
}

/* ------------------------------------------------------------------------ */
/* The inner matrix of a segment                                             */
/* ------------------------------------------------------------------------ */

/* The inner row at which a segment ends when the template is offset bases longer than
 * the read, and the chance of that length. The chance is zero where the inner matrix
 * has no such row. */
typedef struct {
    int32_t jend;
    double  chance;
} template_length;

/* Returns the template length of the segment of row at offset. */
static template_length template_length_at(const outer_row *row, int offset)
{
    template_length length = {
        .jend   = row->len + offset,
        .chance = row->length[offset + PHMM_MAX_LENGTH_ERROR],
    };

    if (length.jend < 0 || length.jend > row->last_inner_row) {
        length.chance = 0.0;
    }

    return length;
}

/* What the cells of an inner row are worth at the end of the segment, in the backward
 * pass: the chance of the template length that ends the segment at the row, divided
 * by the scale of the segment, which is zero where no length ends there; the outer
 * backward cells of the row of the segment, from prefix lo to hi; and what can follow
 * the end. Nothing follows the last segment. */
typedef struct {
    double        chance;
    const double *cells;
    int32_t       lo, hi;
    bool          deletes;   /* a deletion can open after the segment */
    bool          inserts;   /* an insertion left open can end at the row */
} segment_end;

/* Fills end with what inner row j of segment k is worth at the end of the segment.
 * restore undoes the scaling of the segment and is folded into the chance. Every
 * backward value of the segment then carries restore, which the posteriors need and
 * which the outer cells of row k - 1 are defined with. */
static inline void set_segment_end(segment_end *end, const context *ctx, size_t k,
                                   int32_t j, double restore)
{
    const outer_row *row    = row_of(ctx, k);
    int              offset = j - row->len;
    bool             last   = k == ctx->n_segments;

    end->chance  = 0.0;
    end->cells   = backward_at(ctx, k, row->lo);
    end->lo      = row->lo;
    end->hi      = row->hi;
    end->deletes = !last;
    end->inserts = !last && offset <= 0;

    if (offset >= -PHMM_MAX_LENGTH_ERROR && offset <= PHMM_MAX_LENGTH_ERROR) {
        end->chance = template_length_at(row, offset).chance * restore;
    }
}

/* Adds what the paired and the open insertion states of a cell at prefix x are worth
 * at the end of the segment to paired and inserted. The paired state is worth the
 * chance of ending there times what the outer row is worth closed at x, plus what it
 * is worth with a deletion open at x + 1. The open insertion state is worth that
 * chance times what the outer row is worth with an insertion open at x. */
static inline void add_end_worth(const segment_end *end, int32_t x, double *paired,
                                 double *inserted)
{
    const double *cell;
    double        worth;

    if (end->chance <= 0.0 || x < end->lo || x > end->hi) {
        return;
    }

    cell  = end->cells + (size_t)(x - end->lo) * N_OUTER;
    worth = cell[CLOSED];

    if (end->deletes && x < end->hi) {
        worth += cell[N_OUTER + DELETION];
    }

    *paired += worth * end->chance;

    if (end->inserts) {
        *inserted += cell[INSERTION] * end->chance;
    }
}

/* What the step from inner row j to inner row j + 1 reads: the two spans, and whether
 * a pairing from row j is the read's 5'-most, which follows the end state. The
 * backward pass also reads what row j is worth at the end of the segment, and where
 * the backward values of row j go. They go to the slot of the row's parity, or, for
 * row 0, to the outer backward cells of row k - 1. */
typedef struct {
    const inner_span  *span;      /* of row j */
    const inner_span  *next;      /* of row j + 1 */
    double            *line;      /* where the backward values of row j go */
    segment_end        end;
    bool               ends;      /* a pairing from row j is the read's 5'-most */
} inner_step;

/* Fills s with the step from inner row j of segment k, for the forward pass. The
 * struct is filled in place, since a struct returned by value is copied through the
 * stack on every inner row. */
static inline void set_inner_step(inner_step *s, const context *ctx, size_t k, int32_t j)
{
    const outer_row *row = row_of(ctx, k);

    s->span = span_at(ctx, row, j);
    s->next = span_at(ctx, row, j + 1);
    s->ends = k == 1 && j == 0;
}

/* Fills s with the step from inner row j of segment k, for the backward pass. */
static inline void set_backward_step(inner_step *s, const context *ctx, size_t k,
                                     int32_t j, double restore)
{
    const outer_row *above = row_of(ctx, k - 1);

    set_inner_step(s, ctx, k, j);
    s->line = j == 0 ? backward_at(ctx, k - 1, above->lo) : inner_back_row_at(ctx, j);
    set_segment_end(&s->end, ctx, k, j, restore);
}

static inline int32_t max_prefix(int32_t a, int32_t b)
{
    return a > b ? a : b;
}

static inline int32_t min_prefix(int32_t a, int32_t b)
{
    return a < b ? a : b;
}

/* Returns the weight a pairing of base b pays after the paired state of the cell before:
 * nothing more for the read's 5'-most pairing, whose cell before is the termination
 * state, and the weight after no gap otherwise. */
static inline double after_paired(const phmm_base *b, bool first)
{
    return first ? 1.0 : b->after_no_gap;
}

/* Returns the weight of the paths from cell from into the pairing of reference base b,
 * which lands one prefix further along the next row. */
static inline double pairing_into(const double *from, const phmm_base *b, bool first)
{
    return from[PAIRED] * after_paired(b, first) + from[GAP_DELETION] * b->after_deletion
         + from[GAP_INSERTION] * b->after_insertion;
}

/* Fills inner row j + 1 from inner row j, whose read base is uncounted or not, and
 * whose pairings are the read's 5'-most or not. Each cell of the next row is written
 * once: its pairing comes from the cell one prefix back, its insertion from the cell
 * at the same prefix, and its open deletion from the cell before it in its own row. A
 * cell whose source lies outside the band of the row above takes nothing from it. */
ALWAYS_INLINE void forward_cells(const context *ctx, size_t k, const inner_step *s,
                                 bool uncounted, bool first)
{
    const phmm_base        *bases    = ctx->model->bases;
    const inner_span       *span     = s->span;
    const inner_span       *after    = s->next;
    const comparisons      *o        = span->comparisons;
    gap_weights             w        = ctx->gaps;
    const double           *line     = inner_row_at(ctx, k, span);
    double                 *next     = inner_row_at(ctx, k, after);
    int32_t                 from     = span->from;
    int32_t                 to       = span->to;
    int32_t                 next_from = after->from;
    int32_t                 next_to  = after->to;
    bool                    deletes  = after->deletes;
    int32_t                 pair_lo  = max_prefix(next_from, from + 1);
    int32_t                 pair_hi  = min_prefix(next_to, to + 1);
    int32_t                 ins_hi   = !uncounted && span->inserts
                                     ? min_prefix(after->to, to) : -1;
    double                  deletion = 0.0;
    int32_t                 x        = next_from;

    for (; x < pair_lo; x++) {
        const double *same = line + (size_t)(x - from) * N_INNER;
        double       *cell = next + (size_t)(x - next_from) * N_INNER;
        double        inserted = x <= ins_hi
            ? same[PAIRED] * w.open_insertion + same[GAP_INSERTION] * w.extend_insertion
            : 0.0;

        cell[PAIRED]        = 0.0;
        cell[GAP_DELETION]  = deletion;
        cell[GAP_INSERTION] = inserted;
        deletion            = deletion * w.extend_deletion;
    }

    for (; x <= pair_hi; x++) {
        const double    *back   = line + (size_t)(x - 1 - from) * N_INNER;
        const phmm_base *b      = &bases[x - 1];
        double          *cell   = next + (size_t)(x - next_from) * N_INNER;
        double           paired = pairing_into(back, b, first)
                                * (uncounted ? 1.0 : emission_of(o, b));
        double           inserted = x <= ins_hi
            ? back[N_INNER + PAIRED] * w.open_insertion
              + back[N_INNER + GAP_INSERTION] * w.extend_insertion
            : 0.0;

        cell[PAIRED]        = paired;
        cell[GAP_DELETION]  = deletion;
        cell[GAP_INSERTION] = inserted;
        deletion            = deletion * w.extend_deletion + (deletes ? paired : 0.0);
    }

    for (; x <= next_to; x++) {
        double *cell = next + (size_t)(x - next_from) * N_INNER;

        cell[PAIRED]        = 0.0;
        cell[GAP_DELETION]  = deletion;
        cell[GAP_INSERTION] = 0.0;
        deletion            = deletion * w.extend_deletion;
    }
}

/* Fills inner row j + 1 from inner row j. A template base the sequencer did not count
 * emits one and cannot be inserted, so its step is the step without the emission and
 * the insertion. The read's first base is a base of its first segment, so it is
 * counted. */
static void inner_forward_step(const context *ctx, size_t k, const inner_step *s)
{
    if (s->ends) {
        forward_cells(ctx, k, s, false, true);
    } else if (s->span->uncounted) {
        forward_cells(ctx, k, s, true, false);
    } else {
        forward_cells(ctx, k, s, false, false);
    }
}

/* Fills the inner matrix of segment k from its inner row 0. */
static void inner_forward(const context *ctx, size_t k)
{
    const outer_row *row = row_of(ctx, k);

    for (int32_t j = 0; j < row->last_inner_row; j++) {
        inner_step s;

        set_inner_step(&s, ctx, k, j);
        inner_forward_step(ctx, k, &s);
    }
}

/* Returns whether the sequencer can have missed a one-base run after the segment of
 * row k. */
static bool missed_run_possible(const context *ctx, size_t k)
{
    const outer_row *row = row_of(ctx, k);

    return k != ctx->n_segments && row->missed > 0.0 && row->missed_bases != 0;
}

/* Returns whether reference base x is a run of one base that the sequencer can have
 * missed after the segment of row. */
static inline bool missed_run_at(const context *ctx, const outer_row *row, int32_t x)
{
    const phmm_model *model = ctx->model;

    return (size_t)x < model->len && (row->missed_bases & model->bases[x].run_bit) != 0;
}

/* One inner row at which the segment can end, as the exits read it: the chance of the
 * template length that ends there, the cell of the row at the first prefix of the
 * outer row, the last prefix at which the row has a cell, and whether an insertion
 * left open can end there. An insertion can end there where the template is no longer
 * than the read. The band of the row begins at or before the outer row, since the
 * outer row lies around the end of the segment. */
typedef struct {
    double        chance;
    const double *at_lo;
    int32_t       to;
    bool          inserts;
} exit_row;

/* Writes the rows at which the segment of row k can end and reach its outer row, and
 * returns their count. */
static size_t exit_rows_of(const context *ctx, size_t k, exit_row *exits)
{
    const outer_row *row = row_of(ctx, k);
    size_t           n   = 0;

    for (int offset = -PHMM_MAX_LENGTH_ERROR; offset <= PHMM_MAX_LENGTH_ERROR; offset++) {
        template_length   length = template_length_at(row, offset);
        const inner_span *span;

        if (length.chance <= 0.0) {
            continue;
        }

        span = span_at(ctx, row, length.jend);

        if (span->to < row->lo) {
            continue;
        }

        exits[n].chance  = length.chance;
        exits[n].at_lo   = cell_in(inner_row_at(ctx, k, span), span, row->lo);
        exits[n].to      = span->to;
        exits[n].inserts = offset <= 0;
        n++;
    }

    return n;
}

/* Adds what one exit row holds at prefix x, weighted by its chance, to closed and
 * inserted. cell is the offset in doubles of x from the first prefix of the outer
 * row. A row with no cell at x adds nothing. */
static inline void exit_at(const exit_row *exit, int32_t x, size_t cell, double *closed,
                           double *inserted)
{
    const double *at = exit->at_lo + cell;

    if (x > exit->to) {
        return;
    }

    *closed += at[PAIRED] * exit->chance;

    if (exit->inserts) {
        *inserted += at[GAP_INSERTION] * exit->chance;
    }
}

/* Scales the three states of an outer cell by factor. */
static inline void scale_outer_cell(double *cell, double factor)
{
    cell[CLOSED]    *= factor;
    cell[DELETION]  *= factor;
    cell[INSERTION] *= factor;
}

/* Returns the sum of the three states of an outer cell. */
static inline double outer_cell_total(const double *cell)
{
    return cell[CLOSED] + cell[DELETION] + cell[INSERTION];
}

/* Writes the outer cells of row k from the inner rows at which the segment can end,
 * each weighted by its chance, and returns their total. A deletion that opens after
 * the segment and an insertion left open end in their own states, which the last row
 * cannot hold. The deletion open at a prefix comes from the pairings closed at the
 * prefix before, and extends from the deletion open there. Every cell then holds the
 * whole of its open deletion state as the next segment reads it. */
static double fill_outer_row(const context *ctx, size_t k)
{
    const outer_row *row      = row_of(ctx, k);
    double          *cells    = outer_cells(ctx, k);
    bool             last     = k == ctx->n_segments;
    double           extend   = ctx->rates->deletion_to_deletion;
    exit_row         exits[PHMM_LENGTH_OFFSETS];
    size_t           n        = exit_rows_of(ctx, k, exits);
    double           deletion = 0.0;
    double           total    = 0.0;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        double *to       = outer_at(cells, row->lo, x);
        size_t  cell     = (size_t)(x - row->lo) * N_INNER;
        double  closed   = 0.0;
        double  inserted = 0.0;

        for (size_t e = 0; e < n; e++) {
            exit_at(&exits[e], x, cell, &closed, &inserted);
        }

        to[CLOSED]    = closed;
        to[DELETION]  = last ? 0.0 : deletion;
        to[INSERTION] = last ? 0.0 : inserted;
        total        += outer_cell_total(to);
        deletion      = deletion * extend + closed;
    }

    return total;
}

/* Passes over the one-base runs the sequencer can have missed after the row's
 * segment, and returns what the pass adds to the total of the row. A path that skips a
 * run takes the chance of the missed run, and every other path takes the chance that
 * no run was missed. skipped is what the paths that skip the run at x - 1 carry to x. */
static double missed_run_forward(const context *ctx, size_t k)
{
    const outer_row *row     = row_of(ctx, k);
    double          *cells   = outer_cells(ctx, k);
    double           keep    = 1.0 - row->missed;
    double           skipped = 0.0;
    double           added   = 0.0;

    if (!missed_run_possible(ctx, k)) {
        return 0.0;
    }

    for (int32_t x = row->lo; x <= row->hi; x++) {
        double *cell   = outer_at(cells, row->lo, x);
        double  closed = cell[CLOSED] + skipped;
        double  before = outer_cell_total(cell);

        scale_outer_cell(cell, keep);
        cell[CLOSED] += skipped;
        added        += outer_cell_total(cell) - before;
        skipped       = x < row->hi && missed_run_at(ctx, row, x)
                      ? closed * row->missed * ctx->model->bases[x].after_no_gap
                      : 0.0;
    }

    return added;
}

/* Maps the total of a row to a status. A negative, infinite, or NaN total is a fault.
 * A zero or subnormal total means the band admits no path. */
static phmm_status total_status(double total)
{
    if (!(total >= 0.0) || total > DBL_MAX) {
        return PHMM_UNSOUND;
    }
    if (total < DBL_MIN) {
        return PHMM_NO_PATH;
    }

    return PHMM_OK;
}

/* Keeps the scale of row k. The scale is the total of the row where the row is the
 * last or where the total has left the range the rows are kept in, and the row is
 * then divided by it. Elsewhere the scale is one. */
static phmm_status scale_row(const context *ctx, size_t k, double total)
{
    const outer_row *row    = row_of(ctx, k);
    double          *cells  = outer_cells(ctx, k);
    phmm_status      status = total_status(total);
    double           share;

    if (status != PHMM_OK) {
        return status;
    }

    if (k == ctx->n_segments) {
        ctx->scratch->scale[k] = total;
        return PHMM_OK;
    }

    if (total >= SCALE_FLOOR && total <= SCALE_CEILING) {
        ctx->scratch->scale[k] = 1.0;
        return PHMM_OK;
    }

    ctx->scratch->scale[k] = total;
    share                  = 1.0 / total;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        scale_outer_cell(outer_at(cells, row->lo, x), share);
    }

    return PHMM_OK;
}

/* Fills row 0, the end state, which reaches the 5'-most pairing of the read at
 * any base of the first window with the termination rate of that base, and returns its
 * total. */
static double forward_first_row(const context *ctx)
{
    const outer_row *row   = row_of(ctx, 0);
    double          *cells = outer_cells(ctx, 0);
    double           total = 0.0;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        double *cell = outer_at(cells, row->lo, x);

        cell[CLOSED]    = rate_at(ctx->rates->termination, ctx->model->len, x);
        cell[DELETION]  = 0.0;
        cell[INSERTION] = 0.0;
        total          += cell[CLOSED];
    }

    return total;
}

static phmm_status forward(const context *ctx)
{
    phmm_status status = scale_row(ctx, 0, forward_first_row(ctx));

    for (size_t k = 1; k < n_rows(ctx) && status == PHMM_OK; k++) {
        double total;

        inner_forward(ctx, k);
        total  = fill_outer_row(ctx, k);
        total += missed_run_forward(ctx, k);
        status = scale_row(ctx, k, total);
    }

    return status;
}

/* ------------------------------------------------------------------------ */
/* Backward, with the posteriors                                             */
/* ------------------------------------------------------------------------ */

/* Returns the window position of reference base x. */
static phmm_position *landing(const context *ctx, int32_t x)
{
    return ctx->scratch->window + (size_t)(x - ctx->origin);
}

/* Passes back over the one-base runs the sequencer can have missed after the row's
 * segment, as the transpose of the forward pass. A cell then holds what an exit of the
 * segment at its prefix is worth. onward is what a path that skips the run at x is
 * worth at x + 1. */
static void missed_run_backward(const context *ctx, size_t k)
{
    const outer_row *row    = row_of(ctx, k);
    double           keep   = 1.0 - row->missed;
    double           onward = 0.0;

    if (!missed_run_possible(ctx, k)) {
        return;
    }

    for (int32_t x = row->hi; x >= row->lo; x--) {
        double *cell  = backward_at(ctx, k, x);
        double  skips = x < row->hi && missed_run_at(ctx, row, x)
                      ? row->missed * ctx->model->bases[x].after_no_gap * onward
                      : 0.0;

        onward = cell[CLOSED] + skips;
        scale_outer_cell(cell, keep);
        cell[CLOSED] += skips;
    }
}

/* Writes the backward values of the last inner row of segment k. Each cell is worth
 * what it is worth at the end of the segment. */
static void start_inner_backward(const context *ctx, size_t k, double restore)
{
    const outer_row  *row   = row_of(ctx, k);
    const inner_span *span  = span_at(ctx, row, row->last_inner_row);
    double           *cells = inner_back_row_at(ctx, row->last_inner_row);
    segment_end       end;

    set_segment_end(&end, ctx, k, row->last_inner_row, restore);

    for (int32_t x = span->from; x <= span->to; x++) {
        double *cell     = cell_at(cells, span, x);
        double  paired   = 0.0;
        double  inserted = 0.0;

        add_end_worth(&end, x, &paired, &inserted);

        cell[PAIRED]        = paired;
        cell[GAP_DELETION]  = 0.0;
        cell[GAP_INSERTION] = inserted;
    }
}

/* The posterior of one pairing of a reference base: the paths through it, the paths
 * through it on which the template differed, and the paths that close a deletion or an
 * insertion at it. */
typedef struct {
    double covered;
    double modified;
    double deleted;
    double inserted;
} pairing_posterior;

/* Accumulates one pairing of reference base x into the window. The deletion closed at
 * x belongs to the base before, so it is carried to the next cell, which the pass
 * reaches next, and the deletion carried from the cell after lands here. A pairing of
 * a template base the sequencer did not count observed no read base, so it covers
 * nothing and is not counted. Returns the deletion to carry. */
static inline double accumulate_pairing(const context *ctx, int32_t x,
                                        const pairing_posterior *p, bool counted,
                                        bool ends, double carried)
{
    phmm_position *at = landing(ctx, x);

    if (counted) {
        at->coverage   += p->covered;
        at->mismatches += p->modified;
    }

    at->insertions += p->inserted;
    at->deletions  += carried;

    if (ends) {
        at->ends += p->covered - p->deleted - p->inserted;
    }

    return p->deleted;
}

/* Writes a backward cell at prefix x whose read base pairs with nothing. The cell is
 * worth the end of the segment, plus the insertion of its read base into gap, which is
 * the cell of the next row at the same prefix or NULL where that cell does not exist,
 * plus the open deletion of the cell after it in its own row. A deletion opens there
 * with weight deletes, which is one or zero. Returns the open deletion state of the
 * cell. */
static inline double backward_gap_cell(double *cell, const segment_end *end, int32_t x,
                                       const double *gap, const gap_weights *w,
                                       double deletes, double deletion)
{
    double paired   = deletes * deletion;
    double deleted  = deletion * w->extend_deletion;
    double inserted = 0.0;

    add_end_worth(end, x, &paired, &inserted);

    if (gap) {
        paired   += w->open_insertion * gap[GAP_INSERTION];
        inserted += w->extend_insertion * gap[GAP_INSERTION];
    }

    cell[PAIRED]        = paired;
    cell[GAP_DELETION]  = deleted;
    cell[GAP_INSERTION] = inserted;

    return deleted;
}

/* Writes the backward inner row j of segment k, whose read base is uncounted or not.
 * Each cell is worth the end of the segment plus what it reaches through inner row
 * j + 1, and each pairing is accumulated with its posterior. A cell reads the cell one
 * prefix on in the next row, which its pairing lands in, the cell at the same prefix,
 * which its insertion lands in, and the cell after it in its own row, which its open
 * deletion extends into. A cell whose target lies outside the next row takes nothing
 * from it. The posterior of a template difference at a pairing is the posterior of the
 * pairing times the modified part of its emission over the whole emission, so the
 * emission divides out. The backward values carry the inverse of the scale of the
 * segment, so the posteriors need no further factor. */
ALWAYS_INLINE void backward_cells(const context *ctx, size_t k, const inner_step *s,
                                  int32_t j, bool uncounted, bool first)
{
    const phmm_base        *bases    = ctx->model->bases;
    const inner_span       *span     = s->span;
    const inner_span       *after    = s->next;
    const comparisons      *o        = span->comparisons;
    const segment_end       end      = s->end;
    gap_weights             w        = ctx->gaps;
    double                  deletes  = span->deletes ? 1.0 : 0.0;
    const double           *forward  = inner_row_at(ctx, k, span);
    double                 *line     = s->line;
    const double           *next     = inner_back_row_at(ctx, j + 1);
    int32_t                 from     = span->from;
    int32_t                 to       = span->to;
    int32_t                 next_from = after->from;
    int32_t                 pair_lo  = max_prefix(from, next_from - 1);
    int32_t                 pair_hi  = min_prefix(to, after->to - 1);
    int32_t                 ins_lo   = !uncounted && span->inserts
                                     ? max_prefix(from, next_from) : INT32_MAX;
    double                  deletion = 0.0;  /* of the cell to the right */
    double                  carried  = 0.0;  /* deletion posterior of the cell after */
    int32_t                 x        = to;

    for (; x > pair_hi; x--) {
        const double *gap = x >= ins_lo ? next + (size_t)(x - next_from) * N_INNER : NULL;

        deletion = backward_gap_cell(line + (size_t)(x - from) * N_INNER, &end, x, gap,
                                     &w, deletes, deletion);
    }

    for (; x >= pair_lo; x--) {
        const phmm_base   *b        = &bases[x];
        const double      *fwd      = forward + (size_t)(x - from) * N_INNER;
        const double      *on       = next + (size_t)(x + 1 - next_from) * N_INNER;
        double            *cell     = line + (size_t)(x - from) * N_INNER;
        double             emission = uncounted ? 1.0 : emission_of(o, b);
        double             reach    = on[PAIRED] * emission;
        double             no_gap   = fwd[PAIRED] * after_paired(b, first);
        double             deleted  = fwd[GAP_DELETION] * b->after_deletion;
        double             inserted = fwd[GAP_INSERTION] * b->after_insertion;
        double             through  = (no_gap + deleted + inserted) * on[PAIRED];
        pairing_posterior  p        = {
            .covered  = through * emission,
            .modified = uncounted ? 0.0 : through * modified_of(o, b),
            .deleted  = deleted * reach,
            .inserted = inserted * reach,
        };
        double             to_paired   = after_paired(b, first) * reach
                                       + deletes * deletion;
        double             to_deleted  = b->after_deletion * reach
                                       + deletion * w.extend_deletion;
        double             to_inserted = b->after_insertion * reach;

        add_end_worth(&end, x, &to_paired, &to_inserted);

        if (x >= ins_lo) {
            double gap = on[-N_INNER + GAP_INSERTION];

            to_paired   += w.open_insertion * gap;
            to_inserted += w.extend_insertion * gap;
        }

        carried = accumulate_pairing(ctx, x, &p, !uncounted, first, carried);

        cell[PAIRED]        = to_paired;
        cell[GAP_DELETION]  = to_deleted;
        cell[GAP_INSERTION] = to_inserted;
        deletion            = to_deleted;
    }

    landing(ctx, x)->deletions += carried;

    for (; x >= from; x--) {
        const double *gap = x >= ins_lo ? next + (size_t)(x - next_from) * N_INNER : NULL;

        deletion = backward_gap_cell(line + (size_t)(x - from) * N_INNER, &end, x, gap,
                                     &w, deletes, deletion);
    }
}

/* Writes the backward inner row j of segment k. A template base the sequencer did not
 * count emits one, cannot be inserted, and covers nothing, so its step is the step
 * without the emission, the insertion, and the coverage. The read's first base is
 * counted. */
static void inner_backward_step(const context *ctx, size_t k, const inner_step *s,
                                int32_t j)
{
    if (s->ends) {
        backward_cells(ctx, k, s, j, false, true);
    } else if (s->span->uncounted) {
        backward_cells(ctx, k, s, j, true, false);
    } else {
        backward_cells(ctx, k, s, j, false, false);
    }
}

/* Fills the backward inner matrix of segment k row by row from its last. Its row 0
 * is written as the outer cells of row k - 1. */
static void inner_backward(const context *ctx, size_t k)
{
    const outer_row *row     = row_of(ctx, k);
    double           restore = 1.0 / ctx->scratch->scale[k];

    start_inner_backward(ctx, k, restore);

    for (int32_t j = row->last_inner_row - 1; j >= 0; j--) {
        inner_step s;

        set_backward_step(&s, ctx, k, j, restore);
        inner_backward_step(ctx, k, &s, j);
    }
}

static void backward_last_row(const context *ctx)
{
    size_t            k   = ctx->n_segments;
    const outer_row *row = row_of(ctx, k);

    for (int32_t x = row->lo; x <= row->hi; x++) {
        double *cell = backward_at(ctx, k, x);

        cell[CLOSED]    = 1.0;
        cell[DELETION]  = 0.0;
        cell[INSERTION] = 0.0;
    }
}

/* Returns the first row's total of forward times backward, which is one when the two
 * passes describe the same set of paths. Inner row 0 of the first inner matrix holds
 * the forward values of row 0. */
static double agreement(const context *ctx)
{
    const outer_row  *row   = row_of(ctx, 0);
    const inner_span *span  = span_at(ctx, row_of(ctx, 1), 0);
    const double     *cells = inner_row_at(ctx, 1, span);
    double            total = 0.0;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        total += cell_in(cells, span, x)[PAIRED]
               * backward_at(ctx, 0, x)[CLOSED];
    }

    return total;
}

/* A total that is not finite means the backward pass overflowed. Every forward row is
 * scaled by its own total, but a backward value is the probability of the rest of the
 * read from a cell relative to the probability of the read itself, and on a read whose
 * every path is astronomically unlikely, as under rates set to zero and a band too
 * narrow for its gaps, that ratio leaves the range of a double. The read has no path
 * worth scoring, and is rejected as one with none. A finite total away from one means
 * the two passes described different sets of paths, which is a fault. */
static phmm_status passes_agree(const context *ctx)
{
    double total = agreement(ctx);

    if (!isfinite(total)) {
        return PHMM_NO_PATH;
    }

    return fabs(total - 1.0) < NORMALIZATION_TOLERANCE ? PHMM_OK : PHMM_UNSOUND;
}

static phmm_status backward(const context *ctx)
{
    backward_last_row(ctx);

    for (size_t k = ctx->n_segments; k >= 1; k--) {
        missed_run_backward(ctx, k);
        inner_backward(ctx, k);
    }

    return passes_agree(ctx);
}

/* ------------------------------------------------------------------------ */
/* Entry                                                                     */
/* ------------------------------------------------------------------------ */

void phmm_window_bounds(const phmm_window *window, size_t len,
                        size_t *from, size_t *to)
{
    hts_pos_t low  = window->origin < 0 ? -window->origin : 0;
    hts_pos_t high = (hts_pos_t)len - window->origin;

    *from = (size_t)low;
    *to   = high < 0 ? 0 : (size_t)high;

    if (*to > window->len) {
        *to = window->len;
    }
    if (*from > *to) {
        *from = *to;
    }
}

phmm_status phmm_run(const phmm_model *model, const phmm_emissions *emissions,
                     const cm_bam_record *read, const run *runs, size_t n_runs,
                     const phmm_length_errors *errors, phmm_scratch *scratch,
                     phmm_window *out)
{
    context ctx = {
        .model       = model,
        .rates       = model->rates,
        .emissions   = emissions,
        .read        = read,
        .read_runs   = runs,
        .n_read_runs = n_runs,
        .errors      = errors,
        .scratch     = scratch,
        .gaps        = gap_weights_of(model->rates),
    };
    phmm_status status;

    status = prepare(&ctx);

    if (status != PHMM_OK) {
        return status;
    }

    status = forward(&ctx);

    if (status != PHMM_OK) {
        return status;
    }

    status = backward(&ctx);

    if (status != PHMM_OK) {
        return status;
    }

    out->origin = ctx.origin;
    out->len    = ctx.window_len;
    out->at     = scratch->window;

    return PHMM_OK;
}
