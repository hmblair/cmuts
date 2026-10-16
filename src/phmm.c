/* phmm.c -- a banded pair HMM.
 *
 * The read is cut into segments. A read whose sequencer reports run length errors has
 * one segment per run of one base, since the errors enter where a run ends. A read
 * without them has segments of up to SEGMENT_BASES bases, since its cuts carry nothing.
 *
 * The matrix has one outer row per segment, and row 0 before the first. A cell of an
 * outer row is a prefix length of the reference, and holds three states: the segment
 * is complete and the last event was a pairing, a deletion is open after it, or an
 * insertion is open after it.
 *
 * Between two outer rows lies the inner matrix of a segment, which pairs the bases of
 * the segment with the reference one at a time. Its rows are the bases of the segment,
 * and up to two more for template bases that an Ultima sequencer did not count. Each
 * row holds the prefix lengths within the band of its base, and nothing else. Its
 * states are the three of a pair HMM over bases, with the weights of that model: a
 * pairing pays the weight of its base after no gap, or the opening rate of the gap it
 * closes, times the emission of the comparison, and a gap pays its extension rate per
 * base after the first.
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

/* The most bases a segment of a read without run length errors holds. The paths
 * through this many bases stay far above the smallest double. */
#define SEGMENT_BASES 32

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

/* Returns the weight a pairing of base b pays when no gap precedes it. */
static double after_no_gap_at(const phmm_rates *rates, size_t len, hts_pos_t b)
{
    return 1.0 - rate_at(rates->open_insertion, len, b)
               - rate_at(rates->open_deletion, len, b - 1);
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
        .is_one_base_run = false,
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
        if (r->len == 1) {
            model->bases[r->start].is_one_base_run = true;
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
    double  length[PHMM_LENGTH_OFFSETS];  /* chance the template run is s + offset long */
    double  missed;                  /* chance a one-base run after it was not read */
    uint8_t missed_bases;            /* the bases that run can be, as bits by nuc */
    size_t  inner;                   /* offset of the inner matrix in scratch */
    size_t  first_span;              /* index of the span of inner row 0 */
} outer_row;

/* One inner row: the prefixes of its band, from from to to, and where its cells lie,
 * as an offset in doubles from the start of the inner matrix of its segment.
 *
 * The cells run from from to end. Most rows end at to. A row whose band ends where
 * the band of the row before ends, which happens after an inserted read base, has one
 * cell more, at to + 1. The pairing out of the last cell of the row before lands in
 * that cell. The cell is read only where a run ends at this row, one or two bases
 * before its last read base. */
typedef struct {
    int32_t from;
    int32_t to;
    int32_t end;
    size_t  cells;
} inner_span;

struct phmm_scratch {
    aln_place     *places;
    size_t         cap_places;
    outer_row    *rows;
    size_t         cap_rows;
    inner_span    *spans;        /* one per inner row of every inner matrix */
    size_t         cap_spans;
    double        *forward;      /* outer, one row of widest * N_OUTER */
    double        *backward;     /* outer, two rows of widest * N_OUTER */
    size_t         cap_outer;
    double        *scale;        /* one per outer row */
    size_t         cap_scale;
    double        *inner;        /* the inner matrices of every segment */
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

    aln_span span;
    size_t   n_segments;
    int32_t  widest;     /* cells of the widest outer row */
    int32_t  widest_span;  /* cells of the widest inner row */
    size_t   n_spans;    /* inner rows of every inner matrix */
    int32_t  origin;     /* reference position of window value 0 */
    size_t   window_len;
} context;

static size_t n_rows(const context *ctx)
{
    return ctx->n_segments + 1;
}

static inline outer_row *row_of(const context *ctx, size_t k)
{
    return &ctx->scratch->rows[k];
}

/* The forward pass keeps one outer row: scaling a row moves it into inner row 0 of the
 * next inner matrix. */
static inline double *forward_at(const context *ctx, size_t k, int32_t x)
{
    const outer_row *row = row_of(ctx, k);

    return ctx->scratch->forward + (size_t)(x - row->lo) * N_OUTER;
}

/* The backward pass fills row k - 1 from the inner matrix of row k, so it keeps two
 * rows, and row k lies in the slot of its parity. */
static inline double *backward_at(const context *ctx, size_t k, int32_t x)
{
    const outer_row *row = row_of(ctx, k);

    return ctx->scratch->backward
         + ((k & 1) * (size_t)ctx->widest + (size_t)(x - row->lo)) * N_OUTER;
}

/* Returns the span of inner row j of row. */
static inline const inner_span *span_at(const context *ctx, const outer_row *row, int32_t j)
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
static inline const double *cell_in(const double *cells, const inner_span *span, int32_t x)
{
    return cells + (size_t)(x - span->from) * N_INNER;
}

/* ------------------------------------------------------------------------ */
/* One base against the reference                                            */
/* ------------------------------------------------------------------------ */

/* The comparison of a read base with a reference base, as lines in the modification
 * rate m of the reference base. The emission is intercept + slope * m. The part of the
 * emission on which the template differed is share * m. A read base agrees with its
 * template base with chance (1 - m)(1 - e) + m e / 3, and differs from it into one
 * given base with chance (1 - m) e / 3 + m (1 - e / 3) / 3, where e is the chance the
 * base was misread. */
typedef struct {
    double intercept;
    double slope;
    double share;
} comparison;

/* The comparisons of one read base with a reference base that agrees with it, one that
 * differs, and one that is not named. A template base the sequencer did not count was
 * never observed, so it emits one and carries nothing. An unnamed base on either side
 * emits the uninformative chance and carries nothing. */
typedef struct {
    comparison agree;
    comparison differ;
    comparison unnamed;
} comparisons;

static void set_comparison(comparison *o, double intercept, double slope, double share)
{
    o->intercept = intercept;
    o->slope     = slope;
    o->share     = share;
}

/* Writes the comparisons of a read base that is named or not, counted or not, and was
 * misread with chance error. */
static void set_comparisons(comparisons *o, bool named, bool uncounted, double error)
{
    double other = error / OTHER_BASES;
    double shift = 1.0 - error - other;
    double blank = uncounted ? 1.0 : UNINFORMATIVE;

    set_comparison(&o->unnamed, blank, 0.0, 0.0);

    if (uncounted || !named) {
        set_comparison(&o->agree, blank, 0.0, 0.0);
        set_comparison(&o->differ, blank, 0.0, 0.0);
        return;
    }

    set_comparison(&o->agree, 1.0 - error, -shift, other);
    set_comparison(&o->differ, other, shift / OTHER_BASES, (1.0 - other) / OTHER_BASES);
}

/* The comparisons of a read base at each quality score. A read that stores no
 * qualities has bases that are never misread. An unnamed read base and a template base
 * the sequencer did not count each have an entry of their own. */
struct phmm_emissions {
    comparisons at_quality[PHRED_MAX + 1];
    comparisons without_quality;
    comparisons unnamed;
    comparisons uncounted;
};

phmm_emissions *phmm_emissions_create(const phred *quality)
{
    phmm_emissions *emissions = malloc(sizeof *emissions);

    if (!emissions) {
        return NULL;
    }

    for (unsigned q = 0; q <= PHRED_MAX; q++) {
        set_comparisons(&emissions->at_quality[q], true, false,
                        phred_error(quality, (uint8_t)q));
    }

    set_comparisons(&emissions->without_quality, true, false, 0.0);
    set_comparisons(&emissions->unnamed, false, false, 0.0);
    set_comparisons(&emissions->uncounted, false, true, 0.0);

    return emissions;
}

void phmm_emissions_destroy(phmm_emissions *emissions)
{
    free(emissions);
}

/* Returns the comparisons of read base q, whose base is ours. Every base of an Ultima
 * read takes one score, since its quality string describes run lengths. */
static const comparisons *comparisons_of(const context *ctx, nuc ours, int32_t q)
{
    const phmm_emissions *emissions = ctx->emissions;

    if (!nuc_is_base(ours)) {
        return &emissions->unnamed;
    }
    if (ctx->errors) {
        return &emissions->at_quality[(uint8_t)ctx->errors->phred];
    }
    if (!ctx->read->qual) {
        return &emissions->without_quality;
    }

    return &emissions->at_quality[ctx->read->qual[q]];
}

static inline const comparison *comparison_at(const comparisons *o, const phmm_base *base,
                                       nuc ours)
{
    nuc theirs = (nuc)base->code;

    if (!nuc_is_base(theirs)) {
        return &o->unnamed;
    }

    return theirs == ours ? &o->agree : &o->differ;
}

static inline double emission_of(const comparison *o, const phmm_base *base)
{
    return o->intercept + o->slope * base->modification;
}

static inline double modified_of(const comparison *o, const phmm_base *base)
{
    return o->share * base->modification;
}

/* ------------------------------------------------------------------------ */
/* Layout                                                                    */
/* ------------------------------------------------------------------------ */

/* Returns where the CIGAR has reached after j bases of the segment that starts at read
 * base start: the prefix once the base is taken, and the prefix once the deletion
 * straight after it is taken too. */
static const aln_place *cigar_place(const context *ctx, int32_t start, int32_t j)
{
    size_t prefix = (size_t)(start - ctx->span.begin) + (size_t)j;

    return &ctx->scratch->places[prefix];
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

/* Writes the window of a row: the band either side of the prefix the CIGAR has
 * reached after the segment. */
static void lay_out_window(const context *ctx, outer_row *row)
{
    const aln_place *place = cigar_place(ctx, row->start, row->len);

    row->lo = clamp_prefix(ctx, (int32_t)place->first - ctx->model->band);
    row->hi = clamp_prefix(ctx, (int32_t)place->last + ctx->model->band);
}

/* Returns the prefixes inner row j may hold: the band either side of the prefix the
 * CIGAR has reached after j bases of the segment. A row past the bases of the segment
 * takes the band of the last base. */
static inner_span inner_span_of(const context *ctx, const outer_row *row, int32_t j)
{
    int32_t          real  = j < row->len ? j : row->len;
    const aln_place *place = cigar_place(ctx, row->start, real);

    return (inner_span){
        .from = clamp_prefix(ctx, (int32_t)place->first - ctx->model->band),
        .to   = clamp_prefix(ctx, (int32_t)place->last + ctx->model->band),
    };
}

/* Returns how many cells the row holds. */
static int32_t span_width(const inner_span *span)
{
    return span->end - span->from + 1;
}

/* Returns the prefix of the last cell of a row. The row has one cell past its band
 * where its band ends where the band of the row before ends, and that cell lies in
 * the window of the outer row. */
static int32_t span_end_of(const outer_row *row, const inner_span *before,
                           const inner_span *span)
{
    if (before && before->to == span->to && span->to < row->hi) {
        return span->to + 1;
    }

    return span->to;
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
        row->length[PHMM_MAX_LENGTH_ERROR] = 1.0;
        return;
    }

    memcpy(row->length, ctx->errors->of_run[k - 1].length, sizeof row->length);
    row->last_inner_row = row->len + uncounted_rows(row->length);
}

/* Writes the chance of each template length of segment k. A read without the tags has
 * the read length. */
static void lay_out_lengths(const context *ctx, outer_row *row, size_t k)
{
    memset(row->length, 0, sizeof row->length);
    row->last_inner_row = row->len;

    if (ctx->errors) {
        lay_out_run_lengths(ctx, row, k);
    } else {
        row->length[PHMM_MAX_LENGTH_ERROR] = 1.0;
    }
}

/* Writes the chance that a one-base run is missing after the run that is segment k,
 * where that run is not the last of the read. */
static void lay_out_missed(const context *ctx, outer_row *row, size_t k)
{
    row->missed         = 0.0;
    row->missed_bases = 0;

    if (ctx->errors && k < ctx->n_segments) {
        row->missed         = ctx->errors->of_run[k - 1].missed;
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
    outer_row       *rows   = ctx->scratch->rows;
    const aln_place *start;
    int32_t          widest = 0;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        lay_out_segment(ctx, &rows[k], k);
        lay_out_window(ctx, &rows[k]);
        lay_out_lengths(ctx, &rows[k], k);
        lay_out_missed(ctx, &rows[k], k);
    }

    start = cigar_place(ctx, rows[1].start, 0);

    rows[0].lo             = clamp_prefix(ctx, (int32_t)start->first - ctx->model->band);
    rows[0].hi             = clamp_prefix(ctx, (int32_t)start->last + ctx->model->band);
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
            inner_span *span   = &ctx->scratch->spans[row->first_span + (size_t)j];
            inner_span *before = j > 0 ? span - 1 : NULL;

            *span       = inner_span_of(ctx, row, j);
            span->end   = span_end_of(row, before, span);
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
static size_t outer_cells(const context *ctx)
{
    return (size_t)ctx->widest * N_OUTER;
}

/* Grows every buffer the passes need. The outer buffers are a pair, so they are grown
 * together. */
static int grow_buffers(context *ctx, size_t inner)
{
    phmm_scratch *s     = ctx->scratch;
    size_t          outer = 2 * outer_cells(ctx);
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
             sizeof *s->places) < 0) {
        return PHMM_NO_MEMORY;
    }

    ctx->span = aln_places(read, s->places);

    if (ctx->span.begin >= ctx->span.end) {
        return PHMM_NO_PATH;
    }

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

/* The weights of the gaps, which hold along every reference. A deletion pays its
 * closing once, with the opening rate of the base it closes at. An insertion pays its
 * closing once, as the first inserted base, and the opening rate at the base it closes
 * at. Each inserted base pays the chance of being the base it is. */
typedef struct {
    double extend_deletion;
    double close_deletion;
    double extend_insertion;  /* per inserted base after the first */
    double open_insertion;    /* the first inserted base */
} gap_weights;

static gap_weights gap_weights_of(const context *ctx)
{
    const phmm_rates *rates = ctx->rates;

    return (gap_weights){
        .extend_deletion  = rates->deletion_to_deletion,
        .close_deletion   = rates->match_to_deletion,
        .extend_insertion = rates->insertion_to_insertion * UNINFORMATIVE,
        .open_insertion   = rates->match_to_insertion * UNINFORMATIVE,
    };
}

/* Whether read base j of a segment can be inserted. The first base of the read cannot
 * be, nor can a template base the sequencer did not count. */
static bool may_insert(const outer_row *row, size_t k, int32_t j)
{
    return j < row->len && (j >= 1 || k >= 2);
}

/* Whether a deletion can open inside the segment after inner row j. */
static bool may_delete_inside(const outer_row *row, int32_t j)
{
    return j >= 1 && j < row->last_inner_row;
}

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
 * pass: the chance of the template length that ends the segment at the row, which is
 * zero where none does; the outer backward cells of the row of the segment, from
 * prefix lo to hi; and what can follow the end. Nothing follows the last segment. */
typedef struct {
    double        chance;
    const double *cells;
    int32_t       lo, hi;
    bool          deletes;   /* a deletion can open after the segment */
    bool          inserts;   /* an insertion left open can end at the row */
} segment_end;

/* Fills end with what inner row j of segment k is worth at the end of the segment. */
static inline void set_segment_end(segment_end *end, const context *ctx, size_t k, int32_t j)
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
        end->chance = template_length_at(row, offset).chance;
    }
}

/* Returns what the paired state of a cell at prefix x is worth at the end of the
 * segment: the chance of ending there times what the outer row is worth closed at x,
 * and with a deletion open at x + 1. */
static inline double end_paired(const segment_end *end, int32_t x)
{
    const double *cell;
    double        worth;

    if (end->chance <= 0.0 || x < end->lo || x > end->hi) {
        return 0.0;
    }

    cell  = end->cells + (size_t)(x - end->lo) * N_OUTER;
    worth = cell[CLOSED];

    if (end->deletes && x < end->hi) {
        worth += cell[N_OUTER + DELETION];
    }

    return worth * end->chance;
}

/* Returns what the open insertion state of a cell at prefix x is worth at the end of
 * the segment. */
static inline double end_inserted(const segment_end *end, int32_t x)
{
    if (end->chance <= 0.0 || !end->inserts || x < end->lo || x > end->hi) {
        return 0.0;
    }

    return end->cells[(size_t)(x - end->lo) * N_OUTER + INSERTION] * end->chance;
}

/* What the step from inner row j to inner row j + 1 reads: the two spans, the
 * comparisons of the read base of row j, what the step may do, and, for the backward
 * pass, the scale of the segment and what row j is worth at the end of the segment. */
typedef struct {
    const inner_span  *span;      /* of row j */
    const inner_span  *next;      /* of row j + 1 */
    const comparisons *comparisons;
    nuc                base;      /* NUC_N on an uncounted row */
    double             restore;   /* undoes the scaling of the segment */
    segment_end        end;
    bool               uncounted; /* row j is a template base the sequencer did not count */
    bool               deletes;   /* a deletion can open inside after row j */
    bool               inserts;   /* the read base of row j can be inserted */
    bool               ends;      /* a pairing from row j is the read's 5'-most */
} inner_step;

/* Fills s with the step from inner row j of segment k, for the forward pass. The
 * struct is filled in place, since a struct returned by value is copied through the
 * stack on every inner row. */
static inline void set_inner_step(inner_step *s, const context *ctx, size_t k, int32_t j)
{
    const outer_row *row = row_of(ctx, k);

    s->span        = span_at(ctx, row, j);
    s->next        = span_at(ctx, row, j + 1);
    s->uncounted   = j >= row->len;
    s->deletes     = may_delete_inside(row, j);
    s->inserts     = may_insert(row, k, j);
    s->ends        = k == 1 && j == 0;
    s->base        = s->uncounted ? NUC_N : nuc_from_read(ctx->read->seq, row->start + j);
    s->comparisons = s->uncounted ? &ctx->emissions->uncounted
                                  : comparisons_of(ctx, s->base, row->start + j);
}

/* Fills s with the step from inner row j of segment k, for the backward pass. */
static inline void set_backward_step(inner_step *s, const context *ctx, size_t k, int32_t j,
                              double restore)
{
    set_inner_step(s, ctx, k, j);
    s->restore = restore;
    set_segment_end(&s->end, ctx, k, j);
}

static inline int32_t max_prefix(int32_t a, int32_t b)
{
    return a > b ? a : b;
}

static inline int32_t min_prefix(int32_t a, int32_t b)
{
    return a < b ? a : b;
}

/* Returns the weight of the paths from cell from to the pairing of reference base x,
 * which lands one prefix further along the next row. */
static inline double pairing_into(const double *from, const phmm_base *b)
{
    return from[PAIRED] * b->after_no_gap + from[GAP_DELETION] * b->after_deletion
         + from[GAP_INSERTION] * b->after_insertion;
}

/* Fills inner row j + 1 from inner row j. Each cell of the next row is written once:
 * its pairing comes from the cell one prefix back, its insertion from the cell at the
 * same prefix, and its open deletion from the cell before it in its own row. A cell
 * whose source lies outside the band of the row above takes nothing from it. */
static void inner_forward_step(const context *ctx, size_t k, const inner_step *s,
                               int32_t j)
{
    const phmm_base  *bases    = ctx->model->bases;
    gap_weights        w        = gap_weights_of(ctx);
    const double      *line     = inner_row_at(ctx, k, s->span);
    double            *next     = inner_row_at(ctx, k, s->next);
    bool               deletes  = may_delete_inside(row_of(ctx, k), j + 1);
    int32_t            pair_lo  = max_prefix(s->next->from, s->span->from + 1);
    int32_t            pair_hi  = min_prefix(s->next->end, s->span->to + 1);
    int32_t            ins_hi   = s->inserts ? min_prefix(s->next->to, s->span->to) : -1;
    double             deletion = 0.0;
    int32_t            x        = s->next->from;

    for (; x < pair_lo; x++) {
        const double *same = cell_in(line, s->span, x);
        double       *cell = cell_at(next, s->next, x);
        double        inserted = x <= ins_hi
            ? same[PAIRED] * w.open_insertion + same[GAP_INSERTION] * w.extend_insertion
            : 0.0;

        cell[PAIRED]        = 0.0;
        cell[GAP_DELETION]  = deletion;
        cell[GAP_INSERTION] = inserted;
        deletion            = deletion * w.extend_deletion;
    }

    for (; x <= pair_hi; x++) {
        const double    *back   = cell_in(line, s->span, x - 1);
        const phmm_base *b      = &bases[x - 1];
        double          *cell   = cell_at(next, s->next, x);
        const comparison *o     = comparison_at(s->comparisons, b, s->base);
        double           paired = pairing_into(back, b) * emission_of(o, b);
        double           inserted = x <= ins_hi
            ? back[N_INNER + PAIRED] * w.open_insertion
              + back[N_INNER + GAP_INSERTION] * w.extend_insertion
            : 0.0;

        cell[PAIRED]        = paired;
        cell[GAP_DELETION]  = deletion;
        cell[GAP_INSERTION] = inserted;
        deletion            = deletion * w.extend_deletion + (deletes ? paired : 0.0);
    }

    for (; x <= s->next->end; x++) {
        double *cell = cell_at(next, s->next, x);

        cell[PAIRED]        = 0.0;
        cell[GAP_DELETION]  = deletion;
        cell[GAP_INSERTION] = 0.0;
        deletion            = deletion * w.extend_deletion;
    }
}

/* Fills the inner matrix of segment k from its inner row 0. */
static void inner_forward(const context *ctx, size_t k)
{
    const outer_row *row = row_of(ctx, k);

    for (int32_t j = 0; j < row->last_inner_row; j++) {
        inner_step s;

        set_inner_step(&s, ctx, k, j);
        inner_forward_step(ctx, k, &s, j);
    }
}

/* Whether the sequencer can have missed a one-base run after the segment of row k. */
static bool missed_run_possible(const context *ctx, size_t k)
{
    const outer_row *row = row_of(ctx, k);

    return k != ctx->n_segments && row->missed > 0.0 && row->missed_bases != 0;
}

/* Whether reference base x is a run of one base that the sequencer can have missed
 * after the row's segment. */
static bool missed_run_at(const context *ctx, const outer_row *row, int32_t x)
{
    const phmm_model *model = ctx->model;
    const phmm_base  *base;

    if ((size_t)x >= model->len) {
        return false;
    }

    base = &model->bases[x];

    return base->is_one_base_run && nuc_is_base((nuc)base->code)
        && (row->missed_bases & (1U << base->code)) != 0;
}

/* One inner row at which the segment can end, as the exits read it: the chance of the
 * template length that ends there, the row's cells, and whether an insertion left
 * open can end there, which it can where the template is no longer than the read. */
typedef struct {
    double            chance;
    const inner_span *span;
    const double     *cells;
    bool              inserts;
} exit_row;

/* Writes the rows at which the segment of row k can end and returns their count. */
static size_t exit_rows_of(const context *ctx, size_t k, exit_row *exits)
{
    const outer_row *row = row_of(ctx, k);
    size_t           n   = 0;

    for (int offset = -PHMM_MAX_LENGTH_ERROR; offset <= PHMM_MAX_LENGTH_ERROR; offset++) {
        template_length length = template_length_at(row, offset);

        if (length.chance <= 0.0) {
            continue;
        }

        exits[n].chance  = length.chance;
        exits[n].span    = span_at(ctx, row, length.jend);
        exits[n].cells   = inner_row_at(ctx, k, exits[n].span);
        exits[n].inserts = offset <= 0;
        n++;
    }

    return n;
}

/* Adds what one exit row holds at prefix x, weighted by its chance, to closed and
 * inserted. A row with no cell at x adds nothing. */
static inline void exit_at(const exit_row *exit, int32_t x, double *closed, double *inserted)
{
    const double *cell;

    if (x < exit->span->from || x > exit->span->end) {
        return;
    }

    cell     = cell_in(exit->cells, exit->span, x);
    *closed += cell[PAIRED] * exit->chance;

    if (exit->inserts) {
        *inserted += cell[GAP_INSERTION] * exit->chance;
    }
}

/* Writes the outer cells of row k from the inner rows at which the segment can end,
 * each weighted by its chance, and returns their total. A deletion that opens after
 * the segment and an insertion left open end in their own states, which the last row
 * cannot hold. The deletion open at a prefix comes from the pairings closed at the
 * prefix before. */
static double fill_outer_row(const context *ctx, size_t k)
{
    const outer_row *row    = row_of(ctx, k);
    bool             last   = k == ctx->n_segments;
    exit_row         exits[PHMM_LENGTH_OFFSETS];
    size_t           n      = exit_rows_of(ctx, k, exits);
    double           before = 0.0;  /* closed at the prefix before */
    double           total  = 0.0;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        double *to       = forward_at(ctx, k, x);
        double  closed   = 0.0;
        double  inserted = 0.0;

        for (size_t e = 0; e < n; e++) {
            exit_at(&exits[e], x, &closed, &inserted);
        }

        to[CLOSED]    = closed;
        to[DELETION]  = last ? 0.0 : before;
        to[INSERTION] = last ? 0.0 : inserted;
        total        += to[CLOSED] + to[DELETION] + to[INSERTION];
        before        = closed;
    }

    return total;
}

/* Passes over the one-base runs the sequencer can have missed after the row's
 * segment, and returns what the pass adds to the total of the row. */
static double missed_run_forward(const context *ctx, size_t k)
{
    const outer_row *row   = row_of(ctx, k);
    double           added = 0.0;

    if (!missed_run_possible(ctx, k)) {
        return 0.0;
    }

    for (int32_t x = row->lo; x < row->hi; x++) {
        if (missed_run_at(ctx, row, x)) {
            double missed = forward_at(ctx, k, x)[CLOSED] * row->missed;

            forward_at(ctx, k, x + 1)[CLOSED] += missed;
            added                             += missed;
        }
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

/* Keeps the total of row k. Unless the row is the last, divides the row by its total
 * into inner row 0 of the next inner matrix. The open deletion of an outer cell
 * extends along inner row 0 as it is moved, so that every cell of that row holds the
 * whole of its open deletion state. */
static phmm_status scale_row(const context *ctx, size_t k, double total)
{
    const outer_row  *row      = row_of(ctx, k);
    double            deletion = 0.0;
    double            share;
    phmm_status       status;
    const inner_span *span;
    double           *first;

    status = total_status(total);

    if (status != PHMM_OK) {
        return status;
    }

    ctx->scratch->scale[k] = total;

    if (k + 1 == n_rows(ctx)) {
        return PHMM_OK;
    }

    share = 1.0 / total;
    span  = span_at(ctx, row_of(ctx, k + 1), 0);
    first = inner_row_at(ctx, k + 1, span);

    for (int32_t x = row->lo; x <= row->hi; x++) {
        double *cell = forward_at(ctx, k, x);
        double *to   = cell_at(first, span, x);

        deletion = deletion * ctx->rates->deletion_to_deletion + cell[DELETION] * share;

        to[PAIRED]        = cell[CLOSED] * share;
        to[GAP_DELETION]  = deletion;
        to[GAP_INSERTION] = cell[INSERTION] * share;
    }

    return PHMM_OK;
}

/* Fills row 0, where the read can start at any base of the first window with that
 * base's weight in ends, and returns its total. */
static double forward_first_row(const context *ctx)
{
    const outer_row *row   = row_of(ctx, 0);
    double           total = 0.0;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        double *cell = forward_at(ctx, 0, x);

        cell[CLOSED]    = rate_at(ctx->rates->ends, ctx->model->len, x);
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

/* The window position of reference base x. */
static phmm_position *landing(const context *ctx, int32_t x)
{
    return ctx->scratch->window + (size_t)(x - ctx->origin);
}

static void missed_run_backward(const context *ctx, size_t k)
{
    const outer_row *row = row_of(ctx, k);

    if (!missed_run_possible(ctx, k)) {
        return;
    }

    for (int32_t x = row->hi - 1; x >= row->lo; x--) {
        if (missed_run_at(ctx, row, x)) {
            backward_at(ctx, k, x)[CLOSED] += row->missed
                                            * backward_at(ctx, k, x + 1)[CLOSED];
        }
    }
}

/* Writes the backward values of the last inner row of segment k. Each cell is worth
 * what it is worth at the end of the segment. */
static void start_inner_backward(const context *ctx, size_t k)
{
    const outer_row  *row   = row_of(ctx, k);
    const inner_span *span  = span_at(ctx, row, row->last_inner_row);
    double           *cells = inner_back_row_at(ctx, row->last_inner_row);
    segment_end       end;

    set_segment_end(&end, ctx, k, row->last_inner_row);

    for (int32_t x = span->from; x <= span->end; x++) {
        double *cell = cell_at(cells, span, x);

        cell[PAIRED]        = end_paired(&end, x);
        cell[GAP_DELETION]  = 0.0;
        cell[GAP_INSERTION] = end_inserted(&end, x);
    }
}

/* Accumulates one pairing of reference base x into the window. modified is the
 * posterior of the paths through the pairing on which the template differed, and real
 * whether a read base was observed. */
static inline void accumulate_pairing(const context *ctx, const inner_step *s, int32_t x,
                               double modified, double paired, double deleted,
                               double inserted)
{
    phmm_position *at = landing(ctx, x);

    if (!s->uncounted) {
        at->coverage   += paired + deleted + inserted;
        at->mismatches += modified;
    }

    landing(ctx, x - 1)->deletions += deleted;
    at->insertions                 += inserted;

    if (s->ends) {
        at->ends += paired;
    }
}

/* Writes a backward cell at prefix x that has no pairing: what it is worth at the end
 * of the segment, through the insertion of its read base into gap, the cell of the
 * next row at its prefix, where that cell exists, and through the open deletion of
 * the cell after it in its own row. Returns the open deletion state of the cell. */
static inline double backward_gap_cell(double *cell, const inner_step *s, int32_t x,
                                const double *gap, const gap_weights *w, double deletion)
{
    double paired   = end_paired(&s->end, x);
    double inserted = end_inserted(&s->end, x);

    if (gap) {
        paired   += w->open_insertion * gap[GAP_INSERTION];
        inserted += w->extend_insertion * gap[GAP_INSERTION];
    }

    cell[PAIRED]        = paired + (s->deletes ? deletion : 0.0);
    cell[GAP_DELETION]  = deletion * w->extend_deletion;
    cell[GAP_INSERTION] = inserted;

    return cell[GAP_DELETION];
}

/* Writes the backward inner row j of segment k: what its cells are worth at the end
 * of the segment and through inner row j + 1, accumulating each pairing with its
 * posterior. Each cell is read from the cell one prefix on in the next row, which its
 * pairing lands in, from the cell at the same prefix, which its insertion lands in,
 * and from the cell after it in its own row, which its open deletion extends into. A
 * cell whose target lies outside the next row takes nothing from it. The cell past
 * the band, where the row has one, is worth only the end of the segment. The posterior
 * of a template difference at a pairing is the posterior of the pairing times the
 * modified part of its emission over the whole emission, so the emission divides out
 * of it. */
static void inner_backward_step(const context *ctx, size_t k, const inner_step *s,
                                int32_t j)
{
    const phmm_base  *bases    = ctx->model->bases;
    gap_weights        w        = gap_weights_of(ctx);
    double             restore  = s->restore;
    const double      *forward  = inner_row_at(ctx, k, s->span);
    double            *line     = inner_back_row_at(ctx, j);
    const double      *next     = inner_back_row_at(ctx, j + 1);
    int32_t            pair_lo  = max_prefix(s->span->from, s->next->from - 1);
    int32_t            pair_hi  = min_prefix(s->span->to, s->next->end - 1);
    int32_t            ins_lo   = s->inserts ? max_prefix(s->span->from, s->next->from)
                                             : INT32_MAX;
    double             deletion = 0.0;  /* of the cell to the right */
    int32_t            x        = s->span->end;

    for (; x > s->span->to; x--) {
        double *cell = cell_at(line, s->span, x);

        cell[PAIRED]        = end_paired(&s->end, x);
        cell[GAP_DELETION]  = 0.0;
        cell[GAP_INSERTION] = end_inserted(&s->end, x);
    }

    for (; x > pair_hi; x--) {
        const double *gap = x >= ins_lo ? cell_in(next, s->next, x) : NULL;

        deletion = backward_gap_cell(cell_at(line, s->span, x), s, x, gap, &w, deletion);
    }

    for (; x >= pair_lo; x--) {
        const phmm_base *b        = &bases[x];
        const double    *fwd      = cell_in(forward, s->span, x);
        const double    *on       = cell_in(next, s->next, x + 1);
        double          *cell     = cell_at(line, s->span, x);
        double           below    = on[PAIRED] * restore;
        double           through  = pairing_into(fwd, b) * below;
        const comparison *o       = comparison_at(s->comparisons, b, s->base);
        double           emission = emission_of(o, b);
        double           reach    = on[PAIRED] * emission;
        double           paired   = end_paired(&s->end, x) + b->after_no_gap * reach
                                  + (s->deletes ? deletion : 0.0);
        double           deleted  = b->after_deletion * reach
                                  + deletion * w.extend_deletion;
        double           inserted = end_inserted(&s->end, x) + b->after_insertion * reach;

        if (x >= ins_lo) {
            double gap = on[-N_INNER + GAP_INSERTION];

            paired   += w.open_insertion * gap;
            inserted += w.extend_insertion * gap;
        }

        accumulate_pairing(ctx, s, x, through * modified_of(o, b),
                           fwd[PAIRED] * b->after_no_gap * emission * below,
                           fwd[GAP_DELETION] * b->after_deletion * emission * below,
                           fwd[GAP_INSERTION] * b->after_insertion * emission * below);

        cell[PAIRED]        = paired;
        cell[GAP_DELETION]  = deleted;
        cell[GAP_INSERTION] = inserted;
        deletion            = deleted;
    }

    for (; x >= s->span->from; x--) {
        const double *gap = x >= ins_lo ? cell_in(next, s->next, x) : NULL;

        deletion = backward_gap_cell(cell_at(line, s->span, x), s, x, gap, &w, deletion);
    }
}

/* Writes the outer backward cells of row k - 1 from inner row 0 of segment k. */
static void leave_inner_backward(const context *ctx, size_t k)
{
    const outer_row  *above   = row_of(ctx, k - 1);
    const inner_span *span    = span_at(ctx, row_of(ctx, k), 0);
    const double     *cells   = inner_back_row_at(ctx, 0);
    double            restore = 1.0 / ctx->scratch->scale[k];

    for (int32_t x = above->lo; x <= above->hi; x++) {
        const double *first = cell_in(cells, span, x);
        double       *to    = backward_at(ctx, k - 1, x);

        to[CLOSED]    = first[PAIRED] * restore;
        to[DELETION]  = first[GAP_DELETION] * restore;
        to[INSERTION] = first[GAP_INSERTION] * restore;
    }
}

/* Fills the backward inner matrix of segment k row by row from its last, and sets the
 * outer cells of row k - 1 from its first. */
static void inner_backward(const context *ctx, size_t k)
{
    const outer_row *row     = row_of(ctx, k);
    double           restore = 1.0 / ctx->scratch->scale[k];

    start_inner_backward(ctx, k);

    for (int32_t j = row->last_inner_row - 1; j >= 0; j--) {
        inner_step s;

        set_backward_step(&s, ctx, k, j, restore);
        inner_backward_step(ctx, k, &s, j);
    }

    leave_inner_backward(ctx, k);
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
