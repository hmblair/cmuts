/* phmm.c -- a banded pair HMM.
 *
 * The matrix has one outer row per read run, and row 0 before the first. A cell of an
 * outer row is a prefix length of the reference, and holds three states: the run is
 * complete and the last event was a pairing, a deletion is open after it, or an
 * insertion is open after it.
 *
 * Between two outer rows lies the inner matrix of a read run, which pairs the bases of
 * the run with the reference one at a time. Its rows are the bases of the run, and up
 * to two more for template bases that an Ultima sequencer did not count. Its columns
 * are the prefix lengths the two outer rows reach. Its states are the three of a pair
 * HMM over bases, with the weights of that model: a pairing pays the weight of its base
 * after no gap, or the opening rate of the gap it closes, times the emission of the
 * comparison, and a gap pays its extension rate per base after the first.
 *
 * A deletion that follows the last pairing of a run leaves the inner matrix as the open
 * deletion state of the outer cell, and extends or closes in the next inner matrix. A
 * deletion inside a run opens, extends, and closes within one inner matrix. The two
 * never describe the same path, since an open deletion extends only in the inner
 * matrix that closes it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "phmm.h"

#include <float.h>
#include <math.h>
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
 * first read run, around where the CIGAR starts. */
typedef struct {
    run     read_run;                /* unused on row 0 */
    int32_t lo, hi;                  /* the prefixes the row's cells cover */
    int32_t xlo, xhi;                /* the prefixes the inner matrix covers */
    int32_t last_inner_row;          /* the run's bases, and the uncounted rows */
    double  length[PHMM_LENGTH_OFFSETS];  /* chance the template run is s + offset long */
    double  missed;                  /* chance a one-base run after it was not read */
    uint8_t missed_bases;            /* the bases that run can be, as bits by nuc */
    size_t  inner;                   /* offset of the inner matrix in scratch */
    size_t  first_span;              /* index of the span of inner row 0 */
} outer_row;

/* The prefixes one inner row may hold. */
typedef struct {
    int32_t from;
    int32_t to;
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
    double        *inner;        /* the inner matrices of every read run */
    size_t         cap_inner;
    double        *inner_back;   /* the inner matrix of one read run, backward */
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
    int32_t  widest;     /* cells of the widest outer row */
    size_t   n_spans;    /* inner rows of every inner matrix */
    int32_t  origin;     /* reference position of window value 0 */
    size_t   window_len;
} context;

static size_t n_rows(const context *ctx)
{
    return ctx->n_read_runs + 1;
}

static outer_row *row_of(const context *ctx, size_t k)
{
    return &ctx->scratch->rows[k];
}

/* The forward pass keeps one outer row: scaling a row moves it into inner row 0 of the
 * next inner matrix. */
static double *forward_at(const context *ctx, size_t k, int32_t x)
{
    const outer_row *row = row_of(ctx, k);

    return ctx->scratch->forward + (size_t)(x - row->lo) * N_OUTER;
}

/* The backward pass fills row k - 1 from the inner matrix of row k, so it keeps two
 * rows, and row k lies in the slot of its parity. */
static double *backward_at(const context *ctx, size_t k, int32_t x)
{
    const outer_row *row = row_of(ctx, k);

    return ctx->scratch->backward
         + ((k & 1) * (size_t)ctx->widest + (size_t)(x - row->lo)) * N_OUTER;
}

static size_t inner_width(const outer_row *row)
{
    return (size_t)row->xhi - (size_t)row->xlo + 1;
}

static size_t inner_cells(const outer_row *row)
{
    return (size_t)(row->last_inner_row + 1) * inner_width(row);
}

static double *inner_at(const context *ctx, size_t k, int32_t j, int32_t x)
{
    const outer_row *row = row_of(ctx, k);

    return ctx->scratch->inner + row->inner
         + ((size_t)j * inner_width(row) + (size_t)(x - row->xlo)) * N_INNER;
}

static double *inner_back_at(const context *ctx, size_t k, int32_t j, int32_t x)
{
    const outer_row *row = row_of(ctx, k);

    return ctx->scratch->inner_back
         + ((size_t)j * inner_width(row) + (size_t)(x - row->xlo)) * N_INNER;
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

static const comparison *comparison_at(const comparisons *o, const phmm_base *base,
                                       nuc ours)
{
    nuc theirs = (nuc)base->code;

    if (!nuc_is_base(theirs)) {
        return &o->unnamed;
    }

    return theirs == ours ? &o->agree : &o->differ;
}

static double emission_of(const comparison *o, const phmm_base *base)
{
    return o->intercept + o->slope * base->modification;
}

static double modified_of(const comparison *o, const phmm_base *base)
{
    return o->share * base->modification;
}

/* ------------------------------------------------------------------------ */
/* Layout                                                                    */
/* ------------------------------------------------------------------------ */

/* Returns where the CIGAR has reached after j bases of a read run: the prefix once the
 * base is taken, and the prefix once the deletion straight after it is taken too. */
static const aln_place *cigar_place(const context *ctx, const run *read_run, int32_t j)
{
    size_t prefix = (size_t)(read_run->start - ctx->span.begin) + (size_t)j;

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
 * reached after the run. */
static void lay_out_window(const context *ctx, outer_row *row)
{
    const aln_place *place = cigar_place(ctx, &row->read_run, row->read_run.len);

    row->lo = clamp_prefix(ctx, (int32_t)place->first - ctx->model->band);
    row->hi = clamp_prefix(ctx, (int32_t)place->last + ctx->model->band);
}

/* Returns the prefixes inner row j may hold: the band either side of the prefix the
 * CIGAR has reached after j bases of the run, within the inner matrix. A row past the
 * bases of the run takes the band of the last base. */
static inner_span inner_span_of(const context *ctx, const outer_row *row, int32_t j)
{
    int32_t          real  = j < row->read_run.len ? j : row->read_run.len;
    const aln_place *place = cigar_place(ctx, &row->read_run, real);
    inner_span       span  = {
        .from = clamp_prefix(ctx, (int32_t)place->first - ctx->model->band),
        .to   = clamp_prefix(ctx, (int32_t)place->last + ctx->model->band),
    };

    if (span.from < row->xlo) {
        span.from = row->xlo;
    }
    if (span.to > row->xhi) {
        span.to = row->xhi;
    }

    return span;
}

static const inner_span *span_at(const context *ctx, const outer_row *row, int32_t j)
{
    return &ctx->scratch->spans[row->first_span + (size_t)j];
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

/* Writes the chance of each template length of a run. A run whose length the read
 * does not give in full, or a read without the tags, has the read length. */
static void lay_out_lengths(const context *ctx, outer_row *row, size_t k)
{
    bool clipped = (k == 1 && clipped_before(ctx, &row->read_run))
                                  || (k == ctx->n_read_runs
                                      && clipped_after(ctx, &row->read_run));

    memset(row->length, 0, sizeof row->length);
    row->last_inner_row = row->read_run.len;

    if (!ctx->errors || clipped || !nuc_is_base(row->read_run.base)) {
        row->length[PHMM_MAX_LENGTH_ERROR] = 1.0;
        return;
    }

    memcpy(row->length, ctx->errors->of_run[k - 1].length, sizeof row->length);
    row->last_inner_row = row->read_run.len + uncounted_rows(row->length);
}

static void lay_out_missed(const context *ctx, outer_row *row, size_t k)
{
    row->missed         = 0.0;
    row->missed_bases = 0;

    if (ctx->errors && k >= 1 && k < ctx->n_read_runs) {
        row->missed         = ctx->errors->of_run[k - 1].missed;
        row->missed_bases = ctx->errors->of_run[k - 1].missed_bases;
    }
}

/* Lays out every outer row and returns the cells the inner matrices need in all. */
static size_t lay_out_rows(context *ctx)
{
    outer_row       *rows   = ctx->scratch->rows;
    const aln_place *start;
    size_t           inner  = 0;
    int32_t          widest = 0;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        rows[k].read_run = ctx->read_runs[k - 1];
        lay_out_window(ctx, &rows[k]);
        lay_out_lengths(ctx, &rows[k], k);
        lay_out_missed(ctx, &rows[k], k);
    }

    start = cigar_place(ctx, &rows[1].read_run, 0);

    rows[0].lo             = clamp_prefix(ctx, (int32_t)start->first - ctx->model->band);
    rows[0].hi             = clamp_prefix(ctx, (int32_t)start->last + ctx->model->band);
    rows[0].last_inner_row = 0;

    ctx->n_spans = 0;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        rows[k].xlo        = rows[k - 1].lo < rows[k].lo ? rows[k - 1].lo : rows[k].lo;
        rows[k].xhi        = rows[k - 1].hi > rows[k].hi ? rows[k - 1].hi : rows[k].hi;
        rows[k].inner      = inner * N_INNER;
        rows[k].first_span = ctx->n_spans;
        inner             += inner_cells(&rows[k]);
        ctx->n_spans      += (size_t)rows[k].last_inner_row + 1;
    }

    for (size_t k = 0; k < n_rows(ctx); k++) {
        int32_t width = rows[k].hi - rows[k].lo + 1;

        if (width > widest) {
            widest = width;
        }
    }

    ctx->widest = widest;

    return inner;
}

/* Writes the span of every inner row of every inner matrix. */
static void lay_out_spans(const context *ctx)
{
    for (size_t k = 1; k < n_rows(ctx); k++) {
        const outer_row *row = row_of(ctx, k);

        for (int32_t j = 0; j <= row->last_inner_row; j++) {
            ctx->scratch->spans[row->first_span + (size_t)j] = inner_span_of(ctx, row, j);
        }
    }
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

/* Grows every buffer the read needs. The outer buffers are a pair, so they are grown
 * together. */
static int grow_buffers(context *ctx, size_t inner)
{
    phmm_scratch *s     = ctx->scratch;
    size_t          outer = 2 * outer_cells(ctx);
    size_t          one   = 0;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        size_t cells = inner_cells(row_of(ctx, k)) * N_INNER;

        if (cells > one) {
            one = cells;
        }
    }

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
        || grow((void **)&s->inner, &s->cap_inner, inner * N_INNER, sizeof *s->inner) < 0
        || grow((void **)&s->inner_back, &s->cap_inner_back, one,
                sizeof *s->inner_back) < 0
        || grow((void **)&s->spans, &s->cap_spans, ctx->n_spans, sizeof *s->spans) < 0
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

    if (ctx->n_read_runs == 0) {
        return PHMM_NO_PATH;
    }

    if (grow((void **)&s->places, &s->cap_places, (size_t)read->l_qseq + 1,
             sizeof *s->places) < 0
        || grow((void **)&s->rows, &s->cap_rows, ctx->n_read_runs + 1,
                sizeof *s->rows) < 0) {
        return PHMM_NO_MEMORY;
    }

    ctx->span = aln_places(read, s->places);

    inner = lay_out_rows(ctx);
    lay_out_output(ctx);

    if (grow_buffers(ctx, inner) < 0) {
        return PHMM_NO_MEMORY;
    }

    lay_out_spans(ctx);

    memset(s->forward, 0, outer_cells(ctx) * sizeof *s->forward);
    memset(s->inner, 0, inner * N_INNER * sizeof *s->inner);
    memset(s->window, 0, ctx->window_len * sizeof *s->window);

    return PHMM_OK;
}

/* ------------------------------------------------------------------------ */
/* The inner matrix of a read run                                            */
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

/* Whether read base j of a run can be inserted. The first base of the read cannot be,
 * nor can a template base the sequencer did not count. */
static bool may_insert(const outer_row *row, size_t k, int32_t j)
{
    return j < row->read_run.len && (j >= 1 || k >= 2);
}

/* Whether a deletion can open inside the run after inner row j. */
static bool may_delete_inside(const outer_row *row, int32_t j)
{
    return j >= 1 && j < row->last_inner_row;
}

/* What the loops over one inner row read: where the row and the next lie, the band of
 * the row, and the comparisons of its read base. */
typedef struct {
    double *line;      /* the row's first cell, at prefix xlo */
    double *next;      /* the next row's */
    size_t  stride;    /* doubles per inner row */
    int32_t from, to;  /* the prefixes the row may hold */
    const comparisons *comparisons;
    bool    pairs;     /* whether a row follows to pair into */
    bool    uncounted; /* a template base the sequencer did not count */
    bool    deletes;   /* whether a deletion can open inside after this row */
    bool    inserts;   /* whether this row's read base can be inserted */
} inner_row;

/* Sets r to inner row j of the inner matrix of row k, which begins at matrix. The
 * struct is filled in place: returned by value, it is copied through the stack on every
 * inner row, and the copy stalls on store forwarding. */
static void set_inner_row(inner_row *r, const context *ctx, size_t k, double *matrix,
                          int32_t j)
{
    const outer_row *row  = row_of(ctx, k);
    const inner_span *span = span_at(ctx, row, j);

    r->stride  = inner_width(row) * N_INNER;
    r->line    = matrix + (size_t)j * r->stride;
    r->next    = r->line + r->stride;
    r->from    = span->from;
    r->to      = span->to;
    r->pairs   = j < row->last_inner_row;
    r->uncounted = j >= row->read_run.len;
    r->deletes = may_delete_inside(row, j);
    r->inserts = may_insert(row, k, j);
    r->comparisons = r->uncounted
                   ? &ctx->emissions->uncounted
                   : comparisons_of(ctx, row->read_run.base, row->read_run.start + j);
}

/* Fills the inner matrix of row k from its inner row 0. The open deletion state passes
 * along a row in a local, so that each cell is read and written once. */
static void inner_forward(const context *ctx, size_t k)
{
    const outer_row  *row   = row_of(ctx, k);
    const phmm_base *bases = ctx->model->bases;
    gap_weights        w     = gap_weights_of(ctx);
    int32_t            len   = (int32_t)ctx->model->len;
    nuc                ours  = row->read_run.base;

    for (int32_t j = 0; j <= row->last_inner_row; j++) {
        inner_row r;
        double    open_deletion = 0.0;  /* from the cell to the left */

        set_inner_row(&r, ctx, k, ctx->scratch->inner + row->inner, j);

        for (int32_t x = r.from; x <= r.to; x++) {
            size_t  at       = (size_t)(x - row->xlo) * N_INNER;
            double *cell     = r.line + at;
            double  paired   = cell[PAIRED];
            double  deleted  = cell[GAP_DELETION] + open_deletion;
            double  inserted = cell[GAP_INSERTION];

            cell[GAP_DELETION] = deleted;
            open_deletion      = deleted * w.extend_deletion + (r.deletes ? paired : 0.0);

            if (r.pairs && x < row->xhi && x < len) {
                const phmm_base *b = &bases[x];
                const comparison     *o = comparison_at(r.comparisons, b, ours);

                r.next[at + N_INNER + PAIRED] += (paired * b->after_no_gap
                                                + deleted * b->after_deletion
                                                + inserted * b->after_insertion)
                                               * emission_of(o, b);
            }

            if (r.inserts) {
                r.next[at + GAP_INSERTION] += paired * w.open_insertion
                                            + inserted * w.extend_insertion;
            }
        }
    }
}

/* Whether the sequencer can have missed a one-base run after the run of row k. */
static bool missed_run_possible(const context *ctx, size_t k)
{
    const outer_row *row = row_of(ctx, k);

    return k != ctx->n_read_runs && row->missed > 0.0 && row->missed_bases != 0;
}

/* Whether reference base x is a run of one base that the sequencer can have missed
 * after the row's run. */
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

/* Adds the inner rows at each template length, weighted by its chance, into the outer
 * cells of row k. A deletion that opens after the run and an insertion left open end in
 * their own states, which the last row cannot hold. */
static void fill_outer_row(const context *ctx, size_t k)
{
    const outer_row *row  = row_of(ctx, k);
    bool              last = k == ctx->n_read_runs;

    for (int offset = -PHMM_MAX_LENGTH_ERROR; offset <= PHMM_MAX_LENGTH_ERROR; offset++) {
        double  chance = row->length[offset + PHMM_MAX_LENGTH_ERROR];
        int32_t jend   = row->read_run.len + offset;

        if (chance <= 0.0 || jend < 0 || jend > row->last_inner_row) {
            continue;
        }

        for (int32_t x = row->lo; x <= row->hi; x++) {
            const double *cell = inner_at(ctx, k, jend, x);
            double       *to   = forward_at(ctx, k, x);

            to[CLOSED] += cell[PAIRED] * chance;

            if (last) {
                continue;
            }

            if (x < row->hi) {
                forward_at(ctx, k, x + 1)[DELETION] += cell[PAIRED] * chance;
            }
            if (offset <= 0) {
                to[INSERTION] += cell[GAP_INSERTION] * chance;
            }
        }
    }
}

/* Passes over the one-base runs the sequencer can have missed after the row's run. */
static void missed_run_forward(const context *ctx, size_t k)
{
    const outer_row *row = row_of(ctx, k);

    if (!missed_run_possible(ctx, k)) {
        return;
    }

    for (int32_t x = row->lo; x < row->hi; x++) {
        if (missed_run_at(ctx, row, x)) {
            forward_at(ctx, k, x + 1)[CLOSED] += forward_at(ctx, k, x)[CLOSED]
                                               * row->missed;
        }
    }
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
 * into inner row 0 of the next inner matrix, and clears it for the row after. */
static phmm_status scale_row(const context *ctx, size_t k)
{
    const outer_row *row   = row_of(ctx, k);
    double            total = 0.0;
    double            share;
    phmm_status       status;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        const double *cell = forward_at(ctx, k, x);

        total += cell[CLOSED] + cell[DELETION] + cell[INSERTION];
    }

    status = total_status(total);

    if (status != PHMM_OK) {
        return status;
    }

    ctx->scratch->scale[k] = total;

    if (k + 1 == n_rows(ctx)) {
        return PHMM_OK;
    }

    share = 1.0 / total;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        double *cell = forward_at(ctx, k, x);
        double *first = inner_at(ctx, k + 1, 0, x);

        first[PAIRED]        = cell[CLOSED] * share;
        first[GAP_DELETION]  = cell[DELETION] * share;
        first[GAP_INSERTION] = cell[INSERTION] * share;

        cell[CLOSED]    = 0.0;
        cell[DELETION]  = 0.0;
        cell[INSERTION] = 0.0;
    }

    return PHMM_OK;
}

/* Fills row 0, where the read can start at any base of the first window with that
 * base's weight in ends. */
static void forward_first_row(const context *ctx)
{
    const outer_row *row = row_of(ctx, 0);

    for (int32_t x = row->lo; x <= row->hi; x++) {
        forward_at(ctx, 0, x)[CLOSED] = rate_at(ctx->rates->ends, ctx->model->len, x);
    }
}

static phmm_status forward(const context *ctx)
{
    phmm_status status;

    forward_first_row(ctx);
    status = scale_row(ctx, 0);

    for (size_t k = 1; k < n_rows(ctx) && status == PHMM_OK; k++) {
        inner_forward(ctx, k);
        fill_outer_row(ctx, k);
        missed_run_forward(ctx, k);
        status = scale_row(ctx, k);
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

/* Starts the backward inner matrix of row k: what each inner cell at a template
 * length is worth, from the outer backward cells of row k. */
static void start_inner_backward(const context *ctx, size_t k)
{
    const outer_row *row  = row_of(ctx, k);
    bool              last = k == ctx->n_read_runs;

    memset(ctx->scratch->inner_back, 0,
           inner_cells(row) * N_INNER * sizeof *ctx->scratch->inner_back);

    for (int offset = -PHMM_MAX_LENGTH_ERROR; offset <= PHMM_MAX_LENGTH_ERROR; offset++) {
        double  chance = row->length[offset + PHMM_MAX_LENGTH_ERROR];
        int32_t jend   = row->read_run.len + offset;

        if (chance <= 0.0 || jend < 0 || jend > row->last_inner_row) {
            continue;
        }

        for (int32_t x = row->lo; x <= row->hi; x++) {
            const double *from = backward_at(ctx, k, x);
            double       *cell = inner_back_at(ctx, k, jend, x);
            double        worth = from[CLOSED];

            if (!last && x < row->hi) {
                worth += backward_at(ctx, k, x + 1)[DELETION];
            }

            cell[PAIRED] += worth * chance;

            if (!last && offset <= 0) {
                cell[GAP_INSERTION] += from[INSERTION] * chance;
            }
        }
    }
}

/* Accumulates one pairing of reference base x into the window. modified is the
 * posterior of the paths through the pairing on which the template differed, and real
 * whether a read base was observed. */
static void accumulate_pairing(const context *ctx, size_t k, int32_t j, int32_t x,
                           bool real, double modified, double paired, double deleted,
                           double inserted)
{
    phmm_position *at = landing(ctx, x);

    if (real) {
        at->coverage   += paired + deleted + inserted;
        at->mismatches += modified;
    }

    landing(ctx, x - 1)->deletions += deleted;
    at->insertions                 += inserted;

    if (k == 1 && j == 0) {
        at->ends += paired;
    }
}

/* Fills the backward inner matrix of row k from its start, accumulating each pairing
 * with its posterior, and sets the outer cells of row k - 1 to what inner row 0 is
 * worth. The open deletion state passes along a row in a local, so that each cell is
 * read and written once. The posterior of a template difference at a pairing is the
 * posterior of the pairing times the modified part of its emission over the whole
 * emission, so the emission divides out of it. */
static void inner_backward(const context *ctx, size_t k)
{
    const outer_row  *row     = row_of(ctx, k);
    const outer_row  *above   = row_of(ctx, k - 1);
    const phmm_base *bases   = ctx->model->bases;
    gap_weights        w       = gap_weights_of(ctx);
    int32_t            len     = (int32_t)ctx->model->len;
    nuc                ours    = row->read_run.base;
    double             restore = 1.0 / ctx->scratch->scale[k];
    const double      *forward = ctx->scratch->inner + row->inner;

    for (int32_t j = row->last_inner_row; j >= 0; j--) {
        inner_row r;
        double    open_deletion = 0.0;  /* from the cell to the right */

        set_inner_row(&r, ctx, k, ctx->scratch->inner_back, j);

        for (int32_t x = r.to; x >= r.from; x--) {
            size_t        at       = (size_t)(x - row->xlo) * N_INNER;
            double       *cell     = r.line + at;
            const double *fwd      = forward + (size_t)j * r.stride + at;
            double        paired   = cell[PAIRED] + (r.deletes ? open_deletion : 0.0);
            double        deleted  = cell[GAP_DELETION]
                                   + open_deletion * w.extend_deletion;
            double        inserted = cell[GAP_INSERTION];

            if (r.pairs && x < row->xhi && x < len) {
                const phmm_base *b        = &bases[x];
                const comparison     *o        = comparison_at(r.comparisons, b, ours);
                double             below    = r.next[at + N_INNER + PAIRED];
                double             next     = below * emission_of(o, b);
                double             through  = (fwd[PAIRED] * b->after_no_gap
                                               + fwd[GAP_DELETION] * b->after_deletion
                                               + fwd[GAP_INSERTION] * b->after_insertion)
                                              * below * restore;
                double             to_paired   = b->after_no_gap * next;
                double             to_deleted  = b->after_deletion * next;
                double             to_inserted = b->after_insertion * next;

                paired   += to_paired;
                deleted  += to_deleted;
                inserted += to_inserted;

                accumulate_pairing(ctx, k, j, x, !r.uncounted,
                                   through * modified_of(o, b),
                               fwd[PAIRED] * to_paired * restore,
                               fwd[GAP_DELETION] * to_deleted * restore,
                               fwd[GAP_INSERTION] * to_inserted * restore);
            }

            if (r.inserts) {
                double below = r.next[at + GAP_INSERTION];

                paired   += w.open_insertion * below;
                inserted += w.extend_insertion * below;
            }

            cell[PAIRED]        = paired;
            cell[GAP_DELETION]  = deleted;
            cell[GAP_INSERTION] = inserted;
            open_deletion       = deleted;
        }
    }

    for (int32_t x = above->lo; x <= above->hi; x++) {
        const double *first = inner_back_at(ctx, k, 0, x);
        double       *to    = backward_at(ctx, k - 1, x);

        to[CLOSED]    = first[PAIRED] * restore;
        to[DELETION]  = first[GAP_DELETION] * restore;
        to[INSERTION] = first[GAP_INSERTION] * restore;
    }
}

static void backward_last_row(const context *ctx)
{
    size_t            k   = ctx->n_read_runs;
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
    const outer_row *row   = row_of(ctx, 0);
    double            total = 0.0;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        total += inner_at(ctx, 1, 0, x)[PAIRED] * backward_at(ctx, 0, x)[CLOSED];
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

    for (size_t k = ctx->n_read_runs; k >= 1; k--) {
        missed_run_backward(ctx, k);
        start_inner_backward(ctx, k);
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
