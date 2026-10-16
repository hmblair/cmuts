/* runhmm.c -- a banded pair HMM over the runs of a read.
 *
 * Rows are read runs, and row 0 stands before the first. A cell of row k is a prefix
 * length x of the reference, and holds three states: the run is complete and the last
 * event was a pairing, a deletion is open after it, or an insertion is open after it.
 *
 * The step from row k - 1 to row k walks the bases of read run k against the
 * reference, one base at a time, in a matrix of its own. Its rows are the bases of the
 * run, and two more for the bases that an Ultima sequencer can fail to count. Its
 * columns are the prefix lengths the two outer rows can reach. Its states are the
 * three of the base model. The weights are the base model's: a pairing pays the keep
 * weight of its base, or the opening rate of the gap it closes, times the emission of
 * the comparison, and a gap pays its extension rate per base after the first.
 *
 * A deletion that follows the last pairing of a run leaves the step as the open
 * deletion state of the outer cell, and extends or closes in the next step. A
 * deletion inside a run is opened, extended, and closed within the step. The two
 * never describe the same path, since an open deletion extends only in the step that
 * closes it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "runhmm.h"

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

/* How many bases a run length can be off by. */
#define MAX_LENGTH_ERROR 2
#define LENGTH_OFFSETS   (2 * MAX_LENGTH_ERROR + 1)

/* Emission of a comparison that carries no information. */
#define UNINFORMATIVE (1.0 / NUC_BASES)

/* A misread or modified base is equally likely to be any base the reference does not
 * hold. */
#define OTHER_BASES ((double)(NUC_BASES - 1))

/* Forward times backward sums to one on the first row within this tolerance. */
#define NORMALIZATION_TOLERANCE 1e-6

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
static double keep_at(const phmm_rates *rates, size_t len, hts_pos_t b)
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

static int grow_model(runhmm_model *model, size_t len)
{
    if (grow((void **)&model->runs, &model->cap_runs, len, sizeof *model->runs) < 0) {
        return -1;
    }

    return grow((void **)&model->bases, &model->cap_bases, len, sizeof *model->bases);
}

static runhmm_base base_at(const phmm_rates *rates, const cm_fasta_record *ref,
                           hts_pos_t b)
{
    size_t len = ref->len;

    return (runhmm_base){
        .keep           = keep_at(rates, len, b),
        .close_deletion = rate_at(rates->open_deletion, len, b - 1)
                        * rates->match_to_deletion,
        .open_insertion = rate_at(rates->open_insertion, len, b),
        .modification   = rates->modification[b],
        .code           = (uint8_t)nuc_from_char(ref->seq[b]),
        .one_base_run   = 0,
    };
}

static void fill_bases(runhmm_model *model, const cm_fasta_record *ref)
{
    for (size_t b = 0; b < ref->len; b++) {
        model->bases[b] = base_at(model->rates, ref, (hts_pos_t)b);
    }

    for (size_t m = 0; m < model->n_runs; m++) {
        const run *r = &model->runs[m];

        if (r->len == 1 && nuc_is_base(r->letter)) {
            model->bases[r->start].one_base_run = (uint8_t)(1U << r->letter);
        }
    }
}

int runhmm_model_prepare(runhmm_model *model, const phmm_rates *rates,
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

void runhmm_model_free(runhmm_model *model)
{
    free(model->runs);
    free(model->bases);
    *model = (runhmm_model){ 0 };
}

/* ------------------------------------------------------------------------ */
/* Scratch                                                                   */
/* ------------------------------------------------------------------------ */

/* Where one outer row lies and what its step needs. Row 0 stands before the first read
 * run, around where the CIGAR starts. */
typedef struct {
    run     piece;                   /* the read run; unused on row 0 */
    int32_t lo, hi;                  /* the prefixes the row's cells cover */
    int32_t xlo, xhi;                /* the prefixes the step into the row covers */
    int32_t depth;                   /* inner rows: the run's bases and the phantoms */
    double  length[LENGTH_OFFSETS];  /* chance the template run is s + offset long */
    double  missed;                  /* chance a one-base run after it was not read */
    uint8_t missed_letters;
    size_t  inner;                   /* offset of the step's matrix in scratch */
    size_t  first_span;              /* index of the span of the step's inner row 0 */
} row_layout;

/* The prefixes one inner row of a step may hold. */
typedef struct {
    int32_t from;
    int32_t to;
} inner_span;

struct runhmm_scratch {
    aln_place     *places;
    size_t         cap_places;
    run           *pieces;
    size_t         cap_pieces;
    row_layout    *rows;
    size_t         cap_rows;
    inner_span    *spans;        /* one per inner row of every step */
    size_t         cap_spans;
    flow_errors    errors;
    double        *error;        /* the misread chance of each read base */
    size_t         cap_error;
    double        *forward;      /* outer, one row of widest * N_OUTER */
    double        *backward;     /* outer, two rows of widest * N_OUTER */
    size_t         cap_outer;
    double        *scale;        /* one per outer row */
    size_t         cap_scale;
    double        *inner;        /* the step matrices of every row */
    size_t         cap_inner;
    double        *inner_back;   /* the step matrix of one row, backward */
    size_t         cap_inner_back;
    phmm_position *window;
    size_t         cap_window;
};

runhmm_scratch *runhmm_scratch_create(void)
{
    return calloc(1, sizeof(runhmm_scratch));
}

void runhmm_scratch_destroy(runhmm_scratch *scratch)
{
    if (!scratch) {
        return;
    }

    free(scratch->places);
    free(scratch->pieces);
    free(scratch->rows);
    free(scratch->spans);
    flow_errors_free(&scratch->errors);
    free(scratch->error);
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
    const runhmm_model    *model;
    const phmm_rates      *rates;
    const phred           *quality;
    const cm_bam_record   *read;
    runhmm_scratch        *scratch;
    const flow_tables     *flow;    /* NULL where no read carries the tags */
    const flow_errors     *errors;  /* NULL for a read without the tags */

    aln_span span;
    size_t   n_pieces;   /* read runs; outer rows are one more */
    int32_t  widest;     /* cells of the widest outer row */
    size_t   n_spans;    /* inner rows of every step */
    int32_t  origin;     /* reference position of window value 0 */
    size_t   window_len;
} context;

static size_t n_rows(const context *ctx)
{
    return ctx->n_pieces + 1;
}

static row_layout *row_of(const context *ctx, size_t k)
{
    return &ctx->scratch->rows[k];
}

/* The forward pass keeps one outer row: scaling a row moves it into the seeds of the
 * next step. */
static double *forward_at(const context *ctx, size_t k, int32_t x)
{
    const row_layout *row = row_of(ctx, k);

    return ctx->scratch->forward + (size_t)(x - row->lo) * N_OUTER;
}

/* The backward pass fills row k - 1 from the step into row k, so it keeps two rows, and
 * row k lies in the slot of its parity. */
static double *backward_at(const context *ctx, size_t k, int32_t x)
{
    const row_layout *row = row_of(ctx, k);

    return ctx->scratch->backward
         + ((k & 1) * (size_t)ctx->widest + (size_t)(x - row->lo)) * N_OUTER;
}

static size_t inner_width(const row_layout *row)
{
    return (size_t)(row->xhi - row->xlo + 1);
}

static size_t inner_cells(const row_layout *row)
{
    return (size_t)(row->depth + 1) * inner_width(row);
}

static double *inner_at(const context *ctx, size_t k, int32_t j, int32_t x)
{
    const row_layout *row = row_of(ctx, k);

    return ctx->scratch->inner + row->inner
         + ((size_t)j * inner_width(row) + (size_t)(x - row->xlo)) * N_INNER;
}

static double *inner_back_at(const context *ctx, size_t k, int32_t j, int32_t x)
{
    const row_layout *row = row_of(ctx, k);

    return ctx->scratch->inner_back
         + ((size_t)j * inner_width(row) + (size_t)(x - row->xlo)) * N_INNER;
}

/* ------------------------------------------------------------------------ */
/* One base against the reference                                            */
/* ------------------------------------------------------------------------ */

/* Returns the chance that read base q was misread. Every base of an Ultima read takes
 * one score, since its quality string describes run lengths. */
static double error_of(const context *ctx, int32_t q)
{
    if (ctx->errors) {
        return phred_error(ctx->quality, (uint8_t)ctx->errors->phred);
    }

    return ctx->read->qual ? phred_error(ctx->quality, ctx->read->qual[q]) : 0.0;
}

/* Writes the misread chance of every placed read base. */
static void fill_errors(const context *ctx)
{
    for (int32_t q = ctx->span.begin; q < ctx->span.end; q++) {
        ctx->scratch->error[q] = error_of(ctx, q);
    }
}

/* What a pairing of a read base against a reference base pays and credits, as lines in
 * the modification rate m of the reference base. The emission is intercept + slope * m.
 * The part of the emission on which the template differed is share * m. Both come
 * from the chances of the base model: a read base agrees with its template base with
 * chance (1 - m)(1 - e) + m e / 3 and differs from it into one given letter with
 * chance (1 - m) e / 3 + m (1 - e / 3) / 3, where e is the chance the base was
 * misread. */
typedef struct {
    double intercept;
    double slope;
    double share;
} outcome;

/* The outcomes of one read base against a reference base that agrees with it, one
 * that differs, and one that is not named. A phantom is a template base the sequencer
 * did not count, whose letter was never observed, so it emits one and carries nothing.
 * An unnamed base on either side emits the uninformative chance and carries nothing. */
typedef struct {
    outcome agree;
    outcome differ;
    outcome unnamed;
} outcomes;

static void set_outcome(outcome *o, double intercept, double slope, double share)
{
    o->intercept = intercept;
    o->slope     = slope;
    o->share     = share;
}

static void set_outcomes(outcomes *o, nuc ours, bool phantom, double error)
{
    double other = error / OTHER_BASES;
    double shift = 1.0 - error - other;
    double blank = phantom ? 1.0 : UNINFORMATIVE;

    set_outcome(&o->unnamed, blank, 0.0, 0.0);

    if (phantom || !nuc_is_base(ours)) {
        set_outcome(&o->agree, blank, 0.0, 0.0);
        set_outcome(&o->differ, blank, 0.0, 0.0);
        return;
    }

    set_outcome(&o->agree, 1.0 - error, -shift, other);
    set_outcome(&o->differ, other, shift / OTHER_BASES, (1.0 - other) / OTHER_BASES);
}

static const outcome *outcome_at(const outcomes *o, const runhmm_base *base, nuc ours)
{
    nuc theirs = (nuc)base->code;

    if (!nuc_is_base(theirs)) {
        return &o->unnamed;
    }

    return theirs == ours ? &o->agree : &o->differ;
}

static double emission_of(const outcome *o, const runhmm_base *base)
{
    return o->intercept + o->slope * base->modification;
}

static double modified_of(const outcome *o, const runhmm_base *base)
{
    return o->share * base->modification;
}

/* ------------------------------------------------------------------------ */
/* Layout                                                                    */
/* ------------------------------------------------------------------------ */

/* Returns where the CIGAR has reached after j bases of a piece: the prefix once the
 * base is taken, and the prefix once the deletion straight after it is taken too. */
static const aln_place *cigar_place(const context *ctx, const run *piece, int32_t j)
{
    size_t prefix = (size_t)(piece->start - ctx->span.begin) + (size_t)j;

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
static void lay_out_window(const context *ctx, row_layout *row)
{
    const aln_place *place = cigar_place(ctx, &row->piece, row->piece.len);

    row->lo = clamp_prefix(ctx, (int32_t)place->first - ctx->model->band);
    row->hi = clamp_prefix(ctx, (int32_t)place->last + ctx->model->band);
}

/* Returns the prefixes inner row j of a step may hold: the band either side of the
 * prefix the CIGAR has reached after j bases of the run, within the step. A phantom
 * row takes the band of the last base. */
static inner_span inner_span_of(const context *ctx, const row_layout *row, int32_t j)
{
    int32_t          real  = j < row->piece.len ? j : row->piece.len;
    const aln_place *place = cigar_place(ctx, &row->piece, real);
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

static const inner_span *span_at(const context *ctx, const row_layout *row, int32_t j)
{
    return &ctx->scratch->spans[row->first_span + (size_t)j];
}

/* Returns whether a piece continues into the clipped bases on one side, which leaves
 * its length unknown. */
static bool cut_before(const context *ctx, const run *piece)
{
    int32_t before = ctx->span.begin - 1;

    return before >= 0 && nuc_is_base(piece->letter)
        && nuc_from_read(ctx->read->seq, before) == piece->letter;
}

static bool cut_after(const context *ctx, const run *piece)
{
    int32_t after = ctx->span.end;

    return after < ctx->read->l_qseq && nuc_is_base(piece->letter)
        && nuc_from_read(ctx->read->seq, after) == piece->letter;
}

/* Returns how many bases past the read run the longest template length with any chance
 * reaches, which is how many phantom rows the step needs. */
static int32_t phantom_rows(const double *length)
{
    for (int offset = MAX_LENGTH_ERROR; offset > 0; offset--) {
        if (length[MAX_LENGTH_ERROR + offset] > 0.0) {
            return offset;
        }
    }

    return 0;
}

/* Writes the chance of each template length of a run. A run whose length the read
 * does not give in full, or a read without the tags, has the read length. */
static void lay_out_lengths(const context *ctx, row_layout *row, size_t k)
{
    const flow_run_errors *errors;
    double                 longer, shorter;
    bool                   cut = (k == 1 && cut_before(ctx, &row->piece))
                              || (k == ctx->n_pieces && cut_after(ctx, &row->piece));

    memset(row->length, 0, sizeof row->length);
    row->depth = row->piece.len;

    if (!ctx->errors || cut || !nuc_is_base(row->piece.letter)) {
        row->length[MAX_LENGTH_ERROR] = 1.0;
        return;
    }

    errors  = &ctx->errors->of_run[k - 1];
    longer  = errors->longer;
    shorter = errors->shorter;

    row->length[MAX_LENGTH_ERROR + 1] = longer;
    row->length[MAX_LENGTH_ERROR - 1] = shorter;
    row->length[MAX_LENGTH_ERROR + 2] = longer * FLOW_TWO_BASE_FRACTION;
    row->length[MAX_LENGTH_ERROR - 2] = shorter * FLOW_TWO_BASE_FRACTION;
    row->length[MAX_LENGTH_ERROR]     = 1.0 - (longer + shorter)
                                            * (1.0 + FLOW_TWO_BASE_FRACTION);
    row->depth = row->piece.len + phantom_rows(row->length);
}

static void lay_out_missed(const context *ctx, row_layout *row, size_t k)
{
    row->missed         = 0.0;
    row->missed_letters = 0;

    if (ctx->errors && k >= 1 && k < ctx->n_pieces) {
        row->missed         = ctx->errors->of_run[k - 1].missed;
        row->missed_letters = ctx->errors->of_run[k - 1].missed_letters;
    }
}

/* Lays out every row and returns the cells the step matrices need in all. */
static size_t lay_out_rows(context *ctx)
{
    row_layout *rows   = ctx->scratch->rows;
    size_t      inner  = 0;
    int32_t     widest = 0;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        rows[k].piece = ctx->scratch->pieces[k - 1];
        lay_out_window(ctx, &rows[k]);
        lay_out_lengths(ctx, &rows[k], k);
        lay_out_missed(ctx, &rows[k], k);
    }

    rows[0].lo    = clamp_prefix(ctx, (int32_t)cigar_place(ctx, &rows[1].piece, 0)->first
                                      - ctx->model->band);
    rows[0].hi    = clamp_prefix(ctx, (int32_t)cigar_place(ctx, &rows[1].piece, 0)->last
                                      + ctx->model->band);
    rows[0].depth = 0;

    ctx->n_spans = 0;

    for (size_t k = 1; k < n_rows(ctx); k++) {
        rows[k].xlo        = rows[k - 1].lo < rows[k].lo ? rows[k - 1].lo : rows[k].lo;
        rows[k].xhi        = rows[k - 1].hi > rows[k].hi ? rows[k - 1].hi : rows[k].hi;
        rows[k].inner      = inner * N_INNER;
        rows[k].first_span = ctx->n_spans;
        inner             += inner_cells(&rows[k]);
        ctx->n_spans      += (size_t)rows[k].depth + 1;
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

/* Writes the span of every inner row of every step. */
static void lay_out_spans(const context *ctx)
{
    for (size_t k = 1; k < n_rows(ctx); k++) {
        const row_layout *row = row_of(ctx, k);

        for (int32_t j = 0; j <= row->depth; j++) {
            ctx->scratch->spans[row->first_span + (size_t)j] = inner_span_of(ctx, row, j);
        }
    }
}

static void lay_out_output(context *ctx)
{
    const row_layout *rows = ctx->scratch->rows;
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

/* Grows every buffer the read needs. The outer buffers are a pair, so they are grown
 * together. */
/* Returns the doubles one outer row holds. */
static size_t outer_cells(const context *ctx)
{
    return (size_t)ctx->widest * N_OUTER;
}

static int grow_buffers(context *ctx, size_t inner)
{
    runhmm_scratch *s     = ctx->scratch;
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
    runhmm_scratch      *s    = ctx->scratch;
    const cm_bam_record *read = ctx->read;
    size_t               inner;

    if (grow((void **)&s->places, &s->cap_places, (size_t)read->l_qseq + 1,
             sizeof *s->places) < 0
        || grow((void **)&s->pieces, &s->cap_pieces, (size_t)read->l_qseq + 1,
                sizeof *s->pieces) < 0
        || grow((void **)&s->rows, &s->cap_rows, (size_t)read->l_qseq + 2,
                sizeof *s->rows) < 0) {
        return PHMM_NO_MEMORY;
    }

    ctx->span     = aln_places(read, s->places);
    ctx->n_pieces = runs_of_read(read->seq, ctx->span.begin, ctx->span.end, s->pieces);

    if (ctx->n_pieces == 0) {
        return PHMM_NO_PATH;
    }

    ctx->errors = NULL;

    if (ctx->flow) {
        if (flow_errors_decode(&s->errors, ctx->flow, read, s->pieces,
                               ctx->n_pieces) < 0) {
            return PHMM_NO_MEMORY;
        }
        ctx->errors = s->errors.tagged ? &s->errors : NULL;
    }

    inner = lay_out_rows(ctx);
    lay_out_output(ctx);

    if (grow_buffers(ctx, inner) < 0
        || grow((void **)&s->error, &s->cap_error, (size_t)read->l_qseq,
                sizeof *s->error) < 0) {
        return PHMM_NO_MEMORY;
    }

    fill_errors(ctx);
    lay_out_spans(ctx);

    memset(s->forward, 0, outer_cells(ctx) * sizeof *s->forward);
    memset(s->inner, 0, inner * N_INNER * sizeof *s->inner);
    memset(s->window, 0, ctx->window_len * sizeof *s->window);

    return PHMM_OK;
}

/* ------------------------------------------------------------------------ */
/* The step into a row                                                       */
/* ------------------------------------------------------------------------ */

/* The weights that hold along every step. A deletion pays its closing once, with the
 * opening rate of the base it closes at. An insertion pays its closing once, as the
 * first inserted base, and the opening rate at the base it closes at. Each inserted
 * base pays the chance of its letter. */
typedef struct {
    double extend_deletion;
    double close_deletion;
    double extend_insertion;  /* per inserted base after the first */
    double open_insertion;    /* the first inserted base */
} step_weights;

static step_weights step_weights_of(const context *ctx)
{
    const phmm_rates *rates = ctx->rates;

    return (step_weights){
        .extend_deletion  = rates->deletion_to_deletion,
        .close_deletion   = rates->match_to_deletion,
        .extend_insertion = rates->insertion_to_insertion * UNINFORMATIVE,
        .open_insertion   = rates->match_to_insertion * UNINFORMATIVE,
    };
}

/* Whether read base j of a piece can be inserted. The first base of the read cannot
 * be, nor can any phantom. */
static bool may_insert(const row_layout *row, size_t k, int32_t j)
{
    return j < row->piece.len && (j >= 1 || k >= 2);
}

/* Whether a deletion can open inside the piece after inner row j. */
static bool may_delete_inside(const row_layout *row, int32_t j)
{
    return j >= 1 && j < row->depth;
}

/* What the loops over one inner row read: where the row and the next lie, and the
 * error of the row's read base. */
typedef struct {
    double *line;      /* the row's first cell, at prefix xlo */
    double *next;      /* the next row's */
    size_t  stride;    /* doubles per inner row */
    int32_t from, to;  /* the prefixes the row may hold */
    outcomes outcomes; /* of the row's read base */
    double  error;
    bool    pairs;     /* whether a row follows to pair into */
    bool    phantom;
    bool    deletes;   /* whether a deletion can open inside after this row */
    bool    inserts;   /* whether this row's read base can be inserted */
} inner_row;

/* Sets r to inner row j of the step into row k, whose matrix begins at matrix. Filled
 * in place: returned by value, the struct is copied through the stack on every inner
 * row, and the copy stalls on store forwarding. */
static void set_inner_row(inner_row *r, const context *ctx, size_t k, double *matrix,
                          int32_t j)
{
    const row_layout *row  = row_of(ctx, k);
    const inner_span *span = span_at(ctx, row, j);

    r->stride  = inner_width(row) * N_INNER;
    r->line    = matrix + (size_t)j * r->stride;
    r->next    = r->line + r->stride;
    r->from    = span->from;
    r->to      = span->to;
    r->pairs   = j < row->depth;
    r->phantom = j >= row->piece.len;
    r->deletes = may_delete_inside(row, j);
    r->inserts = may_insert(row, k, j);
    r->error   = r->phantom ? 0.0 : ctx->scratch->error[row->piece.start + j];
    set_outcomes(&r->outcomes, row->piece.letter, r->phantom, r->error);
}

/* Fills the step matrix of row k from its seeds. The open deletion state passes along
 * a row in a local, so that each cell is read and written once. */
static void step_forward(const context *ctx, size_t k)
{
    const row_layout  *row   = row_of(ctx, k);
    const runhmm_base *bases = ctx->model->bases;
    step_weights       w     = step_weights_of(ctx);
    int32_t            len   = (int32_t)ctx->model->len;
    nuc                ours  = row->piece.letter;

    for (int32_t j = 0; j <= row->depth; j++) {
        inner_row r;
        double    open_deletion = 0.0;  /* what the cell to the left passes on */

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
                const runhmm_base *b = &bases[x];
                const outcome     *o = outcome_at(&r.outcomes, b, ours);

                r.next[at + N_INNER + PAIRED] += (paired * b->keep
                                                + deleted * b->close_deletion
                                                + inserted * b->open_insertion)
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
static bool may_miss(const context *ctx, size_t k)
{
    const row_layout *row = row_of(ctx, k);

    return k != ctx->n_pieces && row->missed > 0.0 && row->missed_letters != 0;
}

/* Whether reference base x is a run of one base that the sequencer can have missed
 * after the row's run. */
static bool missed_at(const context *ctx, const row_layout *row, int32_t x)
{
    const runhmm_model *model = ctx->model;

    return (size_t)x < model->len
        && (model->bases[x].one_base_run & row->missed_letters) != 0;
}

/* Leaves the step at each template length, into the outer cells of row k. A deletion
 * that opens after the run and an insertion left open end in their own states, which
 * the last row cannot hold. */
static void leave_step(const context *ctx, size_t k)
{
    const row_layout *row  = row_of(ctx, k);
    bool              last = k == ctx->n_pieces;

    for (int offset = -MAX_LENGTH_ERROR; offset <= MAX_LENGTH_ERROR; offset++) {
        double  chance = row->length[offset + MAX_LENGTH_ERROR];
        int32_t jend   = row->piece.len + offset;

        if (chance <= 0.0 || jend < 0 || jend > row->depth) {
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

/* Passes over the one-base runs the sequencer missed after the row's run. */
static void pass_missed_forward(const context *ctx, size_t k)
{
    const row_layout *row = row_of(ctx, k);

    if (!may_miss(ctx, k)) {
        return;
    }

    for (int32_t x = row->lo; x < row->hi; x++) {
        if (missed_at(ctx, row, x)) {
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
 * into the seeds of the step into the next row, and clears it for the row after. */
static phmm_status scale_row(const context *ctx, size_t k)
{
    const row_layout *row   = row_of(ctx, k);
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
        double *seed = inner_at(ctx, k + 1, 0, x);

        seed[PAIRED]        = cell[CLOSED] * share;
        seed[GAP_DELETION]  = cell[DELETION] * share;
        seed[GAP_INSERTION] = cell[INSERTION] * share;

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
    const row_layout *row = row_of(ctx, 0);

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
        step_forward(ctx, k);
        leave_step(ctx, k);
        pass_missed_forward(ctx, k);
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

static void pass_missed_backward(const context *ctx, size_t k)
{
    const row_layout *row = row_of(ctx, k);

    if (!may_miss(ctx, k)) {
        return;
    }

    for (int32_t x = row->hi - 1; x >= row->lo; x--) {
        if (missed_at(ctx, row, x)) {
            backward_at(ctx, k, x)[CLOSED] += row->missed
                                            * backward_at(ctx, k, x + 1)[CLOSED];
        }
    }
}

/* Writes the ends of the step matrix of row k: what each leaving cell is worth, from
 * the outer backward cells of row k. */
static void end_step(const context *ctx, size_t k)
{
    const row_layout *row  = row_of(ctx, k);
    bool              last = k == ctx->n_pieces;

    memset(ctx->scratch->inner_back, 0,
           inner_cells(row) * N_INNER * sizeof *ctx->scratch->inner_back);

    for (int offset = -MAX_LENGTH_ERROR; offset <= MAX_LENGTH_ERROR; offset++) {
        double  chance = row->length[offset + MAX_LENGTH_ERROR];
        int32_t jend   = row->piece.len + offset;

        if (chance <= 0.0 || jend < 0 || jend > row->depth) {
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

/* Credits one pairing of reference base x to the window. modified is the posterior
 * of the paths through the pairing on which the template differed, and real whether a
 * read base was observed. */
static void credit_pairing(const context *ctx, size_t k, int32_t j, int32_t x,
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

/* Fills the step matrix of row k backward from its ends, crediting each pairing with
 * its posterior, and sets the outer cells of row k - 1 to what the seeds are worth.
 * The open deletion state passes along a row in a local, so that each cell is read
 * and written once. The posterior of a template difference at a pairing is the
 * posterior of the pairing times the modified part of its emission over the whole
 * emission, so the emission divides out of it. */
static void step_backward(const context *ctx, size_t k)
{
    const row_layout  *row     = row_of(ctx, k);
    const row_layout  *above   = row_of(ctx, k - 1);
    const runhmm_base *bases   = ctx->model->bases;
    step_weights       w       = step_weights_of(ctx);
    int32_t            len     = (int32_t)ctx->model->len;
    nuc                ours    = row->piece.letter;
    double             restore = 1.0 / ctx->scratch->scale[k];
    const double      *forward = ctx->scratch->inner + row->inner;

    for (int32_t j = row->depth; j >= 0; j--) {
        inner_row r;
        double    open_deletion = 0.0;  /* of the cell to the right */

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
                const runhmm_base *b        = &bases[x];
                const outcome     *o        = outcome_at(&r.outcomes, b, ours);
                double             below    = r.next[at + N_INNER + PAIRED];
                double             next     = below * emission_of(o, b);
                double             through  = (fwd[PAIRED] * b->keep
                                               + fwd[GAP_DELETION] * b->close_deletion
                                               + fwd[GAP_INSERTION] * b->open_insertion)
                                              * below * restore;
                double             to_paired   = b->keep * next;
                double             to_deleted  = b->close_deletion * next;
                double             to_inserted = b->open_insertion * next;

                paired   += to_paired;
                deleted  += to_deleted;
                inserted += to_inserted;

                credit_pairing(ctx, k, j, x, !r.phantom, through * modified_of(o, b),
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
        const double *seed = inner_back_at(ctx, k, 0, x);
        double       *to   = backward_at(ctx, k - 1, x);

        to[CLOSED]    = seed[PAIRED] * restore;
        to[DELETION]  = seed[GAP_DELETION] * restore;
        to[INSERTION] = seed[GAP_INSERTION] * restore;
    }
}

static void backward_last_row(const context *ctx)
{
    size_t            k   = ctx->n_pieces;
    const row_layout *row = row_of(ctx, k);

    for (int32_t x = row->lo; x <= row->hi; x++) {
        double *cell = backward_at(ctx, k, x);

        cell[CLOSED]    = 1.0;
        cell[DELETION]  = 0.0;
        cell[INSERTION] = 0.0;
    }
}

/* Returns the first row's total of forward times backward, which is one when the two
 * passes describe the same set of paths. The seeds of the first step hold the forward
 * values of row 0. */
static double agreement(const context *ctx)
{
    const row_layout *row   = row_of(ctx, 0);
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

    for (size_t k = ctx->n_pieces; k >= 1; k--) {
        pass_missed_backward(ctx, k);
        end_step(ctx, k);
        step_backward(ctx, k);
    }

    return passes_agree(ctx);
}

/* ------------------------------------------------------------------------ */
/* Entry                                                                     */
/* ------------------------------------------------------------------------ */

phmm_status runhmm_run(const runhmm_model *model, const phred *quality,
                       const flow_tables *flow, const cm_bam_record *read,
                       runhmm_scratch *scratch, phmm_window *out)
{
    context ctx = {
        .model   = model,
        .rates   = model->rates,
        .quality = quality,
        .read    = read,
        .scratch = scratch,
        .flow    = flow,
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
