/* phmm.c -- a banded pair HMM over the reference and the cDNA.
 *
 * The matrix has one row per node of the sequence graph of the read. A cell of a row
 * is a prefix length of the reference within the band of the row, which is the band
 * around the CIGAR at the read prefix the node stands at, and holds three states: the
 * last event was a pairing, an insertion is open, or a deletion is open.
 *
 * A row is filled from the row before it by one step of the reverse transcriptase over
 * the base edge of its node: a pairing of the base with the next reference base, which
 * pays the weight of that base after no gap or the opening rate of the gap it closes,
 * times the emission of the comparison; an insertion of the base, which pays the
 * opening or the extension of the insertion; and, within the row, a deletion of the
 * next reference base, which pays its extension per base after the first. Every move
 * of the step carries the weight of the base edge. A row is also filled over each skip
 * edge of its node, which carries all three states of the cell at the same prefix of
 * the node the edge leaves, at the weight of the edge. A deletion beside a skip edge
 * is placed before the edge, in the row the edge leaves, so that each path through
 * the graph is one alignment of one cDNA.
 *
 * The reverse transcriptase copies the reference 3' to 5', so the read's 5'-most
 * pairing is where it reached the end state. The passes run 5' to 3', so the end state
 * is where every path starts. It is held as one weight per row, the start, which is
 * one on the source and reaches a later row over skip edges alone; the pairing out of
 * it pays the termination rate of the base before in place of the weight after no gap,
 * since that is the base the enzyme stopped at instead of copying it. Stopping is one
 * arm of the decision at each later pairing, beside the openings of the two gaps.
 *
 * The forward pass keeps every row. Every SCALE_EVERY rows, a row whose total leaves a
 * wide range is divided by it, except a row that a skip edge passes over, so that the
 * rows a later row reads carry one scale. The backward pass keeps only the rows a row
 * can reach. It applies to each row the weights of the edges into it and one over the
 * scale of its forward row, and to the last row one over the total of the read, so
 * that a forward value times a transition times the backward value it reaches is the
 * posterior of that transition, which the pass adds to the window as it goes.
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

/* Emission of a comparison that carries no information. */
#define UNINFORMATIVE (1.0 / NUC_BASES)

/* A misread or modified base is equally likely to be any base the reference does not
 * hold. */
#define OTHER_BASES ((double)(NUC_BASES - 1))

/* The posteriors of the read's 5'-most pairing sum to one within this tolerance. */
#define NORMALIZATION_TOLERANCE 1e-6

/* A row is divided by its total only when the total leaves this range. The range keeps
 * the totals far from the smallest and the largest double. */
#define SCALE_FLOOR   1e-50
#define SCALE_CEILING 1e50

/* How often the rows are checked. The paths through this many rows stay far above the
 * smallest double. */
#define SCALE_EVERY 32

/* The backward pass keeps the rows a row can reach: the row after it and the rows
 * its skip edges land on. */
#define BACKWARD_ROWS (GRAPH_SKIPS + 1)

/* The three states of one cell. */
typedef struct {
    double paired;
    double inserted;
    double deleted;
} cell;

void phmm_rates_set_transitions(phmm_rates *rates, const phmm_uniform_rates *uniform)
{
    rates->match_to_insertion     = 1.0 - uniform->extend_insertion;
    rates->match_to_deletion      = 1.0 - uniform->extend_deletion;
    rates->insertion_to_insertion = uniform->extend_insertion;
    rates->deletion_to_deletion   = uniform->extend_deletion;
}

/* ------------------------------------------------------------------------ */
/* The model of reverse transcription                                        */
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
               - rate_at(rates->termination, len, b - 1);
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

static phmm_base base_at(const phmm_rates *rates, const cm_fasta_record *ref,
                         hts_pos_t b)
{
    size_t len = ref->len;

    return (phmm_base){
        .after_no_gap    = after_no_gap_at(rates, len, b),
        .after_deletion  = rate_at(rates->open_deletion, len, b - 1)
                         * rates->match_to_deletion,
        .after_insertion = rate_at(rates->open_insertion, len, b),
        .termination     = rate_at(rates->termination, len, b - 1),
        .modification    = rates->modification[b],
        .code            = (uint8_t)nuc_from_char(ref->seq[b]),
    };
}

int phmm_model_prepare(phmm_model *model, const phmm_rates *rates,
                       const cm_fasta_record *ref, int band)
{
    if (grow((void **)&model->bases, &model->cap_bases, ref->len,
             sizeof *model->bases) < 0) {
        return -1;
    }

    model->rates = rates;
    model->len   = ref->len;
    model->band  = band;

    for (size_t b = 0; b < ref->len; b++) {
        model->bases[b] = base_at(rates, ref, (hts_pos_t)b);
    }

    return 0;
}

void phmm_model_free(phmm_model *model)
{
    free(model->bases);
    *model = (phmm_model){ 0 };
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
/* The emissions                                                             */
/* ------------------------------------------------------------------------ */

/* The comparison of one base of z, as the read gives it, with each reference base, by
 * the code of the reference base, as lines in the modification rate m of the reference
 * base. The emission is intercept + slope * m. The part of the emission on which the
 * base of z differed from the reference is share * m. The base agrees with the
 * reference base with chance (1 - m)(1 - e) + m e / 3, and differs from it into one
 * given base with chance (1 - m) e / 3 + m (1 - e / 3) / 3, where e is the chance the
 * base was misread. An unnamed base on either side emits the uninformative chance and
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

/* Writes the comparisons of a read base ours that was misread with chance error. */
static void set_comparisons(comparisons *o, nuc ours, double error)
{
    double other = error / OTHER_BASES;
    double shift = 1.0 - error - other;

    for (nuc theirs = NUC_N; theirs < NUC_COUNT; theirs++) {
        set_comparison(o, theirs, UNINFORMATIVE, 0.0, 0.0);
    }

    if (!nuc_is_base(ours)) {
        return;
    }

    for (nuc theirs = NUC_A; theirs < NUC_COUNT; theirs++) {
        set_comparison(o, theirs, other, shift / OTHER_BASES,
                       (1.0 - other) / OTHER_BASES);
    }

    set_comparison(o, ours, 1.0 - error, -shift, other);
}

/* The comparisons of a read base of each base at each quality score, and of such a
 * base in a read that stores no qualities, which is never misread. */
struct phmm_emissions {
    comparisons at_quality[NUC_COUNT][PHRED_MAX + 1];
    comparisons without_quality[NUC_COUNT];
};

phmm_emissions *phmm_emissions_create(const phred *quality)
{
    phmm_emissions *emissions = malloc(sizeof *emissions);

    if (!emissions) {
        return NULL;
    }

    for (nuc ours = NUC_N; ours < NUC_COUNT; ours++) {
        for (unsigned q = 0; q <= PHRED_MAX; q++) {
            set_comparisons(&emissions->at_quality[ours][q], ours,
                            phred_error(quality, (uint8_t)q));
        }

        set_comparisons(&emissions->without_quality[ours], ours, 0.0);
    }

    return emissions;
}

void phmm_emissions_destroy(phmm_emissions *emissions)
{
    free(emissions);
}

/* Returns the comparisons of a base of a node. */
static const comparisons *comparisons_of(const phmm_emissions *emissions,
                                         const graph_node *node)
{
    if (node->quality == GRAPH_NO_QUALITY) {
        return &emissions->without_quality[node->base];
    }

    return &emissions->at_quality[node->base][(uint8_t)node->quality];
}

/* Returns the emission of the comparison with reference base base. */
static inline double emission_of(const comparisons *o, const phmm_base *base)
{
    return o->intercept[base->code] + o->slope[base->code] * base->modification;
}

/* Returns the part of that emission on which the base of z differed. */
static inline double modified_of(const comparisons *o, const phmm_base *base)
{
    return o->share[base->code] * base->modification;
}

/* ------------------------------------------------------------------------ */
/* Scratch                                                                   */
/* ------------------------------------------------------------------------ */

/* One row of the matrix: the comparisons of its node's base, NULL on the source; the
 * index of its first cell in the forward matrix; its band, from lo to hi; the weight
 * of the start state in it; and one over the total it was divided by in the forward
 * pass. */
typedef struct {
    const comparisons *comparisons;
    size_t             first;
    int32_t            lo;
    int32_t            hi;
    double             start;
    double             restore;
} matrix_row;

struct phmm_scratch {
    aln_place     *places;
    size_t         cap_places;
    matrix_row    *rows;         /* one per node */
    size_t         cap_rows;
    cell          *forward;      /* every row */
    size_t         cap_forward;
    cell          *backward;     /* BACKWARD_ROWS rows of the widest row */
    size_t         cap_backward;
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
    free(scratch->forward);
    free(scratch->backward);
    free(scratch->window);
    free(scratch);
}

/* ------------------------------------------------------------------------ */
/* One read                                                                  */
/* ------------------------------------------------------------------------ */

typedef struct {
    const phmm_model     *model;
    const phmm_emissions *emissions;
    const cm_bam_record  *read;
    const graph          *g;
    phmm_scratch         *scratch;

    gap_weights gaps;         /* the weights of the gaps under the rates */
    size_t      n_rows;
    size_t      cells;        /* cells of every row */
    int32_t     widest;       /* cells of the widest row */
    int32_t     lo;           /* the lowest prefix of any row */
    int32_t     hi;           /* and the highest */
    int32_t     origin;       /* reference position of window value 0 */
    size_t      window_len;
    double      ends;         /* the posteriors of the read's 5'-most pairing, summed */
} context;

static inline const graph_node *node_of(const context *ctx, size_t r)
{
    return &ctx->g->nodes[r];
}

static inline matrix_row *row_of(const context *ctx, size_t r)
{
    return &ctx->scratch->rows[r];
}

/* Returns the forward cells of row r, the first at the first prefix of the row. */
static inline cell *forward_cells_of(const context *ctx, size_t r)
{
    return ctx->scratch->forward + row_of(ctx, r)->first;
}

/* Returns the backward cells of row r, the first at the first prefix of the row. The
 * pass keeps BACKWARD_ROWS rows, and a row lies in the slot of its index modulo their
 * count. */
static inline cell *backward_cells_of(const context *ctx, size_t r)
{
    return ctx->scratch->backward + (r % BACKWARD_ROWS) * (size_t)ctx->widest;
}

/* The cells of one row as a step reads them, with the band they cover. */
typedef struct {
    const cell *cells;
    int32_t     lo;
    int32_t     hi;
} row_view;

/* Returns a view that holds no cell. */
static row_view empty_view(void)
{
    return (row_view){ .cells = NULL, .lo = 0, .hi = -1 };
}

static row_view forward_view(const context *ctx, size_t r)
{
    const matrix_row *row = row_of(ctx, r);

    return (row_view){ .cells = forward_cells_of(ctx, r), .lo = row->lo, .hi = row->hi };
}

static row_view backward_view(const context *ctx, size_t r)
{
    const matrix_row *row = row_of(ctx, r);

    return (row_view){ .cells = backward_cells_of(ctx, r), .lo = row->lo, .hi = row->hi };
}

static inline bool view_holds(const row_view *view, int32_t x)
{
    return x >= view->lo && x <= view->hi;
}

static inline const cell *view_cell(const row_view *view, int32_t x)
{
    return view->cells + (x - view->lo);
}

/* ------------------------------------------------------------------------ */
/* Layout                                                                    */
/* ------------------------------------------------------------------------ */

/* The band around one point of the CIGAR path, within the reference. */
typedef struct {
    int32_t from;
    int32_t to;
} band_edges;

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

/* Returns the band around the point of the CIGAR path at a prefix of the placed span.
 * It runs from the band before the prefix reached once the base is taken, to the band
 * after the prefix reached once the deletion straight after the base is taken too. */
static band_edges band_at(const context *ctx, int32_t prefix)
{
    const aln_place *place = &ctx->scratch->places[prefix];

    return (band_edges){
        .from = clamp_prefix(ctx, (int32_t)place->first - ctx->model->band),
        .to   = clamp_prefix(ctx, (int32_t)place->last + ctx->model->band),
    };
}

/* Lays out row r: its band, around the prefix its node stands at; its cells, after
 * those of the rows before it; and its comparisons. */
static void lay_out_row(context *ctx, size_t r)
{
    const graph_node *node  = node_of(ctx, r);
    band_edges        band  = band_at(ctx, node->prefix);
    int32_t           width = band.to - band.from + 1;

    *row_of(ctx, r) = (matrix_row){
        .comparisons = r > 0 ? comparisons_of(ctx->emissions, node) : NULL,
        .first       = ctx->cells,
        .lo          = band.from,
        .hi          = band.to,
        .start       = r == 0 ? 1.0 : 0.0,
        .restore     = 1.0,
    };

    ctx->cells += (size_t)width;
    ctx->widest = width > ctx->widest ? width : ctx->widest;
    ctx->lo     = band.from < ctx->lo ? band.from : ctx->lo;
    ctx->hi     = band.to > ctx->hi ? band.to : ctx->hi;
}

/* Lays out every row, and the window that covers them all. */
static void lay_out_rows(context *ctx)
{
    ctx->cells  = 0;
    ctx->widest = 0;
    ctx->lo     = INT32_MAX;
    ctx->hi     = INT32_MIN;

    for (size_t r = 0; r < ctx->n_rows; r++) {
        lay_out_row(ctx, r);
    }

    ctx->origin     = ctx->lo - 1;
    ctx->window_len = (size_t)(ctx->hi - ctx->origin) + 1;
}

/* Grows the buffers of the passes to the layout. */
static int grow_buffers(context *ctx)
{
    phmm_scratch *s = ctx->scratch;

    if (grow((void **)&s->forward, &s->cap_forward, ctx->cells, sizeof *s->forward) < 0
        || grow((void **)&s->backward, &s->cap_backward,
                BACKWARD_ROWS * (size_t)ctx->widest, sizeof *s->backward) < 0
        || grow((void **)&s->window, &s->cap_window, ctx->window_len,
                sizeof *s->window) < 0) {
        return -1;
    }

    return 0;
}

/* Places the CIGAR, lays out the matrix over the graph, and clears the window. */
static phmm_status prepare(context *ctx)
{
    phmm_scratch        *s    = ctx->scratch;
    const cm_bam_record *read = ctx->read;

    if (grow((void **)&s->places, &s->cap_places, (size_t)read->l_qseq + 1,
             sizeof *s->places) < 0
        || grow((void **)&s->rows, &s->cap_rows, ctx->n_rows, sizeof *s->rows) < 0) {
        return PHMM_NO_MEMORY;
    }

    aln_places(read, s->places);
    lay_out_rows(ctx);

    if (grow_buffers(ctx) < 0) {
        return PHMM_NO_MEMORY;
    }

    memset(s->window, 0, ctx->window_len * sizeof *s->window);
    ctx->ends = 0.0;

    return PHMM_OK;
}

/* ------------------------------------------------------------------------ */
/* The moves of a step                                                       */
/* ------------------------------------------------------------------------ */

/* Returns the weight of the paths from cell from into the pairing of reference base b,
 * where start is the weight of the end state in the row of from. */
static inline double pairing_from(const cell *from, const phmm_base *b, double start)
{
    return from->paired * b->after_no_gap + from->deleted * b->after_deletion
         + from->inserted * b->after_insertion + start * b->termination;
}

/* Returns the weight of the paths from cell from into an insertion. */
static inline double insertion_from(const cell *from, const gap_weights *w)
{
    return from->paired * w->open_insertion + from->inserted * w->extend_insertion;
}

/* Adds the states of cell from to cell to over a skip edge of the given weight. */
static inline void add_skipped(cell *to, const cell *from, double weight)
{
    to->paired   += weight * from->paired;
    to->inserted += weight * from->inserted;
    to->deleted  += weight * from->deleted;
}

/* ------------------------------------------------------------------------ */
/* Forward                                                                   */
/* ------------------------------------------------------------------------ */

/* What one step of the forward pass reads: the reference, the comparisons of the
 * node's base, the gap weights, the row before with its start, the row the skip edge
 * of two nodes back leaves, and the weights of the three edges into the node. */
typedef struct {
    const phmm_base   *bases;
    const comparisons *o;
    gap_weights        w;
    row_view           from;
    row_view           far;
    double             start;
    double             weight;
    double             skip_one;
    double             skip_two;
} forward_step;

/* Returns the cell at prefix x. deletion carries the open deletion of the step along
 * the row: what it is worth at x on the way in, and at x + 1 on the way out, which the
 * pairing of the step opens and the skip edges leave alone. */
static inline cell forward_cell(const forward_step *s, int32_t x, double *deletion)
{
    cell   c      = { .deleted = *deletion };
    double paired = 0.0;

    if (view_holds(&s->from, x - 1)) {
        const phmm_base *b = &s->bases[x - 1];

        paired = s->weight * pairing_from(view_cell(&s->from, x - 1), b, s->start)
               * emission_of(s->o, b);
    }

    *deletion = *deletion * s->w.extend_deletion + paired;
    c.paired  = paired;

    if (view_holds(&s->from, x)) {
        const cell *same = view_cell(&s->from, x);

        c.inserted = s->weight * insertion_from(same, &s->w);
        add_skipped(&c, same, s->skip_one);
    }

    if (view_holds(&s->far, x)) {
        add_skipped(&c, view_cell(&s->far, x), s->skip_two);
    }

    return c;
}

/* Fills row r from the row before it and from the row its further skip edge leaves,
 * and carries the start state over the skip edges. */
static void forward_row(context *ctx, size_t r)
{
    const graph_node *node   = node_of(ctx, r);
    const matrix_row *before = row_of(ctx, r - 1);
    matrix_row       *row    = row_of(ctx, r);
    bool              far    = node->skip[1] > 0.0;
    forward_step      s      = {
        .bases    = ctx->model->bases,
        .o        = row->comparisons,
        .w        = ctx->gaps,
        .from     = forward_view(ctx, r - 1),
        .far      = far ? forward_view(ctx, r - 2) : empty_view(),
        .start    = before->start,
        .weight   = node->weight,
        .skip_one = node->skip[0],
        .skip_two = node->skip[1],
    };
    cell             *to       = forward_cells_of(ctx, r);
    double            deletion = 0.0;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        to[x - row->lo] = forward_cell(&s, x, &deletion);
    }

    row->start = node->skip[0] * before->start
               + (far ? node->skip[1] * row_of(ctx, r - 2)->start : 0.0);
}

/* Returns the total of row r. */
static double row_total(const context *ctx, size_t r)
{
    const matrix_row *row   = row_of(ctx, r);
    const cell       *cells = forward_cells_of(ctx, r);
    double            total = 0.0;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        const cell *c = &cells[x - row->lo];

        total += c->paired + c->inserted + c->deleted;
    }

    return total;
}

/* Maps the total of a row that every path passes through to a status. A negative,
 * infinite, or NaN total is a fault. A zero or subnormal total means the band admits
 * no path. */
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

/* Divides row r and its start by its total, keeping one over it. */
static void divide_row(const context *ctx, size_t r, double total)
{
    matrix_row *row   = row_of(ctx, r);
    cell       *cells = forward_cells_of(ctx, r);
    double      share = 1.0 / total;

    row->restore  = share;
    row->start   *= share;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        cells[x - row->lo].paired   *= share;
        cells[x - row->lo].inserted *= share;
        cells[x - row->lo].deleted  *= share;
    }
}

/* Checks the total of row r, and divides the row by it where it has left the range the
 * rows are kept in. A row that a skip edge passes over is left as it is, since the row
 * the edge lands on reads it beside the row before, and the two must carry one scale. */
static phmm_status scale_row(const context *ctx, size_t r)
{
    double      total;
    phmm_status status;

    if (graph_node_passed_over(ctx->g, r)) {
        return PHMM_OK;
    }

    total  = row_total(ctx, r);
    status = total_status(total);

    if (status == PHMM_OK && (total < SCALE_FLOOR || total > SCALE_CEILING)) {
        divide_row(ctx, r, total);
    }

    return status;
}

/* Returns the probability of the read relative to the scales: the paired state of the
 * last row summed, since every path ends with a pairing. */
static double read_total(const context *ctx)
{
    size_t            r     = ctx->n_rows - 1;
    const matrix_row *row   = row_of(ctx, r);
    const cell       *cells = forward_cells_of(ctx, r);
    double            total = 0.0;

    for (int32_t x = row->lo; x <= row->hi; x++) {
        total += cells[x - row->lo].paired;
    }

    return total;
}

/* Clears the cells of the source, whose start alone holds weight. */
static void forward_source(const context *ctx)
{
    const matrix_row *row = row_of(ctx, 0);

    memset(forward_cells_of(ctx, 0), 0, (size_t)(row->hi - row->lo + 1) * sizeof(cell));
}

/* Fills every row and leaves the probability of the read in total. The rows are
 * checked and scaled every SCALE_EVERY rows, which keeps them far inside the range of
 * a double; a fault or an empty row before a check reaches the check as a total that
 * is not finite or is zero. The last row is checked but never divided. */
static phmm_status forward(context *ctx, double *total)
{
    forward_source(ctx);

    for (size_t r = 1; r < ctx->n_rows; r++) {
        forward_row(ctx, r);

        if (r % SCALE_EVERY == 0 && r + 1 < ctx->n_rows) {
            phmm_status status = scale_row(ctx, r);

            if (status != PHMM_OK) {
                return status;
            }
        }
    }

    *total = read_total(ctx);

    return total_status(*total);
}

/* ------------------------------------------------------------------------ */
/* Backward, with the posteriors                                             */
/* ------------------------------------------------------------------------ */

/* What one step of the backward pass reads, for row a and the row b after it: the
 * reference, the comparisons of the base of b, the gap weights, the forward cells of a
 * with its start, the backward cells of b and of the row the further skip edge into
 * the row after b leaves, the weights of the edges into those rows, each times one
 * over the scale of its row, and the window the posteriors go to. */
typedef struct {
    const phmm_base   *bases;
    const comparisons *o;
    gap_weights        w;
    row_view           forward;
    row_view           next;
    row_view           far;
    double             start;
    double             weight;
    double             skip_one;
    double             skip_two;
    phmm_position     *window;
    int32_t            origin;
} backward_step;

/* What a cell of the next row is worth once the weight of its base edge and the scale
 * of its row are applied and its open deletion is extended from the right. The pass
 * carries the values of the cell one prefix on as it sweeps the cells of a from the
 * right. */
typedef struct {
    double paired;
    double deleted;
} completed;

/* Completes the cell of the next row at prefix x, whose cell one prefix on k holds on
 * the way in. Returns the completed insertion state of the cell. */
static inline double complete_next(const backward_step *s, int32_t x, completed *k)
{
    const cell *c       = view_cell(&s->next, x);
    double      deleted = s->weight * c->deleted + s->w.extend_deletion * k->deleted;

    k->paired  = s->weight * c->paired + k->deleted;
    k->deleted = deleted;

    return s->weight * c->inserted;
}

/* The posterior of one pairing of a reference base: the paths through it, the paths
 * through it on which the base of z differed, the paths that close an insertion or a
 * deletion at it, and the paths on which it is the read's 5'-most. */
typedef struct {
    double covered;
    double modified;
    double inserted;
    double deleted;
    double ended;
} pairing_posterior;

/* Returns the posterior of the pairing of reference base x out of cell fwd of row a,
 * where landing is the completed paired state of the cell the pairing lands in. The
 * posterior of a difference is the posterior of the pairing times the modified part
 * of its emission over the whole emission, so the emission divides out. The posterior
 * of the read ending at the base is multiplied in the order of the coverage, so that
 * the two are the same number where nothing else reaches the pairing. */
static inline pairing_posterior posterior_of(const backward_step *s, int32_t x,
                                             const cell *fwd, double landing)
{
    const phmm_base *b        = &s->bases[x];
    double           emission = emission_of(s->o, b);
    double           reach    = landing * emission;
    double           no_gap   = fwd->paired * b->after_no_gap;
    double           closed   = fwd->deleted * b->after_deletion;
    double           joined   = fwd->inserted * b->after_insertion;
    double           started  = s->start * b->termination;
    double           through  = (no_gap + closed + joined + started) * landing;

    return (pairing_posterior){
        .covered  = through * emission,
        .modified = through * modified_of(s->o, b),
        .inserted = joined * reach,
        .deleted  = closed * reach,
        .ended    = (started * landing) * emission,
    };
}

/* Accumulates the posterior of a pairing of reference base x into the window. The
 * deletion the pairing closes and the stop it makes both belong to the base before,
 * which the enzyme did not copy. */
static inline void accumulate_pairing(const backward_step *s, int32_t x,
                                      const pairing_posterior *p)
{
    phmm_position *at = s->window + (x - s->origin);

    at->coverage     += p->covered;
    at->mismatches   += p->modified;
    at->insertions   += p->inserted;
    at[-1].deletions += p->deleted;
    at[-1].ends      += p->ended;
}

/* Returns the cell of row a at prefix x: what the row reaches through the row after it
 * by the pairing of reference base x into the cell one prefix on, whose completed
 * values k holds, by the insertion into the cell at the same prefix, and over the skip
 * edges that leave the cell at the same prefix. The pairing is accumulated with its
 * posterior. k is left holding the completed cell at x. */
static inline cell backward_cell(const backward_step *s, int32_t x, completed *k,
                                 double *ends)
{
    cell c = { 0 };

    if (view_holds(&s->next, x + 1)) {
        const phmm_base  *b     = &s->bases[x];
        double            reach = k->paired * emission_of(s->o, b);
        pairing_posterior p     = posterior_of(s, x, view_cell(&s->forward, x),
                                               k->paired);

        accumulate_pairing(s, x, &p);
        *ends += p.ended;

        c.paired   = b->after_no_gap * reach;
        c.deleted  = b->after_deletion * reach;
        c.inserted = b->after_insertion * reach;
    }

    if (view_holds(&s->next, x)) {
        const cell *same = view_cell(&s->next, x);
        double      open = complete_next(s, x, k);

        c.paired   += s->w.open_insertion * open;
        c.inserted += s->w.extend_insertion * open;
        add_skipped(&c, same, s->skip_one);
    } else {
        *k = (completed){ 0.0, 0.0 };
    }

    if (view_holds(&s->far, x)) {
        add_skipped(&c, view_cell(&s->far, x), s->skip_two);
    }

    return c;
}

/* Fills the backward cells of row r from those of the rows after it, from the right,
 * and accumulates the posteriors of the pairings out of r. The cells of the next row
 * past the band of r are completed first, since the open deletion of each extends
 * from the one after it. */
static void backward_row(context *ctx, size_t r)
{
    const matrix_row *row  = row_of(ctx, r);
    const matrix_row *next = row_of(ctx, r + 1);
    const graph_node *node = node_of(ctx, r + 1);
    bool              far  = r + 2 < ctx->n_rows && node_of(ctx, r + 2)->skip[1] > 0.0;
    backward_step     s    = {
        .bases    = ctx->model->bases,
        .o        = next->comparisons,
        .w        = ctx->gaps,
        .forward  = forward_view(ctx, r),
        .next     = backward_view(ctx, r + 1),
        .far      = far ? backward_view(ctx, r + 2) : empty_view(),
        .start    = row->start,
        .weight   = node->weight * next->restore,
        .skip_one = node->skip[0] * next->restore,
        .skip_two = far ? node_of(ctx, r + 2)->skip[1] * row_of(ctx, r + 2)->restore
                        : 0.0,
        .window   = ctx->scratch->window,
        .origin   = ctx->origin,
    };
    cell             *to   = backward_cells_of(ctx, r);
    completed         k    = { 0.0, 0.0 };

    for (int32_t x = next->hi; x > row->hi && x >= next->lo; x--) {
        complete_next(&s, x, &k);
    }

    for (int32_t x = row->hi; x >= row->lo; x--) {
        to[x - row->lo] = backward_cell(&s, x, &k, &ctx->ends);
    }
}

/* Fills the backward cells of the last row: every path ends with a pairing, worth one
 * over the probability of the read. */
static void backward_last_row(const context *ctx, double total)
{
    size_t            r     = ctx->n_rows - 1;
    const matrix_row *row   = row_of(ctx, r);
    cell             *cells = backward_cells_of(ctx, r);

    for (int32_t x = row->lo; x <= row->hi; x++) {
        cells[x - row->lo] = (cell){ .paired = 1.0 / total };
    }
}

/* The posteriors of the read's 5'-most pairing sum to one when the two passes describe
 * the same set of paths. A sum that is not finite means the backward pass overflowed:
 * every forward row is kept in range, but a backward value is the probability of the
 * rest of the read from a cell relative to the probability of the read itself, and on
 * a read whose every path is astronomically unlikely, as under rates set to zero and a
 * band too narrow for its gaps, that ratio leaves the range of a double. The read has
 * no path worth scoring, and is rejected as one with none. A finite sum away from one
 * means the two passes described different sets of paths, which is a fault. */
static phmm_status passes_agree(const context *ctx)
{
    if (!isfinite(ctx->ends)) {
        return PHMM_NO_PATH;
    }

    return fabs(ctx->ends - 1.0) < NORMALIZATION_TOLERANCE ? PHMM_OK : PHMM_UNSOUND;
}

static phmm_status backward(context *ctx, double total)
{
    backward_last_row(ctx, total);

    for (size_t r = ctx->n_rows - 1; r-- > 0;) {
        backward_row(ctx, r);
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
                     const cm_bam_record *read, const graph *g, phmm_scratch *scratch,
                     phmm_window *out)
{
    context ctx = {
        .model     = model,
        .emissions = emissions,
        .read      = read,
        .g         = g,
        .scratch   = scratch,
        .gaps      = gap_weights_of(model->rates),
        .n_rows    = g->n_nodes,
    };
    phmm_status status;
    double      total = 0.0;

    if (g->n_nodes < 2) {
        return PHMM_NO_PATH;
    }

    status = prepare(&ctx);

    if (status != PHMM_OK) {
        return status;
    }

    status = forward(&ctx, &total);

    if (status != PHMM_OK) {
        return status;
    }

    status = backward(&ctx, total);

    if (status != PHMM_OK) {
        return status;
    }

    out->origin = ctx.origin;
    out->len    = ctx.window_len;
    out->at     = scratch->window;

    return PHMM_OK;
}
