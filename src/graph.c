/* graph.c -- the sequence graph of the cDNA a read admits.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "graph.h"

#include <stdbool.h>
#include <stdlib.h>

/* The chances of the lengths of a run whose length the read gives in full. */
static const double READ_LENGTH[GRAPH_LENGTH_OFFSETS] = {
    [GRAPH_MAX_LENGTH_ERROR] = 1.0,
};

/* The chances of the three lengths of one run of z. */
typedef struct {
    double shorter;
    double same;
    double longer;
} run_lengths;

typedef struct {
    graph                   *g;
    const cm_bam_record     *read;
    aln_span                 span;
    const run               *runs;
    size_t                   n_runs;
    const graph_read_errors *errors;
} context;

static int grow(graph *g, size_t needed)
{
    graph_node *nodes;

    if (needed <= g->cap) {
        return 0;
    }

    nodes = realloc(g->nodes, needed * sizeof *nodes);

    if (!nodes) {
        return -1;
    }

    g->nodes = nodes;
    g->cap   = needed;

    return 0;
}

/* Appends a node entered by its base edge at weight, and by a skip edge from the node
 * before it and from the node before that at the two further weights. */
static void add_node(graph *g, int32_t prefix, nuc base, int quality, double weight,
                     double skip_one, double skip_two)
{
    g->nodes[g->n_nodes++] = (graph_node){
        .prefix  = prefix,
        .base    = base,
        .quality = (int16_t)quality,
        .weight  = weight,
        .skip    = { skip_one, skip_two },
    };
}

/* Appends a node that z holds for certain. */
static void add_certain_node(graph *g, int32_t prefix, nuc base, int quality)
{
    add_node(g, prefix, base, quality, 1.0, 0.0, 0.0);
}

/* Returns the read prefix the base at offset q of the read stands at. */
static int32_t prefix_of(const context *ctx, int32_t q)
{
    return q + 1 - ctx->span.begin;
}

/* ------------------------------------------------------------------------ */
/* A sequencer that reports every base                                       */
/* ------------------------------------------------------------------------ */

/* Returns the quality the base at offset q of a read without run-length errors was
 * read at. */
static int quality_of(const context *ctx, int32_t q)
{
    if (!ctx->read->qual) {
        return GRAPH_NO_QUALITY;
    }

    return ctx->read->qual[q];
}

/* Lays out the chain of a read whose every base of z the sequencer reported. */
static void lay_out_chain(const context *ctx)
{
    for (int32_t q = ctx->span.begin; q < ctx->span.end; q++) {
        add_certain_node(ctx->g, prefix_of(ctx, q), nuc_from_read(ctx->read->seq, q),
                         quality_of(ctx, q));
    }
}

/* ------------------------------------------------------------------------ */
/* A sequencer that measures run lengths                                     */
/* ------------------------------------------------------------------------ */

/* Returns whether a read run continues into the clipped bases before the placed
 * span, which leaves its length unknown. */
static bool clipped_before(const context *ctx, const run *read_run)
{
    int32_t before = ctx->span.begin - 1;

    return before >= 0 && nuc_is_base(read_run->base)
        && nuc_from_read(ctx->read->seq, before) == read_run->base;
}

/* Returns whether a read run continues into the clipped bases after the placed span. */
static bool clipped_after(const context *ctx, const run *read_run)
{
    int32_t after = ctx->span.end;

    return after < ctx->read->l_qseq && nuc_is_base(read_run->base)
        && nuc_from_read(ctx->read->seq, after) == read_run->base;
}

/* Returns whether run k has the length the read shows for certain: a run cut by the
 * clip at either end of the placed span, or a run of unnamed bases. */
static bool length_is_certain(const context *ctx, size_t k)
{
    const run *read_run = &ctx->runs[k];

    return (k == 0 && clipped_before(ctx, read_run))
        || (k + 1 == ctx->n_runs && clipped_after(ctx, read_run))
        || !nuc_is_base(read_run->base);
}

/* Returns the chances of the lengths of run k. */
static run_lengths lengths_of(const context *ctx, size_t k)
{
    const double *length = length_is_certain(ctx, k) ? READ_LENGTH
                                                     : ctx->errors->of_run[k].length;

    return (run_lengths){
        .shorter = length[GRAPH_MAX_LENGTH_ERROR - 1],
        .same    = length[GRAPH_MAX_LENGTH_ERROR],
        .longer  = length[GRAPH_MAX_LENGTH_ERROR + 1],
    };
}

/* Appends the last node of a run that can be one base longer than the read shows. The
 * node holds the extra base, entered at the chance of that length, and is entered
 * without it from the node that ends the read length and from the node that ends the
 * shorter length. */
static void add_longer_run_end(const context *ctx, const run *read_run,
                               run_lengths lengths)
{
    int32_t after = prefix_of(ctx, read_run->start + read_run->len - 1);

    add_node(ctx->g, after, read_run->base, GRAPH_NO_QUALITY, lengths.longer,
             lengths.same, lengths.shorter);
}

/* Appends the last node of a run that cannot be longer than the read shows. The node
 * holds the last read base, entered at the chance of the read length, and is entered
 * without it from the node that ends the shorter length. */
static void add_run_end(const context *ctx, const run *read_run, run_lengths lengths)
{
    int32_t after = prefix_of(ctx, read_run->start + read_run->len - 1);

    add_node(ctx->g, after, read_run->base, GRAPH_NO_QUALITY, lengths.same,
             lengths.shorter, 0.0);
}

/* Lays out run k: one node per base the run of z can hold, the last of them weighing
 * the lengths. */
static void lay_out_run(const context *ctx, size_t k)
{
    const run  *read_run = &ctx->runs[k];
    run_lengths lengths  = lengths_of(ctx, k);
    int32_t     certain  = lengths.longer > 0.0 ? read_run->len : read_run->len - 1;

    for (int32_t j = 0; j < certain; j++) {
        add_certain_node(ctx->g, prefix_of(ctx, read_run->start + j), read_run->base,
                         GRAPH_NO_QUALITY);
    }

    if (lengths.longer > 0.0) {
        add_longer_run_end(ctx, read_run, lengths);
    } else {
        add_run_end(ctx, read_run, lengths);
    }
}

/* Lays out the bases the sequencer can have missed between run k and the next: one
 * node per base whose flow lies between the two runs, entered with the base at the
 * chance of a missed run and without it at the rest. The nodes stand at the read
 * prefix after the run. */
static void lay_out_missed_bases(const context *ctx, size_t k)
{
    const graph_run_errors *errors   = &ctx->errors->of_run[k];
    const run              *read_run = &ctx->runs[k];
    int32_t                 after    = prefix_of(ctx, read_run->start + read_run->len - 1);

    if (errors->missed <= 0.0) {
        return;
    }

    for (uint8_t m = 0; m < errors->n_missed; m++) {
        add_node(ctx->g, after, errors->missed_base[m], GRAPH_NO_QUALITY,
                 errors->missed, 1.0 - errors->missed, 0.0);
    }
}

static void lay_out_runs(const context *ctx)
{
    for (size_t k = 0; k < ctx->n_runs; k++) {
        lay_out_run(ctx, k);

        if (k + 1 < ctx->n_runs) {
            lay_out_missed_bases(ctx, k);
        }
    }
}

/* Returns the most nodes the read can need: the source, one node per placed base, and
 * per run one node for an extra base and one per base that can be missed after it. */
static size_t most_nodes(const context *ctx)
{
    size_t placed = (size_t)(ctx->span.end - ctx->span.begin);

    return 1 + placed + ctx->n_runs * (1 + GRAPH_MISSED_BASES);
}

/* ------------------------------------------------------------------------ */
/* Entry                                                                     */
/* ------------------------------------------------------------------------ */

int graph_lay_out(graph *g, const cm_bam_record *read, aln_span span, const run *runs,
                  size_t n_runs, const graph_read_errors *errors)
{
    context ctx = {
        .g      = g,
        .read   = read,
        .span   = span,
        .runs   = runs,
        .n_runs = n_runs,
        .errors = errors,
    };

    if (grow(g, most_nodes(&ctx)) < 0) {
        return -1;
    }

    g->nodes[0] = (graph_node){ 0 };
    g->n_nodes  = 1;

    if (errors) {
        lay_out_runs(&ctx);
    } else {
        lay_out_chain(&ctx);
    }

    return 0;
}

bool graph_node_passed_over(const graph *g, size_t n)
{
    return n + 1 < g->n_nodes && g->nodes[n + 1].skip[1] > 0.0;
}

void graph_free(graph *g)
{
    free(g->nodes);
    *g = (graph){ 0 };
}
