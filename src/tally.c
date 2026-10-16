/* tally.c -- one read's contribution to a reference.
 *
 * Runs the marginal over one read and adds the window it returns to the reference's
 * accumulator, clipping to the reference's bounds. The counted quantities are
 * described in phmm.h.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "tally.h"

#include <stdlib.h>

#include "align.h"
#include "filter.h"
#include "runs.h"

/* A value of ends equal at every base, which leaves every placement of a read's 5'-most
 * paired base within its band equally likely. */
#define UNIFORM_ENDS 1.0

/* The substitution score of every base of an Ultima read. */
#define DEFAULT_ULTIMA_PHRED 40

typedef struct {
    const cm_bam_record   *read;
    const cm_fasta_record *ref;
    const phmm_model      *model;
    const tally_tables    *tables;
    accum                 *target;
    pairs                 *target_pairs;   /* NULL where no pairs are counted */
} context;

struct tally_scratch {
    phmm_scratch  *model;
    run           *runs;      /* of the placed span of the read being marginalized */
    size_t         cap_runs;
    ultima_errors  errors;    /* the run-length errors of that read */
};

/* Counts the read in its length bin. Binned by stored length, as the length filters are
 * applied, so inserted and soft-clipped bases count. A read longer than the range the
 * bins cover falls in none of them; the reads total gives how many those were.
 *
 * A supplementary alignment stores the piece it places and not the read, so binning one
 * would record a length no molecule had.
 *
 * The bins begin at length 1, since a read storing no sequence has been refused
 * already, so the guard against zero is against that filter changing and not against
 * anything reachable from here. */
static void add_length(const context *ctx)
{
    double *bins;
    size_t  length;

    if (filter_is_supplementary(ctx->read)) {
        return;
    }

    bins   = accum_data(ctx->target, ACCUM_LENGTHS);
    length = (size_t)ctx->read->l_qseq;

    if (length > 0 && length <= SHAPE_LENGTH_BINS(ctx->target->cap)) {
        bins[length - 1] += 1.0;
    }
}

/* The histogram has a bin for every quality the filter lets through, so the guard below
 * is against that filter changing and not against anything reachable from here. */
_Static_assert(SHAPE_QUALITY_BINS == FILTER_MAPQ_UNAVAILABLE,
               "the quality histogram must have a bin for every quality the filter passes");

/* Counts the read's mapping quality in its bin. Only primary reads are counted, so the
 * histogram sums to the primary reads counted, as the length histogram does. */
static void add_quality(const context *ctx)
{
    double *bins;

    if (filter_is_supplementary(ctx->read)) {
        return;
    }

    bins = accum_data(ctx->target, ACCUM_QUALITIES);

    if (ctx->read->mapq < SHAPE_QUALITY_BINS) {
        bins[ctx->read->mapq] += 1.0;
    }
}

/* Adds the window's stretch on the reference to the per-base fields, clipping once
 * rather than testing each position. The window holds one read, so each event value is
 * the posterior chance this read carries an event of that kind at the position. */
static void add_window(const context *ctx, const phmm_window *window)
{
    double *coverage   = accum_data(ctx->target, ACCUM_COVERAGE);
    double *mismatches = accum_data(ctx->target, ACCUM_MISMATCHES);
    double *insertions = accum_data(ctx->target, ACCUM_INSERTIONS);
    double *deletions  = accum_data(ctx->target, ACCUM_DELETIONS);
    double *ends       = accum_data(ctx->target, ACCUM_ENDS);
    size_t  begin;
    size_t  end;

    phmm_window_bounds(window, ctx->ref->len, &begin, &end);

    for (size_t i = begin; i < end; i++) {
        size_t               pos = (size_t)(window->origin + (hts_pos_t)i);
        const phmm_position *at  = &window->at[i];

        coverage[pos]   += at->coverage;
        mismatches[pos] += at->mismatches;
        insertions[pos] += at->insertions;
        deletions[pos]  += at->deletions;
        ends[pos]       += at->ends;
    }
}

/* Splits the placed span of the read into runs and returns their count, or -1 where
 * the runs cannot be grown. */
static ssize_t split_into_runs(tally_scratch *scratch, const cm_bam_record *read)
{
    aln_span span = aln_placed_span(read);
    run     *runs;

    if ((size_t)read->l_qseq > scratch->cap_runs) {
        runs = realloc(scratch->runs, (size_t)read->l_qseq * sizeof *runs);

        if (!runs) {
            return -1;
        }

        scratch->runs     = runs;
        scratch->cap_runs = (size_t)read->l_qseq;
    }

    return (ssize_t)runs_of_read(read->seq, span.begin, span.end, scratch->runs);
}

/* Splits a tagged read into runs and decodes its run-length errors. Returns the count
 * of runs, or -1 where a buffer cannot be grown. */
static ssize_t decode_errors(const context *ctx, tally_scratch *scratch)
{
    ssize_t n_runs = split_into_runs(scratch, ctx->read);

    if (n_runs < 0
        || ultima_errors_decode(&scratch->errors, &ctx->tables->ultima, ctx->read,
                                scratch->runs, (size_t)n_runs) < 0) {
        return -1;
    }

    return n_runs;
}

static phmm_status marginalize(const context *ctx, tally_scratch *scratch)
{
    phmm_window        window;
    phmm_status        status;
    phmm_length_errors errors;
    bool               tagged = ultima_read_tagged(ctx->read);
    ssize_t            n_runs = tagged ? decode_errors(ctx, scratch) : 0;

    if (n_runs < 0) {
        return PHMM_NO_MEMORY;
    }

    errors = (phmm_length_errors){
        .phred  = scratch->errors.phred,
        .of_run = scratch->errors.of_run,
    };
    status = phmm_run(ctx->model, ctx->tables->emissions, ctx->read, scratch->runs,
                      (size_t)n_runs, tagged ? &errors : NULL, scratch->model, &window);

    if (status == PHMM_OK) {
        add_window(ctx, &window);

        if (ctx->target_pairs) {
            pairs_count(ctx->target_pairs, ctx->ref->len, &window);
        }
    }

    return status;
}

/* ------------------------------------------------------------------------ */
/* Setup                                                                     */
/* ------------------------------------------------------------------------ */

tally_config tally_defaults(void)
{
    return (tally_config){
        .band         = PHMM_DEFAULT_BAND,
        .min_phred    = 0,
        .ultima_phred = DEFAULT_ULTIMA_PHRED,
        .uniform      = phmm_uniform_defaults(),
    };
}

/* Fills one array of rates with a single rate. */
static void fill_uniform(double *values, size_t len, double rate)
{
    for (size_t i = 0; i < len; i++) {
        values[i] = rate;
    }
}

static int add_read_groups(ultima_decoder *ultima, const cm_bam_stream *stream)
{
    for (size_t file = 0; file < cm_bam_stream_count(stream); file++) {
        if (ultima_decoder_add_read_groups(ultima,
                                           cm_bam_stream_reader(stream, file)) < 0) {
            return -1;
        }
    }

    return 0;
}

int tally_tables_build(tally_tables *tables, const tally_config *config,
                       const cm_bam_stream *stream)
{
    phred_build(&tables->quality, config->min_phred);
    ultima_decoder_build(&tables->ultima, config->ultima_phred);
    tables->uniform   = config->uniform;
    tables->band      = config->band;
    tables->emissions = phmm_emissions_create(&tables->quality);

    if (!tables->emissions) {
        return -1;
    }

    return add_read_groups(&tables->ultima, stream);
}

void tally_tables_free(tally_tables *tables)
{
    ultima_decoder_free(&tables->ultima);
    phmm_emissions_destroy(tables->emissions);
    tables->emissions = NULL;
}

/* Points each array of rates at its own len values of storage. */
static void bind_arrays(phmm_rates *rates, const double *storage, size_t len)
{
    rates->modification   = storage;
    rates->open_insertion = storage + len;
    rates->open_deletion  = storage + 2 * len;
    rates->ends           = storage + 3 * len;
}

/* Writes the uniform rates into storage, in the order bind_arrays lays the arrays out. */
static void fill_arrays(double *storage, size_t len, const phmm_uniform_rates *uniform)
{
    fill_uniform(storage, len, uniform->modification);
    fill_uniform(storage + len, len, uniform->open_insertion);
    fill_uniform(storage + 2 * len, len, uniform->open_deletion);
    fill_uniform(storage + 3 * len, len, UNIFORM_ENDS);
}

void tally_rates_fill(phmm_rates *rates, double *storage, size_t len,
                      const tally_tables *tables)
{
    phmm_rates_set_transitions(rates, &tables->uniform);
    bind_arrays(rates, storage, len);
    fill_arrays(storage, len, &tables->uniform);
}

tally_scratch *tally_scratch_create(void)
{
    tally_scratch *scratch = calloc(1, sizeof *scratch);

    if (!scratch) {
        return NULL;
    }

    scratch->model = phmm_scratch_create();

    if (!scratch->model) {
        free(scratch);
        return NULL;
    }

    return scratch;
}

void tally_scratch_destroy(tally_scratch *scratch)
{
    if (!scratch) {
        return;
    }

    phmm_scratch_destroy(scratch->model);
    ultima_errors_free(&scratch->errors);
    free(scratch->runs);
    free(scratch);
}

/* Returns the accumulator that counts this record once it is tallied. A supplementary
 * record places a piece of a read, and is counted apart from the read itself. */
static accum_field_id counted_field(const cm_bam_record *read)
{
    return filter_is_supplementary(read) ? ACCUM_SUPPLEMENTARY_COUNTED
                                         : ACCUM_PRIMARY_COUNTED;
}

/* Returns the accumulator that counts this record once it is turned away. */
static accum_field_id rejected_field(const cm_bam_record *read)
{
    return filter_is_supplementary(read) ? ACCUM_SUPPLEMENTARY_REJECTED
                                         : ACCUM_PRIMARY_REJECTED;
}

phmm_status tally(const cm_bam_record *read, const cm_fasta_record *ref,
                  const phmm_model *model, const tally_tables *tables,
                  tally_scratch *scratch, accum *target, pairs *target_pairs)
{
    context ctx = {
        .read         = read,
        .ref          = ref,
        .model        = model,
        .tables       = tables,
        .target       = target,
        .target_pairs = target_pairs,
    };
    phmm_status status = marginalize(&ctx, scratch);

    /* A record the model gives no path is counted where a record a filter turned away is
     * counted, and adds no other value. No value has reached the target: the window is added
     * only on PHMM_OK, and the failure is seen before either count below. */
    if (status == PHMM_NO_PATH) {
        *accum_data(target, rejected_field(read)) += 1.0;
        return status;
    }

    if (status != PHMM_OK) {
        return status;
    }

    *accum_data(target, counted_field(read)) += 1.0;

    add_length(&ctx);
    add_quality(&ctx);

    return PHMM_OK;
}
