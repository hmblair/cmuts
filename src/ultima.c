/* ultima.c -- decoding the run-length errors of an Ultima read, per run.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "ultima.h"

#include <stdlib.h>
#include <string.h>

/* The tp values of a run one base longer and one base shorter. Any other value gives
 * no error. */
#define LENGTH_LONGER  1
#define LENGTH_SHORTER (-1)

/* ------------------------------------------------------------------------ */
/* The bases between two runs                                              */
/* ------------------------------------------------------------------------ */

static size_t cycle_index(const nuc *cycle, nuc base)
{
    size_t i = 0;

    while (cycle[i] != base) {
        i++;
    }

    return i;
}

/* Writes the bases of the cycle between from and to, in the order of the cycle. Their
 * flows between the two runs measured no base. */
static void bases_between(const nuc *cycle, nuc from, nuc to, ultima_missed_bases *out)
{
    out->n = 0;

    for (size_t i = (cycle_index(cycle, from) + 1) % NUC_BASES; cycle[i] != to;
         i = (i + 1) % NUC_BASES) {
        out->base[out->n++] = cycle[i];
    }
}

/* Writes the complements of the bases of sequenced in reverse order. */
static void reverse_complement(const ultima_missed_bases *sequenced,
                               ultima_missed_bases *out)
{
    out->n = sequenced->n;

    for (uint8_t i = 0; i < sequenced->n; i++) {
        out->base[i] = nuc_complement(sequenced->base[sequenced->n - 1 - i]);
    }
}

/* Writes the bases that can be missing between a run of from and a run of to, in the
 * orientation of the stored read. A reverse-strand read was sequenced as the
 * complement of to followed by the complement of from, so the bases between are the
 * complements of the flows between those, in reverse order. */
static void bases_missed_between(const nuc *cycle, nuc from, nuc to, bool reversed,
                                 ultima_missed_bases *out)
{
    ultima_missed_bases sequenced;

    if (!reversed) {
        bases_between(cycle, from, to, out);
        return;
    }

    bases_between(cycle, nuc_complement(to), nuc_complement(from), &sequenced);
    reverse_complement(&sequenced, out);
}

/* Writes the bases that can be missed between a run of each base and a run of each
 * other base, for each orientation. Two runs of one base are never adjacent, so that
 * entry holds no base. */
static void fill_bases_missed(ultima_read_group *group)
{
    for (int reversed = 0; reversed < 2; reversed++) {
        for (nuc from = NUC_A; from < NUC_COUNT; from++) {
            for (nuc to = NUC_A; to < NUC_COUNT; to++) {
                ultima_missed_bases *missed = &group->missed[reversed][from][to];

                missed->n = 0;

                if (from != to) {
                    bases_missed_between(group->cycle, from, to, reversed != 0, missed);
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Flow orders                                                               */
/* ------------------------------------------------------------------------ */

/* Reads the first four bases of a flow order as a cycle. Returns false unless they
 * are the four bases, each once. */
static bool cycle_of(const char *order, size_t len, nuc *cycle)
{
    bool seen[NUC_COUNT] = { false };

    if (len < NUC_BASES) {
        return false;
    }

    for (size_t i = 0; i < NUC_BASES; i++) {
        nuc base = nuc_from_char(order[i]);

        if (!nuc_is_base(base) || seen[base]) {
            return false;
        }

        seen[base] = true;
        cycle[i]   = base;
    }

    return true;
}

/* Returns whether the whole flow order repeats the cycle. */
static bool repeats_cycle(const char *order, size_t len, const nuc *cycle)
{
    if (len % NUC_BASES != 0) {
        return false;
    }

    for (size_t i = NUC_BASES; i < len; i++) {
        if (nuc_from_char(order[i]) != cycle[i % NUC_BASES]) {
            return false;
        }
    }

    return true;
}

static int read_group_from(ultima_read_group *group, const cm_bam_read_group *declared)
{
    group->id = strndup(declared->id, declared->id_len);

    if (!group->id) {
        return -1;
    }

    group->has_cycle = declared->flow_order
                    && cycle_of(declared->flow_order, declared->flow_order_len,
                                group->cycle)
                    && repeats_cycle(declared->flow_order, declared->flow_order_len,
                                     group->cycle);

    if (group->has_cycle) {
        fill_bases_missed(group);
    }

    return 0;
}

static int add_read_group(ultima_decoder *tables, const cm_bam_read_group *declared)
{
    ultima_read_group *groups = realloc(tables->groups,
                                      (tables->n_groups + 1) * sizeof *groups);

    if (!groups) {
        return -1;
    }

    tables->groups = groups;

    if (read_group_from(&groups[tables->n_groups], declared) < 0) {
        return -1;
    }

    tables->n_groups++;

    return 0;
}

int ultima_decoder_add_read_groups(ultima_decoder *tables, const cm_bam_reader *reader)
{
    cm_bam_rg_cursor  cursor;
    cm_bam_read_group declared;

    cm_bam_rg_open(&cursor, reader);

    while (cm_bam_rg_next(&cursor, &declared)) {
        if (add_read_group(tables, &declared) < 0) {
            return -1;
        }
    }

    return 0;
}

/* A read with no RG tag belongs to the only read group, if there is exactly one. */
static const ultima_read_group *group_of(const ultima_decoder *tables, const char *id)
{
    if (!id) {
        return tables->n_groups == 1 ? &tables->groups[0] : NULL;
    }

    for (size_t g = 0; g < tables->n_groups; g++) {
        if (strcmp(tables->groups[g].id, id) == 0) {
            return &tables->groups[g];
        }
    }

    return NULL;
}

void ultima_decoder_build(ultima_decoder *tables)
{
    for (unsigned quality = 0; quality <= PHRED_MAX; quality++) {
        tables->error[quality] = phred_error_of(quality);
    }

    tables->groups   = NULL;
    tables->n_groups = 0;
}

void ultima_decoder_free(ultima_decoder *tables)
{
    for (size_t g = 0; g < tables->n_groups; g++) {
        free(tables->groups[g].id);
    }

    free(tables->groups);
    tables->groups   = NULL;
    tables->n_groups = 0;
}

/* ------------------------------------------------------------------------ */
/* The errors                                                                */
/* ------------------------------------------------------------------------ */

static int grow_errors(ultima_errors *errors, size_t n_runs)
{
    graph_run_errors *of_run;

    if (n_runs <= errors->cap) {
        return 0;
    }

    of_run = realloc(errors->of_run, n_runs * sizeof *of_run);

    if (!of_run) {
        return -1;
    }

    errors->of_run = of_run;
    errors->cap    = n_runs;

    return 0;
}

void ultima_errors_free(ultima_errors *errors)
{
    free(errors->of_run);
    *errors = (ultima_errors){ 0 };
}

/* The chances that a run is one base longer and one base shorter than called. */
typedef struct {
    double longer;
    double shorter;
} one_base_errors;

/* Sums the chances of the bases of a run. A run of two or more bases carries half of
 * the chance of an alternative length on each of two bases placed symmetrically in
 * the run, so the sum over the bases is the whole chance. A run of one base carries
 * the whole chance on its one base. GATK's FlowBasedRead sums the bases the same way.
 * A run of unnamed bases has no chance. */
static one_base_errors sum_lengths(const ultima_decoder *tables,
                                   const cm_bam_record *read, const run *r)
{
    one_base_errors errors = { 0.0, 0.0 };

    if (!nuc_is_base(r->base)) {
        return errors;
    }

    /* The multiplication by a comparison avoids a branch. */
    for (int32_t j = r->start; j < r->start + r->len; j++) {
        double error = tables->error[read->qual[j]];

        errors.longer  += error * (read->tp[j] == LENGTH_LONGER);
        errors.shorter += error * (read->tp[j] == LENGTH_SHORTER);
    }

    return errors;
}

/* Scales the two chances down so that the called length keeps the minimal chance.
 * GATK's FlowBasedRead raises the called length to the minimum and leaves the two
 * chances as they are, so its column can sum past one. */
static one_base_errors keep_called_length(one_base_errors errors)
{
    double total = errors.longer + errors.shorter;
    double scale;

    if (total <= ULTIMA_MAXIMAL_ERROR_CHANCE) {
        return errors;
    }

    scale = ULTIMA_MAXIMAL_ERROR_CHANCE / total;

    return (one_base_errors){ errors.longer * scale, errors.shorter * scale };
}

/* Writes the chance of each template length from the one-base chances, with the rest
 * on the called length. */
static void write_lengths(graph_run_errors *out, one_base_errors errors)
{
    double *length = out->length;

    length[GRAPH_MAX_LENGTH_ERROR + 1] = errors.longer;
    length[GRAPH_MAX_LENGTH_ERROR - 1] = errors.shorter;
    length[GRAPH_MAX_LENGTH_ERROR]     = 1.0 - errors.longer - errors.shorter;
}

static double missed_chance_at(const ultima_decoder *tables, const cm_bam_record *read,
                               int32_t j)
{
    int score = (unsigned char)read->t0[j] - ULTIMA_PHRED_OFFSET;

    return tables->error[score < 0 ? 0 : score];
}

/* Returns the lesser of chance and the greatest error chance, so that the called gap
 * keeps the minimal chance. */
static double capped_gap_chance(double chance)
{
    return chance < ULTIMA_MAXIMAL_ERROR_CHANCE ? chance : ULTIMA_MAXIMAL_ERROR_CHANCE;
}

/* Writes the chance that a base unlike either run is missing between a run and the
 * next, and the bases that can be. The t0 tag gives a chance at the last base of one
 * run and at the first of the next. The lesser is the chance of the gap, as in GATK's
 * FlowBasedRead, and it holds for each flow between the two runs. The function caps
 * the chance so that the called gap keeps the minimal chance. */
static void decode_gap_after(graph_run_errors *out, const ultima_decoder *tables,
                             const cm_bam_record *read, const ultima_read_group *group,
                             const run *r, const run *next)
{
    const ultima_missed_bases *missed;
    double                     at_before, at_after;
    bool                       reversed = (read->flag & BAM_FREVERSE) != 0;

    out->missed   = 0.0;
    out->n_missed = 0;

    if (!next || !read->t0 || !group || !group->has_cycle) {
        return;
    }

    if (!nuc_is_base(r->base) || !nuc_is_base(next->base)) {
        return;
    }

    at_before = missed_chance_at(tables, read, r->start + r->len - 1);
    at_after  = missed_chance_at(tables, read, next->start);

    missed = &group->missed[reversed][r->base][next->base];

    out->missed   = capped_gap_chance(at_before < at_after ? at_before : at_after);
    out->n_missed = missed->n;

    for (uint8_t i = 0; i < missed->n; i++) {
        out->missed_base[i] = missed->base[i];
    }
}

bool ultima_read_tagged(const cm_bam_record *read)
{
    return read->tp && read->qual && read->seq;
}

int ultima_errors_decode(ultima_errors *errors, const ultima_decoder *tables,
                       const cm_bam_record *read, const run *runs, size_t n_runs)
{
    const ultima_read_group *group;

    if (grow_errors(errors, n_runs) < 0) {
        return -1;
    }

    group = group_of(tables, read->rg);

    for (size_t k = 0; k < n_runs; k++) {
        graph_run_errors *out  = &errors->of_run[k];
        const run       *next = k + 1 < n_runs ? &runs[k + 1] : NULL;

        write_lengths(out, keep_called_length(sum_lengths(tables, read, &runs[k])));
        decode_gap_after(out, tables, read, group, &runs[k], next);
    }

    return 0;
}
