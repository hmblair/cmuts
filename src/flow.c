/* flow.c -- decoding the run length errors of an Ultima read, per run.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "flow.h"

#include <stdlib.h>
#include <string.h>

/* The tp values of a run one base longer and one base shorter. Any other value gives
 * no error. */
#define LENGTH_LONGER  1
#define LENGTH_SHORTER (-1)

/* The offset that the t0 string adds to a PHRED score to make it printable. */
#define PHRED_OFFSET 33

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

static int read_group_from(flow_read_group *group, const cm_bam_read_group *declared)
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

    return 0;
}

static int add_read_group(flow_tables *tables, const cm_bam_read_group *declared)
{
    flow_read_group *groups = realloc(tables->groups,
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

int flow_tables_add_read_groups(flow_tables *tables, const cm_bam_reader *reader)
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
static const flow_read_group *group_of(const flow_tables *tables, const char *id)
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

void flow_tables_build(flow_tables *tables, int phred)
{
    for (unsigned quality = 0; quality <= PHRED_MAX; quality++) {
        tables->error[quality] = phred_error_of(quality);
    }

    tables->phred    = phred;
    tables->groups   = NULL;
    tables->n_groups = 0;
}

void flow_tables_free(flow_tables *tables)
{
    for (size_t g = 0; g < tables->n_groups; g++) {
        free(tables->groups[g].id);
    }

    free(tables->groups);
    tables->groups   = NULL;
    tables->n_groups = 0;
}

/* ------------------------------------------------------------------------ */
/* The letters between two runs                                              */
/* ------------------------------------------------------------------------ */

static size_t cycle_index(const nuc *cycle, nuc base)
{
    size_t i = 0;

    while (cycle[i] != base) {
        i++;
    }

    return i;
}

/* Returns the letters of the cycle between from and to, as bits indexed by nuc. Their
 * flows between the two runs measured zero bases. */
static uint8_t letters_between(const nuc *cycle, nuc from, nuc to)
{
    uint8_t letters = 0;

    for (size_t i = (cycle_index(cycle, from) + 1) % NUC_BASES; cycle[i] != to;
         i = (i + 1) % NUC_BASES) {
        letters |= (uint8_t)(1U << cycle[i]);
    }

    return letters;
}

static uint8_t complement_letters(uint8_t letters)
{
    uint8_t complemented = 0;

    for (nuc base = NUC_A; base < NUC_COUNT; base++) {
        if (letters & (1U << base)) {
            complemented |= (uint8_t)(1U << nuc_complement(base));
        }
    }

    return complemented;
}

/* Returns the letters that can be missing between a run of from and a run of to, in
 * the orientation of the stored read. A reverse-strand read was sequenced as the
 * complement of to followed by the complement of from. */
static uint8_t letters_missed_between(const nuc *cycle, nuc from, nuc to,
                                      bool reversed)
{
    if (!reversed) {
        return letters_between(cycle, from, to);
    }

    return complement_letters(letters_between(cycle, nuc_complement(to),
                                              nuc_complement(from)));
}

/* ------------------------------------------------------------------------ */
/* The errors                                                                */
/* ------------------------------------------------------------------------ */

static int grow_errors(flow_errors *errors, size_t n_runs)
{
    flow_run_errors *of_run;

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

void flow_errors_free(flow_errors *errors)
{
    free(errors->of_run);
    *errors = (flow_errors){ 0 };
}

/* Sums the chances of the bases of a run. A run of unnamed bases has none. */
static void sum_lengths(flow_run_errors *out, const flow_tables *tables,
                        const cm_bam_record *read, const run *r)
{
    out->longer  = 0.0;
    out->shorter = 0.0;

    if (!nuc_is_base(r->letter)) {
        return;
    }

    for (int32_t j = r->start; j < r->start + r->len; j++) {
        double error = tables->error[read->qual[j]];

        if (read->tp[j] == LENGTH_LONGER) {
            out->longer += error;
        } else if (read->tp[j] == LENGTH_SHORTER) {
            out->shorter += error;
        }
    }
}

/* Scales the two chances down so that the called length keeps the minimal chance. */
static void keep_called_length(flow_run_errors *out)
{
    double total = (out->longer + out->shorter) * (1.0 + FLOW_TWO_BASE_FRACTION);
    double scale;

    if (total <= 1.0 - FLOW_MINIMAL_CALL_CHANCE) {
        return;
    }

    scale         = (1.0 - FLOW_MINIMAL_CALL_CHANCE) / total;
    out->longer  *= scale;
    out->shorter *= scale;
}

static double missed_chance_at(const flow_tables *tables, const cm_bam_record *read,
                               int32_t j)
{
    int score = (unsigned char)read->t0[j] - PHRED_OFFSET;

    return tables->error[score < 0 ? 0 : score];
}

/* Writes the chance that a base of another letter is missing between a run and the
 * next. The t0 tag gives a chance at the last base of one run and at the first of the
 * next. The lesser is the chance of the gap, as in the reference decoder. */
static void decode_gap_after(flow_run_errors *out, const flow_tables *tables,
                             const cm_bam_record *read, const flow_read_group *group,
                             const run *r, const run *next)
{
    double at_before, at_after;
    bool   reversed = (read->flag & BAM_FREVERSE) != 0;

    out->missed         = 0.0;
    out->missed_letters = 0;

    if (!next || !read->t0 || !group || !group->has_cycle) {
        return;
    }

    if (!nuc_is_base(r->letter) || !nuc_is_base(next->letter)) {
        return;
    }

    at_before = missed_chance_at(tables, read, r->start + r->len - 1);
    at_after  = missed_chance_at(tables, read, next->start);

    out->missed         = at_before < at_after ? at_before : at_after;
    out->missed_letters = letters_missed_between(group->cycle, r->letter,
                                                 next->letter, reversed);
}

int flow_errors_decode(flow_errors *errors, const flow_tables *tables,
                       const cm_bam_record *read, const run *runs, size_t n_runs)
{
    const flow_read_group *group;

    errors->tagged = read->tp && read->qual && read->seq;

    if (!errors->tagged) {
        return 0;
    }

    if (grow_errors(errors, n_runs) < 0) {
        return -1;
    }

    errors->phred = tables->phred;
    group         = group_of(tables, read->rg);

    for (size_t k = 0; k < n_runs; k++) {
        flow_run_errors *out  = &errors->of_run[k];
        const run       *next = k + 1 < n_runs ? &runs[k + 1] : NULL;

        sum_lengths(out, tables, read, &runs[k]);
        keep_called_length(out);
        decode_gap_after(out, tables, read, group, &runs[k], next);
    }

    return 0;
}
