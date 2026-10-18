/* ultima.h -- the run-length errors of an Ultima read, per run.
 *
 * An Ultima sequencer measures the length of each run of one base. Its errors are
 * runs called one base too long or too short. The tp tag gives, per base, the
 * alternative length that the base's quality score is the chance of: one longer, one
 * shorter, or none. The chance of an alternative length is the sum over the bases of
 * the run that refer to it. The t0 tag gives, per base, the chance that a base unlike
 * either run was missed between the run and the next run.
 *
 * The decoding follows the class FlowBasedRead of GATK, which reads the same tags, and
 * the description of the base calling in Almogy et al., bioRxiv 2022.05.29.493900. Each
 * flow between two runs that measured no base can have held one base, and the t0 tag
 * gives the chance of each such miss, so the bases that can be missed between two runs
 * are the bases of the cycle between theirs, in the order of the cycle.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "bam.h"
#include "nuc.h"
#include "graph.h"
#include "phred.h"
#include "runs.h"

/* The least chance of the called run length. The tags of a run can add up past one,
 * and GATK's FlowBasedRead gives the called length at least this much, as its
 * MINIMAL_CALL_PROB, which has the same value. */
#define ULTIMA_MINIMAL_CALL_CHANCE 0.1

/* The greatest chance of an error at a run or at the gap after it, which leaves the
 * called length or the called gap the minimal chance. */
#define ULTIMA_MAXIMAL_ERROR_CHANCE (1.0 - ULTIMA_MINIMAL_CALL_CHANCE)

/* The offset that the t0 string adds to a Phred score to make it printable. */
#define ULTIMA_PHRED_OFFSET 33

/* The bases that can be missed between a run of one base and a run of the next, in the
 * order the reverse transcriptase wrote them. */
typedef struct {
    uint8_t n;
    nuc     base[GRAPH_MISSED_BASES];
} ultima_missed_bases;

/* The flow order of one read group: the cycle of four bases the sequencer repeats, and
 * the bases that can be missed between a run of one base and a run of the next, for a
 * forward read and for a reversed one. */
typedef struct {
    char                *id;
    nuc                  cycle[NUC_BASES];
    bool                 has_cycle;  /* false if the header declares an order of another form */
    ultima_missed_bases  missed[2][NUC_COUNT][NUC_COUNT];
} ultima_read_group;

/* What decoding the tags needs: the error chance of each quality score, the
 * substitution score of every base of an Ultima read, and the flow order of each read
 * group. Built once and shared by every worker. */
typedef struct {
    double           error[PHRED_MAX + 1];
    int              phred;
    ultima_read_group *groups;
    size_t           n_groups;
} ultima_decoder;

/* Builds the tables with no read groups. phred is the substitution score of every base
 * of an Ultima read, whose quality string describes run lengths only. */
void ultima_decoder_build(ultima_decoder *tables, int phred);

/* Adds the read groups of a header. Returns -1 if the table cannot be grown. */
int ultima_decoder_add_read_groups(ultima_decoder *tables, const cm_bam_reader *reader);

void ultima_decoder_free(ultima_decoder *tables);

/* The length errors of one read, decoded from its tags, one entry per run. */
typedef struct {
    int               phred;   /* the substitution score of every base */
    graph_run_errors *of_run;
    size_t            cap;     /* runs of_run is sized for */
} ultima_errors;

/* Returns whether a read carries the tp tag and the base qualities that the decoding
 * reads. */
bool ultima_read_tagged(const cm_bam_record *read);

/* Decodes the errors of the given runs of a tagged read. Returns -1 if the entries
 * cannot be grown. */
int ultima_errors_decode(ultima_errors *errors, const ultima_decoder *tables,
                       const cm_bam_record *read, const run *runs, size_t n_runs);

void ultima_errors_free(ultima_errors *errors);
