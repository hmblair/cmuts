/* flow.h -- the run length errors of an Ultima read, per run.
 *
 * An Ultima sequencer measures the length of each run of one base. Its errors are
 * runs called one base too long or too short. The tp tag gives, per base, the
 * alternative length that the base's quality score is the chance of: one longer, one
 * shorter, or none. The chance of an alternative length is the sum over the bases of
 * the run that refer to it. The t0 tag gives, per base, the chance that a base of
 * base unlike either run was missed between the run and the next run.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "bam.h"
#include "nuc.h"
#include "phred.h"
#include "runs.h"

/* The least chance of the called run length, as in the reference decoder. */
#define FLOW_MINIMAL_CALL_CHANCE 0.1

/* The chance of a length error of two bases, as a fraction of the chance of the
 * one-base error in the same direction. Measured on the 2A3 sample. */
#define FLOW_TWO_BASE_FRACTION 0.014

/* The flow order of one read group: the cycle of four bases the sequencer repeats, and
 * the bases that can be missed between a run of one base and a run of the next,
 * as bits indexed by nuc, for a forward read and for a reversed one. */
typedef struct {
    char    *id;
    nuc      cycle[NUC_BASES];
    bool     has_cycle;  /* false if the header declares an order of another form */
    uint8_t  bases_missed[2][NUC_COUNT][NUC_COUNT];
} flow_read_group;

/* The error chance of each quality score, the substitution score of every base of an
 * Ultima read, and the flow order of each read group. Built once and shared by every
 * worker. */
typedef struct {
    double           error[PHRED_MAX + 1];
    int              phred;
    flow_read_group *groups;
    size_t           n_groups;
} flow_tables;

/* Builds the tables with no read groups. phred is the substitution score of every base
 * of an Ultima read, whose quality string describes run lengths only. */
void flow_tables_build(flow_tables *tables, int phred);

/* Adds the read groups of a header. Returns -1 if the table cannot be grown. */
int flow_tables_add_read_groups(flow_tables *tables, const cm_bam_reader *reader);

void flow_tables_free(flow_tables *tables);

/* The length errors of one read run. longer is the chance that the template run is
 * one base longer than the read run, and shorter that it is one base shorter. missed
 * is the chance that a base unlike either run is missing between this run and the
 * next, and missed_bases holds the bases it can be, as bits indexed by nuc. */
typedef struct {
    double  longer;
    double  shorter;
    double  missed;
    uint8_t missed_bases;
} flow_run_errors;

/* The length errors of one read, one entry per run. */
typedef struct {
    bool             tagged;  /* whether the read carries the tags; the entries are
                                 valid only if so */
    int              phred;   /* the substitution score of every base */
    flow_run_errors *of_run;
    size_t           cap;     /* runs of_run is sized for */
} flow_errors;

/* Decodes the errors of the given runs of a read. Sets tagged to false if the read has
 * no tp tag or no base qualities. Returns -1 if the entries cannot be grown. */
int flow_errors_decode(flow_errors *errors, const flow_tables *tables,
                       const cm_bam_record *read, const run *runs, size_t n_runs);

void flow_errors_free(flow_errors *errors);
