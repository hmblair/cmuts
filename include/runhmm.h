/* runhmm.h -- a banded pair HMM over the runs of a read.
 *
 * The read is scored one run at a time. A cell pairs one read run with a prefix of the
 * reference, and the band keeps every cell within a fixed number of bases of the CIGAR
 * path, as the base model does. Inside a cell, the bases of the run are paired one at
 * a time as the base model does, so the rates and the quality of each base enter as
 * they do there. A sequencer that measures run lengths, as Ultima does, reports its
 * length errors per run, and those enter where a run ends.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "bam.h"
#include "fasta.h"
#include "flow.h"
#include "phmm.h"
#include "phred.h"
#include "runs.h"

/* How far either side of the CIGAR path the model may look, in reference bases. */
#define RUNHMM_DEFAULT_BAND 2

/* One reference base as every pairing of it reads it: what the pairing pays after no
 * gap, after a deletion whose last base is the one before, and after an insertion;
 * the modification rate; the base as nuc names it; and whether the base is, by
 * itself, a run of one base. */
typedef struct {
    double  keep;
    double  close_deletion;
    double  open_insertion;
    double  modification;
    uint8_t code;
    bool    is_one_base_run;
} runhmm_base;

/* One reference prepared for the marginal: its runs, and each base as the pairings
 * read it. Built once per reference and read by every read against it, from any
 * thread. */
typedef struct {
    const phmm_rates *rates;
    run              *runs;      /* owned */
    runhmm_base      *bases;     /* owned; one per base */
    size_t            n_runs;
    size_t            len;
    size_t            cap_runs;
    size_t            cap_bases;
    int               band;
} runhmm_model;

/* Prepares ref under rates, whose arrays must hold one value per base of ref and must
 * outlive every use of the model. Returns -1 where the arrays cannot be grown. */
int runhmm_model_prepare(runhmm_model *model, const phmm_rates *rates,
                         const cm_fasta_record *ref, int band);

void runhmm_model_free(runhmm_model *model);

/* Buffers one thread reuses across reads. */
typedef struct runhmm_scratch runhmm_scratch;

runhmm_scratch *runhmm_scratch_create(void);
void            runhmm_scratch_destroy(runhmm_scratch *scratch);

/* Marginalizes one read over the alignments its band admits, leaving the result in
 * out, which borrows from scratch and is valid until the next call on it. flow gives
 * the tables for an Ultima read; a read without the tags is scored by the rates and
 * its base qualities alone. The model must be prepared for the read's reference. The
 * read must store a sequence and place at least one base. */
phmm_status runhmm_run(const runhmm_model *model, const phred *quality,
                       const flow_tables *flow, const cm_bam_record *read,
                       runhmm_scratch *scratch, phmm_window *out);
