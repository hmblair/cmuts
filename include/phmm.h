/* phmm.h -- a banded pair HMM.
 *
 * The model aligns a read to its reference through pairings, deletions, and insertions,
 * at the rates of the reverse transcriptase, and reads the posterior of each event at
 * each reference base off the forward and backward passes.
 *
 * The read is scored one run of one base at a time. The bases of a run are paired one
 * at a time, as a pair HMM over bases does, within a band around the CIGAR path. A
 * sequencer that measures run lengths, as Ultima does, reports its errors per run, and
 * they enter where a run ends: the template run can be longer or shorter than the read
 * run by one or two bases, at the chances the read's tags give, and a run of one base
 * can be missing from the read between two runs. A read without the tags is scored
 * exactly as a pair HMM over bases would score it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "bam.h"
#include "fasta.h"
#include "params.h"
#include "phred.h"
#include "runs.h"

/* How far either side of the CIGAR path the model may look. */
#define PHMM_DEFAULT_BAND 2

/* Every parameter of the model for one reference.
 *
 * The transition weights hold along the whole reference: the extending of each run and
 * the closing that leaves what the extension does not take. Insertion to deletion and
 * deletion to insertion are absent: an inserted base followed by a deleted one is a
 * substitution.
 *
 * The arrays hold one value per reference base, indexed as the counted events are.
 * modification[b] is the chance base b was modified before it was read.
 * open_insertion[b] opens the insertion counted at base b, and open_deletion[b] the
 * deletion run whose 3'-most deleted base is b; the two are the arms of the decision
 * made at the pairing 3' of each run. ends[b] is proportional to the number of reads
 * whose 5'-most paired base is b, which is the termination rate at b times the coverage
 * at b; only the ratios between bases in one read's band matter. A cell outside the
 * reference reads the nearest value, which reaches only terms a zero emission clears. */
typedef struct {
    double        match_to_insertion;
    double        match_to_deletion;
    double        insertion_to_insertion;
    double        deletion_to_deletion;
    const double *modification;
    const double *open_insertion;
    const double *open_deletion;
    const double *ends;
} phmm_rates;

/* The arrays of phmm_rates. */
#define PHMM_RATE_ARRAYS 4

/* Sets the transition weights from the extension rates in uniform. The arrays are left
 * as they are. */
void phmm_rates_set_transitions(phmm_rates *rates, const phmm_uniform_rates *uniform);

/* One reference base as every pairing of it reads it: what the pairing pays after no
 * gap, after a deletion whose last base is the one before, and after an insertion; the
 * modification rate; the base as nuc names it; and whether the base is, by itself, a
 * run of one base. */
typedef struct {
    double  after_no_gap;
    double  after_deletion;
    double  after_insertion;
    double  modification;
    uint8_t code;
    bool    is_one_base_run;
} phmm_base;

/* One reference prepared for the marginal: its runs and each base as the pairings read
 * it. Built once per reference and read by every read against it, from any thread. */
typedef struct {
    const phmm_rates *rates;
    run              *runs;      /* owned */
    phmm_base        *bases;     /* owned; one per base */
    size_t            n_runs;
    size_t            len;
    size_t            cap_runs;
    size_t            cap_bases;
    int               band;
} phmm_model;

/* Prepares ref under rates, whose arrays must hold one value per base of ref and must
 * outlive every use of the model. The buffers are grown to fit and never shrink, so a
 * model can be prepared again for another reference. Returns -1 where they cannot be
 * grown. */
int phmm_model_prepare(phmm_model *model, const phmm_rates *rates,
                       const cm_fasta_record *ref, int band);

void phmm_model_free(phmm_model *model);

/* Buffers one thread reuses across reads, grown to fit the longest read seen. */
typedef struct phmm_scratch phmm_scratch;

phmm_scratch *phmm_scratch_create(void);
void          phmm_scratch_destroy(phmm_scratch *scratch);

/* The computed values for one read at one reference position. */
typedef struct {
    double coverage;     /* base read there */
    double mismatches;   /* template differences under the read bases */
    double insertions;   /* insertions opened after this base */
    double deletions;    /* deletions opened at this base */
    double ends;         /* the read's 5'-most pairing lands here */
} phmm_position;

/* Computed values for one read, in a window around the reference. */
typedef struct {
    hts_pos_t            origin;  /* reference position of value 0 */
    size_t               len;
    const phmm_position *at;      /* one per position of the window */
} phmm_window;

/* Gives the window indices that fall inside a reference of len bases: from the first to
 * one past the last, an empty range where none do. */
void phmm_window_bounds(const phmm_window *window, size_t len,
                        size_t *from, size_t *to);

/* How a marginalization ended. PHMM_NO_PATH skips the read, while the other
 * two failures end the run. */
typedef enum {
    PHMM_OK,          /* out holds the read's contribution */
    PHMM_NO_MEMORY,   /* the matrix could not be allocated */
    PHMM_NO_PATH,     /* no alignment the band admits has any probability */
    PHMM_UNSOUND,     /* the two passes did not agree */
} phmm_status;

/* How many bases a run length can be off by. */
#define PHMM_MAX_LENGTH_ERROR 2
#define PHMM_LENGTH_OFFSETS   (2 * PHMM_MAX_LENGTH_ERROR + 1)

/* What a sequencer that measures run lengths reports on one read run: the chance of
 * each template length, from PHMM_MAX_LENGTH_ERROR bases shorter than the read run to
 * as many longer, indexed by the offset plus PHMM_MAX_LENGTH_ERROR; the chance that a
 * run of one base is missing between this run and the next; and the bases that run
 * can be, as bits indexed by nuc. */
typedef struct {
    double  length[PHMM_LENGTH_OFFSETS];
    double  missed;
    uint8_t missed_bases;
} phmm_run_errors;

/* The run length errors of one read: one entry per run of its placed span, and the
 * substitution score every base takes, since the quality string of such a read
 * describes run lengths. */
typedef struct {
    int                    phred;
    const phmm_run_errors *of_run;
} phmm_length_errors;

/* Marginalizes one read over the alignments its band admits, leaving the result in
 * out, which borrows from scratch and is valid until the next call on it.
 *
 * runs holds the n_runs runs of the placed span of the read, as runs_of_read gives
 * them. errors holds the run length errors of the read, and is NULL for a read whose
 * sequencer reports none, which is then scored by the rates and its base qualities
 * alone. The model must be prepared for the read's reference. The read must store a
 * sequence and place at least one base, which the filter enforces. */
phmm_status phmm_run(const phmm_model *model, const phred *quality,
                     const cm_bam_record *read, const run *runs, size_t n_runs,
                     const phmm_length_errors *errors, phmm_scratch *scratch,
                     phmm_window *out);
