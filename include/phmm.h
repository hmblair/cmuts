/* phmm.h -- a banded pair HMM over the reference and the cDNA.
 *
 * The read x is produced from the reference y in two stages. The reverse transcriptase
 * walks along y and writes the cDNA z, pairing, deleting, and inserting bases at its
 * rates, and stopping at its termination rate; the path of the walk is the alignment a.
 * The sequencer then reads z and reports x, with errors of its own. The joint
 * probability is P(a | y) P(z | a, y) P(x | z). The first two factors are the model of
 * reverse transcription, a pair HMM over y and z, which this module holds. The third
 * is the model of the sequencer, which graph.h holds as the sequence graph of the
 * cDNAs the read admits.
 *
 * Both a and z are hidden. The posterior of every event of reverse transcription at
 * each reference base, given x and y, is read off a forward and a backward pass over
 * the product of the pair HMM and the graph, within a band around the CIGAR path.
 * Every base of z is a step of the reverse transcriptase and enters the counts,
 * whether or not the sequencer reported it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "bam.h"
#include "fasta.h"
#include "graph.h"
#include "params.h"
#include "phred.h"

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
 * deletion run whose 3'-most deleted base is b. termination[b] is the chance that the
 * reverse transcriptase stops once it has copied base b: the transition from the
 * pairing of b to the end state, which emits nothing and has no transitions out of
 * it. Each pairing makes one decision with four arms: it opens an insertion, opens a
 * deletion, stops, or pairs the next base with what the other three leave. A cell
 * outside the reference reads the nearest value, which reaches only terms a zero
 * emission clears. */
typedef struct {
    double        match_to_insertion;
    double        match_to_deletion;
    double        insertion_to_insertion;
    double        deletion_to_deletion;
    const double *modification;
    const double *open_insertion;
    const double *open_deletion;
    const double *termination;
} phmm_rates;

/* The arrays of phmm_rates. */
#define PHMM_RATE_ARRAYS 4

/* Sets the transition weights from the extension rates in uniform. The arrays are left
 * as they are. */
void phmm_rates_set_transitions(phmm_rates *rates, const phmm_uniform_rates *uniform);

/* One reference base as every pairing of it reads it: what the pairing pays after no
 * gap, after a deletion whose last base is the one before, and after an insertion; the
 * termination rate, which the first pairing of a read pays in place of the three; the
 * modification rate; and the base as nuc names it. */
typedef struct {
    double  after_no_gap;
    double  after_deletion;
    double  after_insertion;
    double  termination;
    double  modification;
    uint8_t code;
} phmm_base;

/* One reference prepared for the marginal: each base as the pairings read it. Built
 * once per reference and read by every read against it, from any thread. */
typedef struct {
    const phmm_rates *rates;
    phmm_base        *bases;     /* owned; one per base */
    size_t            len;
    size_t            cap_bases;
    int               band;
} phmm_model;

/* Prepares ref under rates, whose arrays must hold one value per base of ref and must
 * outlive every use of the model. The buffer is grown to fit and never shrinks, so a
 * model can be prepared again for another reference. Returns -1 where it cannot be
 * grown. */
int phmm_model_prepare(phmm_model *model, const phmm_rates *rates,
                       const cm_fasta_record *ref, int band);

void phmm_model_free(phmm_model *model);

/* The comparisons of a base of z with the reference at each quality score. Built once
 * from the quality table and read by every thread. */
typedef struct phmm_emissions phmm_emissions;

/* Returns NULL where the table cannot be allocated. */
phmm_emissions *phmm_emissions_create(const phred *quality);
void            phmm_emissions_destroy(phmm_emissions *emissions);

/* Buffers one thread reuses across reads, grown to fit the longest read seen. */
typedef struct phmm_scratch phmm_scratch;

phmm_scratch *phmm_scratch_create(void);
void          phmm_scratch_destroy(phmm_scratch *scratch);

/* The computed values for one read at one reference position. */
typedef struct {
    double coverage;     /* a base of z pairs with it */
    double mismatches;   /* the base of z differs from it, under those pairings */
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

/* Marginalizes one read over the cDNAs its graph admits and the alignments its band
 * admits, leaving the result in out, which borrows from scratch and is valid until the
 * next call on it.
 *
 * g is the graph of the read, laid out over its placed span. The band of each node is
 * taken from the CIGAR of the read. The model must be prepared for the read's
 * reference. The read must store a sequence, which the filter enforces. */
phmm_status phmm_run(const phmm_model *model, const phmm_emissions *emissions,
                     const cm_bam_record *read, const graph *g, phmm_scratch *scratch,
                     phmm_window *out);
