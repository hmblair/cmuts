/* graph.h -- the sequence graph of the cDNA a read admits.
 *
 * The reverse transcriptase writes the cDNA z, and the sequencer reads z and reports
 * the read x with errors of its own, so z is hidden. The read determines a sequence
 * graph: a directed acyclic graph in which every node past the source holds one base
 * that z can contain, every edge has a weight, and every path from the source spells
 * one candidate z with the weight the sequencer model gives it.
 *
 * The nodes lie in the order the reverse transcriptase wrote them. A node is entered
 * from the node before it by its base edge, which carries the base of the node at the
 * weight of the edge, or from one of the nodes before it by a skip edge, which carries
 * no base: z lacks the base of the node and of every node the edge passes over. A
 * sequencer that reports every base, as Illumina does, gives a chain: one node per
 * placed base of x, with the quality score of the base, and no skip edge.
 *
 * A sequencer that measures the length of each run of one base, as Ultima does,
 * reports the chance that each run is one base longer or shorter than the read shows,
 * and the chance that a run of one base was missed between two runs. A run of such a
 * read gives one node per base the run can hold. The last of them is entered by its
 * base edge at the chance of the longest length, and by skip edges from the nodes
 * that end the shorter lengths at theirs. Each base whose flow lies between two runs
 * gives one node, entered by its base edge at the chance of a missed run and by a
 * skip edge at the rest. No node of such a read has a quality, since the quality
 * string of the read describes run lengths and every sequencing error of the read is a
 * length error.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "align.h"
#include "bam.h"
#include "nuc.h"
#include "runs.h"

/* How many bases a run length can be off by. */
#define GRAPH_MAX_LENGTH_ERROR 1
#define GRAPH_LENGTH_OFFSETS   (2 * GRAPH_MAX_LENGTH_ERROR + 1)

/* The most bases that can be missed between two runs: the named bases that are neither
 * run's base. */
#define GRAPH_MISSED_BASES (NUC_BASES - 2)

/* What a sequencer that measures run lengths reports on one run of a read: the chance
 * of each length of the run of z, from GRAPH_MAX_LENGTH_ERROR bases shorter than the
 * read run to as many longer, indexed by the offset plus GRAPH_MAX_LENGTH_ERROR; the
 * chance that a run of one base is missing between this run and the next, which holds
 * for each base that can be missing; and those bases, in the order the reverse
 * transcriptase wrote them. */
typedef struct {
    double  length[GRAPH_LENGTH_OFFSETS];
    double  missed;
    uint8_t n_missed;
    nuc     missed_base[GRAPH_MISSED_BASES];
} graph_run_errors;

/* The run-length errors of one read: one entry per run of its placed span. */
typedef struct {
    const graph_run_errors *of_run;
} graph_read_errors;

/* The quality of a base that is never misread: a base of a read that stores no quality
 * string, and every base of a read with run-length errors. */
#define GRAPH_NO_QUALITY (-1)

/* The skip edges into a node come from the node before it and from the node before
 * that: a run can end two nodes before the node that holds its longest length. */
#define GRAPH_SKIPS 2

/* One node: the weight of its base edge; the weight of the skip edge from each node
 * before it, zero where there is none; the read prefix it stands at, which places its
 * band on the reference; and its base, as nuc names it, with the quality the base was
 * read at. The source holds no base, and nothing reads its other fields. */
typedef struct {
    double  weight;
    double  skip[GRAPH_SKIPS];
    int32_t prefix;
    int16_t quality;
    uint8_t base;
} graph_node;

/* The graph of one read, grown to fit and reused across reads. Node 0 is the source. */
typedef struct {
    graph_node *nodes;
    size_t      n_nodes;
    size_t      cap;
} graph;

/* Lays out the graph of a read whose placed span is span. errors holds the run-length
 * errors of the read, and is NULL for a read whose sequencer reports none. With errors,
 * runs holds the n_runs runs of the placed span, as runs_of_read gives them; without,
 * runs is not read. Returns -1 where the nodes cannot be grown. */
int graph_lay_out(graph *g, const cm_bam_record *read, aln_span span, const run *runs,
                  size_t n_runs, const graph_read_errors *errors);

/* Returns whether a skip edge passes over node n without landing on it: the edge into
 * node n + 1 from node n - 1. */
static inline bool graph_node_passed_over(const graph *g, size_t n)
{
    return n + 1 < g->n_nodes && g->nodes[n + 1].skip[1] > 0.0;
}

/* Returns whether z holds the base of node n for certain: the node has no skip edge
 * and a base edge of weight one. */
static inline bool graph_node_certain(const graph *g, size_t n)
{
    const graph_node *node = &g->nodes[n];

    /* The certain nodes of a read with run-length errors follow no regular pattern,
     * so a branch on them is often mispredicted. The bitwise operators avoid one. */
    return (node->weight == 1.0) & (node->skip[0] == 0.0) & (node->skip[1] == 0.0);
}

void graph_free(graph *g);
