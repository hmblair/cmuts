/* bam.h -- sequential iteration over the alignments of a BAM file.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <htslib/sam.h>

#include "iter.h"

/* A single alignment record.
 *
 * Every pointer borrows memory owned by the reader that produced the record, and is
 * invalidated by the next cm_bam_next() or cm_bam_close() on that reader. Copy anything
 * that must outlive the current iteration step.
 *
 * No filtering is applied, so unmapped reads reach the caller, carrying
 * tid == CM_BAM_NO_REF, pos == -1 and n_cigar == 0. Test flag against BAM_FUNMAP and not
 * tid: a paired read may be unmapped and still carry its mate's reference. */
typedef struct {
    int32_t         tid;      /* reference index into the file header */
    hts_pos_t       pos;      /* 0-based leftmost coordinate on that reference */
    uint16_t        flag;     /* SAM flags, for BAM_FUNMAP and the rest */
    uint8_t         mapq;     /* mapping quality */
    const uint32_t *cigar;    /* raw ops; bam_cigar_op / bam_cigar_oplen decode */
    uint32_t        n_cigar;  /* number of ops in cigar */
    const uint8_t  *seq;      /* bases packed two to a byte; nuc_from_read reads
                                 one. NULL where the record stores no sequence */
    const uint8_t  *qual;     /* raw PHRED scores, no +33 offset; NULL if absent */
    int32_t         l_qseq;   /* number of bases, and length of seq and qual */
    /* The run length error tags of an Ultima read, one value per base; see flow.h.
     * NULL if the record has no such tag, or the tag has another length. */
    const int8_t   *tp;       /* the alternative run length of each base's quality,
                                 relative to the called length */
    const char     *t0;       /* PHRED scores with the +33 offset, of a base unlike
                                 either run missing between two runs */
    const char     *rg;       /* the read group id, or NULL if the record has none */
} cm_bam_record;

/* Reference index of an unmapped read, and of any tid with no header entry. */
#define CM_BAM_NO_REF (-1)

typedef struct cm_bam_reader cm_bam_reader;

/* Opens path and reads its header. No .bai index is required or consulted, iteration being
 * a linear pass in file order. Returns NULL on failure, leaving the reason in why. */
cm_bam_reader *cm_bam_open(const char *path, const char **why);

/* Points CRAM decoding at the given reference. Without this, htslib locates one on its
 * own, from the header's UR path or M5 checksum, which need not be the reference the
 * caller works against. Harmless for SAM and BAM, which carry their own sequence. */
int cm_bam_set_reference(cm_bam_reader *reader, const char *fasta_path);

/* Fills out with the next alignment. Returns a cm_iter_status; on
 * CM_ITER_ERROR the cause is available from cm_bam_error(). */
int cm_bam_next(cm_bam_reader *reader, cm_bam_record *out);

void cm_bam_close(cm_bam_reader *reader);

/* Returns a description of the reader's failure, or NULL if it has not failed. */
const char *cm_bam_error(const cm_bam_reader *reader);

/* Give how far the reader has read, and how far it has to go, in compressed bytes of the
 * alignments only, excluding the header. Span is 0 where the size is not known. */
uint64_t cm_bam_position(const cm_bam_reader *reader);
uint64_t cm_bam_span(const cm_bam_reader *reader);

/* Draws the reader's decompression threads from a pool, parallelizing inflation only:
 * reading and record parsing stay sequential. Several readers may draw from one pool,
 * which must outlive every reader given it. */
int cm_bam_use_pool(cm_bam_reader *reader, htsThreadPool *pool);

/* ------------------------------------------------------------------------ */
/* The record underlying a view                                              */
/* ------------------------------------------------------------------------ */

/* Gives the htslib record behind the most recent cm_bam_next(), for callers that must
 * retain a read past the next advance. Copy it with bam_copy1. */
const bam1_t *cm_bam_raw(const cm_bam_reader *reader);

/* Presents an htslib record as a cm_bam_record, with the tags left NULL. The view
 * borrows from record, so it stays valid exactly as long as record does. */
void cm_bam_record_view(const bam1_t *record, cm_bam_record *out);

/* Fills the tags of a view, which searches the record's auxiliary data. */
void cm_bam_record_view_tags(const bam1_t *record, cm_bam_record *out);

/* ------------------------------------------------------------------------ */
/* Header queries                                                            */
/* ------------------------------------------------------------------------ */

/* Gives the name of the reference with the given index, or NULL if there is no such
 * reference. The returned string is owned by the reader's header. */
const char *cm_bam_refname(const cm_bam_reader *reader, int32_t tid);

/* Gives the length of the reference with the given index, or -1 if there is no such
 * reference. */
hts_pos_t cm_bam_reflen(const cm_bam_reader *reader, int32_t tid);

/* Gives the length of the longest reference in the header. */
hts_pos_t cm_bam_max_reflen(const cm_bam_reader *reader);

/* Gives the number of references declared in the header. */
int32_t cm_bam_nref(const cm_bam_reader *reader);

/* Returns whether the header declares SO:coordinate. */
bool cm_bam_is_coordinate_sorted(const cm_bam_reader *reader);

/* ------------------------------------------------------------------------ */
/* Reference declarations                                                    */
/* ------------------------------------------------------------------------ */

/* A forward-only walk over the header's @SQ lines. The cursor holds the state of the walk
 * and belongs to the caller. */
typedef struct {
    const char *line;  /* the @SQ line describing tid, or NULL past the last */
    int32_t     tid;
} cm_bam_sq_cursor;

void cm_bam_sq_open(cm_bam_sq_cursor *cursor, const cm_bam_reader *reader);

/* Gives the M5 checksum declared for a reference, its length in len, or NULL where the
 * header declares none. The result points into the header text and lives as long as the
 * reader does. References must be requested in non-decreasing order. */
const char *cm_bam_sq_checksum(cm_bam_sq_cursor *cursor, int32_t tid, size_t *len);

/* ------------------------------------------------------------------------ */
/* Read groups                                                               */
/* ------------------------------------------------------------------------ */

/* One @RG line of the header. Both strings point into the header text, have no
 * terminator, and live as long as the reader does. */
typedef struct {
    const char *id;
    size_t      id_len;
    const char *flow_order;      /* the FO field, or NULL where the line has none */
    size_t      flow_order_len;
} cm_bam_read_group;

/* A forward-only walk over the header's @RG lines. The cursor holds the state of the walk
 * and belongs to the caller. */
typedef struct {
    const char *line;  /* the next @RG line to give, or NULL past the last */
} cm_bam_rg_cursor;

void cm_bam_rg_open(cm_bam_rg_cursor *cursor, const cm_bam_reader *reader);

/* Returns false after the last read group. A line without an ID is passed over. */
bool cm_bam_rg_next(cm_bam_rg_cursor *cursor, cm_bam_read_group *out);
