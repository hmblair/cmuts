/* runs.h -- the runs of one base in a sequence.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

#include "nuc.h"

/* A maximal stretch of one base. An unnamed base is a run of its own. */
typedef struct {
    nuc     base;
    int32_t start;  /* offset of the first base */
    int32_t len;
} run;

/* Writes the runs of bases [from, to) of a packed read sequence and returns their count.
 * out must hold to - from entries. */
size_t runs_of_read(const uint8_t *seq, int32_t from, int32_t to, run *out);

/* Writes the runs of a text sequence of len bases and returns their count. out must
 * hold len entries. */
size_t runs_of_text(const char *seq, size_t len, run *out);
