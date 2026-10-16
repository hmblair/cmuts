/* runs.c -- the runs of one base in a sequence.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "runs.h"

/* Appends a base to the last run if it is the same base, and opens a run
 * otherwise. Returns the new count. An unnamed base never joins a run. */
static size_t take_base(run *out, size_t count, nuc base, int32_t offset)
{
    if (count > 0 && nuc_is_base(base) && out[count - 1].base == base) {
        out[count - 1].len++;
        return count;
    }

    out[count] = (run){ .base = base, .start = offset, .len = 1 };

    return count + 1;
}

size_t runs_of_read(const uint8_t *seq, int32_t from, int32_t to, run *out)
{
    size_t count = 0;

    for (int32_t j = from; j < to; j++) {
        count = take_base(out, count, nuc_from_read(seq, j), j);
    }

    return count;
}

size_t runs_of_text(const char *seq, size_t len, run *out)
{
    size_t count = 0;

    for (size_t b = 0; b < len; b++) {
        count = take_base(out, count, nuc_from_char(seq[b]), (int32_t)b);
    }

    return count;
}
