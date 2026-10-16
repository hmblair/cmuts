/* phred.c -- the table that gives each score as a probability.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "phred.h"

#include <math.h>

#include "nuc.h"

/* The error of a base that gives no information about the template. At this error a
 * comparison emits the same value whether the two bases agree or differ, and the chance
 * of a modification stays at the rate the model uses.
 *
 * The table holds no larger error. A larger error reverses the comparison, so a base
 * that agrees with the reference raises the chance of a modification. */
#define UNINFORMATIVE_ERROR ((double)(NUC_BASES - 1) / NUC_BASES)

double phred_error_of(unsigned quality)
{
    return pow(10.0, -(double)quality / 10.0);
}

/* A base takes the uninformative error when its score is below the minimum, or when
 * its score gives a larger error. */
void phred_build(phred *table, int minimum)
{
    for (unsigned quality = 0; quality <= PHRED_MAX; quality++) {
        double called = phred_error_of(quality);

        table->error[quality] = (int)quality < minimum || called > UNINFORMATIVE_ERROR
                              ? UNINFORMATIVE_ERROR
                              : called;
    }
}
