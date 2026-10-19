/* rates.c -- the events of each kind at a position over that channel's denominator.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "rates.h"

#include <math.h>

rate_config rate_defaults(void)
{
    return (rate_config){ .min_depth = 1 };
}

/* The accumulated events each channel counts. */
static const accum_field_id EVENTS_OF[FMT_N_CHANNELS] = {
    [FMT_CHANNEL_MISMATCHES]   = ACCUM_MISMATCHES,
    [FMT_CHANNEL_INSERTIONS]   = ACCUM_INSERTIONS,
    [FMT_CHANNEL_DELETIONS]    = ACCUM_DELETIONS,
    [FMT_CHANNEL_TERMINATIONS] = ACCUM_ENDS,
};

/* How many bases 3' of the position each channel's trials sit. A mismatch is tried on
 * every pairing of the position. An insertion, a deletion, and a termination are arms
 * of the one decision that each pairing makes, so their trials are the pairings that
 * make it. That pairing sits at the position itself for an insertion, and one base 3'
 * of it for a deletion and for a termination, which are both counted at the first base
 * the enzyme did not copy. */
static const size_t TRIALS_SHIFT[FMT_N_CHANNELS] = {
    [FMT_CHANNEL_MISMATCHES]   = 0,
    [FMT_CHANNEL_INSERTIONS]   = 0,
    [FMT_CHANNEL_DELETIONS]    = 1,
    [FMT_CHANNEL_TERMINATIONS] = 1,
};

/* Returns a channel's trials at one position, or zero where the pairing that tries it
 * falls past the reference. */
static double depth_at(fmt_channel channel, const double *coverage, size_t i, size_t len)
{
    size_t at = i + TRIALS_SHIFT[channel];

    return at < len ? coverage[at] : 0.0;
}

/* Returns whether a position carries enough depth to report on. Some is always required:
 * a depth of zero does not report on a position with none. */
static bool meets_min_depth(double wanted, double depth)
{
    return depth > 0.0 && depth >= wanted;
}

/* Returns whether a position falls within a masked end of the reference. */
static bool masked_at(const rate_config *cfg, size_t i, size_t len)
{
    return i < cfg->nan_5p || len - i <= cfg->nan_3p;
}

/* Returns the events at a position over its depth, held to one. Every event is a share
 * of the depth it is tried against, so the ratio cannot exceed one except by rounding. */
static double ratio_of(double events, double depth)
{
    double rate = depth > 0.0 ? events / depth : 0.0;

    return rate > 1.0 ? 1.0 : rate;
}

/* Returns whether a position reports a rate: it carries enough depth and lies outside
 * the masked ends. */
static bool reported_at(const rate_config *cfg, double depth, size_t i, size_t len)
{
    return meets_min_depth(cfg->min_depth, depth) && !masked_at(cfg, i, len);
}

/* Returns one position's rate, or NaN where none is reported. */
static double rate_at(const rate_config *cfg, double events, double depth,
                      size_t i, size_t len)
{
    if (!reported_at(cfg, depth, i, len)) {
        return (double)NAN;
    }

    return ratio_of(events, depth);
}

/* Returns a value as the output stores it, in single precision. */
static double as_stored(double value)
{
    return (double)(float)value;
}

/* Returns the binomial standard error of one position's rate, or NaN where none is
 * reported. The rate and the depth enter at the precision the output stores them, so
 * that the error is a function of the stored values alone. */
static double error_at(const rate_config *cfg, double events, double depth,
                       size_t i, size_t len)
{
    double rate;
    double stored_depth;

    if (!reported_at(cfg, depth, i, len)) {
        return (double)NAN;
    }

    rate         = as_stored(ratio_of(events, depth));
    stored_depth = as_stored(depth);

    return sqrt(rate * (1.0 - rate) / stored_depth);
}

void rate_of(const rate_config *cfg, const accum *acc, size_t len,
             fmt_channel channel, double *restrict out)
{
    const double *events   = accum_const_data(acc, EVENTS_OF[channel]);
    const double *coverage = accum_const_data(acc, ACCUM_COVERAGE);

    for (size_t i = 0; i < len; i++) {
        double depth = depth_at(channel, coverage, i, len);

        out[i] = rate_at(cfg, events[i], depth, i, len);
    }
}

void rate_error_of(const rate_config *cfg, const accum *acc, size_t len,
                   fmt_channel channel, double *restrict out)
{
    const double *events   = accum_const_data(acc, EVENTS_OF[channel]);
    const double *coverage = accum_const_data(acc, ACCUM_COVERAGE);

    for (size_t i = 0; i < len; i++) {
        double depth = depth_at(channel, coverage, i, len);

        out[i] = error_at(cfg, events[i], depth, i, len);
    }
}
