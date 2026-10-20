/* version.c -- the version with the commit it was built from.
 *
 * The Makefile defines CMUTS_COMMIT when the source directory is a git
 * checkout. The value ends in -dirty when the checkout had uncommitted changes.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "version.h"

#ifdef CMUTS_COMMIT
#define CMUTS_VERSION_STAMP CMUTS_VERSION " (" CMUTS_COMMIT ")"
#else
#define CMUTS_VERSION_STAMP CMUTS_VERSION
#endif

const char *cmuts_version(void)
{
    return CMUTS_VERSION_STAMP;
}
