/* version.h -- the version every program reports.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#define CMUTS_VERSION "2.1.0"

/* Returns the version, followed in brackets by the commit the binary was built
 * from. Returns the version alone when CMUTS_COMMIT was not defined at compile
 * time. */
const char *cmuts_version(void);
