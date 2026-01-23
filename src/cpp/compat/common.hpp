#ifndef _CMUTS_COMPAT_COMMON_HEADER
#define _CMUTS_COMPAT_COMMON_HEADER

#pragma message("Warning: #include \"common.hpp\" is deprecated. Use specific headers from io/common/")

// Re-export all types, streams, alignments, and file handling from the new locations
#include "io/common/types.hpp"
#include "io/common/stream.hpp"
#include "io/common/alignment.hpp"
#include "io/common/file.hpp"

// Include infrastructure dependencies that were previously bundled
#include "infra/utils.hpp"
#include "infra/mutex.hpp"

#endif
