#pragma once

#ifdef HAVE_NLOHMANN_JSON
#undef HAVE_NLOHMANN_JSON
#endif

#if __has_include("NCPA/extern/nlohmann/include/nlohmann/json.hpp")
#include "NCPA/extern/nlohmann/include/nlohmann/json.hpp"
#define HAVE_NLOHMANN_JSON true
#else
#define HAVE_NLOHMANN_JSON false
#endif
