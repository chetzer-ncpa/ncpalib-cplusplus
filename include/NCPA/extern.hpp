#pragma once

// undefine any existing defines.  This file is definitive.
#ifdef NCPA_HAVE_NLOHMANN_JSON
#  undef NCPA_HAVE_NLOHMANN_JSON
#endif

#ifdef NCPA_HAVE_GSL
#  undef NCPA_HAVE_GSL
#endif

// check for existence and usage of nlohmann::json
#if __has_include( "nlohmann/json.hpp" )

#  define NCPA_HAVE_NLOHMANN_JSON 1
#  ifndef NCPA_USE_NLOHMANN_JSON
#    define NCPA_USE_NLOHMANN_JSON 1
#  endif
#else
#  define NCPA_HAVE_NLOHMANN_JSON_HPP 0
#  ifdef NCPA_USE_NLOHMANN_JSON
#    undef NCPA_USE_NLOHMANN_JSON
#  endif
#  define NCPA_USE_NLOHMANN_JSON 0
#endif

// check for existence and usage of GSL
#if __has_include( "gsl/gsl_spline.h" )
#  define NCPA_HAVE_GSL 1
#  ifndef NCPA_USE_GSL
#    define NCPA_USE_GSL 1
#  endif
#else
#  define NCPA_USE_GSL 0
#endif

// check for existence and usage of NetCDF library
