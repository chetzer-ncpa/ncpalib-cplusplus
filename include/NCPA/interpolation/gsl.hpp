#pragma once

#include "NCPA/extern.hpp"

#if NCPA_HAVE_GSL && NCPA_USE_GSL
#  ifndef NCPA_USE_GSL_INTERPOLATION
#    define NCPA_USE_GSL_INTERPOLATION 1
#  endif
#  if NCPA_USE_GSL_INTERPOLATION
#    include "NCPA/interpolation/gsl/declarations.hpp"
#    include "NCPA/interpolation/gsl/gsl_interpolators.hpp"
#  endif
#endif
