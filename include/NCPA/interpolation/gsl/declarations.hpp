#pragma once

#include "NCPA/interpolation/exceptions.hpp"
#include "NCPA/extern.hpp"

#if NCPA_HAVE_GSL && NCPA_USE_GSL
#  ifndef NCPA_USE_GSL_INTERPOLATION
#    define NCPA_USE_GSL_INTERPOLATION 1
#  endif
#  if NCPA_USE_GSL_INTERPOLATION
#    include "gsl/gsl_spline.h"
#    if GSL_MAJOR_VERSION >= 2
#      define NCPA_INTERPOLATION_HAVE_GSL_STEFFEN_SPLINE 1
#    else
#      define NCPA_INTERPOLATION_HAVE_GSL_STEFFEN_SPLINE 0
#    endif
namespace NCPA {
    namespace interpolation {
        namespace GSL {
            // template<typename INDEPTYPE, typename DEPTYPE,
            //          typename PARAMTYPE = const gsl_interp_type *,
            //          typename = void, typename = void>
            // class gsl_spline_1d;
            template<typename FloatType>
            class GSLSpline;
        }  // namespace GSL
    }
}

#  else
#    define NCPA_USE_GSL_INTERPOLATION                 0
#    define NCPA_INTERPOLATION_HAVE_GSL_STEFFEN_SPLINE 0
#  endif
#endif