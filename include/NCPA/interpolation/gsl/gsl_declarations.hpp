#pragma once

#ifndef NCPA_INTERPOLATION_IGNORE_GSL
#  include "NCPA/interpolation/defines.hpp"

#  ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE
#    include "gsl/gsl_interp.h"
#    include "gsl/gsl_spline.h"
#    include "NCPA/interpolation/types.hpp"

namespace NCPA {
    namespace interpolation {
        namespace GSL {
            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE_WITH_PARAM(
                gsl_spline_1d, NCPA::interpolation::_spline_1d,
                const gsl_interp_type * );
        }
    }  // namespace interpolation
}  // namespace NCPA
#  endif
#endif
