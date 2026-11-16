#pragma once


#include "NCPA/interpolation/declarations.hpp"
#include "NCPA/interpolation/gsl/declarations.hpp"

#if NCPA_USE_GSL_INTERPOLATION
#  include "gsl/gsl_spline.h"

namespace NCPA {
    namespace interpolation {
        namespace GSL {
            static interpolator_type_t interpolator_type(
                const gsl_interp_type *gsl_type ) {
                std::string gslname = gsl_type->name;
                if (gslname == "linear") {
                    return interpolator_type_t::GSL_LINEAR;
                } else if (gslname == "cspline") {
                    return interpolator_type_t::GSL_CUBIC;
                } else if (gslname == "akima") {
                    return interpolator_type_t::GSL_AKIMA;
                } else if (gslname == "akima-periodic") {
                    return interpolator_type_t::GSL_AKIMA_PERIODIC;
                } else if (gslname == "polynomial") {
                    return interpolator_type_t::GSL_POLYNOMIAL;
                } else if (gslname == "cspline-periodic") {
                    return interpolator_type_t::GSL_CUBIC_PERIODIC;
#  if NCPA_INTERPOLATION_HAVE_GSL_STEFFEN_SPLINE
                } else if (gslname == "steffen") {
                    return NCPA_INTERPOLATION_HAVE_GSL_STEFFEN_SPLINE::
                        GSL_STEFFEN;
#  endif
                } else {
                    throw UnsupportedInterpolator(
                        "Unrecognized GSL interpolator" );
                }
            }

            static const gsl_interp_type *gsl_interpolator_type(
                interpolator_type_t inttype ) {
                switch (inttype) {
                    case interpolator_type_t::GSL_LINEAR:
                        return gsl_interp_linear;
                    case interpolator_type_t::GSL_CUBIC:
                        return gsl_interp_cspline;
                    case interpolator_type_t::GSL_AKIMA:
                        return gsl_interp_akima;
                    case interpolator_type_t::GSL_AKIMA_PERIODIC:
                        return gsl_interp_akima_periodic;
                    case interpolator_type_t::GSL_POLYNOMIAL:
                        return gsl_interp_polynomial;
                    case interpolator_type_t::GSL_CUBIC_PERIODIC:
                        return gsl_interp_cspline_periodic;
#  if NCPA_INTERPOLATION_HAVE_GSL_STEFFEN_SPLINE
                    case interpolator_type_t::GSL_STEFFEN:
                        return gsl_interp_steffen;
#  endif
                    default:
                        throw UnsupportedInterpolator(
                            "Unrecognized GSL interpolator" );
                }
            }
        }  // namespace GSL
    }  // namespace interpolation
}  // namespace NCPA
#endif
