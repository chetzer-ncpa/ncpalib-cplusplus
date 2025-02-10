#pragma once

#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/types.hpp"

namespace NCPA {
    namespace interpolation {
        namespace LANL {
            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
                _lanl_spline_1d, NCPA::interpolation::_spline_1d );
            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( linear_spline_1d,
                                                   _lanl_spline_1d );
            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( natural_cubic_spline_1d,
                                                   _lanl_spline_1d );

            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
                bicubic_spline_2d, NCPA::interpolation::_spline_2d );
            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
                natural_spline_2d, NCPA::interpolation::_spline_2d );

            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
                hybrid_spline_3d, NCPA::interpolation::_spline_3d );
        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA
