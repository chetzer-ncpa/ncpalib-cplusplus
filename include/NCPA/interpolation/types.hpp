#pragma once

namespace NCPA {
    namespace interpolation {

        enum class interpolator_type_t {
            NEAREST_NEIGHBOR,
            LANL_LINEAR,
            LANL_CUBIC,
            GSL_LINEAR,
            GSL_POLYNOMIAL,
            GSL_CUBIC,
            GSL_AKIMA,
            GSL_STEFFEN,
            GSL_CUBIC_PERIODIC,
            GSL_AKIMA_PERIODIC
        };

        enum class extrapolator_type_t {
            FORBIDDEN,
            CONSTANT,
            LINEAR,
            PERIODIC
        };
    }
}  // namespace NCPA
