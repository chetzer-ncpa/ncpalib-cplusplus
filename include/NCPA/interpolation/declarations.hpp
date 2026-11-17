#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/extern.hpp"
#include "NCPA/interpolation/exceptions.hpp"
#include "NCPA/interpolation/gsl/declarations.hpp"
#include "NCPA/interpolation/lanl/declarations.hpp"

namespace NCPA {
    namespace interpolation {
        enum class interpolator_type_t {
            NEAREST_NEIGHBOR,
            STRATIFIED_2D,
            NEAREST_NEIGHBOR_2D,
            LANL_LINEAR,
            LANL_CUBIC,
            LANL_NATURAL_2D,
            LANL_BICUBIC_2D,
            LANL_HYBRID_3D,
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
            ZERO,
            LINEAR,
            PERIODIC
        };

        template<typename FloatType = double>
        class abstract_spline;

        template<typename FloatType = double,
                 typename AxesType  = std::vector<FloatType>,
                 typename FType     = NCPA::arrays::SmartVector<FloatType>>
        class typed_spline;

        template<typename FloatType>
        using spline_ptr_t = std::unique_ptr<abstract_spline<FloatType>>;

        template<typename FloatType = double>
        class NearestNeighborSpline;

        template<typename FloatType = double, typename Enabler = void>
        class Interpolator;

        class InterpolatorFactory;
    }  // namespace interpolation
}  // namespace NCPA
