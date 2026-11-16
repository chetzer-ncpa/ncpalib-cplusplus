#pragma once

#include "NCPA/exceptions.hpp"

#include <stdexcept>

namespace NCPA {
    namespace interpolation {
        NCPA_DEFINE_EXCEPTION( AmbiguousCoordinates, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( BadInput, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( DimensionsDoNotMatch, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( RequestedPointOutOfBounds, std::out_of_range )
        NCPA_DEFINE_EXCEPTION( InputSizesDoNotMatch, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( InterpolatorNotReady, std::logic_error )
        NCPA_DEFINE_EXCEPTION( InvalidDimensions, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( UnsupportedDerivative, std::out_of_range )
        NCPA_DEFINE_EXCEPTION( UnsupportedInterpolator, std::out_of_range )
    }  // namespace interpolation
}  // namespace NCPA
