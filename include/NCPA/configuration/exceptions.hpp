#pragma once

#include <stdexcept>
#include "NCPA/exceptions.hpp"

namespace NCPA {
    namespace config {
        NCPA_DEFINE_EXCEPTION( ArgumentNotDefined, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( BadMarker, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( ExpectedValueMissing, std::runtime_error )
        NCPA_DEFINE_EXCEPTION( NoConversionDefined, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( NoUnitsAssigned, std::logic_error )
        NCPA_DEFINE_EXCEPTION( SubsectionNotFound, std::out_of_range )
        NCPA_DEFINE_EXCEPTION( TooManyValues, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( UnknownMarker, std::out_of_range )
        NCPA_DEFINE_EXCEPTION( UnknownMarkerType, std::out_of_range )
        NCPA_DEFINE_EXCEPTION( ValidationFailed, std::out_of_range )
    }  // namespace config
}  // namespace NCPA
