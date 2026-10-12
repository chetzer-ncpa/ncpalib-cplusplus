#pragma once

#include <stdexcept>

#define NCPA_DEFINE_EXCEPTION( _EXCEPTION_NAME_, _EXCEPTION_BASE_ )       \
    class _EXCEPTION_NAME_ : public _EXCEPTION_BASE_ {                    \
        public:                                                           \
            _EXCEPTION_NAME_( const std::string& msg = "" ) :             \
                _EXCEPTION_BASE_( std::string( #_EXCEPTION_NAME_ ) + ": " \
                                  + msg ) {}                              \
    };

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
    }  // namespace config
}  // namespace NCPA
