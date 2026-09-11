#pragma once

#include <stdexcept>
#include <string>
#include <unordered_map>

#define DECLARE_CONFIGURATION_EXCEPTION( _CLASSNAME_, _ERROR_CODE_ ) \
    class _CLASSNAME_ : public ConfigurationException {              \
        public:                                                      \
            static constexpr int error_code = _ERROR_CODE_;          \
            _CLASSNAME_( const std::string& message ) :              \
                ConfigurationException( #_CLASSNAME_, message,       \
                                        _ERROR_CODE_ ) {}            \
    };

namespace NCPA {
    namespace config {

        class ConfigurationException : public std::runtime_error {
            private:
                int _error_code;

            public:
                ConfigurationException( const std::string& exception_name,
                                        const std::string& message,
                                        int code ) :
                    std::runtime_error( exception_name + ": " + message ),
                    _error_code { code } {}

                virtual int code() const noexcept { return _error_code; }
        };

        // Argument-related exceptions, starting at 300
        DECLARE_CONFIGURATION_EXCEPTION( ExpectedValueNotFound, 301 )
        DECLARE_CONFIGURATION_EXCEPTION( ValueRequestedFromFlag, 302 )
        DECLARE_CONFIGURATION_EXCEPTION( TooManyValuesAssigned, 303 )
        DECLARE_CONFIGURATION_EXCEPTION( NoDefaultValue, 304 )
        DECLARE_CONFIGURATION_EXCEPTION( ValueIndexOutOfRange, 305 )
        DECLARE_CONFIGURATION_EXCEPTION( ArgumentNotFound, 306 )
    }  // namespace config
}  // namespace NCPA
