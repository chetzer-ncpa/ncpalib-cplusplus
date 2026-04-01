#pragma once

#include "NCPA/defines.hpp"

#include <memory>

namespace NCPA {
    namespace config {
        enum class test_status_t { NONE, PENDING, FAILED, PASSED };

        enum class parameter_form_t { UNDEF, SCALAR, VECTOR };

        inline std::string to_string( parameter_form_t t ) {
            switch (t) {
                case parameter_form_t::SCALAR:
                    return "scalar";
                    break;
                case parameter_form_t::VECTOR:
                    return "vector";
                    break;
                default:
                    return "undefined";
            }
        }

        enum class parameter_type_t {
            UNDEF,
            INTEGER,
            UNSIGNED_INTEGER,
            FLOAT,
            STRING,
            BOOLEAN,
            ENUM,
            COMPLEX,
            OTHER
        };

        inline std::string to_string( parameter_type_t t ) {
            switch (t) {
                case parameter_type_t::FLOAT:
                    return "float";
                    break;
                case parameter_type_t::INTEGER:
                    return "signed integer";
                    break;
                case parameter_type_t::UNSIGNED_INTEGER:
                    return "unsigned integer";
                    break;
                case parameter_type_t::STRING:
                    return "string";
                    break;
                case parameter_type_t::BOOLEAN:
                    return "boolean";
                    break;
                case parameter_type_t::COMPLEX:
                    return "complex";
                    break;
                case parameter_type_t::ENUM:
                    return "enum";
                    break;
                case parameter_type_t::OTHER:
                    return "other";
                    break;
                default:
                    return "undefined";
            }
        }

        template<typename KEYTYPE>
        class Configurable;
        // class ConfigurationParameter;
        // template<typename T>
        // class TypedParameter;


    }  // namespace config
}  // namespace NCPA

// class _configuration_parameter;
DECLARE_CLASS_AND_SWAP_2NAMESPACE( Parameter, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( TypedParameter, NCPA, config )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_2NAMESPACE( ScalarParameter, NCPA,
                                                  config )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_2NAMESPACE( VectorParameter, NCPA,
                                                  config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( ConfigurationMap, NCPA, config )

template<typename DERIVEDTYPE, typename KEYTYPE>
void swap( NCPA::config::Configurable<KEYTYPE>& a,
           NCPA::config::Configurable<KEYTYPE>& b ) noexcept;

namespace NCPA {
    namespace config {
        typedef std::unique_ptr<Parameter> param_ptr_t;
        
        template<typename KEYTYPE>
        using param_pair_t = std::pair<KEYTYPE, param_ptr_t>;

        // floating point
        template<typename PARAMTYPE,
                 typename std::enable_if<
                     std::is_floating_point<PARAMTYPE>::value, int>::type = 0>
        parameter_type_t parameter_type() {
            return parameter_type_t::FLOAT;
        }

        // signed integer
        template<typename PARAMTYPE,
                 typename std::enable_if<
                     ( std::is_integral<PARAMTYPE>::value
                       && !( std::is_same<PARAMTYPE, bool>::value )
                       && std::is_signed<PARAMTYPE>::value ),
                     int>::type = 0>
        parameter_type_t parameter_type() {
            return parameter_type_t::INTEGER;
        }

        // unsigned integer
        template<typename PARAMTYPE,
                 typename std::enable_if<
                     ( std::is_integral<PARAMTYPE>::value
                       && !( std::is_same<PARAMTYPE, bool>::value )
                       && std::is_unsigned<PARAMTYPE>::value ),
                     int>::type = 0>
        parameter_type_t parameter_type() {
            return parameter_type_t::UNSIGNED_INTEGER;
        }

        // boolean
        template<typename PARAMTYPE,
                 typename std::enable_if<std::is_same<PARAMTYPE, bool>::value,
                                         int>::type = 0>
        parameter_type_t parameter_type() {
            return parameter_type_t::BOOLEAN;
        }

        // string
        template<typename PARAMTYPE,
                 typename std::enable_if<
                     ( !( std::is_arithmetic<PARAMTYPE>::value )
                       && std::is_convertible<PARAMTYPE, std::string>::value ),
                     int>::type = 0>
        parameter_type_t parameter_type() {
            return parameter_type_t::STRING;
        }

        // complex
        template<typename PARAMTYPE,
                 typename std::enable_if<
                     ( !( std::is_scalar<PARAMTYPE>::value )
                       && NCPA::types::is_complex<PARAMTYPE>::value ),
                     int>::type = 0>
        parameter_type_t parameter_type() {
            return parameter_type_t::COMPLEX;
        }

        // enum
        template<typename PARAMTYPE,
                 typename std::enable_if<std::is_enum<PARAMTYPE>::value,
                                         int>::type = 0>
        parameter_type_t parameter_type() {
            return parameter_type_t::ENUM;
        }

        // everything else
        template<
            typename PARAMTYPE,
            typename std::enable_if<
                ( !( std::is_enum<PARAMTYPE>::value
                     || std::is_arithmetic<PARAMTYPE>::value
                     || std::is_convertible<PARAMTYPE, std::string>::value
                     || NCPA::types::is_complex<PARAMTYPE>::value
                     || ( !( std::is_scalar<PARAMTYPE>::value )
                          && NCPA::types::is_complex<PARAMTYPE>::value ) ) ),
                int>::type = 0>
        parameter_type_t parameter_type() {
            return parameter_type_t::OTHER;
        }

    }  // namespace config
}  // namespace NCPA
