#pragma once
#include "NCPA/strings.hpp"
#include "NCPA/types.hpp"

#include <complex>
#include <memory>
#include <string>
#include <type_traits>

namespace NCPA {
    namespace config {
        enum class parameter_form_t { UNDEF, SCALAR, VECTOR };

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


        class BaseParameter;
        template<typename PARAMTYPE>
        class TypedParameter;

        template<typename PARAMTYPE, typename Enable = void>
        class ScalarParameter;

        template<typename PARAMTYPE, typename Enable = void>
        class VectorParameter;

        template<typename PARAMTYPE, typename Enable = void>
        class ScalarParameterWithUnits;

        template<typename PARAMTYPE, typename Enable = void>
        class VectorParameterWithUnits;


        // builder class
        class Parameter;

        typedef std::unique_ptr<BaseParameter> param_ptr_t;
        typedef ScalarParameter<double> DoubleParameter;
        typedef ScalarParameter<int> IntegerParameter;
        typedef ScalarParameter<std::string> StringParameter;
        typedef ScalarParameter<bool> BooleanParameter;
        typedef ScalarParameter<std::complex<double>> ComplexParameter;
        typedef VectorParameter<double> DoubleVectorParameter;

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
    }  // namespace config
}  // namespace NCPA

inline void from_string( const std::string& s,
                         NCPA::config::parameter_form_t& t ) {
    std::string lcs = NCPA::strings::to_lower( s );
    if (lcs == "scalar") {
        t = NCPA::config::parameter_form_t::SCALAR;
    } else if (lcs == "vector") {
        t = NCPA::config::parameter_form_t::VECTOR;
    } else {
        t = NCPA::config::parameter_form_t::UNDEF;
    }
}

inline void from_string( const std::string& s,
                         NCPA::config::parameter_type_t& t ) {
    using NCPA::strings::equals_up_to;

    std::string lcs = NCPA::strings::to_lower( s );
    if (equals_up_to( lcs, "float", 5 ) || equals_up_to( lcs, "double", 6 )) {
        t = NCPA::config::parameter_type_t::FLOAT;
    } else if (equals_up_to( lcs, "int", 3 )) {
        t = NCPA::config::parameter_type_t::INTEGER;
    } else if (lcs == "unsigned integer" || lcs == "size_t"
               || equals_up_to( lcs, "uint", 4 )) {
        t = NCPA::config::parameter_type_t::UNSIGNED_INTEGER;
    } else if (lcs == "string" || equals_up_to( lcs, "char", 4 )
               || lcs == "std::string" || lcs == "text") {
        t = NCPA::config::parameter_type_t::STRING;
    } else if (equals_up_to( lcs, "bool", 4 )) {
        t = NCPA::config::parameter_type_t::BOOLEAN;
    } else if (equals_up_to( lcs, "enum", 4 )) {
        t = NCPA::config::parameter_type_t::ENUM;
    } else if (equals_up_to( lcs, "complex", 7 )
               || equals_up_to( lcs, "std::complex", 12 )) {
        t = NCPA::config::parameter_type_t::COMPLEX;
    } else {
        t = NCPA::config::parameter_type_t::OTHER;
    }
}
