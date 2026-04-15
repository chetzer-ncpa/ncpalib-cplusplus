#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/types.hpp"

#include <string>

namespace NCPA {
    namespace config {
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

        template<typename PARAMTYPE>
        TypedParameter<PARAMTYPE> *make_typed( Parameter& p ) {
            return dynamic_cast<TypedParameter<PARAMTYPE> *>( &p );
        }

        template<typename PARAMTYPE>
        ScalarParameter<PARAMTYPE> *make_scalar( Parameter& p ) {
            return dynamic_cast<ScalarParameter<PARAMTYPE> *>( &p );
        }

        template<typename PARAMTYPE>
        VectorParameter<PARAMTYPE> *make_vector( Parameter& p ) {
            return dynamic_cast<VectorParameter<PARAMTYPE> *>( &p );
        }

        // // floating point
        // template<typename PARAMTYPE,
        //          typename std::enable_if<
        //              std::is_floating_point<PARAMTYPE>::value, int>::type = 0>
        // parameter_type_t parameter_type() {
        //     return parameter_type_t::FLOAT;
        // }

        // // signed integer
        // template<typename PARAMTYPE,
        //          typename std::enable_if<
        //              ( std::is_integral<PARAMTYPE>::value
        //                && !( std::is_same<PARAMTYPE, bool>::value )
        //                && std::is_signed<PARAMTYPE>::value ),
        //              int>::type = 0>
        // parameter_type_t parameter_type() {
        //     return parameter_type_t::INTEGER;
        // }

        // // unsigned integer
        // template<typename PARAMTYPE,
        //          typename std::enable_if<
        //              ( std::is_integral<PARAMTYPE>::value
        //                && !( std::is_same<PARAMTYPE, bool>::value )
        //                && std::is_unsigned<PARAMTYPE>::value ),
        //              int>::type = 0>
        // parameter_type_t parameter_type() {
        //     return parameter_type_t::UNSIGNED_INTEGER;
        // }

        // // boolean
        // template<typename PARAMTYPE,
        //          typename std::enable_if<std::is_same<PARAMTYPE, bool>::value,
        //                                  int>::type = 0>
        // parameter_type_t parameter_type() {
        //     return parameter_type_t::BOOLEAN;
        // }

        // // string
        // template<typename PARAMTYPE,
        //          typename std::enable_if<
        //              ( !( std::is_arithmetic<PARAMTYPE>::value )
        //                && std::is_convertible<PARAMTYPE, std::string>::value ),
        //              int>::type = 0>
        // parameter_type_t parameter_type() {
        //     return parameter_type_t::STRING;
        // }

        // // complex
        // template<typename PARAMTYPE,
        //          typename std::enable_if<
        //              ( !( std::is_scalar<PARAMTYPE>::value )
        //                && NCPA::types::is_complex<PARAMTYPE>::value ),
        //              int>::type = 0>
        // parameter_type_t parameter_type() {
        //     return parameter_type_t::COMPLEX;
        // }

        // // enum
        // template<typename PARAMTYPE,
        //          typename std::enable_if<std::is_enum<PARAMTYPE>::value,
        //                                  int>::type = 0>
        // parameter_type_t parameter_type() {
        //     return parameter_type_t::ENUM;
        // }

        // // everything else
        // template<
        //     typename PARAMTYPE,
        //     typename std::enable_if<
        //         ( !( std::is_enum<PARAMTYPE>::value
        //              || std::is_arithmetic<PARAMTYPE>::value
        //              || std::is_convertible<PARAMTYPE, std::string>::value
        //              || NCPA::types::is_complex<PARAMTYPE>::value
        //              || ( !( std::is_scalar<PARAMTYPE>::value )
        //                   && NCPA::types::is_complex<PARAMTYPE>::value ) ) ),
        //         int>::type = 0>
        // parameter_type_t parameter_type() {
        //     return parameter_type_t::OTHER;
        // }

        // convenience functions
        template<typename T>
        test_ptr_t IsLessThan( const T& val ) {
            return test_ptr_t( new IsLessThanTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsLessThan( const T& val, bool equalOK ) {
            return equalOK
                     ? test_ptr_t( new IsLessThanOrEqualToTest<T>( val ) )
                     : test_ptr_t( new IsLessThanTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsGreaterThan( const T& val ) {
            return test_ptr_t( new IsGreaterThanTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsGreaterThan( const T& val, bool equalOK ) {
            return equalOK
                     ? test_ptr_t( new IsGreaterThanOrEqualToTest<T>( val ) )
                     : test_ptr_t( new IsGreaterThanTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsEqualTo( const T& val ) {
            return test_ptr_t( new IsEqualToTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsNotEqualTo( const T& val ) {
            return test_ptr_t( new IsNotEqualToTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsZero() {
            return IsEqualTo<T>( (T)0 );
        }

        template<typename T>
        test_ptr_t IsNotZero() {
            return IsNotEqualTo<T>( (T)0 );
        }

        template<typename T>
        test_ptr_t IsNegative() {
            return IsLessThan<T>( (T)0 );
        }

        template<typename T>
        test_ptr_t IsPositive() {
            return IsGreaterThan<T>( (T)0 );
        }

        inline test_ptr_t IsNotEmptyString() {
            return IsNotEqualTo<std::string>( "" );
        }

        inline test_ptr_t IsEmptyString() {
            return IsEqualTo<std::string>( "" );
        }
    }  // namespace config
}  // namespace NCPA
