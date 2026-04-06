#pragma once

// #include "NCPA/defines.hpp"
#include "NCPA/types.hpp"

#include <memory>

namespace NCPA {
    namespace config {
        enum class test_status_t { NONE, PENDING, FAILED, PASSED };

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

        class ValidationTest;
        class ValidationTestSuite;
        class NullaryValidationTest;
        class Parameter;

        template<typename T>
        class TypedParameter;

        template<typename T, typename Enable = void>
        class ScalarParameter;

        template<typename T, typename Enable = void>
        class VectorParameter;

        template<typename T>
        class ConfigurationMap;

        template<typename KEYTYPE>
        class Configurable;


        template<typename T>
        class TypedValidationTest;
        template<typename T>
        class UnaryValidationTest;
        template<typename T>
        class BinaryValidationTest;
        template<typename T>
        class ListValidationTest;
        template<typename T>
        class IsEqualToTest;
        template<typename T>
        class IsNotEqualToTest;
        template<typename T>
        class IsNotOneOfTest;
        template<typename T>
        class IsOneOfTest;
        class WasSetTest;
        template<typename T, typename Enable = void>
        class IsBetweenTest;
        template<typename T, typename Enable = void>
        class IsGreaterThanTest;
        template<typename T, typename Enable = void>
        class IsGreaterThanOrEqualToTest;
        template<typename T, typename Enable = void>
        class IsLessThanTest;
        template<typename T, typename Enable = void>
        class IsLessThanOrEqualToTest;

        typedef std::unique_ptr<Parameter> param_ptr_t;
        typedef ScalarParameter<double> DoubleParameter;
        typedef ScalarParameter<int> IntegerParameter;
        typedef ScalarParameter<std::string> StringParameter;
        typedef ScalarParameter<bool> BooleanParameter;
        typedef std::unique_ptr<ValidationTest> test_ptr_t;

        template<typename KEYTYPE>
        using param_pair_t = std::pair<KEYTYPE, param_ptr_t>;

    }  // namespace config
}  // namespace NCPA

void swap( NCPA::config::Parameter& a, NCPA::config::Parameter& b ) noexcept;

template<typename T>
void swap( NCPA::config::TypedParameter<T>& a,
           NCPA::config::TypedParameter<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::ScalarParameter<T>& a,
           NCPA::config::ScalarParameter<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::VectorParameter<T>& a,
           NCPA::config::VectorParameter<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::ConfigurationMap<T>& a,
           NCPA::config::ConfigurationMap<T>& b ) noexcept;

template<typename DERIVEDTYPE, typename KEYTYPE>
void swap( NCPA::config::Configurable<KEYTYPE>& a,
           NCPA::config::Configurable<KEYTYPE>& b ) noexcept;


void swap( NCPA::config::ValidationTest& a,
           NCPA::config::ValidationTest& b ) noexcept;


void swap( NCPA::config::ValidationTestSuite& a,
           NCPA::config::ValidationTestSuite& b ) noexcept;


void swap( NCPA::config::NullaryValidationTest& a,
           NCPA::config::NullaryValidationTest& b ) noexcept;


template<typename T>
void swap( NCPA::config::TypedValidationTest<T>& a,
           NCPA::config::TypedValidationTest<T>& b ) noexcept;
template<typename T>
void swap( NCPA::config::UnaryValidationTest<T>& a,
           NCPA::config::UnaryValidationTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::BinaryValidationTest<T>& a,
           NCPA::config::BinaryValidationTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::ListValidationTest<T>& a,
           NCPA::config::ListValidationTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::IsEqualToTest<T>& a,
           NCPA::config::IsEqualToTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::IsNotEqualToTest<T>& a,
           NCPA::config::IsNotEqualToTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::IsNotOneOfTest<T>& a,
           NCPA::config::IsNotOneOfTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::IsOneOfTest<T>& a,
           NCPA::config::IsOneOfTest<T>& b ) noexcept;

void swap( NCPA::config::WasSetTest& a, NCPA::config::WasSetTest& b ) noexcept;

template<typename T>
void swap( NCPA::config::IsBetweenTest<T>& a,
           NCPA::config::IsBetweenTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::IsGreaterThanTest<T>& a,
           NCPA::config::IsGreaterThanTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::IsGreaterThanOrEqualToTest<T>& a,
           NCPA::config::IsGreaterThanOrEqualToTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::IsLessThanTest<T>& a,
           NCPA::config::IsLessThanTest<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::IsLessThanOrEqualToTest<T>& a,
           NCPA::config::IsLessThanOrEqualToTest<T>& b ) noexcept;

namespace NCPA {
    namespace config {
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
