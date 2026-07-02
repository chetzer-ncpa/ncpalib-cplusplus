#pragma once

// #include "NCPA/defines.hpp"
#include "NCPA/configuration/types/parameter_form_t.hpp"
#include "NCPA/configuration/types/parameter_type_t.hpp"
#include "NCPA/configuration/types/parse_result_t.hpp"
#include "NCPA/configuration/types/test_status_t.hpp"
#include "NCPA/types.hpp"
#include "NCPA/exceptions.hpp"

#include <complex>
#include <memory>
#include <string>

#ifndef NCPA_CONFIG_DEFAULT_NEWLINE_MARKER
#  define NCPA_CONFIG_DEFAULT_NEWLINE_MARKER "<br>"
#endif

namespace NCPA {
    namespace config {
        struct help_text_formatter_options_t {
                size_t indent_spaces   = 4;
                size_t max_width       = 80;
                // std::string newline_marker = "<br>";
                std::string word_regex = "([^\\s]+)";
        };

        struct help_text_section_formatter_options_t {
                size_t first_line_indent        = 0;
                size_t hanging_indent           = 0;
                size_t title_indent             = 0;
                bool indent_subsections         = true;
                bool reset_indent_after_newline = false;
                bool newline_before_title       = true;
                bool newline_after_title        = false;
        };

        // static help_text_formatter_options_t HELP_TEXT_FORMATTER_OPTIONS;
        static std::string NEWLINE_MARKER = NCPA_CONFIG_DEFAULT_NEWLINE_MARKER;

        // struct HelpTextSectionFormattingOptions {
        //         size_t first_line_indent = 0;
        //         size_t hanging_indent    = 0;
        //         size_t title_indent      = 0;
        //         bool indent_subsections  = true;
        // };

        class ValidationTest;
        class ValidationTestSuite;
        class NullaryValidationTest;
        class BaseParameter;
        class ArgumentSet;
        class Parser;

        class HelpTextSection;
        class HelpTextParagraphSection;
        class HelpTextOrganizerSection;
        class HelpTextArgumentSection;
        class HelpTextFormatter;

        template<typename T>
        class TypedParameter;

        template<typename T, typename Enable = void>
        class ScalarParameter;

        template<typename T, typename Enable = void>
        class VectorParameter;

        template<typename T, typename Enable = void>
        class ScalarParameterWithUnits;

        template<typename T, typename Enable = void>
        class VectorParameterWithUnits;

        // builder class
        class Parameter;

        template<typename KEYTYPE = std::string>
        class ConfigurationMap;

        template<typename KEYTYPE = std::string>
        class Configurable;

        class Argument;

        template<typename T>
        class TypedArgument;

        template<typename INTYPE, typename KEYTYPE = std::string>
        class Mapping;

        template<typename INTYPE, typename KEYTYPE = std::string>
        using mapping_ptr_t = std::unique_ptr<Mapping<INTYPE, KEYTYPE>>;

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

        typedef std::unique_ptr<BaseParameter> param_ptr_t;
        typedef ScalarParameter<double> DoubleParameter;
        typedef ScalarParameter<int> IntegerParameter;
        typedef ScalarParameter<std::string> StringParameter;
        typedef ScalarParameter<bool> BooleanParameter;
        typedef ScalarParameter<std::complex<double>> ComplexParameter;
        typedef VectorParameter<double> DoubleVectorParameter;
        typedef std::unique_ptr<ValidationTest> test_ptr_t;

        template<typename KEYTYPE>
        using param_pair_t = std::pair<KEYTYPE, param_ptr_t>;

    }  // namespace config
}  // namespace NCPA

inline void swap( NCPA::config::BaseParameter& a,
                  NCPA::config::BaseParameter& b ) noexcept;

template<typename T>
void swap( NCPA::config::TypedParameter<T>& a,
           NCPA::config::TypedParameter<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::ScalarParameter<T>& a,
           NCPA::config::ScalarParameter<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::ScalarParameterWithUnits<T>& a,
           NCPA::config::ScalarParameterWithUnits<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::VectorParameter<T>& a,
           NCPA::config::VectorParameter<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::VectorParameterWithUnits<T>& a,
           NCPA::config::VectorParameterWithUnits<T>& b ) noexcept;

template<typename T>
void swap( NCPA::config::ConfigurationMap<T>& a,
           NCPA::config::ConfigurationMap<T>& b ) noexcept;

template<typename KEYTYPE>
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

        template<typename T>
        class can_use_from_string {
                template<typename U>
                static auto test( int )
                    -> decltype( from_string(
                                     std::declval<const std::string&>(),
                                     std::declval<U&>() ),
                                 std::true_type() );
                template<typename>
                static std::false_type test( ... );

            public:
                static constexpr bool value = decltype( test<T>( 0 ) )::value;
        };

        template<typename T>
        typename std::enable_if<( std::is_integral<T>::value
                                  && !( std::is_enum<T>::value
                                        || std::is_same<T, bool>::value ) ),
                                void>::type
            parse_string( const std::string& str, T& val ) {
            val = static_cast<T>( std::stoi( str ) );
        }

        template<typename T>
        typename std::enable_if<std::is_floating_point<T>::value, void>::type
            parse_string( const std::string& str, T& val ) {
            val = static_cast<T>( std::stod( str ) );
        }

        template<typename T>
        typename std::enable_if<( std::is_same<T, bool>::value ), void>::type
            parse_string( const std::string& str, T& val ) {
            val = ( str.size() > 0 && ( str[ 0 ] == 't' || str[ 0 ] == 'T' ) );
        }

        template<typename T>
        typename std::enable_if<
            ( !( std::is_arithmetic<T>::value )
              && std::is_convertible<std::string, T>::value ),
            void>::type
            parse_string( const std::string& str, T& val ) {
            val = str;
        }

        template<typename T>
        typename std::enable_if<
            ( !( std::is_arithmetic<T>::value
                 || std::is_convertible<std::string, T>::value )
              && can_use_from_string<T>::value ),
            void>::type
            parse_string( const std::string& str, T& val ) {
            from_string( str, val );
        }

        template<typename T>
        typename std::enable_if<
            ( !( std::is_arithmetic<T>::value
                 || std::is_convertible<std::string, T>::value
                 || can_use_from_string<T>::value ) ),
            void>::type
            parse_string( const std::string& str, T& val ) {
            throw NCPA::NotImplementedError( "parse_string() not implemented for this type!" );
        }
    }  // namespace config
}  // namespace NCPA
