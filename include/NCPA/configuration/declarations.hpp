#pragma once

// #include "NCPA/defines.hpp"
#include "NCPA/configuration/types/parameter_form_t.hpp"
#include "NCPA/configuration/types/parameter_type_t.hpp"
#include "NCPA/configuration/types/parse_result_t.hpp"
#include "NCPA/configuration/types/test_status_t.hpp"
#include "NCPA/types.hpp"
#include "NCPA/exceptions.hpp"

#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <utility>

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

        enum class test_result_t {
            NONE = 0,
            PENDING,
            FAILED,
            PASSED
        };
        struct validation_status_t {
                test_result_t result;
                std::string message;
        };

        template<typename T>
        using validation_function_t = std::function<validation_status_t( T )>;

        class AbstractValidation;

        template<typename T>
        class TypedValidation;

        template<typename T>
        class PredefinedValidation;

        class Validation;

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

        typedef BaseParameter parameter_t;
        typedef std::unique_ptr<parameter_t> param_ptr_t;
        typedef ScalarParameter<double> DoubleParameter;
        typedef ScalarParameter<int> IntegerParameter;
        typedef ScalarParameter<std::string> StringParameter;
        typedef ScalarParameter<bool> BooleanParameter;
        typedef ScalarParameter<std::complex<double>> ComplexParameter;
        typedef VectorParameter<double> DoubleVectorParameter;
        typedef std::unique_ptr<ValidationTest> test_ptr_t;

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
