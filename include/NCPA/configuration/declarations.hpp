#pragma once

// #include "NCPA/defines.hpp"
#include "NCPA/configuration/arguments/declarations.hpp"
#include "NCPA/configuration/exceptions.hpp"
#include "NCPA/configuration/help/declarations.hpp"
#include "NCPA/configuration/parameters/declarations.hpp"
#include "NCPA/configuration/types/declarations.hpp"
#include "NCPA/configuration/validation/declarations.hpp"
#include "NCPA/exceptions.hpp"
#include "NCPA/types.hpp"

#include <complex>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace NCPA {
    namespace config {
        
        class Parser;

        template<typename KEYTYPE = std::string>
        class ConfigurationMap;

        template<typename KEYTYPE = std::string>
        class Configurable;

        template<typename INTYPE, typename KEYTYPE = std::string>
        class Mapping;

        template<typename INTYPE, typename KEYTYPE = std::string>
        using mapping_ptr_t = std::unique_ptr<Mapping<INTYPE, KEYTYPE>>;

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
            throw NCPA::NotImplementedError(
                "parse_string() not implemented for this type!" );
        }
    }  // namespace config
}  // namespace NCPA
