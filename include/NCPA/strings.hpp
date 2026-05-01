#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <sstream>
#include <string>
#include <vector>

namespace NCPA {
    namespace strings {

        template<typename T>
        std::string to_string( const std::vector<T>& v,
                               const std::string& sep = ", " ) {
            std::ostringstream oss;
            for (size_t i = 0; i < v.size(); ++i) {
                if (i > 0) {
                    oss << sep;
                }
                oss << v[ i ];
            }
            return oss.str();
        }

        template<typename T>
        std::string to_string_indexed( const std::vector<T>& v ) {
            std::ostringstream oss;
            for (size_t i = 0; i < v.size(); ++i) {
                oss << i << ": " << v[ i ] << std::endl;
            }
            oss << std::endl;
            return oss.str();
        }

        /**
        Removes leading and trailing whitespace from a string, using
        specified whitespace characters.
        @brief Removes leading and trailing whitespace from a string.
        @param orig The original string to deblank.
        @param whitespace A string containing all characters to treat as
        whitespace
        @returns The original string with leading and trailing whitespace
        removed.
        */
        static inline std::string deblank( const std::string& str,
                                           const std::string& whitespace
                                           = " \t\n\r" ) {
            const size_t strBegin = str.find_first_not_of( whitespace );
            if (strBegin == std::string::npos) {
                return "";
            }

            const size_t strEnd   = str.find_last_not_of( whitespace );
            const size_t strRange = strEnd - strBegin + 1;

            return str.substr( strBegin, strRange );
        }

        /**
        Splits a string into tokens using supplied delimiter characters.
        Similar to the Perl split() function.
        @brief Breaks a string into tokens using delimiter characters.
        @param input The string to split.
        @param delimiters The delimiter characters.  Defaults to whitespace.
        @returns A vector of the string tokens.
         */
        static inline std::vector<std::string> split( std::string input,
                                                      std::string delimiters
                                                      = " \t\n\r" ) {
            std::vector<std::string> tokens;
            tokens.clear();
            size_t firstind, lastind;
            firstind = input.find_first_not_of( delimiters );
            lastind  = input.find_first_of( delimiters, firstind );
            while (firstind != std::string::npos) {
                tokens.push_back(
                    input.substr( firstind, lastind - firstind ) );
                firstind = input.find_first_not_of( delimiters, lastind );
                lastind  = input.find_first_of( delimiters, firstind );
            }
            return tokens;
        }

        /**
        Converts a numerical time to its equivalent time string in the
        format yyyy-mm-dd hh:mm:ss.mmm
        @brief Returns a human-readable time from an epoch time.
        @param d The time to print, as an epoch time_t value.
        @returns The string representation of the time.
        */
        static inline std::string time_as_string( double d ) {
            time_t temptime = (time_t)d;
            tm *uttime      = std::gmtime( &temptime );
            double ipart, fpart;
            fpart = std::modf( d, &ipart );
            char holder[ 25 ];

            std::snprintf( holder, 25, "%4d-%02d-%02dT%02d:%02d:%06.3fZ",
                           (unsigned char)( uttime->tm_year ) + 1900,
                           (unsigned char)( uttime->tm_mon + 1 ),
                           (unsigned char)( uttime->tm_mday ),
                           (unsigned char)( uttime->tm_hour ),
                           (unsigned char)( uttime->tm_min ),
                           (double)( uttime->tm_sec ) + fpart );
            std::string s = holder;
            return s;
        }

        // Taken from https://en.cppreference.com/w/cpp/string/byte/tolower
        static inline std::string to_lower( std::string s ) {
            std::transform(
                s.begin(), s.end(), s.begin(),
                []( unsigned char c ) { return std::tolower( c ); } );
            return s;
        }

        // Taken from https://en.cppreference.com/w/cpp/string/byte/tolower
        static inline std::string to_upper( std::string s ) {
            std::transform(
                s.begin(), s.end(), s.begin(),
                []( unsigned char c ) { return std::toupper( c ); } );
            return s;
        }

        static inline bool equals_up_to( const std::string& s1,
                                         const std::string& s2, size_t n ) {
            if (s1.size() < n || s2.size() < n) {
                return false;
            }
            if (s1.substr( n ) == s2.substr( n )) {
                return true;
            }
            return false;
        }
    }  // namespace strings
}  // namespace NCPA
