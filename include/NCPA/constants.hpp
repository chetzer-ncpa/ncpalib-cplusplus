#pragma once

#include "defines.hpp"
#include <complex>
#include <cmath>

namespace NCPA {
    namespace constants {
        // double precision to enable DOUBLE_EQUAL unit tests
        constexpr double PI
            = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;

        constexpr std::complex<double> I = std::complex<double>( 0.0, 1.0 );

        /**
         * Returns zero as the specified type.
         */
        template<typename T, ENABLE_FUNCTION_IF_ARITHMETIC( T )>
        constexpr T zero() {
            // T z = 0;
            return (T)( 0 );
        }

        template<typename T, ENABLE_FUNCTION_IF_COMPLEX( T )>
        constexpr T zero() {
            // T z( 0, 0 );
            return T( 0, 0 );
        }

        /**
         * Returns unity as the specified type.
         */
        template<typename T, ENABLE_FUNCTION_IF_ARITHMETIC( T )>
        constexpr T one() {
            // T z = 1;
            return (T)( 1 );
        }

        template<typename T, ENABLE_FUNCTION_IF_COMPLEX( T )>
        constexpr T one() {
            // T z( 1, 0 );
            return T( 1, 0 );
        }

        template<typename T, ENABLE_FUNCTION_IF_ARITHMETIC( T )>
        bool is_zero( T val ) {
            return ( std::fpclassify( val ) == FP_ZERO );
        }

        template<typename T, ENABLE_FUNCTION_IF_COMPLEX( T )>
        bool is_zero( T val ) {
            return is_zero( val.real() ) && is_zero( val.imag() );
        }
    }
}