/**
 * @file constants.hpp
 * @author Claus Hetzer
 * @date April 16, 2026
 * @version 1.0.0
 * @brief Header file providing templates and expressions for useful constants
 */

#pragma once

#include "defines.hpp"

#include <cmath>
#include <complex>

namespace NCPA {
    namespace constants {
        /**
         * @brief Pi to double precision.
         *
         * Useful for DOUBLE_EQ unit tests and other uses requiring double
         * precision
         */
        constexpr double PI
            = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;

        /**
         * @brief Euler-Mascheroni constant to double precision.
         *
         * Useful for DOUBLE_EQ unit tests and other uses requiring double
         * precision
         */
        constexpr double EGAMMA
            = 0.57721566490153286060651209008240243104215933593992;

        /**
         * @brief Convenience constant for I.
         */
        constexpr std::complex<double> I = std::complex<double>( 0.0, 1.0 );

        /**
         * @brief Returns zero as the specified type.
         * @tparam T The type to return.  Must be arithmetic or complex.
         * @return The value 0, static_cast to the correct type (if arithmetic)
         * or constructed (if complex).
         */
        template<typename T, ENABLE_FUNCTION_IF_ARITHMETIC( T )>
        constexpr T zero() {
            return static_cast<T>( 0 );
        }

        template<typename T, ENABLE_FUNCTION_IF_COMPLEX( T )>
        constexpr T zero() {
            return T( 0, 0 );
        }

        /**
         * @brief Returns zero as the specified type.
         * @tparam T The type to return.  Must be arithmetic or complex.
         * @return The value 1, static_cast to the correct type (if arithmetic)
         *      or constructed (if complex).
         */
        template<typename T, ENABLE_FUNCTION_IF_ARITHMETIC( T )>
        constexpr T one() {
            return static_cast<T>( 1 );
        }

        template<typename T, ENABLE_FUNCTION_IF_COMPLEX( T )>
        constexpr T one() {
            return T( 1, 0 );
        }

        /**
         * @brief Determines if a value evaluates to zero.
         * @tparam T The type to check.  Must be arithmetic or complex.
         * @param val The value to evaluate against zero.
         * @return true if the value evaluates to zero, false otherwise.
         */
        template<typename T, ENABLE_FUNCTION_IF_ARITHMETIC( T )>
        bool is_zero( T val ) {
            return ( std::fpclassify( val ) == FP_ZERO );
        }

        template<typename T, ENABLE_FUNCTION_IF_COMPLEX( T )>
        bool is_zero( T val ) {
            return is_zero( val.real() ) && is_zero( val.imag() );
        }

        /**
         * @brief Sets a variable to zero if its magnitude is less than a
         * tolerance value.
         *
         * This function zeros out floating-point values if their magnitude
         * is below a specified tolerance.  This is useful e.g. to deal with
         * calculations that may return insignificantly small, but still
         * nonzero, results that should be treated as zero.  Real and imaginary
         * parts of complex values are evaluated, and zeroed out, separately.
         *
         * @tparam T U The types of the variable and the tolerance.  If T is
         *      arithmetic, U defaults to T.  If T is complex, U defaults to
         *      double.
         * @param val A reference to the value to evaluate and potentially zero
         *      out.
         * @param tol The tolerance for the evaluation.
         * @return true if the value was zeroed out, false otherwise
         */
        template<typename T, typename U = T,
                 ENABLE_FUNCTION_IF_ARITHMETIC( T )>
        bool zero_out( T& val, const U& tol ) {
            if (std::abs( val ) < std::abs( tol )) {
                val = zero<T>();
                return true;
            } else {
                return false;
            }
        }

        template<typename T, typename U = double,
                 ENABLE_FUNCTION_IF_COMPLEX( T )>
        bool zero_out( T& val, const U& tol ) {
            bool nonzero = false;
            if (std::abs( val.real() ) < std::abs( tol )) {
                val.real( zero<T>().real() );
            } else {
                nonzero = true;
            }
            if (std::abs( val.imag() ) < std::abs( tol )) {
                val.imag( zero<T>().imag() );
            } else {
                nonzero = true;
            }
            return !nonzero;
        }
    }  // namespace constants
}  // namespace NCPA
