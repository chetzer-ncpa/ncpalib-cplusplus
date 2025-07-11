#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/constants.hpp"
#include "NCPA/defines.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstring>
#include <functional>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

namespace NCPA {
    namespace math {

        using NCPA::constants::I;
        using NCPA::constants::is_zero;
        using NCPA::constants::one;
        using NCPA::constants::PI;
        using NCPA::constants::zero;

        // template<typename T, ENABLE_FUNCTION_IF_ARITHMETIC( T )>
        // bool is_zero( T val ) {
        // return ( std::fpclassify( val ) == FP_ZERO );
        // }

        // template<typename T, ENABLE_FUNCTION_IF_COMPLEX( T )>
        // bool is_zero( T val ) {
        //     return is_zero( val.real() ) && is_zero( val.imag() );
        // }


        template<typename T>
        bool equals( T x, T y, size_t n = 1,
                     ENABLE_FUNCTION_IF_ARITHMETIC( T ) ) {
            if (std::numeric_limits<T>::is_exact) {
                return x == y;
            } else {
                // Since `epsilon()` is the gap size (ULP, unit in the last
                // place) of floating-point numbers in interval [1, 2), we can
                // scale it to the gap size in interval [2^e, 2^{e+1}), where
                // `e` is the exponent of `x` and `y`.

                // If `x` and `y` have different gap sizes (which means they
                // have different exponents), we take the smaller one. Taking
                // the bigger one is also reasonable, I guess.
                const T m = std::min( std::abs( x ), std::abs( y ) );

                // Subnormal numbers have fixed exponent, which is
                // `min_exponent - 1`.
                const int exp = m < std::numeric_limits<T>::min()
                                  ? std::numeric_limits<T>::min_exponent - 1
                                  : std::ilogb( m );

                // We consider `x` and `y` equal if the difference between them
                // is within `n` ULPs.
                return std::abs( x - y )
                    <= n
                           * std::ldexp( std::numeric_limits<T>::epsilon(),
                                         exp );
            }
        }

        template<typename T>
        bool is_even( T x, ENABLE_FUNCTION_IF_INTEGRAL( T ) ) {
            return ( ( x & 1 ) == 0 );
        }

        template<typename T>
        bool is_odd( T x, ENABLE_FUNCTION_IF_INTEGRAL( T ) ) {
            return !is_even( x );
        }

        template<typename T>
        bool equals( T x, T y, size_t n = 1,
                     ENABLE_FUNCTION_IF_COMPLEX( T ) ) {
            return equals( x.real(), y.real() )
                && equals( x.imag(), y.imag() );
        }

        /**
         * Returns -1, 0, or 1 depending on the sign of the argument.
         * @brief Returns -1, 0, or 1 depending on the sign of the argument.
         * @param n The value to test for sign
         * @returns -1 if n<0, 1 if n>0, 0 otherwise
         */
        template<typename T>
        int sign( T n ) {
            if (n > 0) {
                return 1;
            } else if (n < 0) {
                return -1;
            } else {
                return 0;
            }
        }

        /**
         * Computes the simple square with a multiply instead of an
         * exponentiation.
         */
        template<typename T>
        T square( T n ) {
            return n * n;
        }

        /**
         * Finds and returns the indices of the array elements before and after
         * the supplied value.
         */
        template<typename T>
        bool find_interval_inclusive( const T *z, size_t NZ, T val,
                                      size_t& bottom, size_t& top ) {
            const double *it = std::lower_bound( z, z + NZ, val );
            size_t diff      = it - z;
            if (diff == 0) {
                bottom = 0;
                if (val < z[ 0 ]) {
                    top = 0;
                    return false;
                } else {
                    top = 1;
                    return true;
                }
            } else if (diff == NZ) {
                top = NZ;
                if (val > z[ NZ - 1 ]) {
                    bottom = NZ;
                    return false;
                } else {
                    bottom = NZ - 1;
                    return true;
                }
            } else {
                bottom = diff - 1;
                top    = diff;
                return true;
            }
        }

        template<typename T>
        bool find_interval_inclusive( const std::vector<T>& z, T val,
                                      size_t& bottom, size_t& top ) {
            return find_interval_inclusive( &z[ 0 ], z.size(), val, bottom,
                                            top );
        }

        /**
        Finds and returns the index of the array element that is closest
        in value to a supplied reference value.  Minimizes abs( z - zref ).
        @brief Finds the index of the closest array element to a value.
        @param z The array to check.
        @param NZ The size of the array.
        @param zs The reference value.
        @returns The index of the array element closest to the reference value.
        */
        template<typename T>
        size_t find_closest_index( T *z, size_t NZ, T zs ) {
            if (NZ == 1) {
                return 0;
            }
            double diff = 0.0, mindiff = std::numeric_limits<T>::max();
            size_t tmpind = 0;

            for (size_t i = 0; i < NZ; i++) {
                diff = std::fabs( ( (double)z[ i ] ) - ( (double)zs ) );
                if (diff < mindiff) {
                    tmpind  = i;
                    mindiff = diff;
                }
            }

            return tmpind;
        }

        /**
        Finds and returns the index of the vector element that is closest
        in value to a supplied reference value.  Minimizes abs( z - zref ).
        @brief Finds the index of the closest vector element to a value.
        @param z The vector to check.
        @param zs The reference value.
        @returns The index of the vector element closest to the reference
        value.
        */
        template<typename T>
        size_t find_closest_index( std::vector<T> z, T zs ) {
            if (z.size() == 1) {
                return 0;
            }
            T diff = 0.0, mindiff = std::numeric_limits<T>::max();
            size_t tmpind = 0, i;
            for (i = 0; i < z.size(); i++) {
                diff = std::abs( z[ i ] - zs );
                if (diff < mindiff) {
                    tmpind  = i;
                    mindiff = diff;
                }
            }
            return tmpind;
        }

        /**
        @brief Converts from Cartesian to polar coordinates.
        @param x The x Cartesian input coordinate.
        @param y The y Cartesian input coordinate.
        @param r The r polar output coordinate.
        @param theta_rad The theta output polar coordinate, in radians.  Will
        return 0.0 if x==0 and y==0
        */
        template<typename T>
        void cart2pol( T x, T y, T& r, T& theta_rad ) {
            r = std::sqrt( x * x + y * y );
            try {
                theta_rad = std::atan2( y, x );
            } catch (std::domain_error e) {
                theta_rad = 0.0;
            }
        }

        /**
        @brief Converts from polar to Cartesian coordinates.
        @param r The r polar input coordinate.
        @param theta_rad The theta input polar coordinate, in radians.
        @param x The x Cartesian output coordinate.
        @param y The y Cartesian output coordinate.
        */
        template<typename T>
        void pol2cart( T r, T theta_rad, T& x, T& y ) {
            x = r * std::cos( theta_rad );
            y = r * std::sin( theta_rad );
        }

        /**
        @brief Converts from degrees to radians.
        @param deg_in The input value in degrees.
        @returns The output value in radians.
        */
        template<typename T>
        T deg2rad( T deg_in ) {
            return (T)( ( (double)deg_in ) * PI / 180.0 );
        }

        /**
        @brief Converts from radians to degrees.
        @param rad_in The input value in radians.
        @returns The output value in degrees.
        */
        template<typename T>
        T rad2deg( T rad_in ) {
            return (T)( ( (double)rad_in ) * 180.0 / PI );
        }

        /**
         * @brief Converts from math angles to geographic azimuth
         * @param deg_in The input value in degrees CCW from the x-axis
         * @returns The same vale in degrees CW from north
         */
        template<typename T>
        T math2az( T deg_in ) {
            T deg_out = ( (T)90.0 ) - deg_in;
            while (deg_out < 0.0) {
                deg_out += 360.0;
            }
            return deg_out;
        }

        /**
         * @brief Converts from geographic azimuth to math angle
         * @param deg_in The input value in degrees CW from north
         * @returns The same vale in degrees CCW from the x-axis
         */
        template<typename T>
        T az2math( T deg_in ) {
            T deg_out = ( (T)90.0 ) - deg_in;
            while (deg_out < 0.0) {
                deg_out += 360.0;
            }
            return deg_out;
        }

        /**
        @brief Returns the maximum value from a vector.
        @param vals The vector holding the values to check.
        @returns The maximum value found in the vector.
        */
        template<typename T>
        T max( const std::vector<T>& vals ) {
            T maxval = vals.front();
            for (auto cit = vals.cbegin(); cit != vals.cend(); ++cit) {
                maxval = std::max( *cit, maxval );
            }
            return maxval;
        }

        /**
        @brief Returns the index of the maximum value from an array.
        @param vals The array holding the values to check.
        @param size The size of the array.
        @returns The index of the maximum value found in the array.
        */
        template<typename T>
        size_t index_of_max( const T *vals, size_t size ) {
            T maxval   = vals[ 0 ];
            size_t ind = 0;
            for (size_t i = 1; i < size; i++) {
                if (vals[ i ] > maxval) {
                    maxval = vals[ i ];
                    ind    = i;
                }
            }
            return ind;
        }

        /**
        @brief Returns the maximum value from an array.
        @param vals The array holding the values to check.
        @param size The size of the array.
        @returns The maximum value found in the array.
        */
        template<typename T>
        T max( const T *vals, size_t size ) {
            T maxval = vals[ 0 ];
            for (size_t i = 1; i < size; i++) {
                maxval = std::max( vals[ i ], maxval );
            }
            return maxval;
        }

        /**
        @brief Returns the index of the maximum value from a vector.
        @param vals The vector holding the values to check.
        @returns The index of the maximum value found in the vector.
        */
        template<typename T>
        size_t index_of_max( const std::vector<T>& vals ) {
            T maxval   = vals[ 0 ];
            size_t ind = 0;
            for (size_t i = 1; i < vals.size(); i++) {
                if (vals[ i ] > maxval) {
                    maxval = vals[ i ];
                    ind    = i;
                }
            }
            return ind;
        }

        /**
        @brief Returns the minimum value from an array.
        @param vals The array holding the values to check.
        @param size The size of the array.
        @returns The minimum value found in the array.
        */
        template<typename T>
        T min( const T *vals, size_t size ) {
            T minval = vals[ 0 ];
            for (size_t i = 1; i < size; i++) {
                minval = std::min( vals[ i ], minval );
            }
            return minval;
        }

        /**
        @brief Returns the index of the minimum value from an array.
        @param vals The array holding the values to check.
        @param size The size of the array.
        @returns The index of the minimum value found in the array.
        */
        template<typename T>
        size_t index_of_min( const T *vals, size_t size ) {
            T minval   = vals[ 0 ];
            size_t ind = 0;
            for (size_t i = 1; i < size; i++) {
                if (vals[ i ] < minval) {
                    minval = vals[ i ];
                    ind    = i;
                }
            }
            return ind;
        }

        /**
        @brief Returns the minimum value from a vector.
        @param vals The vector holding the values to check.
        @returns The minimum value found in the vector.
        */
        template<typename T>
        T min( const std::vector<T>& vals ) {
            T minval = vals.front();
            for (auto cit = vals.cbegin(); cit != vals.cend(); ++cit) {
                minval = std::min( *cit, minval );
            }
            return minval;
        }

        /**
        @brief Returns the index of the minimum value from a vector.
        @param vals The vector holding the values to check.
        @returns The index of the minimum value found in the vector.
        */
        template<typename T>
        size_t index_of_min( const std::vector<T>& vals ) {
            T minval   = vals.front();
            size_t ind = 0;
            for (size_t i = 1; i < vals.size(); i++) {
                if (vals[ i ] < minval) {
                    minval = vals[ i ];
                    ind    = i;
                }
            }
            return ind;
        }

        /**
        Finds and returns the indices of a grid point nearest to a supplied
        point.
        @brief Finds the index of the closest grid point to a coordinate pair.
        @param _x1 The first dimension grid points
        @param _x2 The second dimension grid points
        @param x1 The first coordinate to check
        @param x2 The second coordinate to check
        @returns The coordinates of the closest grid point
        */
        template<typename T>
        std::pair<size_t, size_t> find_closest_point( std::vector<T> _x1,
                                                      std::vector<T> _x2, T x1,
                                                      T x2 ) {
            size_t minind1, maxind1, minind2, maxind2, ind1, ind2;
            if (!NCPA::math::find_interval_inclusive( _x1, x1, minind1,
                                                      maxind1 )
                && minind1 > 0) {
                --minind1;
                --maxind1;
            }

            if (!NCPA::math::find_interval_inclusive( _x2, x2, minind2,
                                                      maxind2 )
                && minind2 > 0) {
                --minind2;
                --maxind2;
            }

            std::vector<T> distance_to_corners
                = { std::sqrt( NCPA::math::square( x1 - _x1[ minind1 ] )
                               + NCPA::math::square( x2 - _x2[ minind2 ] ) ),
                    std::sqrt( NCPA::math::square( x1 - _x1[ maxind1 ] )
                               + NCPA::math::square( x2 - _x2[ minind2 ] ) ),
                    std::sqrt( NCPA::math::square( x1 - _x1[ minind1 ] )
                               + NCPA::math::square( x2 - _x2[ maxind2 ] ) ),
                    std::sqrt( NCPA::math::square( x1 - _x1[ maxind1 ] )
                               + NCPA::math::square( x2 - _x2[ maxind2 ] ) ) };
            switch (NCPA::math::index_of_min( distance_to_corners )) {
                case 0:
                    return { minind1, minind2 };
                    break;
                case 1:
                    return { maxind1, minind2 };
                    break;
                case 2:
                    return { minind1, maxind2 };
                    break;
                case 3:
                    return { maxind1, maxind2 };
                    break;
                default:
                    throw std::logic_error(
                        "eval_f: something went wrong, this should "
                        "never happen" );
            }
        }

        /**
         * Converts a vector of reals into the equivalent vector of
         * complex<>s with the imaginary value set to zero.
         * @brief Converts a real vector to a complex<> vector
         * @param real The real vector to convert
         * @returns The complex<> vector to return
         */
        template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
        std::vector<std::complex<T>> real2complex( const std::vector<T>& in ) {
            std::vector<std::complex<T>> out;
            out.reserve( in.size() );
            for (auto it = in.cbegin(); it != in.cend(); ++it) {
                out.push_back( std::complex<T>( *it, 0.0 ) );
            }
            return out;
        }

        /**
         * Converts an array of reals into the equivalent array of
         * complex<>s with the imaginary value set to zero.
         * @brief Converts a real array to a complex<> array
         * @param n The size of the arrays
         * @param real The real array to convert
         * @param out The complex<> array to return
         */
        template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
        void real2complex( size_t n, const T *in, std::complex<T> *out ) {
            for (size_t i = 0; i < n; i++) {
                out[ i ] = std::complex<T>( in[ i ], 0.0 );
            }
        }

        /**
         * Converts two vectors of real values into the equivalent vector of
         * complex<>s where one vector represents the real parts and the
         * other the imaginary parts.
         * @brief Converts two real vectors to a complex<> vector
         * @param real The real vector to convert as the real parts
         * @param imag The real vector to convert as the imaginary parts
         */
        template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
        std::vector<std::complex<T>> real2complex(
            const std::vector<T>& real, const std::vector<T>& imag ) {
            size_t Nr = real.size();
            size_t Ni = imag.size();
            size_t N  = std::max( Nr, Ni );
            std::vector<std::complex<T>> out( N );
            for (size_t i = 0; i < N; i++) {
                if (i < Nr) {
                    out[ i ].real( real[ i ] );
                }
                if (i < Ni) {
                    out[ i ].imag( imag[ i ] );
                }
            }
            return out;
        }

        /**
         * Converts two arrays of reals into the equivalent array of
         * complex<>s where one array represents the real parts and the
         * other the imaginary parts.
         * @brief Converts a real array to a complex<> array
         * @param n The size of the arrays
         * @param real The real array to convert as the real parts
         * @param imag The real array to convert as the real parts
         * @param out The complex<> array to return
         */
        template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
        void real2complex( size_t n, const T *real, const T *imag,
                           std::complex<T> *out ) {
            for (size_t i = 0; i < n; i++) {
                out[ i ] = std::complex<T>( real[ i ], imag[ i ] );
            }
        }

        /**
         * Converts a vector of complex<>s into the equivalent vectors of
         * reals where one vector represents the real parts and the other the
         * imaginary parts.
         * @brief Converts a complex<> vector to two real vectors
         * @param in The complex array to decompose
         * @param real The real parts
         * @param imag The imaginary parts
         */
        template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
        void complex2real( const std::vector<std::complex<T>>& in,
                           std::vector<T>& real, std::vector<T>& imag ) {
            real.resize( in.size() );
            imag.resize( in.size() );
            for (size_t i = 0; i < in.size(); i++) {
                real[ i ] = in[ i ].real();
                imag[ i ] = in[ i ].imag();
            }
        }

        /**
         * Converts an array of complex<>s into the equivalent arrays of
         * reals where one array represents the real parts and the other the
         * imaginary parts.
         * @brief Converts a complex<> array to two real arrays
         * @param n The size of the arrays
         * @param in The complex array to decompose
         * @param real The real parts
         * @param imag The imaginary parts
         */
        template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
        void complex2real( size_t n, const std::complex<T> *in, T *real,
                           T *imag ) {
            for (size_t i = 0; i < n; i++) {
                real[ i ] = in[ i ].real();
                imag[ i ] = in[ i ].imag();
            }
        }

        /**
        Generates a vector of N random real numbers in a given range.
        @brief Generates a vector of random numbers.
        @param N Number of random values to generate.
        @param minrange The lower end of the range.
        @param maxrange The upper end of the range.
        @returns A vector of random numbers.
        */
        template<typename T>
        std::vector<T> random_numbers( size_t N, T minrange, T maxrange,
                                       ENABLE_FUNCTION_IF_INTEGRAL( T ) ) {
            std::vector<T> randn;
            randn.reserve( N );
            std::random_device rd;
            std::mt19937 generator( rd() );
            std::uniform_int_distribution<T> distribution( minrange,
                                                           maxrange );
            for (size_t i = 0; i < N; i++) {
                randn.push_back( distribution( generator ) );
            }
            return randn;
        }

        template<typename T>
        std::vector<T> random_numbers( size_t N, T minrange, T maxrange,
                                       ENABLE_FUNCTION_IF_REAL( T ) ) {
            std::vector<T> randn;
            randn.reserve( N );
            std::random_device rd;
            std::mt19937 generator( rd() );
            std::uniform_real_distribution<T> distribution( minrange,
                                                            maxrange );
            for (size_t i = 0; i < N; i++) {
                randn.push_back( distribution( generator ) );
            }
            return randn;
        }

        template<typename T>
        std::vector<T> random_numbers( size_t N,
                                       typename T::value_type minrange,
                                       typename T::value_type maxrange,
                                       ENABLE_FUNCTION_IF_COMPLEX( T ) ) {
            std::vector<typename T::value_type> randr
                = random_numbers<typename T::value_type>( N, minrange,
                                                          maxrange ),
                randi = random_numbers<typename T::value_type>( N, minrange,
                                                                maxrange );
            return real2complex( randr, randi );
        }

        // /**
        // Generates a vector of N random complex numbers in a given range.
        // @brief Generates a vector of random complex numbers.
        // @param N Number of random values to generate.
        // @param minreal The lower end of the range for the real parts.
        // @param maxreal The upper end of the range for the real parts.
        // @param minimag The lower end of the range for the imaginary parts.
        // @param maximag The upper end of the range for the imaginary parts.
        // @returns A vector of random numbers.
        // */
        // template<typename T>
        // std::vector<std::complex<T>> random_numbers( size_t N, T minreal,
        //                                              T maxreal, T minimag,
        //                                              T maximag ) {
        //     std::vector<T> randr = random_numbers<T>( N, minreal, maxreal ),
        //                    randi = random_numbers<T>( N, minimag, maximag );
        //     return real2complex( randr, randi );
        // }

        /**
         * Returns a single random real number in a given range.
         * @param minrange The lower end of the range.
         * @param maxrange The upper end of the range.
         * @returns A single random real number.
         */
        template<typename T, typename U = typename T::value_type>
        T random_number( U minrange, U maxrange ) {
            return random_numbers<T>( 1, minrange, maxrange )[ 0 ];
        }

        // /**
        //  * Returns a single random complex number in a given range.
        //  * @param minreal The lower end of the range for the real part.
        //  * @param maxreal The upper end of the range for the real part.
        //  * @param minimag The lower end of the range for the imaginary part.
        //  * @param maximag The upper end of the range for the imaginary part.
        //  * @returns A single random complex number.
        //  */
        // template<typename T>
        // std::complex<T> random_number( T minreal, T maxreal, T minimag,
        //                                T maximag ) {
        //     return random_numbers<std::complex<T>>( 1, minreal, maxreal,
        //     minimag,
        //                               maximag )[ 0 ];
        // }

        /**
        @brief Evaluate a polynomial from its coefficients.
        @param coeffs The coefficients in increasing order starting with power
        0 (i.e. the constant).
        @param x The value at which to evaluate the polynomial.
        @returns The calculated value of the polynomial at x.
        */
        template<typename T>
        T evalpoly( const std::vector<T>& coeffs, T x ) {
            T val = 0.0;
            for (size_t i = 0; i < coeffs.size(); i++) {
                val += coeffs[ i ] * std::pow( x, (double)i );
            }
            return val;
        }

        /**
        @brief Evaluate a polynomial from its coefficients.
        @param N the number of coefficients (not its degree).
        @param coeffs The coefficients in increasing order starting with power
        0 (i.e. the constant).
        @param x The value at which to evaluate the polynomial.
        @returns The calculated value of the polynomial at x.
        */
        template<typename T>
        T evalpoly( size_t N, const T *coeffs, T x ) {
            std::vector<T> v( coeffs, coeffs + N );
            return evalpoly<T>( v, x );
        }

        /**
        Performs a simple linear interpolation.  Returns the value at x
        between (x1,y1) and (x2,y2).
        @brief Performs a simple linear interpolation.
        @param x1 The x coordinate of the start point.
        @param y1 The y coordinate of the start point.
        @param x2 The x coordinate of the end point.
        @param y2 The y coordinate of the end point.
        @param x The x value at which to perform the interpolation.
        @returns The value of y at x.
        */
        template<typename T>
        T linear_interpolation( T x1, T y1, T x2, T y2, T x ) {
            return y1 + ( x - x1 ) * ( y2 - y1 ) / ( x2 - x1 );
        }

        /**
        Returns a vector of N values, linearly spaced between two supplied end
        points.
        @brief Returns a new linearly spaced vector of values.
        @param firstval The starting value of the array.
        @param lastval The ending value of the array.
        @param N The number of values to return.
        @returns The vector containing the values.
        */
        template<typename T>
        std::vector<T> linspace( T firstval, T lastval, size_t N ) {
            T stepsize = ( lastval - firstval ) / ( (T)( N - 1 ) );
            std::vector<T> vec( N );
            for (size_t i = 0; i < N; i++) {
                vec[ i ] = firstval + (T)i * stepsize;
            }
            return vec;
        }

        /**
        Returns an array of N values, linearly spaced between two
        supplied end points, placed in a pre-allocated array.
        @brief Populates a linearly spaced array of values.
        @param firstval The starting value of the array.
        @param lastval The ending value of the array.
        @param N The number of values to return.
        @param vec The array to place the values in.  If nullptr, will be
        dynamically allocated.
        */
        template<typename T>
        void linspace( T firstval, T lastval, size_t N, T *& vec ) {
            if (vec == nullptr) {
                vec = NCPA::arrays::zeros<T>( N );
            }
            std::vector<T> tmpvec = linspace<T>( firstval, lastval, N );
            std::copy( tmpvec.begin(), tmpvec.end(), vec );
        }

        /**
        Returns a vector of N values, linearly spaced between two supplied end
        points.
        @brief Returns a new linearly spaced vector of values.
        @param firstval The starting value of the array.
        @param lastval The ending value of the array.
        @param N The number of values to return.
        @returns The vector containing the values.
        */
        template<typename T>
        // requires std::floating_point<T>
        std::vector<T> logspace( T firstval, T lastval, size_t N,
                                 ENABLE_FUNCTION_IF_REAL( T ) ) {
            // static_assert( std::floating_point<T>, "logspace() only supports
            // floating-point types" );
            T la                = std::log10( firstval );
            T lb                = std::log10( lastval );
            std::vector<T> logs = NCPA::math::linspace( la, lb, N );
            for (size_t i = 0; i < N; i++) {
                logs[ i ] = std::pow( (T)( 10.0 ), logs[ i ] );
            }
            return logs;
        }

        // logarithmically spaced vector from a to b (NOT from 10^a to 10^b)
        /**
        Returns an array of N values, logarithmically spaced between two
        supplied end points a and b, placed in a pre-allocated array.  Note
        that this returns the values, not the exponents of the values as in
        some other logspace functions; in other words, it returns the array
        (a, ..., b) and NOT (10^a, ..., 10^b).
        @brief Populates a linearly spaced array of values.
        @param a The starting value of the array.
        @param b The ending value of the array.
        @param N The number of values to return.
        @param ls The array to place the values in.
        */
        template<typename T>
        // requires std::floating_point<T>
        void logspace( T a, T b, size_t N, T *& ls,
                       ENABLE_FUNCTION_IF_REAL( T ) ) {
            // static_assert( std::floating_point<T>, "logspace() only supports
            // floating-point types" );
            if (ls == nullptr) {
                ls = NCPA::arrays::zeros<T>( N );
            }
            std::vector<T> tmpvec = logspace<T>( a, b, N );
            std::copy( tmpvec.begin(), tmpvec.end(), ls );
        }

        /**
        Computes the mean of a vector.
        @brief Computes the mean of a vector.
        @param d The vector to average.
        @returns The computed mean.
        */
        template<typename T>
        T mean( const std::vector<T>& d ) {
            T m  = 0;
            T Tn = (T)( d.size() );
            for (T val : d) {
                m += val / Tn;
            }
            return m;
        }

        /**
        For input v, returns the lowest positive integer value n such that
        2^n >= v.
        @brief Returns the next integer power of 2.
        @param v The input number to start at.
        @returns The next integer power of 2.
        */
        template<typename T>
        T nextpow2( T v ) {
            return (T)( std::ceil( std::log2( v ) ) );
        }

        /**
        @brief Peforms a simple trapezoidal integration.
        @param N The number of points to integrate over.
        @param xvec The x values.
        @param yvec The y values.
        @returns The computed integral.
        */
        template<typename T>
        T trapz( size_t N, T *xvec, T *yvec ) {
            assert( N > 1 );
            T sum = 0.0;
            for (size_t i = 1; i < N; i++) {
                sum += ( yvec[ i ] + yvec[ i - 1 ] ) * 0.5
                     * ( xvec[ i ] - xvec[ i - 1 ] );
            }
            return sum;
        }
    }  // namespace math
}  // namespace NCPA
