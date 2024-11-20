#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/types.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstring>
#include <functional>
#include <limits>
#include <random>
#include <type_traits>
#include <vector>

namespace NCPA {
    namespace math {

        // double precision to enable DOUBLE_EQUAL unit tests
        constexpr double PI
            = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;

        /**
         * Returns zero as the specified type.
         */
        template<typename T, ENABLE_IF( std::is_arithmetic<T> )>
        const T zero() {
            T z = 0;
            return z;
        }

        template<typename T, ENABLE_IF( NCPA::types::is_complex<T> )>
        const T zero() {
            T z( 0, 0 );
            return z;
        }

        /**
         * Returns unity as the specified type.
         *
         */
        template<typename T, ENABLE_IF( std::is_arithmetic<T> )>
        const T one() {
            T z = 1;
            return z;
        }

        template<typename T, ENABLE_IF( NCPA::types::is_complex<T> )>
        const T one() {
            T z( 1, 0 );
            return z;
        }

        template<typename T>
        bool equals( T x, T y, size_t n = 1, ENABLE_IF( std::is_arithmetic<T> ) ) {
            if ( std::numeric_limits<T>::is_exact ) {
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
        bool equals( T x, T y, size_t n = 1, ENABLE_IF( NCPA::types::is_complex<T> ) ) {
            return equals( x.real(), y.real() ) && equals( x.imag(), y.imag() );
        }

        /**
         * Returns -1, 0, or 1 depending on the sign of the argument.
         * @brief Returns -1, 0, or 1 depending on the sign of the argument.
         * @param n The value to test for sign
         * @returns -1 if n<0, 1 if n>0, 0 otherwise
         */
        template<typename T>
        int sign( T n ) {
            if ( n > 0 ) {
                return 1;
            } else if ( n < 0 ) {
                return -1;
            } else {
                return 0;
            }
        }

        /**
         * Finds and returns the indices of the array elements before and after
         * the supplied value.
         */
        template<typename T>
        bool find_interval_inclusive( T *z, size_t NZ, T val, size_t& bottom,
                                      size_t& top ) {
            double *it  = std::lower_bound( z, z + NZ, val );
            size_t diff = it - z;
            if ( diff == 0 ) {
                bottom = 0;
                if ( val < z[ 0 ] ) {
                    top = 0;
                    return false;
                } else {
                    top = 1;
                    return true;
                }
            } else if ( diff == NZ ) {
                top = NZ;
                if ( val > z[ NZ - 1 ] ) {
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
            if ( NZ == 1 ) {
                return 0;
            }
            double diff = 0.0, mindiff = std::numeric_limits<T>::max();
            size_t tmpind = 0;

            for ( size_t i = 0; i < NZ; i++ ) {
                diff = std::fabs( ( (double)z[ i ] ) - ( (double)zs ) );
                if ( diff < mindiff ) {
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
            if ( z.size() == 1 ) {
                return 0;
            }
            T diff = 0.0, mindiff = std::numeric_limits<T>::max();
            size_t tmpind = 0, i;
            for ( i = 0; i < z.size(); i++ ) {
                diff = std::abs( z[ i ] - zs );
                if ( diff < mindiff ) {
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
            } catch ( std::domain_error e ) {
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
        @brief Returns the maximum value from a vector.
        @param vals The vector holding the values to check.
        @returns The maximum value found in the vector.
        */
        template<typename T>
        T max( std::vector<T> vals ) {
            T maxval = vals.front();
            for ( auto cit = vals.cbegin(); cit != vals.cend(); ++cit ) {
                maxval = std::max( *cit, maxval );
            }
            return maxval;
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
            for ( size_t i = 1; i < size; i++ ) {
                maxval = std::max( vals[ i ], maxval );
            }
            return maxval;
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
            for ( size_t i = 1; i < size; i++ ) {
                minval = std::min( vals[ i ], minval );
            }
            return minval;
        }

        /**
        @brief Returns the maximum value from a vector.
        @param vals The vector holding the values to check.
        @returns The maximum value found in the vector.
        */
        template<typename T>
        T min( std::vector<T> vals ) {
            T minval = vals.front();
            for ( auto cit = vals.cbegin(); cit != vals.cend(); ++cit ) {
                minval = std::min( *cit, minval );
            }
            return minval;
        }

        /**
         * Converts a vector of reals into the equivalent vector of
         * complex<>s with the imaginary value set to zero.
         * @brief Converts a real vector to a complex<> vector
         * @param real The real vector to convert
         * @returns The complex<> vector to return
         */
        template<typename T = double, ENABLE_IF( std::is_floating_point<T> )>
        std::vector<std::complex<T>> real2complex( const std::vector<T>& in ) {
            std::vector<std::complex<T>> out;
            out.reserve( in.size() );
            for ( auto it = in.cbegin(); it != in.cend(); ++it ) {
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
        template<typename T = double, ENABLE_IF( std::is_floating_point<T> )>
        void real2complex( size_t n, const T *in, std::complex<T> *out ) {
            for ( size_t i = 0; i < n; i++ ) {
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
        template<typename T = double, ENABLE_IF( std::is_floating_point<T> )>
        std::vector<std::complex<T>> real2complex(
            const std::vector<T>& real, const std::vector<T>& imag ) {
            size_t Nr = real.size();
            size_t Ni = imag.size();
            size_t N  = std::max( Nr, Ni );
            std::vector<std::complex<T>> out( N );
            for ( size_t i = 0; i < N; i++ ) {
                if ( i < Nr ) {
                    out[ i ].real( real[ i ] );
                }
                if ( i < Ni ) {
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
        template<typename T = double, ENABLE_IF( std::is_floating_point<T> )>
        void real2complex( size_t n, const T *real, const T *imag,
                           std::complex<T> *out ) {
            for ( size_t i = 0; i < n; i++ ) {
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
        template<typename T = double, ENABLE_IF( std::is_floating_point<T> )>
        void complex2real( const std::vector<std::complex<T>>& in,
                           std::vector<T>& real, std::vector<T>& imag ) {
            real.resize( in.size() );
            imag.resize( in.size() );
            for ( size_t i = 0; i < in.size(); i++ ) {
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
        template<typename T = double, ENABLE_IF( std::is_floating_point<T> )>
        void complex2real( size_t n, const std::complex<T> *in, T *real,
                           T *imag ) {
            for ( size_t i = 0; i < n; i++ ) {
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
                                       ENABLE_IF( std::is_integral<T> ) ) {
            std::vector<T> randn;
            randn.reserve( N );
            std::random_device rd;
            std::mt19937 generator( rd() );
            std::uniform_int_distribution<T> distribution( minrange,
                                                           maxrange );
            for ( size_t i = 0; i < N; i++ ) {
                randn.push_back( distribution( generator ) );
            }
            return randn;
        }

        template<typename T>
        std::vector<T> random_numbers(
            size_t N, T minrange, T maxrange,
            ENABLE_IF( std::is_floating_point<T> ) ) {
            // static_assert( std::floating_point<T>, "This template is for
            // integer types" );
            std::vector<T> randn;
            randn.reserve( N );
            std::random_device rd;
            std::mt19937 generator( rd() );
            std::uniform_real_distribution<T> distribution( minrange,
                                                            maxrange );
            for ( size_t i = 0; i < N; i++ ) {
                randn.push_back( distribution( generator ) );
            }
            return randn;
        }

        template<typename T>
        std::vector<T> random_numbers(
            size_t N, typename T::value_type minrange,
            typename T::value_type maxrange,
            ENABLE_IF( NCPA::types::is_complex<T> ) ) {
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
            for ( size_t i = 0; i < coeffs.size(); i++ ) {
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
            for ( size_t i = 0; i < N; i++ ) {
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
            if ( vec == nullptr ) {
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
                                 ENABLE_IF( std::is_floating_point<T> ) ) {
            // static_assert( std::floating_point<T>, "logspace() only supports
            // floating-point types" );
            T la                = std::log10( firstval );
            T lb                = std::log10( lastval );
            std::vector<T> logs = NCPA::math::linspace( la, lb, N );
            for ( size_t i = 0; i < N; i++ ) {
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
                       ENABLE_IF( std::is_floating_point<T> ) ) {
            // static_assert( std::floating_point<T>, "logspace() only supports
            // floating-point types" );
            if ( ls == nullptr ) {
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
            for ( T val : d ) {
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
            for ( size_t i = 1; i < N; i++ ) {
                sum += ( yvec[ i ] + yvec[ i - 1 ] ) * 0.5
                     * ( xvec[ i ] - xvec[ i - 1 ] );
            }
            return sum;
        }

        /**
        @brief Performs element-wise array addition.
        @param v1 The first array to add.
        @param v2 The second array to add.
        @returns The array holding the sum.
        */
        template<typename T>
        T add_vectors( T& v1, T& v2,
                       ENABLE_IF( NCPA::types::is_iterable<T> ) ) {
            size_t N = std::min<size_t>( v1.size(), v2.size() );
            T v3     = v1.size() >= v2.size() ? v1 : v2;
            std::transform( v1.cbegin(), v1.cbegin() + N, v2.cbegin(),
                            v3.begin(), std::plus<typename T::value_type> {} );
            return v3;
        }

        /**
        @brief Performs element-wise array addition.
        @param N The number of points in the array.
        @param v1 The first array to add.
        @param v2 The second array to add.
        @param v12 The array to hold the sum.  Can be the same as either input
        array, in which case the values are replaced.
        */
        template<typename T>
        void add_arrays( size_t N, T *v1, T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for ( size_t i = 0; i < N; i++ ) {
                tempvec[ i ] = v1[ i ] + v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        /**
        Performs element-wise array division.  If arrays are different lengths,
        all values beyond the length of the shorter vector will be zero.
        @brief Performs element-wise array division.
        @param v1 The array to divide.
        @param v2 The array to divide by.
        @returns The array holding the quotient.
        */
        template<typename T>
        T divide_vectors( T& v1, T& v2,
                          ENABLE_IF( NCPA::types::is_iterable<T> ) ) {
            T v3 = T( std::max<size_t>( v1.size(), v2.size() ), 0.0 );
            std::transform( v1.cbegin(),
                            v1.cbegin()
                                + std::min<size_t>( v1.size(), v2.size() ),
                            v2.cbegin(), v3.begin(),
                            std::divides<typename T::value_type> {} );
            return v3;
        }

        /**
        Divides one array by another element-wise, returning the
        quotient in a supplied array.
        @brief Performs element-wise array division.
        @param N The number of points in the array.
        @param v1 The array to divide.
        @param v2 The array to divide by.
        @param v12 The array to hold the quotient.  Can be the same as either
        input array, in which case the values are replaced.
        */
        template<typename T>
        void divide_arrays( size_t N, T *v1, T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for ( size_t i = 0; i < N; i++ ) {
                tempvec[ i ] = v1[ i ] / v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        /**
        Performs element-wise array multiplication.  If arrays are different
        lengths, all values beyond the length of the shorter vector will be
        zero.
        @brief Performs element-wise array multiplication.
        @param v1 The first array to multiply.
        @param v2 The second array to multiply.
        @returns The array holding the product.
        */
        template<typename T>
        T multiply_vectors( T& v1, T& v2,
                            ENABLE_IF( NCPA::types::is_iterable<T> ) ) {
            T v3 = T( std::max<size_t>( v1.size(), v2.size() ), 0.0 );
            std::transform( v1.cbegin(),
                            v1.cbegin()
                                + std::min<size_t>( v1.size(), v2.size() ),
                            v2.cbegin(), v3.begin(),
                            std::multiplies<typename T::value_type> {} );
            return v3;
        }

        /**
        Multiplies two arrays together element-wise, returning the
        product in a supplied array.
        @brief Performs element-wise array multiplication.
        @param N The number of points in the array.
        @param v1 The first array to multiply.
        @param v2 The second array to multiply.
        @param v12 The array to hold the product.  Can be the same as either
        input array, in which case the values are replaced.
        */
        template<typename T>
        void multiply_arrays( size_t N, T *v1, T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for ( size_t i = 0; i < N; i++ ) {
                tempvec[ i ] = v1[ i ] * v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        /**
        Scales an array of values by a constant value, returning the
        result in a dynamically-allocated array.
        @brief Scales an array by a constant value.
        @param N The number of points in the array.
        @param in The array to scale.
        @param factor The factor to scale by.
        @param out The new, dynamically-allocated scaled array.
        */
        template<typename T, typename U>
        void scale_array( size_t N, U *in, T factor, U *& out ) {
            U *tempvec = NCPA::arrays::zeros<U>( N );
            for ( size_t i = 0; i < N; i++ ) {
                tempvec[ i ] = in[ i ] * (U)factor;
            }
            std::memcpy( out, tempvec, N * sizeof( U ) );
            delete[] tempvec;
        }

        /**
        Scales an array.
        @brief Performs array scaling.
        @param v1 The array to multiply.
        @param scalar The scalar to multiply by.
        @returns The scaled array.
        */
        template<typename T, typename U>
        T scale_vector(
            T& v1, U scalar,
            typename std::enable_if<NCPA::types::is_iterable_of<T, U>::value,
                                    int>::type ENABLER
            = 0 ) {
            T v3 = v1;
            std::transform( v3.begin(), v3.end(), v3.begin(),
                            [ scalar ]( U num ) { return num * scalar; } );
            return v3;
        }

        /**
        Scales an array of values in-place by a constant value
        @brief Scales an array by a constant value in place.
        @param N The number of points in the array.
        @param in The array to scale.
        @param factor The factor to scale by.
        */
        template<typename T, typename U>
        void scale_array( size_t N, U *in, T factor ) {
            for ( size_t i = 0; i < N; i++ ) {
                in[ i ] *= factor;
            }
        }

        /**
        @brief Performs element-wise array subtraction.
        @param v1 The first array to add.
        @param v2 The second array to add.
        @returns The array holding the sum.
        */
        template<typename T>
        T subtract_vectors(
            T& v1, T& v2,
            typename std::enable_if<NCPA::types::is_iterable<T>::value,
                                    int>::type ENABLER
            = 0 ) {
            T v3 = scale_vector( v2, -1.0 );
            return add_vectors( v1, v3 );
        }

        /**
        Subtracts one array from another element-wise, returning the
        difference in a supplied array.
        @brief Performs element-wise array subtraction.
        @param N The number of points in the array.
        @param v1 The array to subtract from.
        @param v2 The array to subtract.
        @param v12 The array to hold the difference.  Can be the same as either
        input array, in which case the values are replaced.
        */
        template<typename T>
        void subtract_arrays( size_t N, T *v1, T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for ( size_t i = 0; i < N; i++ ) {
                tempvec[ i ] = v1[ i ] - v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        template<typename T>
        bool is_zero( T val ) {
            return ( std::fpclassify( val ) == FP_ZERO );
        }

        template<typename T>
        bool is_zero( std::complex<T> val ) {
            return is_zero( val.real() ) && is_zero( val.imag() );
        }

        // /**
        // @brief Compares two doubles to a given number of decimal places
        // @param val1 First value
        // @param val2 Second value
        // @param precision Decimal places to compare
        // @returns true if the numbers are equal to that many decimal places,
        // false otherwise
        //  */
        // bool within( double val1, double val2, size_t precision );

        // /**
        // Shifts a point in Cartesian coordinates to a new position
        // relative to a different origin.
        // @brief Returns coordinates relative to a new origin.
        // @param old_x The x coordinate relative to (0,0).
        // @param old_y The y coordinate relative to (0,0).
        // @param x_new_origin The x coordinate of the new origin.
        // @param y_new_origin The y coordinate of the new origin.
        // @param new_x The new x coordinate.
        // @param new_y The new y coordiante.
        // */
        // template<typename T>
        // void move_origin( T old_x, T old_y, T x_new_origin, T y_new_origin,
        //                   T &new_x, T &new_y ) {
        //     new_x = old_x - x_new_origin;
        //     new_y = old_y - y_new_origin;
        // }

        // /**
        // Shifts an array of Cartesian points from an origin at (0,0) to
        // a new origin, and returns the coordinates relative to the new
        // origin.
        // @brief Returns coordinates relative to a new origin.
        // @param npts The number of points to move.
        // @param old_x The x coordinates relative to (0,0)
        // @param old_y The y coordinates relative to (0,0)
        // @param x_new_origin The x coordinate of the new origin
        // @param y_new_origin The y coordinate of the new origin
        // @param new_x The array to hold the new x coordinates.
        // @param new_y The array to hold the new y coordinates.
        // */
        // template<typename T>
        // void move_origin( size_t npts, const T *old_x, const T *old_y,
        //                   T x_new_origin, T y_new_origin, T *new_x, T *new_y
        //                   ) {
        //     for ( size_t i = 0; i < npts; i++ ) {
        //         move_origin<T>( old_x[ i ], old_y[ i ], x_new_origin,
        //         y_new_origin,
        //                         new_x[ i ], new_y[ i ] );
        //     }
        // }


    }  // namespace math
}  // namespace NCPA
