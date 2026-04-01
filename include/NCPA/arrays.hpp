#pragma once

#include "NCPA/constants.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/ndvector.hpp"
#include "NCPA/types.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <map>
#include <numeric>
#include <type_traits>
#include <vector>

// predeclare classes and typedefs
namespace NCPA {
    namespace arrays {
        template<typename T>
        class vector2d_t;

        template<typename T>
        class vector3d_t;

        template<typename T, size_t N>
        class _abstract_arraylike;

        template<typename T, size_t N>
        class ArrayLike;

        template<typename T, size_t N>
        class ArrayView;

        template<typename T, size_t N>
        class NDimensionalArray;

        template<typename T>
        class TwoDimensionalArray;

        template<typename T>
        class ThreeDimensionalArray;

        template<typename T>
        class TwoDimensionalArray;

        template<typename T, size_t N>
        using view_t = ArrayView<T, N - 1>;

        template<typename T, size_t N>
        using local_registry_t
            = std::map<std::pair<size_t, size_t>, view_t<T, N>>;

        template<typename T, size_t N>
        using global_registry_entry_t
            = std::map<std::pair<size_t, size_t>, view_t<T, N>>;

        template<typename T, size_t N>
        using global_registry_t
            = std::map<const ArrayLike<T, N> *, global_registry_entry_t<T, N>>;
    }  // namespace arrays
}  // namespace NCPA

// predeclare swap functions
template<typename T, size_t N>
static void swap( NCPA::arrays::_abstract_arraylike<T, N>& a,
                  NCPA::arrays::_abstract_arraylike<T, N>& b ) noexcept;

template<typename T, size_t N>
static void swap( NCPA::arrays::ArrayLike<T, N>& a,
                  NCPA::arrays::ArrayLike<T, N>& b ) noexcept;

template<typename T>
static void swap( NCPA::arrays::ArrayLike<T, 1>& a,
                  NCPA::arrays::ArrayLike<T, 1>& b ) noexcept;


template<typename T, size_t N>
static void swap( NCPA::arrays::ArrayView<T, N>& a,
                  NCPA::arrays::ArrayView<T, N>& b ) noexcept;

template<typename T>
static void swap( NCPA::arrays::ArrayView<T, 0>& a,
                  NCPA::arrays::ArrayView<T, 0>& b ) noexcept;

template<typename T, size_t N>
static void swap( NCPA::arrays::NDimensionalArray<T, N>& a,
                  NCPA::arrays::NDimensionalArray<T, N>& b ) noexcept;

template<typename T>
static void swap( NCPA::arrays::ThreeDimensionalArray<T>& a,
                  NCPA::arrays::ThreeDimensionalArray<T>& b ) noexcept;

template<typename T>
static void swap( NCPA::arrays::TwoDimensionalArray<T>& a,
                  NCPA::arrays::TwoDimensionalArray<T>& b ) noexcept;

// resizeable 2-d and 3-d vectors
namespace NCPA {
    namespace arrays {
        template<typename T>
        class vector2d_t : public ndvector<2, T> {
            public:
                using ndvector<2, T>::set;

                vector2d_t() : ndvector<2, T>() {}

                vector2d_t( size_t nx1, size_t nx2, const T& val = (T)0 ) :
                    ndvector<2, T>() {
                    this->resize2d( nx1, nx2, val );
                }

                vector2d_t( size_t nx1, size_t nx2, const T **vals ) :
                    vector2d_t( nx1, nx2 ) {
                    for (size_t i = 0; i < nx1; ++i) {
                        this->at( i ).assign( vals[ i ], vals[ i ] + nx2 );
                    }
                }

                virtual void resize2d( size_t nx1, size_t nx2,
                                       const T& val = (T)0 ) {
                    // this->resize( nx1, std::vector<T>( nx2, val ) );
                    this->reshape( { nx1, nx2 }, val );
                }

                virtual void size2d( size_t& nx1, size_t& nx2 ) const {
                    auto dims = this->shape();
                    nx1       = dims[ 0 ];
                    nx2       = dims[ 1 ];
                }

                virtual size_t dim( size_t dimnum ) const {
                    return this->shape()[ dimnum ];
                }

                virtual void fill( T val ) { this->set( val ); }
        };

        template<typename T>
        class vector3d_t : public ndvector<3, T> {
            public:
                vector3d_t() : ndvector<3, T>() {}

                vector3d_t( size_t nx1, size_t nx2, size_t nx3,
                            const T& val = (T)0.0 ) :
                    ndvector<3, T>() {
                    this->resize3d( nx1, nx2, nx3, val );
                }

                vector3d_t( size_t nx1, size_t nx2, size_t nx3,
                            const T ***vals ) :
                    vector3d_t( nx1, nx2, nx3 ) {
                    for (size_t i = 0; i < nx1; ++i) {
                        for (size_t j = 0; j < nx2; ++j) {
                            this->at( i ).at( j ).assign(
                                vals[ i ][ j ], vals[ i ][ j ] + nx3 );
                        }
                    }
                }

                vector3d_t( const vector3d_t<T>& other ) :
                    ndvector<3, T>( other ) {}

                vector3d_t( const ndvector<3, T>& other ) :
                    ndvector<3, T>( other ) {}

                virtual void resize3d( size_t nx1, size_t nx2, size_t nx3,
                                       const T& val = (T)0.0 ) {
                    this->reshape( { nx1, nx2, nx3 }, val );
                }

                virtual void size3d( size_t& nx1, size_t& nx2,
                                     size_t& nx3 ) const {
                    auto dims = this->shape();
                    nx1       = dims[ 0 ];
                    nx2       = dims[ 1 ];
                    nx3       = dims[ 2 ];
                }

                virtual size_t dim( size_t dimnum ) const {
                    return this->shape()[ dimnum ];
                }
        };

        template<typename T>
        T *as_array( std::vector<T>& in ) {
            return &in[ 0 ];
        }

        template<typename T>
        const T *as_array( const std::vector<T>& in ) {
            return &in[ 0 ];
        }

        /**
        Dynamically allocates a new array and sets all elements to zero
        before returning it.
        @brief Returns a new array of all zeros.
        @param n The size of the array.
        @returns A pointer to the newly-allocated array.
        */
        template<typename T>
        T *zeros( size_t n ) {
            T *out = new T[ n ]();
            return out;
        }

        /**
        @brief Dynamically allocates and returns a two-dimensional array.
        @param nr The first dimension of the array.
        @param nc The second dimension of the array.
        @returns A pointer to the newly-allocated array.
        */
        template<typename T>
        T **zeros( size_t nr, size_t nc ) {
            T **v;
            v = new T *[ nr ];
            for (size_t i = 0; i < nr; i++) {
                v[ i ] = new T[ nc ]();
            }
            return v;
        }

        /**
        @brief Dynamically allocates and returns a three-dimensional array.
        @param nr The first dimension of the array.
        @param nc The second dimension of the array.
        @param nd The third dimension of the array
        @returns A pointer to the newly-allocated array.
        */
        template<typename T>
        T ***zeros( size_t nr, size_t nc, size_t nd ) {
            T ***v;
            v = new T **[ nr ];
            for (size_t i = 0; i < nr; i++) {
                v[ i ] = new T *[ nc ];
                for (size_t j = 0; j < nc; j++) {
                    v[ i ][ j ] = new T[ nd ]();
                }
            }
            return v;
        }

        /**
        @brief Frees a dynamically allocated 1-D array.  Most provided
           for consistency of interface.
        @param v The array to free.
        @param nr The dimension of the array.
        */
        template<typename T>
        void free_array( T *& v, size_t nr,
                         ENABLE_FUNCTION_IF_NOT_DELETEABLE( T ) ) {
            delete[] v;
            v = nullptr;
        }

        template<typename T>
        void free_array( T *& v, size_t nr,
                         ENABLE_FUNCTION_IF_DELETEABLE( T ) ) {
            for (size_t i = 0; i < nr; i++) {
                delete v[ i ];
            }
            delete[] v;
            v = nullptr;
        }

        /**
        @brief Frees a dynamically allocated 2-D array.
        @param v The array to free.
        @param nr The first dimension of the array.
        @param nc The second dimension of the array.
        */
        template<typename T>
        void free_array( T **& v, size_t nr, size_t nc,
                         ENABLE_FUNCTION_IF_NOT_DELETEABLE( T ) ) {
            for (size_t i = 0; i < nr; i++) {
                delete[] v[ i ];
            }
            delete[] v;
            v = nullptr;
        }

        template<typename T>
        void free_array( T **& v, size_t nr, size_t nc,
                         ENABLE_FUNCTION_IF_DELETEABLE( T ) ) {
            for (size_t i = 0; i < nr; i++) {
                for (size_t j = 0; j < nc; j++) {
                    delete v[ i ][ j ];
                }
                delete[] v[ i ];
            }
            delete[] v;
            v = nullptr;
        }

        /**
        Frees a dynamically-allocated three-dimensional array, and calls
        delete on each element as it does so.
        @brief Frees a dynamically-allocated three-dimensional array and its
        contents.
        @param data The array to free.
        @param nd1 The first dimension of the array.
        @param nd2 The second dimension of the array.
        @param nd3 The third dimension of the array.
        */
        template<typename T>
        void free_array( T ***& data, size_t nd1, size_t nd2, size_t nd3,
                         ENABLE_FUNCTION_IF_NOT_DELETEABLE( T ) ) {
            size_t i, j, k;
            for (i = 0; i < nd1; ++i) {
                if (data[ i ] != NULL) {
                    for (j = 0; j < nd2; ++j) {
                        delete[] data[ i ][ j ];
                    }
                    delete[] data[ i ];
                }
            }
            delete[] data;
            data = nullptr;
        }

        template<typename T>
        void free_array( T ***& data, size_t nd1, size_t nd2, size_t nd3,
                         ENABLE_FUNCTION_IF_DELETEABLE( T ) ) {
            size_t i, j, k;
            for (i = 0; i < nd1; ++i) {
                if (data[ i ] != NULL) {
                    for (j = 0; j < nd2; ++j) {
                        for (k = 0; i < nd3; k++) {
                            delete data[ i ][ j ][ k ];
                        }
                        delete[] data[ i ][ j ];
                    }
                    delete[] data[ i ];
                }
            }
            delete[] data;
            data = nullptr;
        }

        /**
        Circularly shifts the elements in an array X by K positions. If K is
        positive, then the values of X are circularly shifted from the
        beginning to the end. If K is negative, they are shifted from the
        end to the beginning.  The resultant shifted array is returned in
        a new dynamically-allocated array.
        @brief Circularly shifts array elements.
        @param X The array whose elements are to be shifted.
        @param N The size of the array.
        @param K The number of positions to shift.
        @param out The shifted array.  Can be the same as either input to
        perform in-place.
        */
        template<typename T>
        void circshift( T *X, size_t N, int K, T *& out ) {
            while (std::abs( K ) > N) {
                K -= ( (int)N ) * ( K < 0 ? -1 : 1 );
            }
            size_t i;
            T *tempvec = zeros<T>( N );
            if (K < 0) {
                // move from the back to the front
                size_t negshift = (size_t)( -K );
                // first, move the last -K values to the front
                for (i = 0; i < negshift; i++) {
                    tempvec[ i ] = X[ N - negshift + i ];
                }
                for (i = negshift; i < N; i++) {
                    tempvec[ i ] = X[ i - negshift ];
                }
            } else if (K > 0) {
                // move from the front to the back
                // first, move the first K values to the back
                for (i = 0; i < K; i++) {
                    tempvec[ N - K + i ] = X[ i ];
                }
                // Now, move the next N-K values to the front
                for (i = K; i < N; i++) {
                    tempvec[ i - K ] = X[ i ];
                }
            } else {
                std::memcpy( tempvec, X, N * sizeof( T ) );
            }

            // copy over so it can be done in place if desired
            if (out == nullptr) {
                out = zeros<T>( N );
            }
            std::memcpy( out, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        /**
        @brief Sets each element of an array to a constant value.
        @param v The array to fill.
        @param nr The dimension of the array.
        @param val The value to set each element to.
        */
        template<typename T>
        void fill( T *v, size_t nr, const T& val ) {
            std::fill( v, v + nr, val );
        }

        /**
        @brief Sets each element of a 2-D array to a constant value.
        @param v The array to fill.
        @param nr The first dimension of the array.
        @param nc The second dimension of the array.
        @param val The value to set each element to.
        */
        template<typename T>
        void fill( T **v, size_t nr, size_t nc, const T& val ) {
            for (size_t i = 0; i < nr; i++) {
                fill( v[ i ], nc, val );
            }
        }

        /**
        @brief Sets each element of a three-dimensional array to a constant
        value.
        @param v The array to fill.
        @param nr The first dimension of the array.
        @param nc The second dimension of the array.
        @param nz The third dimension of the array.
        @param val The value to set each element to.
        */
        template<typename T>
        void fill( T ***v, size_t nr, size_t nc, size_t nz, const T& val ) {
            for (size_t i = 0; i < nr; i++) {
                fill( v[ i ], nc, nz, val );
            }
        }

        /**
        @brief Copies one array to another
        @param from The array to copy.
        @param n The number of values to copy.
        @param to The array to copy into.
        */
        template<typename T>
        void copy( const T *from, size_t n, T *to ) {
            if (from != nullptr) {
                std::copy( from, from + n, to );
            } else if (to != nullptr) {
                free_array( to, n );
            }
        }

        /**
        @brief Copies one 2-D array to another
        @param from The array to copy.
        @param n1 The first dimension of the array.
        @param n2 The second dimension of the array.
        @param to The array to copy into.
        */
        template<typename T>
        void copy( const T **from, size_t n1, size_t n2, T **to ) {
            if (from != nullptr) {
                for (size_t i = 0; i < n1; i++) {
                    copy( from[ i ], n2, to[ i ] );
                }
            } else if (to != nullptr) {
                free_array( to, n1, n2 );
            }
        }

        /**
        @brief Copies one 3-D array to another
        @param from The array to copy.
        @param n1 The first dimension of the array.
        @param n2 The second dimension of the array.
        @param n3 The third dimension of the array
        @param to The array to copy into.
        */
        template<typename T>
        void copy( const T ***from, size_t n1, size_t n2, size_t n3,
                   T ***to ) {
            if (from != nullptr) {
                for (size_t i = 0; i < n1; i++) {
                    copy( from[ i ], n2, n3, to[ i ] );
                }
            } else if (to != nullptr) {
                free_array( to, n1, n2, n3 );
            }
        }

        /**
        Dynamically allocates a new array and sets the elements to
        (a, a+1, a+2, ..., a+n-1).
        @brief Returns a new array of index values.
        @param n The size of the array.
        @param a The offset from zero to use for each value.  Defaults to 0.
        @returns A pointer to the newly-allocated array.
        */
        template<typename T>
        std::vector<T> index_vector( size_t n, T a = 0 ) {
            std::vector<T> inds( n );
            for (size_t i = 0; i < n; i++) {
                inds[ i ] = (T)i + a;
            }
            return inds;
        }

        /**
        @brief Reverses the order of the first N elements of an array.
        @param in The input array.
        @param N The number of elements of the array to reverse.
        @param out Holds the output array. Can be the same as the input.
        */
        template<typename T>
        void reverse( T *in, size_t N, T *& out ) {
            T *tempvec = zeros<T>( N );
            size_t j;
            for (j = 0; j < N; j++) {
                tempvec[ N - 1 - j ] = in[ j ];
            }
            std::memcpy( out, tempvec, N * sizeof( T ) );

            delete[] tempvec;
        }

        // use like:
        // auto p = sort_permutation(vectorA,
        //      [](T const& a, T const& b){ /*some comparison*/ });
        //
        // From https://stackoverflow.com/a/17074810
        // @todo generalize for any container
        template<typename T, typename Compare>
        std::vector<std::size_t> sort_permutation( const std::vector<T>& vec,
                                                   Compare compare ) {
            std::vector<std::size_t> p( vec.size() );
            std::iota( p.begin(), p.end(), 0 );
            std::sort( p.begin(), p.end(),
                       [ & ]( std::size_t i, std::size_t j ) {
                           return compare( vec[ i ], vec[ j ] );
                       } );
            return p;
        }

        template<typename T, typename Compare>
        std::vector<std::size_t> sort_permutation( const std::vector<T>& vec1,
                                                   const std::vector<T>& vec2,
                                                   Compare compare ) {
            if (vec1.size() != vec2.size()) {
                throw std::invalid_argument(
                    "Vectors must be the same size!" );
            }
            std::vector<std::size_t> p( vec1.size() );
            std::iota( p.begin(), p.end(), 0 );
            std::sort( p.begin(), p.end(),
                       [ & ]( std::size_t i, std::size_t j ) {
                           return compare( vec1[ i ], vec1[ j ] )
                               || ( vec1[ i ] == vec1[ j ]
                                    && compare( vec2[ i ], vec2[ j ] ) );
                       } );
            return p;
        }

        template<typename T>
        std::vector<std::size_t> sort_permutation_increasing(
            const std::vector<T>& vec ) {
            return sort_permutation<T>(
                vec, []( T const& a, T const& b ) { return a < b; } );
        }

        template<typename T>
        std::vector<std::size_t> sort_permutation_increasing(
            const std::vector<T>& vec1, const std::vector<T>& vec2 ) {
            return sort_permutation<T>(
                vec1, vec2, []( T const& a, T const& b ) { return a < b; } );
        }

        template<typename T>
        std::vector<std::size_t> sort_permutation_decreasing(
            const std::vector<T>& vec ) {
            return sort_permutation<T>(
                vec, []( T const& a, T const& b ) { return a > b; } );
        }

        template<typename T>
        std::vector<std::size_t> sort_permutation_decreasing(
            const std::vector<T>& vec1, const std::vector<T>& vec2 ) {
            return sort_permutation<T>(
                vec1, vec2, []( T const& a, T const& b ) { return a > b; } );
        }

        // From https://stackoverflow.com/a/17074810
        template<typename T>
        std::vector<T> apply_permutation( const std::vector<T>& vec,
                                          const std::vector<std::size_t>& p ) {
            std::vector<T> sorted_vec( vec.size() );
            std::transform( p.begin(), p.end(), sorted_vec.begin(),
                            [ & ]( std::size_t i ) { return vec[ i ]; } );
            return sorted_vec;
        }

        // From https://stackoverflow.com/a/17074810
        template<typename T>
        void apply_permutation_in_place( std::vector<T>& vec,
                                         const std::vector<std::size_t>& p ) {
            std::vector<bool> done( vec.size() );
            for (std::size_t i = 0; i < vec.size(); ++i) {
                if (done[ i ]) {
                    continue;
                }
                done[ i ]          = true;
                std::size_t prev_j = i;
                std::size_t j      = p[ i ];
                while (i != j) {
                    std::swap( vec[ prev_j ], vec[ j ] );
                    done[ j ] = true;
                    prev_j    = j;
                    j         = p[ j ];
                }
            }
        }

        template<typename T>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::string& separator = ", " ) {
            for (size_t i = 0; i < basevec.size(); ++i) {
                if (i > 0) {
                    os << separator;
                }
                os << basevec[ i ];
            }
            // os << std::endl;
        }

        template<typename T, typename U>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::vector<std::vector<U>>& items,
                    const std::string& separator = " " ) {
            for (size_t i = 0; i < basevec.size(); i++) {
                os << basevec[ i ];
                for (size_t j = 0; j < items.size(); j++) {
                    os << separator << items[ j ][ i ];
                }
                os << std::endl;
            }
        }

        template<typename T, typename U>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::vector<U> depvec,
                    const std::string& separator = " ",
                    ENABLE_FUNCTION_IF_ARITHMETIC( U ) ) {
            for (size_t i = 0; i < basevec.size(); i++) {
                os << basevec[ i ] << separator << depvec[ i ] << std::endl;
            }
        }

        template<typename T, typename U>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::vector<U> depvec,
                    const std::string& separator = " ",
                    ENABLE_FUNCTION_IF_COMPLEX( U ) ) {
            for (size_t i = 0; i < basevec.size(); i++) {
                os << basevec[ i ] << separator << depvec[ i ].real()
                   << separator << depvec[ i ].imag() << std::endl;
            }
        }

        /**
        @brief Performs element-wise vector addition.
        @param v1 The first vector to add.
        @param v2 The second vector to add.
        @returns The vector holding the sum.
        */
        template<typename T>
        T add_vectors( const T& v1, const T& v2,
                       ENABLE_FUNCTION_IF_ITERABLE( T ) ) {
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
        void add_arrays( size_t N, const T *v1, const T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for (size_t i = 0; i < N; i++) {
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
        T divide_vectors( const T& v1, const T& v2,
                          ENABLE_FUNCTION_IF_ITERABLE( T ) ) {
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
        void divide_arrays( size_t N, const T *v1, const T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for (size_t i = 0; i < N; i++) {
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
        T multiply_vectors( const T& v1, const T& v2,
                            ENABLE_FUNCTION_IF_ITERABLE( T ) ) {
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
        void multiply_arrays( size_t N, const T *v1, const T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for (size_t i = 0; i < N; i++) {
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
        void scale_array( size_t N, const U *in, T factor, U *& out ) {
            U *tempvec = NCPA::arrays::zeros<U>( N );
            for (size_t i = 0; i < N; i++) {
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
            const T& v1, U scalar,
            typename std::enable_if<NCPA::types::is_iterable_of<T, U>::value,
                                    int>::type ENABLER
            = 0 ) {
            if (scalar == NCPA::constants::one<U>()) {
                return v1;
            }
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
            if (factor != NCPA::constants::one<T>()) {
                for (size_t i = 0; i < N; i++) {
                    in[ i ] *= factor;
                }
            }
        }

        /**
        Offsets an array.
        @brief Performs array offset.
        @param v1 The array to add to.
        @param scalar The scalar to add.
        @returns The scaled array.
        */
        template<typename T, typename U>
        T offset_vector(
            const T& v1, U scalar,
            typename std::enable_if<NCPA::types::is_iterable_of<T, U>::value,
                                    int>::type ENABLER
            = 0 ) {
            if (NCPA::constants::is_zero<U>( scalar )) {
                return v1;
            }
            T v3 = v1;
            for (auto i = 0; i < v1.size(); i++) {
                v3[ i ] += scalar;
            }
            return v3;
        }

        /**
        Offsets an array of values in-place by a constant value
        @brief Offsets an array by a constant value in place.
        @param N The number of points in the array.
        @param in The array to offset.
        @param factor The factor to offset by.
        */
        template<typename T, typename U>
        void offset_array( size_t N, U *in, T factor ) {
            if (!NCPA::constants::is_zero<T>( factor )) {
                for (size_t i = 0; i < N; i++) {
                    in[ i ] += factor;
                }
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
            const T& v1, const T& v2,
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
        void subtract_arrays( size_t N, const T *v1, const T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for (size_t i = 0; i < N; i++) {
                tempvec[ i ] = v1[ i ] - v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        template<typename FROMTYPE,typename TOTYPE>
        std::vector<TOTYPE> cast_vector( const std::vector<FROMTYPE>& vec ) {
            std::vector<TOTYPE> newvec( vec.size() );
            for (size_t i = 0; i < vec.size(); ++i) {
                newvec[i] = static_cast<TOTYPE>( vec.at(i) );
            }
            return newvec;
        }
    }  // namespace arrays
}  // namespace NCPA

// operator overloads for above functions
template<typename T>
std::vector<T> operator+( const std::vector<T>& a, const std::vector<T>& b ) {
    return NCPA::arrays::add_vectors<std::vector<T>>( a, b );
}

template<typename T>
std::vector<T> operator+( const std::vector<T>& a, const T& b ) {
    return NCPA::arrays::offset_vector<std::vector<T>>( a, b );
}

template<typename T>
std::vector<T> operator-( const std::vector<T>& a, const std::vector<T>& b ) {
    return NCPA::arrays::subtract_vectors<std::vector<T>>( a, b );
}

template<typename T>
std::vector<T> operator-( const std::vector<T>& a, const T& b ) {
    return NCPA::arrays::offset_vector<std::vector<T>>( a, -b );
}

template<typename T>
std::vector<T> operator*( const std::vector<T>& a, const std::vector<T>& b ) {
    return NCPA::arrays::multiply_vectors<std::vector<T>>( a, b );
}

template<typename T>
std::vector<T> operator*( const std::vector<T>& a, const T& b ) {
    return NCPA::arrays::scale_vector<std::vector<T>>( a, b );
}

template<typename T>
std::vector<T> operator/( const std::vector<T>& a, const T& b ) {
    return NCPA::arrays::scale_vector<std::vector<T>>( a,
                                                       std::pow( b, -1.0 ) );
}

template<typename T>
std::vector<T> operator-( const std::vector<T>& a ) {
    return NCPA::arrays::scale_vector<std::vector<T>>(
        a, -NCPA::constants::one<T>() );
}

namespace NCPA {
    namespace arrays {

        template<typename T, size_t N>
        class _abstract_arraylike {
            public:
                // virtual T& at( const std::array<size_t, N>& coords ) = 0;
                // virtual const T& at(
                //     const std::array<size_t, N>& coords ) const
                //     = 0;

                virtual ~_abstract_arraylike() {}

                friend void ::swap<>( _abstract_arraylike<T, N>& a,
                                      _abstract_arraylike<T, N>& b ) noexcept;

                virtual size_t buffer( T *b ) const            = 0;
                virtual size_t debuffer( const T *b )          = 0;
                virtual size_t dimension( size_t d ) const     = 0;
                virtual std::istream& read( std::istream& is ) = 0;
                virtual void redimension( const std::array<size_t, N>& dims )
                    = 0;
                virtual size_t size() const                             = 0;
                virtual view_t<T, N>& view( size_t dim, size_t dimind ) = 0;
                virtual const view_t<T, N>& view( size_t dim,
                                                  size_t dimind ) const
                    = 0;
                virtual std::ostream& write( std::ostream& os ) const = 0;

                virtual std::array<size_t, N>& dimensions() = 0;

                virtual const std::array<size_t, N>& dimensions() const = 0;

                virtual view_t<T, N>& operator[]( size_t ind ) {
                    return this->view( 0, ind );
                }

                virtual const view_t<T, N>& operator[]( size_t ind ) const {
                    return this->view( 0, ind );
                }

            protected:
                // std::array<size_t, N> _dimensions;

                virtual void check_dimensions( size_t dim ) const {
                    if (dim >= N) {
                        std::ostringstream oss;
                        oss << "ArrayLike: Requested dimension " << dim
                            << " too large for " << N << "-dimensional array";
                        throw std::range_error( oss.str() );
                    }
                }

                virtual void check_dimension_size( size_t dim,
                                                   size_t dimind ) const {
                    if (dimind >= this->dimensions().at( dim )) {
                        std::ostringstream oss;
                        oss << "ArrayLike: Requested element " << dimind
                            << " too large for " << dim << "th dimension of "
                            << this->dimensions_string() << " array";
                        throw std::range_error( oss.str() );
                    }
                }

                virtual std::string dimensions_string() const {
                    std::ostringstream dimstr;
                    dimstr << this->dimension( 0 );
                    for (size_t i = 1; i < N; ++i) {
                        dimstr << "x" << this->dimension( i );
                    }
                    return dimstr.str();
                }
        };

        /*
        Generalized interface for an N-D array that can be pointed back to a
        1-D array or vector
        */
        template<typename T, size_t N>
        class ArrayLike : public _abstract_arraylike<T, N> {
                static_assert( N > 0, "Dimension must be nonzero." );

            public:
                ArrayLike() {}

                ArrayLike( const std::array<size_t, N>& dims ) :
                    _dimensions { dims } {}

                virtual ~ArrayLike() { deregister_parent( this ); }

                ArrayLike( const ArrayLike<T, N>& other ) : ArrayLike<T, N>() {
                    _dimensions  = other._dimensions;
                    _local_views = other._local_views;
                }

                ArrayLike( ArrayLike<T, N>&& other ) noexcept :
                    ArrayLike<T, N>() {
                    ::swap( *this, other );
                }

                ArrayLike<T, N>& operator=( ArrayLike<T, N> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ArrayLike<T, N>& a,
                                      ArrayLike<T, N>& b ) noexcept;

                virtual T& at( const std::array<size_t, N>& coords ) = 0;
                virtual const T& at(
                    const std::array<size_t, N>& coords ) const
                    = 0;

                virtual size_t buffer( T *b ) const override {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to buffer()" );
                    }
                    size_t counter = 0;
                    for (size_t i = 0; i < this->_dimensions.at( 0 ); ++i) {
                        counter += this->view( 0, i ).buffer( b + counter );
                    }
                    return counter;
                }

                template<typename BUFFERTYPE>
                size_t buffer_as( BUFFERTYPE *b ) const {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to buffer()" );
                    }
                    size_t counter = 0;
                    for (size_t i = 0; i < this->_dimensions.at( 0 ); ++i) {
                        counter += this->view( 0, i ).buffer_as( b + counter );
                    }
                    return counter;
                }

                virtual size_t debuffer( const T *b ) override {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to debuffer()" );
                    }
                    size_t counter = 0;
                    for (size_t i = 0; i < this->_dimensions.at( 0 ); ++i) {
                        counter += this->view( 0, i ).debuffer( b + counter );
                    }
                    return counter;
                }

                template<typename BUFFERTYPE>
                size_t debuffer_as( const BUFFERTYPE *b ) {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to debuffer()" );
                    }
                    size_t counter = 0;
                    for (size_t i = 0; i < this->_dimensions.at( 0 ); ++i) {
                        // auto v   = this->view( 0, i );
                        // counter += v.debuffer_as( b + counter );
                        counter
                            += this->view( 0, i ).debuffer_as( b + counter );
                    }
                    return counter;
                }

                virtual size_t dimension( size_t d ) const override {
                    this->check_dimensions( d );
                    return this->dimensions().at( d );
                }

                // virtual size_t& dimension( size_t d ) {
                //     return ( d >= this->dimensions.size()
                //                  ? 1
                //                  : this->dimensions.at( d ) );
                // }

                // virtual const size_t dimensions() const { return N; }

                virtual std::istream& read( std::istream& is ) override {
                    T *b = new T[ this->size() ];
                    is.read( reinterpret_cast<char *>( b ),
                             this->size() * sizeof( T ) );
                    if (is.gcount() != this->size() * sizeof( T )) {
                        throw std::runtime_error(
                            "ArrayLike.read(): Not enough samples read!" );
                    }
                    size_t elements = this->debuffer( b );
                    if (elements != this->size()) {
                        throw std::runtime_error(
                            "ArrayLike.read(): Mismatch between internal "
                            "array size and size of buffer input" );
                    }
                    delete[] b;
                    return is;
                }

                template<typename BUFFERTYPE>
                std::istream& read_as( std::istream& is ) {
                    BUFFERTYPE *b = new BUFFERTYPE[ this->size() ];
                    is.read( reinterpret_cast<char *>( b ),
                             this->size() * sizeof( BUFFERTYPE ) );
                    if (is.gcount() != this->size() * sizeof( T )) {
                        throw std::runtime_error(
                            "ArrayLike.read(): Not enough samples read!" );
                    }
                    size_t elements = this->debuffer_as<BUFFERTYPE>( b );
                    if (elements != this->size()) {
                        throw std::runtime_error(
                            "ArrayLike.read(): Mismatch between internal "
                            "array size and size of buffer input" );
                    }
                    delete[] b;
                    return is;
                }

                virtual void redimension(
                    const std::array<size_t, N>& dims ) override {
                    this->_dimensions = dims;
                }

                virtual size_t size() const override {
                    size_t count = 1;
                    for (auto it = this->_dimensions.begin();
                         it != this->_dimensions.end(); ++it) {
                        count *= ( *it );
                    }
                    return count;
                }

                virtual view_t<T, N>& view( size_t dim,
                                            size_t dimind ) override {
                    this->check_dimensions( dim );
                    this->check_dimension_size( dim, dimind );
                    std::pair<size_t, size_t> key { dim, dimind };
                    auto keyview = _local_views.find( key );
                    if (keyview == _local_views.end()) {
                        _local_views[ key ]
                            = view_t<T, N>( this, dim, dimind );
                        return _local_views[ key ];
                    } else {
                        return keyview->second;
                    }
                }

                virtual const view_t<T, N>& view(
                    size_t dim, size_t dimind ) const override {
                    this->check_dimensions( dim );
                    this->check_dimension_size( dim, dimind );
                    std::pair<size_t, size_t> key { dim, dimind };
                    auto global_registry = _global_views.find( this );
                    if (global_registry == _global_views.end()) {
                        global_registry
                            = _global_views
                                  .emplace( this,
                                            global_registry_entry_t<T, N> {} )
                                  .first;
                    }
                    global_registry_entry_t<T, N>& entry
                        = global_registry->second;
                    entry[ key ] = view_t<T, N>( this, dim, dimind );
                    return entry[ key ];
                }

                virtual std::ostream& write(
                    std::ostream& os ) const override {
                    T *b            = new T[ this->size() ];
                    size_t elements = this->buffer( b );
                    if (elements != this->size()) {
                        throw std::runtime_error(
                            "ArrayLike.write(): Mismatch between internal "
                            "array size and size of buffer output" );
                    }
                    os.write( reinterpret_cast<const char *>( b ),
                              elements * sizeof( T ) );
                    delete[] b;
                    return os;
                }

                template<typename OUTTYPE>
                std::ostream& write_as( std::ostream& os ) const {
                    OUTTYPE *b      = new OUTTYPE[ this->size() ];
                    size_t elements = this->buffer_as<OUTTYPE>( b );
                    if (elements != this->size()) {
                        throw std::runtime_error(
                            "ArrayLike.write_as(): Mismatch between internal "
                            "array size and size of buffer output" );
                    }
                    os.write( reinterpret_cast<const char *>( b ),
                              elements * sizeof( OUTTYPE ) );
                    delete[] b;
                    return os;
                }

                static bool global_registry_contains( ArrayLike<T, N> *key );

                static global_registry_entry_t<T, N>& register_parent(
                    ArrayLike<T, N> *key );

                static void deregister_parent( ArrayLike<T, N> *key );

                static void register_view( ArrayLike<T, N> *key,
                                           const ArrayView<T, N - 1>& v );

                virtual std::array<size_t, N>& dimensions() override {
                    return _dimensions;
                }

                virtual const std::array<size_t, N>& dimensions()
                    const override {
                    return _dimensions;
                }


            private:
                std::array<size_t, N> _dimensions;
                local_registry_t<T, N> _local_views;

                static global_registry_t<T, N> _global_views;
        };

        template<typename T, size_t N>
        global_registry_t<T, N> ArrayLike<T, N>::_global_views;

        template<typename T, size_t N>
        void ArrayLike<T, N>::deregister_parent( ArrayLike<T, N> *key ) {
            if (!ArrayLike<T, N>::global_registry_contains( key )) {
                ArrayLike<T, N>::_global_views.erase( key );
            }
        }

        template<typename T, size_t N>
        void ArrayLike<T, N>::register_view( ArrayLike<T, N> *key,
                                             const ArrayView<T, N - 1>& v ) {
            ArrayLike<T, N>::register_parent( key );
            ArrayLike<T, N>::_global_views[ key ].emplace( v.key(), v );
        }

        template<typename T, size_t N>
        global_registry_entry_t<T, N>& ArrayLike<T, N>::register_parent(
            ArrayLike<T, N> *key ) {
            if (!ArrayLike<T, N>::global_registry_contains( key )) {
                return ArrayLike<T, N>::_global_views
                    .emplace( { key, global_registry_entry_t<T, N> {} } )
                    ->first;
            } else {
                return ArrayLike<T, N>::_global_views.find( key )->second;
            }
        }

        template<typename T, size_t N>
        bool ArrayLike<T, N>::global_registry_contains(
            ArrayLike<T, N> *key ) {
            return ( ArrayLike<T, N>::_global_views.find( key )
                     == ArrayLike<T, N>::_global_views.end() );
        }

        template<typename T>
        class ArrayLike<T, 1> : public _abstract_arraylike<T, 1> {
            public:
                ArrayLike() {}

                ArrayLike( const std::array<size_t, 1>& dims ) :
                    _dimensions { dims } {}

                ArrayLike( size_t newdim ) :
                    _abstract_arraylike<T, 1>(
                        std::array<size_t, 1>( { newdim } ) ) {}

                virtual ~ArrayLike() {}

                friend void ::swap<>( ArrayLike<T, 1>& a,
                                      ArrayLike<T, 1>& b ) noexcept;

                ArrayLike( const ArrayLike<T, 1>& other ) :
                    ArrayLike<T, 1>( other._dimensions ) {
                    this->_dimensions = other._dimensions;
                    _local_views      = other._local_views;
                }

                ArrayLike( ArrayLike<T, 1>&& other ) noexcept :
                    ArrayLike<T, 1>() {
                    ::swap( *this, other );
                }

                ArrayLike<T, 1>& operator=( ArrayLike<T, 1> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual T& at( const std::array<size_t, 1>& coords ) = 0;
                virtual const T& at(
                    const std::array<size_t, 1>& coords ) const
                    = 0;

                virtual T& at( size_t coord ) {
                    return this->at( std::array<size_t, 1> { coord } );
                }

                virtual const T& at( size_t coord ) const {
                    return this->at( std::array<size_t, 1> { coord } );
                }

                virtual size_t buffer( T *b ) const override {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to buffer()" );
                    }
                    size_t counter = 0;
                    for (size_t i = 0; i < this->_dimensions.at( 0 ); ++i) {
                        b[ i ] = this->at( i );
                    }
                    return this->_dimensions.at( 0 );
                }

                template<typename BUFFERTYPE>
                size_t buffer_as( BUFFERTYPE *b ) const {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to buffer()" );
                    }
                    size_t counter = 0;
                    for (size_t i = 0; i < this->_dimensions.at( 0 ); ++i) {
                        b[ i ] = (BUFFERTYPE)( this->at( i ) );
                    }
                    return this->_dimensions.at( 0 );
                }

                virtual size_t debuffer( const T *b ) override {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to debuffer()" );
                    }
                    for (size_t i = 0; i < this->size(); ++i) {
                        this->at( i ) = b[ i ];
                    }
                    return this->size();
                }

                template<typename BUFFERTYPE>
                size_t debuffer_as( const BUFFERTYPE *b ) {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to debuffer()" );
                    }
                    for (size_t i = 0; i < this->size(); ++i) {
                        this->at( i ) = (T)( b[ i ] );
                    }
                    return this->size();
                }

                virtual size_t dimension( size_t d ) const override {
                    this->check_dimensions( d );
                    return this->dimensions().at( d );
                }

                virtual std::array<size_t, 1>& dimensions() override {
                    return _dimensions;
                }

                virtual const std::array<size_t, 1>& dimensions()
                    const override {
                    return _dimensions;
                }

                // virtual const size_t dimensions() const { return 1; }

                virtual std::istream& read( std::istream& is ) override {
                    T *b = new T[ this->size() ];
                    is.read( reinterpret_cast<char *>( b ),
                             this->size() * sizeof( T ) );
                    if (is.gcount() != this->size() * sizeof( T )) {
                        throw std::runtime_error(
                            "ArrayLike.read(): Not enough samples read!" );
                    }
                    size_t elements = this->debuffer( b );
                    if (elements != this->size()) {
                        throw std::runtime_error(
                            "ArrayLike.read(): Mismatch between internal "
                            "array size and size of buffer input" );
                    }
                    delete[] b;
                    return is;
                }

                template<typename BUFFERTYPE>
                std::istream& read_as( std::istream& is ) {
                    BUFFERTYPE *b = new BUFFERTYPE[ this->size() ];
                    is.read( reinterpret_cast<char *>( b ),
                             this->size() * sizeof( BUFFERTYPE ) );
                    if (is.gcount() != this->size() * sizeof( T )) {
                        throw std::runtime_error(
                            "ArrayLike.read(): Not enough samples read!" );
                    }
                    size_t elements = this->debuffer_as<BUFFERTYPE>( b );
                    if (elements != this->size()) {
                        throw std::runtime_error(
                            "ArrayLike.read(): Mismatch between internal "
                            "array size and size of buffer input" );
                    }
                    delete[] b;
                    return is;
                }

                virtual void redimension(
                    const std::array<size_t, 1>& dims ) override {
                    this->_dimensions = dims;
                }

                virtual void redimension( size_t dim ) {
                    this->_dimensions.at( 0 ) = dim;
                }

                virtual size_t size() const override {
                    return this->_dimensions.at( 0 );
                }

                virtual view_t<T, 1>& view( size_t dim,
                                            size_t dimind ) override {
                    this->check_dimensions( dim );
                    this->check_dimension_size( dim, dimind );
                    if (dimind >= this->dimension( dim )) {}
                    std::pair<size_t, size_t> key { dim, dimind };
                    auto keyview = _local_views.find( key );
                    if (keyview == _local_views.end()) {
                        _local_views[ key ]
                            = view_t<T, 1>( this, dim, dimind );
                        return _local_views[ key ];
                    } else {
                        return keyview->second;
                    }
                }

                virtual const view_t<T, 1>& view(
                    size_t dim, size_t dimind ) const override {
                    this->check_dimensions( dim );
                    this->check_dimension_size( dim, dimind );
                    std::pair<size_t, size_t> key { dim, dimind };
                    auto global_registry = _global_views.find( this );
                    if (global_registry == _global_views.end()) {
                        global_registry
                            = _global_views
                                  .emplace( this,
                                            global_registry_entry_t<T, 1>() )
                                  .first;
                    }

                    auto view_it = global_registry->second.find( key );
                    if (view_it == global_registry->second.end()) {
                        view_it = global_registry->second
                                      .emplace( key, view_t<T, 1>( this, dim,
                                                                   dimind ) )
                                      .first;
                    }
                    return view_it->second;
                }

                virtual std::ostream& write(
                    std::ostream& os ) const override {
                    T *b            = new T[ this->size() ];
                    size_t elements = this->buffer( b );
                    if (elements != this->size()) {
                        throw std::runtime_error(
                            "ArrayLike.write(): Mismatch between internal "
                            "array size and size of buffer output" );
                    }
                    os.write( reinterpret_cast<const char *>( b ),
                              elements * sizeof( T ) );
                    delete[] b;
                    return os;
                }

                template<typename OUTTYPE>
                std::ostream& write_as( std::ostream& os ) const {
                    OUTTYPE *b      = new OUTTYPE[ this->size() ];
                    size_t elements = this->buffer_as<OUTTYPE>( b );
                    if (elements != this->size()) {
                        throw std::runtime_error(
                            "ArrayLike.write_as(): Mismatch between internal "
                            "array size and size of buffer output" );
                    }
                    os.write( reinterpret_cast<const char *>( b ),
                              elements * sizeof( OUTTYPE ) );
                    delete[] b;
                    return os;
                }

                // virtual view_t<T, 1>& operator[]( size_t ind ) {
                //     return this->view( 0, ind );
                // }

                // virtual const view_t<T, 1>& operator[]( size_t ind ) const {
                //     return this->view( 0, ind );
                // }

                static bool global_registry_contains(
                    const ArrayLike<T, 1> *key );


                static global_registry_entry_t<T, 1>& register_parent(
                    const ArrayLike<T, 1> *key );

                static void deregister_parent( const ArrayLike<T, 1> *key );


                static void register_view( const ArrayLike<T, 1> *key,
                                           const ArrayView<T, 0>& v );


            private:
                local_registry_t<T, 1> _local_views;
                static global_registry_t<T, 1> _global_views;
                std::array<size_t, 1> _dimensions;
        };

        template<typename T>
        global_registry_t<T, 1> ArrayLike<T, 1>::_global_views;

        template<typename T, size_t N>
        class ArrayView : public ArrayLike<T, N> {
            public:
                ArrayView() :
                    ArrayLike<T, N>(),
                    _parent { nullptr },
                    _cparent { nullptr },
                    _const_dim { 0 } {}

                ArrayView( ArrayLike<T, N + 1> *parent, size_t const_dim,
                           size_t const_dim_ind ) :
                    ArrayLike<T, N>(),
                    _parent { parent },
                    _cparent { parent },
                    _const_dim { const_dim },
                    _key { const_dim, const_dim_ind } {
                    size_t counter = 0;
                    std::array<size_t, N> newdims;
                    for (size_t i = 0; i <= N; ++i) {
                        if (i != const_dim) {
                            _var_dims[ counter ] = i;
                            newdims[ counter ]   = parent->dimension( i );
                            counter++;
                        }
                    }
                    _mapped_indices[ _const_dim ] = const_dim_ind;

                    this->_dimensions = newdims;
                }

                ArrayView( const ArrayLike<T, N + 1> *parent, size_t const_dim,
                           size_t const_dim_ind ) :
                    ArrayLike<T, N>(),
                    _parent { nullptr },
                    _cparent { parent },
                    _const_dim { const_dim },
                    _key { const_dim, const_dim_ind } {
                    size_t counter = 0;
                    std::array<size_t, N> newdims;
                    for (size_t i = 0; i <= N; ++i) {
                        if (i != const_dim) {
                            _var_dims[ counter ] = i;
                            newdims[ counter ]   = parent->dimension( i );
                            counter++;
                        }
                    }
                    _mapped_indices[ _const_dim ] = const_dim_ind;
                    this->_dimensions             = newdims;
                }

                virtual ~ArrayView() {}

                ArrayView( const ArrayView<T, N>& other ) :
                    ArrayLike<T, N>(),
                    _parent { other._parent },
                    _const_dim { other._const_dim },
                    _key { other._key },
                    _var_dims { other._var_dims },
                    _mapped_indices { other._mapped_indices },
                    _dimensions { other._dimensions } {}

                ArrayView( ArrayView<T, N>&& other ) noexcept :
                    ArrayLike<T, N>() {
                    ::swap( *this, other );
                }

                ArrayView<T, N>& operator=( ArrayView<T, N> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ArrayView<T, N>& a,
                                      ArrayView<T, N>& b ) noexcept;

                virtual T& at( const std::array<size_t, N>& coords ) override {
                    if (_parent == nullptr) {
                        throw std::logic_error(
                            "Attempted non-const access to const view" );
                    }
                    for (size_t i = 0; i < N; ++i) {
                        _mapped_indices[ _var_dims[ i ] ] = coords[ i ];
                    }
                    return _parent->at( _mapped_indices );
                }

                virtual const T& at(
                    const std::array<size_t, N>& coords ) const override {
                    std::array<size_t, N + 1> mapped = _mapped_indices;
                    for (size_t i = 0; i < N; ++i) {
                        mapped[ _var_dims[ i ] ] = coords[ i ];
                    }
                    return _cparent->at( mapped );
                }

                virtual std::array<size_t, N>& dimensions() override {
                    return _dimensions;
                }

                virtual const std::array<size_t, N>& dimensions()
                    const override {
                    return _dimensions;
                }

                virtual std::pair<size_t, size_t> key() const { return _key; }

                virtual ArrayLike<T, N + 1>& parent() const {
                    return *_parent;
                }

                virtual void redimension(
                    const std::array<size_t, N>& dims ) override {
                    throw std::logic_error(
                        "Cannot redimension an array view!" );
                }

            private:
                ArrayLike<T, N + 1> *_parent;
                const ArrayLike<T, N + 1> *_cparent;
                size_t _const_dim;
                std::pair<size_t, size_t> _key;
                std::array<size_t, N> _var_dims;
                std::array<size_t, N + 1> _mapped_indices;
                std::array<size_t, N> _dimensions;
        };

        template<typename T>
        class ArrayView<T, 1> : public ArrayLike<T, 1> {
            public:
                ArrayView() :
                    ArrayLike<T, 1>(),
                    _parent { nullptr },
                    _cparent { nullptr },
                    _const_dim { 0 } {}

                ArrayView( ArrayLike<T, 2> *parent, size_t const_dim,
                           size_t const_dim_ind ) :
                    ArrayLike<T, 1>(),
                    _parent { parent },
                    _cparent { parent },
                    _const_dim { const_dim },
                    _key { const_dim, const_dim_ind } {
                    size_t counter = 0;
                    std::array<size_t, 1> newdims;
                    for (size_t i = 0; i <= 1; ++i) {
                        if (i != const_dim) {
                            _var_dims[ counter ] = i;
                            newdims[ counter ]   = parent->dimension( i );
                            counter++;
                        }
                    }
                    _mapped_indices[ _const_dim ] = const_dim_ind;
                    this->_dimensions             = newdims;
                }

                ArrayView( const ArrayLike<T, 2> *parent, size_t const_dim,
                           size_t const_dim_ind ) :
                    ArrayLike<T, 1>(),
                    _parent { nullptr },
                    _cparent { parent },
                    _const_dim { const_dim },
                    _key { const_dim, const_dim_ind } {
                    size_t counter = 0;
                    std::array<size_t, 1> newdims;
                    for (size_t i = 0; i <= 1; ++i) {
                        if (i != const_dim) {
                            _var_dims[ counter ] = i;
                            newdims[ counter ]   = parent->dimension( i );
                            counter++;
                        }
                    }
                    _mapped_indices[ _const_dim ] = const_dim_ind;
                    this->_dimensions             = newdims;
                }

                virtual ~ArrayView() {}

                ArrayView( const ArrayView<T, 1>& other ) :
                    ArrayLike<T, 1>(),
                    _parent { other._parent },
                    _const_dim { other._const_dim },
                    _key { other._key },
                    _var_dims { other._var_dims },
                    _mapped_indices { other._mapped_indices },
                    _dimensions { other._dimensions } {}

                ArrayView( ArrayView<T, 1>&& other ) noexcept :
                    ArrayLike<T, 1>() {
                    ::swap( *this, other );
                }

                ArrayView<T, 1>& operator=( ArrayView<T, 1> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ArrayView<T, 1>& a,
                                      ArrayView<T, 1>& b ) noexcept;

                virtual T& at( const std::array<size_t, 1>& coords ) override {
                    if (_parent == nullptr) {
                        throw std::logic_error(
                            "Attempted non-const access to const view" );
                    }
                    _mapped_indices[ _var_dims[ 0 ] ] = coords[ 0 ];
                    return _parent->at( _mapped_indices );
                }

                virtual const T& at(
                    const std::array<size_t, 1>& coords ) const override {
                    std::array<size_t, 2> mapped = _mapped_indices;
                    mapped[ _var_dims[ 0 ] ]     = coords[ 0 ];
                    return _cparent->at( mapped );
                }

                virtual std::array<size_t, 1>& dimensions() override {
                    return _dimensions;
                }

                virtual const std::array<size_t, 1>& dimensions()
                    const override {
                    return _dimensions;
                }

                virtual std::pair<size_t, size_t> key() const { return _key; }

                virtual ArrayLike<T, 2>& parent() const { return *_parent; }

                virtual void redimension(
                    const std::array<size_t, 1>& dims ) override {
                    throw std::logic_error(
                        "Cannot redimension an array view!" );
                }

            private:
                ArrayLike<T, 2> *_parent;
                const ArrayLike<T, 2> *_cparent;
                size_t _const_dim;
                std::pair<size_t, size_t> _key;
                std::array<size_t, 1> _var_dims;
                std::array<size_t, 2> _mapped_indices;
                std::array<size_t, 1> _dimensions;
        };

        template<typename T>
        class ArrayView<T, 0> {
            public:
                ArrayView() : _parent { nullptr }, _parent_ind { 0 } {}

                ArrayView( ArrayLike<T, 1> *parent, size_t const_dim,
                           size_t const_dim_ind ) :
                    _parent { parent },
                    _cparent { parent },
                    _parent_ind { const_dim_ind },
                    _key { const_dim, const_dim_ind } {
                    if (const_dim != 0) {
                        throw std::range_error(
                            "Invalid dimension specified!" );
                    }
                }

                ArrayView( const ArrayLike<T, 1> *parent, size_t const_dim,
                           size_t const_dim_ind ) :
                    _parent { nullptr },
                    _cparent { parent },
                    _parent_ind { const_dim_ind },
                    _key { const_dim, const_dim_ind } {
                    if (const_dim != 0) {
                        throw std::range_error(
                            "Invalid dimension specified!" );
                    }
                }

                ArrayView( const ArrayView<T, 0>& other ) :
                    _parent { other._parent },
                    _cparent { other._cparent },
                    _parent_ind { other._parent_ind },
                    _key { other._key } {}

                ArrayView( ArrayView<T, 0>&& other ) noexcept {
                    ::swap( *this, other );
                }

                ArrayView<T, 0>& operator=( ArrayView<T, 0> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                ArrayView<T, 0>& operator=( const T& other ) {
                    this->at() = other;
                    return *this;
                }

                virtual ~ArrayView() {}

                friend void ::swap<>( ArrayView<T, 0>& a,
                                      ArrayView<T, 0>& b ) noexcept;

                // template<typename OTHERTYPE,
                //          typename std::enable_if<
                //              std::is_convertible<OTHERTYPE, T>::value>::type
                //          = true>
                // operator T() const {
                //     return (OTHERTYPE)this->at();
                // }

                virtual T& at() {
                    if (_parent != nullptr) {
                        return _parent->at( { _parent_ind } );
                    } else {
                        throw std::logic_error(
                            "Requested non-const reference from const view!" );
                    }
                }

                virtual const T& at() const {
                    return _cparent->at( { _parent_ind } );
                }

                virtual T& at( const std::array<T, 0>& ) { return this->at(); }

                virtual const T& at( const std::array<T, 0>& ) const {
                    return this->at();
                }

                virtual T& at( size_t ind ) {
                    if (ind > 0) {
                        throw std::range_error(
                            "0-D array view has only one element" );
                    }
                    return this->at();
                }

                virtual const T& at( size_t ind ) const {
                    if (ind > 0) {
                        throw std::range_error(
                            "0-D array view has only one element" );
                    }
                    return this->at();
                }

                virtual size_t buffer( T *b ) const {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to buffer()" );
                    }
                    b[ 0 ] = this->at();
                    return 1;
                }

                virtual size_t debuffer( T *b ) {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to debuffer()" );
                    }
                    this->at() = b[ 0 ];
                    return 1;
                }

                // virtual size_t dimension( size_t n ) const { return 1; }

                virtual std::array<size_t, 0>& dimensions() {
                    return _empty_dim_array;
                }

                virtual const std::array<size_t, 0>& dimensions() const {
                    return _empty_dim_array;
                }

                virtual std::pair<size_t, size_t> key() const { return _key; }

                virtual ArrayLike<T, 1>& parent() const { return *_parent; }

                virtual std::istream& read( std::istream& is ) {
                    is >> this->at();
                    return is;
                }

                // virtual void redimension( const std::array<size_t, 1>& dims
                // ) {
                //     throw std::logic_error(
                //         "Cannot redimension an array view!" );
                // }


                virtual size_t size() const { return 1; }

                // virtual view_t<T, 0>& view( size_t dim, size_t dimind )
                // override {
                //     throw std::logic_error( "Can't get a view of a
                //     0-dimensional array" );
                // }
                // virtual const view_t<T, 0>& view( size_t dim,
                //                                   size_t dimind ) const
                //     override {
                //     throw std::logic_error( "Can't get a view of a
                //     0-dimensional array" );
                // }
                virtual std::ostream& write( std::ostream& os ) const {
                    os << this->at();
                    return os;
                }

                // protected:
                // virtual std::array<size_t, 0>& dimensions override {
                //     return _empty_dim_array;
                // }
                // virtual const std::array<size_t, 0>& dimensions
                // const override {
                //     return _empty_dim_array;
                // }

            private:
                ArrayLike<T, 1> *_parent;
                const ArrayLike<T, 1> *_cparent;
                size_t _parent_ind;
                std::pair<size_t, size_t> _key;
                std::array<size_t, 0> _empty_dim_array;
        };

        template<typename T, size_t N>
        class NDimensionalArray : public ArrayLike<T, N> {
            public:
                NDimensionalArray( const std::array<size_t, N>& dims ) :
                    ArrayLike<T, N>( dims ) {
                    this->redimension( dims );
                }

                virtual ~NDimensionalArray() {}

                NDimensionalArray( NDimensionalArray<T, N>&& other ) noexcept {
                    ::swap( *this, other );
                }

                NDimensionalArray<T, N>& operator=(
                    NDimensionalArray<T, N> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( NDimensionalArray<T, N>& a,
                                      NDimensionalArray<T, N>& b ) noexcept;

                NDimensionalArray( const NDimensionalArray<T, N>& other ) :
                    ArrayLike<T, N>( other ), _internal { other._internal } {}

                virtual T& at( const std::array<size_t, N>& coords ) override {
                    return this->internal().at( this->coords2index( coords ) );
                }

                virtual const T& at(
                    const std::array<size_t, N>& coords ) const override {
                    return this->internal().at( this->coords2index( coords ) );
                }

                virtual void redimension(
                    const std::array<size_t, N>& dims ) override {
                    size_t intsize = 1;
                    for (auto it = dims.begin(); it != dims.end(); ++it) {
                        intsize *= ( *it );
                    }
                    _internal.clear();
                    _internal.resize( intsize );
                }

                virtual size_t size() const override {
                    return _internal.size();
                }

                virtual std::ostream& write( std::ostream& os ) const {
                    os.write(
                        reinterpret_cast<const char *>( _internal.data() ),
                        _internal.size() * sizeof( T ) );
                    return os;
                }

                virtual size_t buffer( T *b ) const {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to buffer()" );
                    }
                    std::memcpy( b, _internal.data(),
                                 _internal.size() * sizeof( T ) );
                    return _internal.size();
                }

                template<typename BUFFERTYPE>
                size_t buffer_as( BUFFERTYPE *b ) const {
                    if (b == nullptr) {
                        throw std::invalid_argument(
                            "Null pointer passed to buffer()" );
                    }
                    for (size_t i = 0; i < _internal.size(); ++i) {
                        b[ i ] = (BUFFERTYPE)( _internal.at( i ) );
                    }
                    return _internal.size();
                }

            protected:
                size_t coords2index(
                    const std::array<size_t, N>& coords ) const {
                    if (N == 1) {
                        return coords[ 0 ];
                    }
                    size_t index   = coords[ N - 1 ];
                    size_t scale   = 1;
                    size_t counter = N - 1;
                    while (counter > 0) {
                        scale *= this->dimension( counter );
                        index += scale * coords[ --counter ];
                    }
                    return index;
                }

                std::vector<T>& internal() { return _internal; }

                const std::vector<T>& internal() const { return _internal; }

            private:
                std::vector<T> _internal;
        };

        template<typename T>
        class TwoDimensionalArray : public NDimensionalArray<T, 2> {
            public:
                TwoDimensionalArray( const std::array<size_t, 2>& dims ) :
                    NDimensionalArray<T, 2>( dims ) {}

                TwoDimensionalArray( size_t N1, size_t N2 ) :
                    TwoDimensionalArray<T>(
                        std::array<size_t, 2> { N1, N2 } ) {}

                TwoDimensionalArray( const TwoDimensionalArray<T>& other ) :
                    NDimensionalArray<T, 2>( other ) {}

                TwoDimensionalArray(
                    TwoDimensionalArray<T>&& other ) noexcept {
                    ::swap( *this, other );
                }

                virtual ~TwoDimensionalArray() {}

                TwoDimensionalArray<T>& operator=(
                    TwoDimensionalArray<T> other ) {
                    ::swap( *this, other );
                    return *this;
                }
        };

        template<typename T>
        class ThreeDimensionalArray : public NDimensionalArray<T, 3> {
            public:
                ThreeDimensionalArray( const std::array<size_t, 3>& dims ) :
                    NDimensionalArray<T, 3>( dims ) {}

                ThreeDimensionalArray( size_t N1, size_t N2, size_t N3 ) :
                    ThreeDimensionalArray<T>(
                        std::array<size_t, 3> { N1, N2, N3 } ) {}

                ThreeDimensionalArray(
                    const ThreeDimensionalArray<T>& other ) :
                    NDimensionalArray<T, 3>( other ) {}

                ThreeDimensionalArray(
                    ThreeDimensionalArray<T>&& other ) noexcept {
                    ::swap( *this, other );
                }

                virtual ~ThreeDimensionalArray() {}

                ThreeDimensionalArray<T>& operator=(
                    ThreeDimensionalArray<T> other ) {
                    ::swap( *this, other );
                    return *this;
                }
        };
    }  // namespace arrays
}  // namespace NCPA

template<typename T, size_t N>
static void swap( NCPA::arrays::_abstract_arraylike<T, N>& a,
                  NCPA::arrays::_abstract_arraylike<T, N>& b ) noexcept {
    using std::swap;
}

template<typename T, size_t N>
static void swap( NCPA::arrays::ArrayLike<T, N>& a,
                  NCPA::arrays::ArrayLike<T, N>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::arrays::_abstract_arraylike<T, N>&>( a ),
            static_cast<NCPA::arrays::_abstract_arraylike<T, N>&>( b ) );
    swap( a._dimensions, b._dimensions );
    swap( a._local_views, b._local_views );
}

template<typename T>
static void swap( NCPA::arrays::ArrayLike<T, 1>& a,
                  NCPA::arrays::ArrayLike<T, 1>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::arrays::_abstract_arraylike<T, 1>&>( a ),
            static_cast<NCPA::arrays::_abstract_arraylike<T, 1>&>( b ) );
    swap( a._dimensions, b._dimensions );
    swap( a._local_views, b._local_views );
}

template<typename T, size_t N>
static void swap( NCPA::arrays::ArrayView<T, N>& a,
                  NCPA::arrays::ArrayView<T, N>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::arrays::ArrayLike<T, N>&>( a ),
            static_cast<NCPA::arrays::ArrayLike<T, N>&>( b ) );
    swap( a._parent, b._parent );
    swap( a._cparent, b._cparent );
    swap( a._const_dim, b._const_dim );
    swap( a._key, b._key );
    swap( a._var_dims, b._var_dims );
    swap( a._mapped_indices, b._mapped_indices );
    swap( a._dimensions, b._dimensions );
}

template<typename T, size_t N>
static void swap( NCPA::arrays::ArrayView<T, 1>& a,
                  NCPA::arrays::ArrayView<T, 1>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::arrays::ArrayLike<T, 1>&>( a ),
            static_cast<NCPA::arrays::ArrayLike<T, 1>&>( b ) );
    swap( a._parent, b._parent );
    swap( a._cparent, b._cparent );
    swap( a._const_dim, b._const_dim );
    swap( a._key, b._key );
    swap( a._var_dims, b._var_dims );
    swap( a._mapped_indices, b._mapped_indices );
    swap( a._dimensions, b._dimensions );
}

template<typename T>
static void swap( NCPA::arrays::ArrayView<T, 0>& a,
                  NCPA::arrays::ArrayView<T, 0>& b ) noexcept {
    using std::swap;
    swap( a._parent, b._parent );
    swap( a._cparent, b._cparent );
    swap( a._parent_ind, b._parent_ind );
    swap( a._key, b._key );
}

template<typename T, size_t N>
static void swap( NCPA::arrays::NDimensionalArray<T, N>& a,
                  NCPA::arrays::NDimensionalArray<T, N>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::arrays::ArrayLike<T, N>&>( a ),
            static_cast<NCPA::arrays::ArrayLike<T, N>&>( b ) );
    swap( a._internal, b._internal );
}

template<typename T>
static void swap( NCPA::arrays::TwoDimensionalArray<T>& a,
                  NCPA::arrays::TwoDimensionalArray<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::arrays::NDimensionalArray<T, 2>&>( a ),
            static_cast<NCPA::arrays::NDimensionalArray<T, 2>&>( b ) );
}

template<typename T>
static void swap( NCPA::arrays::ThreeDimensionalArray<T>& a,
                  NCPA::arrays::ThreeDimensionalArray<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::arrays::NDimensionalArray<T, 3>&>( a ),
            static_cast<NCPA::arrays::NDimensionalArray<T, 3>&>( b ) );
}

