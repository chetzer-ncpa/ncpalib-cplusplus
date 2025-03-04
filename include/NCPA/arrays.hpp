#pragma once

#include "NCPA/constants.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/ndvector.hpp"
#include "NCPA/types.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <numeric>
#include <vector>

namespace NCPA {
    namespace arrays {

        template<typename T>
        // class vector2d_t : public std::vector<std::vector<T>> {
        class vector2d_t : public ndvector<2, T> {
            public:
                using ndvector<2,T>::set;
                vector2d_t() : ndvector<2, T>() {}

                vector2d_t( size_t nx1, size_t nx2, const T& val = (T)0 ) :
                    ndvector<2, T>() {
                    this->resize2d( nx1, nx2, val );
                }

                vector2d_t( size_t nx1, size_t nx2, const T **vals ) :
                    vector2d_t( nx1, nx2 ) {
                    for ( size_t i = 0; i < nx1; ++i ) {
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
                    // nx1 = this->size();
                    // if ( nx1 > 0 ) {
                    //     nx2 = this->at( 0 ).size();
                    // } else {
                    //     nx2 = 0;
                    // }
                }

                virtual size_t dim( size_t dimnum ) const {
                    return this->shape()[ dimnum ];
                    // switch ( dimnum ) {
                    //     case 0:
                    //         return this->size();
                    //         break;
                    //     case 1:
                    //         return ( this->size() > 0 ? this->at( 0 ).size()
                    //                                   : 0 );
                    //         break;
                    //     default:
                    //         throw std::range_error(
                    //             "vector2d_t.dim(): Invalid dimension "
                    //             "requested" );
                    // }
                }

                virtual void fill( T val ) {
                    this->set( val );
                    // for ( auto it = this->begin(); it != this->end(); ++it )
                    // {
                    //     it->assign( this->dim( 1 ), val );
                    // }
                }
        };

        template<typename T>
        // class vector3d_t : public std::vector<vector2d_t<T>> {
        class vector3d_t : public ndvector<3, T> {
            public:
                vector3d_t() : ndvector<3, T>() {}

                vector3d_t( size_t nx1, size_t nx2, size_t nx3,
                            const T& val = (T)0 ) :
                            ndvector<3, T>() {
                    this->resize3d( nx1, nx2, nx3, val );
                }

                vector3d_t( size_t nx1, size_t nx2, size_t nx3,
                            const T ***vals ) :
                    vector3d_t( nx1, nx2, nx3 ) {
                    for ( size_t i = 0; i < nx1; ++i ) {
                        for ( size_t j = 0; j < nx2; ++j ) {
                            this->at( i ).at( j ).assign(
                                vals[ i ][ j ], vals[ i ][ j ] + nx3 );
                        }
                    }
                }

                virtual void resize3d( size_t nx1, size_t nx2, size_t nx3,
                                       const T& val = (T)0 ) {
                    this->reshape( { nx1, nx2, nx3 }, val );
                    // this->resize( nx1, vector2d_t<T>( nx2, nx3, val ) );
                }

                virtual void size3d( size_t& nx1, size_t& nx2,
                                     size_t& nx3 ) const {
                    auto dims = this->shape();
                    nx1 = dims[0];
                    nx2 = dims[1];
                    nx3 = dims[2];
                    // nx1 = this->size();
                    // if ( nx1 > 0 ) {
                    //     this->at( 0 ).size2d( nx2, nx3 );
                    // } else {
                    //     nx2 = 0;
                    //     nx3 = 0;
                    // }
                }

                virtual size_t dim( size_t dimnum ) const {
                    return this->shape()[ dimnum ];
                    // switch ( dimnum ) {
                    //     case 0:
                    //         return this->size();
                    //         break;
                    //     case 1:
                    //         return ( this->size() > 0 ? this->at( 0 ).size()
                    //                                   : 0 );
                    //         break;
                    //     case 2:
                    //         return ( this->size() > 0
                    //                          && this->at( 0 ).size() > 0
                    //                      ? this->at( 0 ).at( 0 ).size()
                    //                      : 0 );
                    //         break;
                    //     default:
                    //         throw std::range_error(
                    //             "vector3d_t.dim(): Invalid dimension "
                    //             "requested" );
                    // }
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
            for ( size_t i = 0; i < nr; i++ ) {
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
            for ( size_t i = 0; i < nr; i++ ) {
                v[ i ] = new T *[ nc ];
                for ( size_t j = 0; j < nc; j++ ) {
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
            for ( size_t i = 0; i < nr; i++ ) {
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
            for ( size_t i = 0; i < nr; i++ ) {
                delete[] v[ i ];
            }
            delete[] v;
            v = nullptr;
        }

        template<typename T>
        void free_array( T **& v, size_t nr, size_t nc,
                         ENABLE_FUNCTION_IF_DELETEABLE( T ) ) {
            for ( size_t i = 0; i < nr; i++ ) {
                for ( size_t j = 0; j < nc; j++ ) {
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
            for ( i = 0; i < nd1; ++i ) {
                if ( data[ i ] != NULL ) {
                    for ( j = 0; j < nd2; ++j ) {
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
            for ( i = 0; i < nd1; ++i ) {
                if ( data[ i ] != NULL ) {
                    for ( j = 0; j < nd2; ++j ) {
                        for ( k = 0; i < nd3; k++ ) {
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
            while ( std::abs( K ) > N ) {
                K -= ( (int)N ) * ( K < 0 ? -1 : 1 );
            }
            size_t i;
            T *tempvec = zeros<T>( N );
            if ( K < 0 ) {
                // move from the back to the front
                size_t negshift = (size_t)( -K );
                // first, move the last -K values to the front
                for ( i = 0; i < negshift; i++ ) {
                    tempvec[ i ] = X[ N - negshift + i ];
                }
                for ( i = negshift; i < N; i++ ) {
                    tempvec[ i ] = X[ i - negshift ];
                }
            } else if ( K > 0 ) {
                // move from the front to the back
                // first, move the first K values to the back
                for ( i = 0; i < K; i++ ) {
                    tempvec[ N - K + i ] = X[ i ];
                }
                // Now, move the next N-K values to the front
                for ( i = K; i < N; i++ ) {
                    tempvec[ i - K ] = X[ i ];
                }
            } else {
                std::memcpy( tempvec, X, N * sizeof( T ) );
            }

            // copy over so it can be done in place if desired
            if ( out == nullptr ) {
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
            for ( size_t i = 0; i < nr; i++ ) {
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
            for ( size_t i = 0; i < nr; i++ ) {
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
            if ( from != nullptr ) {
                std::copy( from, from + n, to );
            } else if ( to != nullptr ) {
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
            if ( from != nullptr ) {
                for ( size_t i = 0; i < n1; i++ ) {
                    copy( from[ i ], n2, to[ i ] );
                }
            } else if ( to != nullptr ) {
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
            if ( from != nullptr ) {
                for ( size_t i = 0; i < n1; i++ ) {
                    copy( from[ i ], n2, n3, to[ i ] );
                }
            } else if ( to != nullptr ) {
                free_array( to, n1, n2, n3 );
            }
        }

        // /**
        // Dynamically allocates a new array and sets the elements to
        // (0, 1, 2, ..., n-1).
        // @brief Returns a new array of index values.
        // @param n The size of the array.
        // @returns A pointer to the newly-allocated array.
        // */
        // template<typename T>
        // T *index_vector( size_t n ) {
        //     T *ivec = zeros<T>( n );
        //     for ( size_t i = 0; i < n; i++ ) {
        //         ivec[ i ] = (T)i;
        //     }
        //     return ivec;
        // }

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
            for ( size_t i = 0; i < n; i++ ) {
                inds[ i ] = (T)i + a;
            }
            return inds;
        }

        // template<typename T>
        // T *index_vector( size_t n, T a = 0 ) {
        //     T *ivec = zeros<T>( n );
        //     for ( size_t i = 0; i < n; i++ ) {
        //         ivec[ i ] = (T)i + a;
        //     }
        //     return ivec;
        // }

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
            for ( j = 0; j < N; j++ ) {
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
            if ( vec1.size() != vec2.size() ) {
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
            for ( std::size_t i = 0; i < vec.size(); ++i ) {
                if ( done[ i ] ) {
                    continue;
                }
                done[ i ]          = true;
                std::size_t prev_j = i;
                std::size_t j      = p[ i ];
                while ( i != j ) {
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
            for ( size_t i = 0; i < basevec.size(); ++i ) {
                if ( i > 0 ) {
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
            for ( size_t i = 0; i < basevec.size(); i++ ) {
                os << basevec[ i ];
                for ( size_t j = 0; j < items.size(); j++ ) {
                    os << separator << items[ j ][ i ];
                }
                os << std::endl;
            }
        }

        // template<typename T, typename U>
        // void write( std::ostream& os, const std::vector<T> basevec,
        //             const std::vector<U> depvec,
        //             const std::string& separator = " " ) {
        //     for ( size_t i = 0; i < basevec.size(); i++ ) {
        //         os << basevec[ i ] << separator << depvec[ i ] << std::endl;
        //     }
        // }

        template<typename T, typename U>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::vector<U> depvec,
                    const std::string& separator = " ",
                    ENABLE_FUNCTION_IF_ARITHMETIC( U ) ) {
            for ( size_t i = 0; i < basevec.size(); i++ ) {
                os << basevec[ i ] << separator << depvec[ i ] << std::endl;
            }
        }

        template<typename T, typename U>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::vector<U> depvec,
                    const std::string& separator = " ",
                    ENABLE_FUNCTION_IF_COMPLEX( U ) ) {
            for ( size_t i = 0; i < basevec.size(); i++ ) {
                os << basevec[ i ] << separator << depvec[ i ].real()
                   << separator << depvec[ i ].imag() << std::endl;
            }
        }

        /**
        @brief Performs element-wise array addition.
        @param v1 The first array to add.
        @param v2 The second array to add.
        @returns The array holding the sum.
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
        void scale_array( size_t N, const U *in, T factor, U *& out ) {
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
            const T& v1, U scalar,
            typename std::enable_if<NCPA::types::is_iterable_of<T, U>::value,
                                    int>::type ENABLER
            = 0 ) {
            if ( scalar == NCPA::constants::one<U>() ) {
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
            if ( factor != NCPA::constants::one<T>() ) {
                for ( size_t i = 0; i < N; i++ ) {
                    in[ i ] *= factor;
                }
            }
        }

        /**
        Scales an array.
        @brief Performs array scaling.
        @param v1 The array to multiply.
        @param scalar The scalar to multiply by.
        @returns The scaled array.
        */
        template<typename T, typename U>
        T offset_vector(
            const T& v1, U scalar,
            typename std::enable_if<NCPA::types::is_iterable_of<T, U>::value,
                                    int>::type ENABLER
            = 0 ) {
            if ( NCPA::constants::is_zero<U>( scalar ) ) {
                return v1;
            }
            T v3 = v1;
            for ( auto i = 0; i < v1.size(); i++ ) {
                v3[ i ] += scalar;
            }
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
        void offset_array( size_t N, U *in, T factor ) {
            if ( !NCPA::constants::is_zero<T>( factor ) ) {
                for ( size_t i = 0; i < N; i++ ) {
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
            for ( size_t i = 0; i < N; i++ ) {
                tempvec[ i ] = v1[ i ] - v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

    }  // namespace arrays
}  // namespace NCPA

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
