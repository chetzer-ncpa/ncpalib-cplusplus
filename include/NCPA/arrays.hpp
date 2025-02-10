#pragma once

#include "NCPA/defines.hpp"
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
        class vector2d_t : public std::vector<std::vector<T>> {
            public:
                vector2d_t() : std::vector<std::vector<T>>() {}

                vector2d_t( size_t nx1, size_t nx2, const T& val ) :
                    std::vector<std::vector<T>>() {
                    this->resize2d( nx1, nx2, val );
                }

                virtual void resize2d( size_t nx1, size_t nx2,
                                       const T& val = (T)0 ) {
                    this->resize( nx1, std::vector<T>( nx2, val ) );
                }

                virtual void size2d( size_t& nx1, size_t& nx2 ) const {
                    nx1 = this->size();
                    if ( nx1 > 0 ) {
                        nx2 = this->at( 0 ).size();
                    } else {
                        nx2 = 0;
                    }
                }
        };

        template<typename T>
        class vector3d_t : public std::vector<vector2d_t<T>> {
            public:
                vector3d_t() : std::vector<vector2d_t<T>>() {}

                vector3d_t( size_t nx1, size_t nx2, size_t nx3,
                            const T& val ) :
                    std::vector<vector2d_t<T>>() {
                    this->resize3d( nx1, nx2, val );
                }

                virtual void resize3d( size_t nx1, size_t nx2, size_t nx3,
                                       const T& val = (T)0 ) {
                    this->resize( nx1, vector2d_t<T>( nx2, nx3, val ) );
                }

                virtual void size3d( size_t& nx1, size_t& nx2,
                                     size_t& nx3 ) const {
                    nx1 = this->size();
                    if ( nx1 > 0 ) {
                        this->at( 0 ).size2d( nx2, nx3 );
                    } else {
                        nx2 = 0;
                        nx3 = 0;
                    }
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
            std::copy( from, from + n, to );
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
            for ( size_t i = 0; i < n1; i++ ) {
                copy( from[ i ], n2, to[ i ] );
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
            for ( size_t i = 0; i < n1; i++ ) {
                copy( from[ i ], n2, n3, to[ i ] );
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
                os << basevec[i];
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


    }  // namespace arrays
}  // namespace NCPA
