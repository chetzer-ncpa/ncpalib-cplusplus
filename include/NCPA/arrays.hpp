#pragma once

#include "NCPA/types.hpp"

#include <cstddef>
#include <cstring>

namespace NCPA {
    namespace arrays {
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
                         ENABLE_IF( !NCPA::types::is_deleteable<T> ) ) {
            delete[] v;
            v = nullptr;
        }

        template<typename T>
        void free_array( T *& v, size_t nr,
                         ENABLE_IF( NCPA::types::is_deleteable<T> ) ) {
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
                         ENABLE_IF( !NCPA::types::is_deleteable<T> ) ) {
            for ( size_t i = 0; i < nr; i++ ) {
                delete[] v[ i ];
            }
            delete[] v;
            v = nullptr;
        }

        template<typename T>
        void free_array( T **& v, size_t nr, size_t nc,
                         ENABLE_IF( NCPA::types::is_deleteable<T> ) ) {
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
                         ENABLE_IF( !NCPA::types::is_deleteable<T> ) ) {
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
                         ENABLE_IF( NCPA::types::is_deleteable<T> ) ) {
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
                K -= ( (int)N ) * (K < 0 ? -1 : 1);
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
        @brief Sets each element of a 2-D array to a constant value.
        @param v The array to fill.
        @param nr The dimension of the array.
        @param val The value to set each element to.
        */
        template<typename T>
        void fill( T *v, size_t nr, const T& val ) {
            for ( size_t i = 0; i < nr; i++ ) {
                v[ i ] = val;
            }
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
                // for ( size_t j = 0; j < nc; j++ ) {
                //     v[ i ][ j ] = val;
                // }
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
                // for ( size_t j = 0; j < nc; j++ ) {
                //     for ( size_t k = 0; k < nz; k++ ) {
                //         v[ i ][ j ][ k ] = val;
                //     }
                // }
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
        T *index_vector( size_t n, T a = 0 ) {
            T *ivec = zeros<T>( n );
            for ( size_t i = 0; i < n; i++ ) {
                ivec[ i ] = (T)i + a;
            }
            return ivec;
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
            for ( j = 0; j < N; j++ ) {
                tempvec[ N - 1 - j ] = in[ j ];
            }
            std::memcpy( out, tempvec, N * sizeof( T ) );

            delete[] tempvec;
        }


    }  // namespace arrays
}  // namespace NCPA