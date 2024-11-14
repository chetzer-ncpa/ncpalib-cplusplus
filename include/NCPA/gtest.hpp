#pragma once

#if __has_include( "gtest/gtest.h" )
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#define EXPECT_NONFATAL_FAILURES( statement_, n_ ) \
    EXPECT_NONFATAL_FAILURE(                       \
        { EXPECT_NONFATAL_FAILURE( statement_, "" ); }, "Actual: " #n_ )

#define ASSERT_ARRAY_EQ( N, AR1, AR2 )             \
    {                                              \
        for ( auto i = 0; i < N; i++ ) {           \
            ASSERT_EQ( AR1[ i ], AR2[ i ] )        \
                << "Arrays differ at index " << i; \
        }                                          \
    }

#define ASSERT_ARRAY_DOUBLE_EQ( N, AR1, AR2 )      \
    {                                              \
        for ( auto i = 0; i < N; i++ ) {           \
            ASSERT_DOUBLE_EQ( AR1[ i ], AR2[ i ] ) \
                << "Arrays differ at index " << i; \
        }                                          \
    }

#define ASSERT_ARRAY2D_EQ( M, N, AR1, AR2 )                        \
    {                                                              \
        for ( auto i = 0; i < M; i++ ) {                           \
            for ( auto j = 0; j < N; j++ ) {                       \
                ASSERT_EQ( AR1[ i ][ j ], AR2[ i ][ j ] )          \
                    << "Arrays differ at index " << i << "," << j; \
            }                                                      \
        }                                                          \
    }

#define ASSERT_ARRAY2D_DOUBLE_EQ( M, N, AR1, AR2 )                 \
    {                                                              \
        for ( auto i = 0; i < M; i++ ) {                           \
            for ( auto j = 0; j < N; j++ ) {                       \
                ASSERT_DOUBLE_EQ( AR1[ i ][ j ], AR2[ i ][ j ] )   \
                    << "Arrays differ at index " << i << "," << j; \
            }                                                      \
        }                                                          \
    }

#define ASSERT_ARRAY3D_EQ( M, N, P, AR1, AR2 )                               \
    {                                                                        \
        for ( auto i = 0; i < M; i++ ) {                                     \
            for ( auto j = 0; j < N; j++ ) {                                 \
                for ( auto k = 0; k < P; k++ ) {                             \
                    ASSERT_EQ( AR1[ i ][ j ][ k ], AR2[ i ][ j ][ k ] )      \
                        << "Arrays differ at index " << i << "," << j << "," \
                        << k;                                                \
                }                                                            \
            }                                                                \
        }                                                                    \
    }

#define ASSERT_ARRAY3D_DOUBLE_EQ( M, N, P, AR1, AR2 )                        \
    {                                                                        \
        for ( auto i = 0; i < M; i++ ) {                                     \
            for ( auto j = 0; j < N; j++ ) {                                 \
                for ( auto k = 0; k < P; k++ ) {                             \
                    ASSERT_DOUBLE_EQ( AR1[ i ][ j ][ k ],                    \
                                      AR2[ i ][ j ][ k ] )                   \
                        << "Arrays differ at index " << i << "," << j << "," \
                        << k;                                                \
                }                                                            \
            }                                                                \
        }                                                                    \
    }

#define ASSERT_COMPLEX_DOUBLE_EQ( X, Y )        \
    {                                           \
        ASSERT_DOUBLE_EQ( X.real(), Y.real() ); \
        ASSERT_DOUBLE_EQ( X.imag(), Y.imag() ); \
    }

#define ASSERT_ARRAY_COMPLEX_DOUBLE_EQ( N, AR1, AR2 )                \
    {                                                                \
        for ( auto i = 0; i < N; i++ ) {                             \
            ASSERT_DOUBLE_EQ( AR1[ i ].real(), AR2[ i ].real() )     \
                << "Arrays' real values differ at index " << i;      \
            ASSERT_DOUBLE_EQ( AR1[ i ].imag(), AR2[ i ].imag() )     \
                << "Arrays' imaginary values differ at index " << i; \
        }                                                            \
    }

#define ASSERT_ARRAY2D_COMPLEX_DOUBLE_EQ( M, N, AR1, AR2 )                \
    {                                                                     \
        for ( auto i = 0; i < M; i++ ) {                                  \
            for ( auto j = 0; j < N; j++ ) {                              \
                ASSERT_DOUBLE_EQ( AR1[ i ][ j ].real(),                   \
                                  AR2[ i ][ j ].real() )                  \
                    << "Arrays' real values differ at index " << i << "," \
                    << j;                                                 \
                ASSERT_DOUBLE_EQ( AR1[ i ][ j ].imag(),                   \
                                  AR2[ i ][ j ].imag() )                  \
                    << "Arrays' imaginary values differ at index " << i   \
                    << "," << j;                                          \
            }                                                             \
        }                                                                 \
    }

#define ASSERT_ARRAY3D_COMPLEX_DOUBLE_EQ( M, N, P, AR1, AR2 )                 \
    {                                                                         \
        for ( auto i = 0; i < M; i++ ) {                                      \
            for ( auto j = 0; j < N; j++ ) {                                  \
                for ( auto k = 0; k < P; k++ ) {                              \
                    ASSERT_DOUBLE_EQ( AR1[ i ][ j ][ k ].real(),              \
                                      AR2[ i ][ j ][ k ].real() )             \
                        << "Arrays' real values differ at index " << i << "," \
                        << j << "," << k;                                     \
                    ASSERT_DOUBLE_EQ( AR1[ i ][ j ][ k ].imag(),              \
                                      AR2[ i ][ j ][ k ].imag() )             \
                        << "Arrays' imaginary values differ at index " << i   \
                        << "," << j << "," << k;                              \
                }                                                             \
            }                                                                 \
        }                                                                     \
    }

#define EXPECT_ARRAY_EQ( N, AR1, AR2 )             \
    {                                              \
        for ( auto i = 0; i < N; i++ ) {           \
            EXPECT_EQ( AR1[ i ], AR2[ i ] )        \
                << "Arrays differ at index " << i; \
        }                                          \
    }

#define EXPECT_ARRAY_DOUBLE_EQ( N, AR1, AR2 )      \
    {                                              \
        for ( auto i = 0; i < N; i++ ) {           \
            EXPECT_DOUBLE_EQ( AR1[ i ], AR2[ i ] ) \
                << "Arrays differ at index " << i; \
        }                                          \
    }

#define EXPECT_ARRAY2D_EQ( M, N, AR1, AR2 )                        \
    {                                                              \
        for ( auto i = 0; i < M; i++ ) {                           \
            for ( auto j = 0; j < N; j++ ) {                       \
                EXPECT_EQ( AR1[ i ][ j ], AR2[ i ][ j ] )          \
                    << "Arrays differ at index " << i << "," << j; \
            }                                                      \
        }                                                          \
    }

#define EXPECT_ARRAY2D_DOUBLE_EQ( M, N, AR1, AR2 )                 \
    {                                                              \
        for ( auto i = 0; i < M; i++ ) {                           \
            for ( auto j = 0; j < N; j++ ) {                       \
                EXPECT_DOUBLE_EQ( AR1[ i ][ j ], AR2[ i ][ j ] )   \
                    << "Arrays differ at index " << i << "," << j; \
            }                                                      \
        }                                                          \
    }

#define EXPECT_ARRAY3D_EQ( M, N, P, AR1, AR2 )                               \
    {                                                                        \
        for ( auto i = 0; i < M; i++ ) {                                     \
            for ( auto j = 0; j < N; j++ ) {                                 \
                for ( auto k = 0; k < P; k++ ) {                             \
                    EXPECT_EQ( AR1[ i ][ j ][ k ], AR2[ i ][ j ][ k ] )      \
                        << "Arrays differ at index " << i << "," << j << "," \
                        << k;                                                \
                }                                                            \
            }                                                                \
        }                                                                    \
    }

#define EXPECT_ARRAY3D_DOUBLE_EQ( M, N, P, AR1, AR2 )                        \
    {                                                                        \
        for ( auto i = 0; i < M; i++ ) {                                     \
            for ( auto j = 0; j < N; j++ ) {                                 \
                for ( auto k = 0; k < P; k++ ) {                             \
                    EXPECT_DOUBLE_EQ( AR1[ i ][ j ][ k ],                    \
                                      AR2[ i ][ j ][ k ] )                   \
                        << "Arrays differ at index " << i << "," << j << "," \
                        << k;                                                \
                }                                                            \
            }                                                                \
        }                                                                    \
    }
#define EXPECT_COMPLEX_DOUBLE_EQ( X, Y )        \
    {                                           \
        EXPECT_DOUBLE_EQ( X.real(), Y.real() ); \
        EXPECT_DOUBLE_EQ( X.imag(), Y.imag() ); \
    }

#define EXPECT_ARRAY_COMPLEX_DOUBLE_EQ( N, AR1, AR2 )                \
    {                                                                \
        for ( auto i = 0; i < N; i++ ) {                             \
            EXPECT_DOUBLE_EQ( AR1[ i ].real(), AR2[ i ].real() )     \
                << "Arrays' real values differ at index " << i;      \
            EXPECT_DOUBLE_EQ( AR1[ i ].imag(), AR2[ i ].imag() )     \
                << "Arrays' imaginary values differ at index " << i; \
        }                                                            \
    }

#define EXPECT_ARRAY2D_COMPLEX_DOUBLE_EQ( M, N, AR1, AR2 )                \
    {                                                                     \
        for ( auto i = 0; i < M; i++ ) {                                  \
            for ( auto j = 0; j < N; j++ ) {                              \
                EXPECT_DOUBLE_EQ( AR1[ i ][ j ].real(),                   \
                                  AR2[ i ][ j ].real() )                  \
                    << "Arrays' real values differ at index " << i << "," \
                    << j;                                                 \
                EXPECT_DOUBLE_EQ( AR1[ i ][ j ].imag(),                   \
                                  AR2[ i ][ j ].imag() )                  \
                    << "Arrays' imaginary values differ at index " << i   \
                    << "," << j;                                          \
            }                                                             \
        }                                                                 \
    }

#define EXPECT_ARRAY3D_COMPLEX_DOUBLE_EQ( M, N, P, AR1, AR2 )                 \
    {                                                                         \
        for ( auto i = 0; i < M; i++ ) {                                      \
            for ( auto j = 0; j < N; j++ ) {                                  \
                for ( auto k = 0; k < P; k++ ) {                              \
                    EXPECT_DOUBLE_EQ( AR1[ i ][ j ][ k ].real(),              \
                                      AR2[ i ][ j ][ k ].real() )             \
                        << "Arrays' real values differ at index " << i << "," \
                        << j << "," << k;                                     \
                    EXPECT_DOUBLE_EQ( AR1[ i ][ j ][ k ].imag(),              \
                                      AR2[ i ][ j ][ k ].imag() )             \
                        << "Arrays' imaginary values differ at index " << i   \
                        << "," << j << "," << k;                              \
                }                                                             \
            }                                                                 \
        }                                                                     \
    }

#endif