#pragma once

#include <gtest/gtest.h>
#include <gmock/gmock.h>


#define ASSERT_ARRAY_EQ( N, AR1, AR2 ) { \
    for (auto i = 0; i < N; i++) { \
        ASSERT_EQ( AR1[i], AR2[i] ) << "Arrays differ at index " << i; \
    } \
}

#define ASSERT_DOUBLE_ARRAY_EQ( N, AR1, AR2 ) { \
    for (auto i = 0; i < N; i++) { \
        ASSERT_DOUBLE_EQ( AR1[i], AR2[i] ) << "Arrays differ at index " << i; \
    } \
}

#define ASSERT_ARRAY2D_EQ( M, N, AR1, AR2 ) { \
    for (auto i = 0; i < M; i++) { \
        for (auto j = 0; j < N; j++) { \
            ASSERT_EQ( AR1[i][j], AR2[i][j] ) << "Arrays differ at index " << i << "," << j; \
        } \
    } \
}

#define ASSERT_ARRAY2D_DOUBLE_EQ( M, N, AR1, AR2 ) { \
    for (auto i = 0; i < M; i++) { \
        for (auto j = 0; j < N; j++) { \
            ASSERT_DOUBLE_EQ( AR1[i][j], AR2[i][j] ) << "Arrays differ at index " << i << "," << j; \
        } \
    } \
}

#define ASSERT_ARRAY3D_EQ( M, N, P, AR1, AR2 ) { \
    for (auto i = 0; i < M; i++) { \
        for (auto j = 0; j < N; j++) { \
            for (auto k = 0; k < P; k++) { \
                ASSERT_EQ( AR1[i][j][k], AR2[i][j][k] ) << "Arrays differ at index " << i << "," << j << "," << k; \
            } \
        } \
    } \
}

#define ASSERT_ARRAY3D_DOUBLE_EQ( M, N, P, AR1, AR2 ) { \
    for (auto i = 0; i < M; i++) { \
        for (auto j = 0; j < N; j++) { \
            for (auto k = 0; k < P; k++) { \
                ASSERT_DOUBLE_EQ( AR1[i][j][k], AR2[i][j][k] ) << "Arrays differ at index " << i << "," << j << "," << k; \
            } \
        } \
    } \
}