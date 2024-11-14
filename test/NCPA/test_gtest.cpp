// #include "NCPA/math.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <complex>
#include <stdexcept>
#include <vector>

#include <gtest/gtest-spi.h>

using namespace std;
using namespace testing;

#define TESTARRAY_SIZE_DIM1 5
#define TESTARRAY_1D        { 0, 1, 2, 3, 4 }
#define TESTARRAY_SIZE_DIM2 4
#define TESTARRAY_2D  \
    {                 \
        TESTARRAY_1D, \
        TESTARRAY_1D, \
        TESTARRAY_1D, \
        TESTARRAY_1D, \
    }
#define TESTARRAY_SIZE_DIM3 3
#define TESTARRAY_3D  \
    {                 \
        TESTARRAY_2D, \
        TESTARRAY_2D, \
        TESTARRAY_2D, \
    }
#define TESTCOMPLEX1 complex<double>( 1.0, -1.0 )
#define TESTCOMPLEX2 complex<double>( -1.0, 1.0 )
#define TESTARRAY_COMPLEX_1D \
    { TESTCOMPLEX1, TESTCOMPLEX1, TESTCOMPLEX2, TESTCOMPLEX2, TESTCOMPLEX2 }
#define TESTARRAY_COMPLEX_2D                                            \
    { TESTARRAY_COMPLEX_1D, TESTARRAY_COMPLEX_1D, TESTARRAY_COMPLEX_1D, \
      TESTARRAY_COMPLEX_1D }
#define TESTARRAY_COMPLEX_3D \
    { TESTARRAY_COMPLEX_2D, TESTARRAY_COMPLEX_2D, TESTARRAY_COMPLEX_2D }

#define MAKE_TESTARRAY_INT_1D( NAME ) \
    int NAME[ TESTARRAY_SIZE_DIM1 ] = TESTARRAY_1D
#define MAKE_TESTARRAY_DOUBLE_1D( NAME ) \
    double NAME[ TESTARRAY_SIZE_DIM1 ] = TESTARRAY_1D
#define MAKE_TESTARRAY_INT_2D( NAME ) \
    int NAME[ TESTARRAY_SIZE_DIM2 ][ TESTARRAY_SIZE_DIM1 ] = TESTARRAY_2D
#define MAKE_TESTARRAY_DOUBLE_2D( NAME ) \
    double NAME[ TESTARRAY_SIZE_DIM2 ][ TESTARRAY_SIZE_DIM1 ] = TESTARRAY_2D
#define MAKE_TESTARRAY_INT_3D( NAME )                      \
    int NAME[ TESTARRAY_SIZE_DIM3 ][ TESTARRAY_SIZE_DIM2 ] \
            [ TESTARRAY_SIZE_DIM1 ]                        \
        = TESTARRAY_3D
#define MAKE_TESTARRAY_DOUBLE_3D( NAME )                      \
    double NAME[ TESTARRAY_SIZE_DIM3 ][ TESTARRAY_SIZE_DIM2 ] \
               [ TESTARRAY_SIZE_DIM1 ]                        \
        = TESTARRAY_3D

#define MAKE_COMPLEX_DOUBLE1( NAME ) complex<double> NAME = TESTCOMPLEX1
#define MAKE_COMPLEX_DOUBLE2( NAME ) complex<double> NAME = TESTCOMPLEX2
#define MAKE_TESTARRAY_COMPLEX_1D( NAME ) \
    complex<double> NAME[ TESTARRAY_SIZE_DIM1 ] = TESTARRAY_COMPLEX_1D
#define MAKE_TESTARRAY_COMPLEX_2D( NAME )                              \
    complex<double> NAME[ TESTARRAY_SIZE_DIM2 ][ TESTARRAY_SIZE_DIM1 ] \
        = TESTARRAY_COMPLEX_2D
#define MAKE_TESTARRAY_COMPLEX_3D( NAME )                              \
    complex<double> NAME[ TESTARRAY_SIZE_DIM3 ][ TESTARRAY_SIZE_DIM2 ] \
                        [ TESTARRAY_SIZE_DIM1 ]                        \
        = TESTARRAY_COMPLEX_3D

// functions to perform the tests.  We separate out the functions so that the
// EXPECT_FATAL_FAILURE and EXPECT_NONFATAL_FAILURE macros work correctly.

void multifail( size_t n ) {
    double d1 = 1.0;
    double d2 = d1 + 2.0;
    for (auto i = 0; i < n; i++) {
      EXPECT_DOUBLE_EQ( d1, d2 );
    }
}

void test_assert_complex_eq_pass() {
    MAKE_COMPLEX_DOUBLE1( c1 );
    MAKE_COMPLEX_DOUBLE1( c2 );
    ASSERT_COMPLEX_DOUBLE_EQ( c1, c2 );
}

void test_assert_complex_eq_fail() {
    MAKE_COMPLEX_DOUBLE1( c1 );
    MAKE_COMPLEX_DOUBLE2( c2 );
    ASSERT_COMPLEX_DOUBLE_EQ( c1, c2 );
}

void test_assert_array_eq_int_arrays_pass() {
    MAKE_TESTARRAY_INT_1D( a1 );
    MAKE_TESTARRAY_INT_1D( a2 );
    ASSERT_ARRAY_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array_eq_double_arrays_pass() {
    MAKE_TESTARRAY_DOUBLE_1D( a1 );
    MAKE_TESTARRAY_DOUBLE_1D( a2 );
    ASSERT_ARRAY_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array_eq_int_arrays_fail() {
    MAKE_TESTARRAY_INT_1D( a1 );
    MAKE_TESTARRAY_INT_1D( a2 );
    a2[ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    ASSERT_ARRAY_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array_eq_double_arrays_fail() {
    MAKE_TESTARRAY_DOUBLE_1D( a1 );
    MAKE_TESTARRAY_DOUBLE_1D( a2 );
    a2[ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    ASSERT_ARRAY_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array2d_eq_int_arrays_pass() {
    MAKE_TESTARRAY_INT_2D( a1 );
    MAKE_TESTARRAY_INT_2D( a2 );
    ASSERT_ARRAY2D_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array2d_eq_double_arrays_pass() {
    MAKE_TESTARRAY_DOUBLE_2D( a1 );
    MAKE_TESTARRAY_DOUBLE_2D( a2 );
    ASSERT_ARRAY2D_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1,
                              a2 );
}

void test_assert_array2d_eq_int_arrays_fail() {
    MAKE_TESTARRAY_INT_2D( a1 );
    MAKE_TESTARRAY_INT_2D( a2 );
    a2[ TESTARRAY_SIZE_DIM2 - 1 ][ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    ASSERT_ARRAY2D_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array2d_eq_double_arrays_fail() {
    MAKE_TESTARRAY_DOUBLE_2D( a1 );
    MAKE_TESTARRAY_DOUBLE_2D( a2 );
    a2[ TESTARRAY_SIZE_DIM2 - 1 ][ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    ASSERT_ARRAY2D_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1,
                              a2 );
}

void test_assert_array3d_eq_int_arrays_pass() {
    MAKE_TESTARRAY_INT_3D( a1 );
    MAKE_TESTARRAY_INT_3D( a2 );
    ASSERT_ARRAY3D_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                       TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array3d_eq_double_arrays_pass() {
    MAKE_TESTARRAY_DOUBLE_3D( a1 );
    MAKE_TESTARRAY_DOUBLE_3D( a2 );
    ASSERT_ARRAY3D_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                              TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array3d_eq_int_arrays_fail() {
    MAKE_TESTARRAY_INT_3D( a1 );
    MAKE_TESTARRAY_INT_3D( a2 );
    a2[ TESTARRAY_SIZE_DIM3 - 1 ][ TESTARRAY_SIZE_DIM2 - 1 ]
      [ TESTARRAY_SIZE_DIM1 - 1 ]
        += 3;
    ASSERT_ARRAY3D_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                       TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array3d_eq_double_arrays_fail() {
    MAKE_TESTARRAY_DOUBLE_3D( a1 );
    MAKE_TESTARRAY_DOUBLE_3D( a2 );
    a2[ TESTARRAY_SIZE_DIM3 - 1 ][ TESTARRAY_SIZE_DIM2 - 1 ]
      [ TESTARRAY_SIZE_DIM1 - 1 ]
        += 3;
    ASSERT_ARRAY3D_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                              TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array_eq_complex_arrays_pass() {
    MAKE_TESTARRAY_COMPLEX_1D( a1 );
    MAKE_TESTARRAY_COMPLEX_1D( a2 );
    ASSERT_ARRAY_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array_eq_complex_arrays_fail() {
    MAKE_TESTARRAY_COMPLEX_1D( a1 );
    MAKE_TESTARRAY_COMPLEX_1D( a2 );
    a2[ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    ASSERT_ARRAY_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array2d_eq_complex_arrays_pass() {
    MAKE_TESTARRAY_COMPLEX_2D( a1 );
    MAKE_TESTARRAY_COMPLEX_2D( a2 );
    ASSERT_ARRAY2D_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1,
                                      a1, a2 );
}

void test_assert_array2d_eq_complex_arrays_fail() {
    MAKE_TESTARRAY_COMPLEX_2D( a1 );
    MAKE_TESTARRAY_COMPLEX_2D( a2 );
    a2[ TESTARRAY_SIZE_DIM2 - 1 ][ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    ASSERT_ARRAY2D_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1,
                                      a1, a2 );
}

void test_assert_array3d_eq_complex_arrays_pass() {
    MAKE_TESTARRAY_COMPLEX_3D( a1 );
    MAKE_TESTARRAY_COMPLEX_3D( a2 );
    ASSERT_ARRAY3D_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                                      TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array3d_eq_complex_arrays_fail() {
    MAKE_TESTARRAY_COMPLEX_3D( a1 );
    MAKE_TESTARRAY_COMPLEX_3D( a2 );
    a2[ TESTARRAY_SIZE_DIM3 - 1 ][ TESTARRAY_SIZE_DIM2 - 1 ]
      [ TESTARRAY_SIZE_DIM1 - 1 ]
        += 3;
    ASSERT_ARRAY3D_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                                      TESTARRAY_SIZE_DIM1, a1, a2 );
}

// tests that fail nonfatally using EXPECT
void test_expect_complex_eq_pass() {
    MAKE_COMPLEX_DOUBLE1( c1 );
    MAKE_COMPLEX_DOUBLE1( c2 );
    EXPECT_COMPLEX_DOUBLE_EQ( c1, c2 );
}

void test_expect_complex_eq_fail() {
    MAKE_COMPLEX_DOUBLE1( c1 );
    MAKE_COMPLEX_DOUBLE2( c2 );
    EXPECT_COMPLEX_DOUBLE_EQ( c1, c2 );
}

void test_expect_array_eq_int_arrays_pass() {
    MAKE_TESTARRAY_INT_1D( a1 );
    MAKE_TESTARRAY_INT_1D( a2 );
    EXPECT_ARRAY_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array_eq_double_arrays_pass() {
    MAKE_TESTARRAY_DOUBLE_1D( a1 );
    MAKE_TESTARRAY_DOUBLE_1D( a2 );
    EXPECT_ARRAY_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array_eq_int_arrays_fail() {
    MAKE_TESTARRAY_INT_1D( a1 );
    MAKE_TESTARRAY_INT_1D( a2 );
    a2[ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    EXPECT_ARRAY_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array_eq_double_arrays_fail() {
    MAKE_TESTARRAY_DOUBLE_1D( a1 );
    MAKE_TESTARRAY_DOUBLE_1D( a2 );
    a2[ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    EXPECT_ARRAY_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array2d_eq_int_arrays_pass() {
    MAKE_TESTARRAY_INT_2D( a1 );
    MAKE_TESTARRAY_INT_2D( a2 );
    EXPECT_ARRAY2D_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array2d_eq_double_arrays_pass() {
    MAKE_TESTARRAY_DOUBLE_2D( a1 );
    MAKE_TESTARRAY_DOUBLE_2D( a2 );
    EXPECT_ARRAY2D_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1,
                              a2 );
}

void test_expect_array2d_eq_int_arrays_fail() {
    MAKE_TESTARRAY_INT_2D( a1 );
    MAKE_TESTARRAY_INT_2D( a2 );
    a2[ TESTARRAY_SIZE_DIM2 - 1 ][ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    EXPECT_ARRAY2D_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array2d_eq_double_arrays_fail() {
    MAKE_TESTARRAY_DOUBLE_2D( a1 );
    MAKE_TESTARRAY_DOUBLE_2D( a2 );
    a2[ TESTARRAY_SIZE_DIM2 - 1 ][ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    EXPECT_ARRAY2D_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1,
                              a2 );
}

void test_expect_array3d_eq_int_arrays_pass() {
    MAKE_TESTARRAY_INT_3D( a1 );
    MAKE_TESTARRAY_INT_3D( a2 );
    EXPECT_ARRAY3D_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                       TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array3d_eq_double_arrays_pass() {
    MAKE_TESTARRAY_DOUBLE_3D( a1 );
    MAKE_TESTARRAY_DOUBLE_3D( a2 );
    EXPECT_ARRAY3D_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                              TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array3d_eq_int_arrays_fail() {
    MAKE_TESTARRAY_INT_3D( a1 );
    MAKE_TESTARRAY_INT_3D( a2 );
    a2[ TESTARRAY_SIZE_DIM3 - 1 ][ TESTARRAY_SIZE_DIM2 - 1 ]
      [ TESTARRAY_SIZE_DIM1 - 1 ]
        += 3;
    EXPECT_ARRAY3D_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                       TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array3d_eq_double_arrays_fail() {
    MAKE_TESTARRAY_DOUBLE_3D( a1 );
    MAKE_TESTARRAY_DOUBLE_3D( a2 );
    a2[ TESTARRAY_SIZE_DIM3 - 1 ][ TESTARRAY_SIZE_DIM2 - 1 ]
      [ TESTARRAY_SIZE_DIM1 - 1 ]
        += 3;
    EXPECT_ARRAY3D_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                              TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array_eq_complex_arrays_pass() {
    MAKE_TESTARRAY_COMPLEX_1D( a1 );
    MAKE_TESTARRAY_COMPLEX_1D( a2 );
    for (auto i = 0; i < TESTARRAY_SIZE_DIM1; i++) {
      EXPECT_COMPLEX_DOUBLE_EQ( a1[i], a2[i] );
    }
    EXPECT_ARRAY_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array_eq_complex_arrays_fail() {
    MAKE_TESTARRAY_COMPLEX_1D( a1 );
    MAKE_TESTARRAY_COMPLEX_1D( a2 );
    a2[ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    EXPECT_ARRAY_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array2d_eq_complex_arrays_pass() {
    MAKE_TESTARRAY_COMPLEX_2D( a1 );
    MAKE_TESTARRAY_COMPLEX_2D( a2 );
    EXPECT_ARRAY2D_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1,
                                      a1, a2 );
}

void test_expect_array2d_eq_complex_arrays_fail() {
    MAKE_TESTARRAY_COMPLEX_2D( a1 );
    MAKE_TESTARRAY_COMPLEX_2D( a2 );
    a2[ TESTARRAY_SIZE_DIM2 - 1 ][ TESTARRAY_SIZE_DIM1 - 1 ] += 3;
    EXPECT_ARRAY2D_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1,
                                      a1, a2 );
}

void test_expect_array3d_eq_complex_arrays_pass() {
    MAKE_TESTARRAY_COMPLEX_3D( a1 );
    MAKE_TESTARRAY_COMPLEX_3D( a2 );
    EXPECT_ARRAY3D_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                                      TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array3d_eq_complex_arrays_fail() {
    MAKE_TESTARRAY_COMPLEX_3D( a1 );
    MAKE_TESTARRAY_COMPLEX_3D( a2 );
    a2[ TESTARRAY_SIZE_DIM3 - 1 ][ TESTARRAY_SIZE_DIM2 - 1 ]
      [ TESTARRAY_SIZE_DIM1 - 1 ]
        += 3;
    EXPECT_ARRAY3D_COMPLEX_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2,
                                      TESTARRAY_SIZE_DIM1, a1, a2 );
}

// Tests begin here
TEST( NCPAGtestLibraryTest, ExpectNonfatalFailuresWorksCorrectlyForNFailures ) {
  double d1 = 1.0, d2 = 2.0;
  EXPECT_NONFATAL_FAILURES( { multifail(3); },3);
}

TEST( NCPAGtestLibraryTest, AssertArrayEqPassesOnEqual ) {
    test_assert_array_eq_int_arrays_pass();
    test_assert_array_eq_double_arrays_pass();
}

TEST( NCPAGtestLibraryTest, AssertArrayEqFailsOnUnequal ) {
    EXPECT_FATAL_FAILURE(
        { test_assert_array_eq_int_arrays_fail(); }, "Arrays differ" );
    EXPECT_FATAL_FAILURE(
        { test_assert_array_eq_double_arrays_fail(); }, "Arrays differ" );
}

TEST( NCPAGtestLibraryTest, AssertArray2DEqPassesOnEqual ) {
    test_assert_array2d_eq_int_arrays_pass();
    test_assert_array2d_eq_double_arrays_pass();
}

TEST( NCPAGtestLibraryTest, AssertArray2DEqFailsOnUnequal ) {
    EXPECT_FATAL_FAILURE(
        { test_assert_array2d_eq_int_arrays_fail(); }, "Arrays differ" );
    EXPECT_FATAL_FAILURE(
        { test_assert_array2d_eq_double_arrays_fail(); }, "Arrays differ" );
}

TEST( NCPAGtestLibraryTest, AssertArray3DEqPassesOnEqual ) {
    test_assert_array3d_eq_int_arrays_pass();
    test_assert_array3d_eq_double_arrays_pass();
}

TEST( NCPAGtestLibraryTest, AssertArray3DEqFailsOnUnequal ) {
    EXPECT_FATAL_FAILURE(
        { test_assert_array3d_eq_int_arrays_fail(); }, "Arrays differ" );
    EXPECT_FATAL_FAILURE(
        { test_assert_array3d_eq_double_arrays_fail(); }, "Arrays differ" );
}

TEST( NCPAGtestLibraryTest, AssertComplexEqPassesOnEqual ) {
    test_assert_complex_eq_pass();
}

TEST( NCPAGtestLibraryTest, AssertComplexEqFailsOnUnequal ) {
    EXPECT_FATAL_FAILURE(
        { test_assert_complex_eq_fail(); }, "Expected equality" );
}

TEST( NCPAGtestLibraryTest, AssertComplexArrayEqPassesOnEqual ) {
    test_assert_array_eq_complex_arrays_pass();
}

TEST( NCPAGtestLibraryTest, AssertComplexArrayEqFailsOnUnequal ) {
    EXPECT_FATAL_FAILURE(
        { test_assert_array_eq_complex_arrays_fail(); }, "values differ" );
}

TEST( NCPAGtestLibraryTest, AssertComplex2DArrayEqPassesOnEqual ) {
    test_assert_array2d_eq_complex_arrays_pass();
}

TEST( NCPAGtestLibraryTest, AssertComplex2DArrayEqFailsOnUnequal ) {
    EXPECT_FATAL_FAILURE(
        { test_assert_array2d_eq_complex_arrays_fail(); }, "values differ" );
}

TEST( NCPAGtestLibraryTest, AssertComplex3DArrayEqPassesOnEqual ) {
    test_assert_array3d_eq_complex_arrays_pass();
}

TEST( NCPAGtestLibraryTest, AssertComplex3DArrayEqFailsOnUnequal ) {
    EXPECT_FATAL_FAILURE(
        { test_assert_array3d_eq_complex_arrays_fail(); }, "values differ" );
}




TEST( NCPAGtestLibraryTest, ExpectArrayEqPassesOnEqual ) {
    test_expect_array_eq_int_arrays_pass();
    test_expect_array_eq_double_arrays_pass();
}

TEST( NCPAGtestLibraryTest, ExpectArrayEqFailsOnUnequal ) {
    EXPECT_NONFATAL_FAILURE(
        { test_expect_array_eq_int_arrays_fail(); }, "Arrays differ" );
    EXPECT_NONFATAL_FAILURE(
        { test_expect_array_eq_double_arrays_fail(); }, "Arrays differ" );
}

TEST( NCPAGtestLibraryTest, ExpectArray2DEqPassesOnEqual ) {
    test_expect_array2d_eq_int_arrays_pass();
    test_expect_array2d_eq_double_arrays_pass();
}

TEST( NCPAGtestLibraryTest, ExpectArray2DEqFailsOnUnequal ) {
    EXPECT_NONFATAL_FAILURE(
        { test_expect_array2d_eq_int_arrays_fail(); }, "Arrays differ" );
    EXPECT_NONFATAL_FAILURE(
        { test_expect_array2d_eq_double_arrays_fail(); }, "Arrays differ" );
}

TEST( NCPAGtestLibraryTest, ExpectArray3DEqPassesOnEqual ) {
    test_expect_array3d_eq_int_arrays_pass();
    test_expect_array3d_eq_double_arrays_pass();
}

TEST( NCPAGtestLibraryTest, ExpectArray3DEqFailsOnUnequal ) {
    EXPECT_NONFATAL_FAILURE(
        { test_expect_array3d_eq_int_arrays_fail(); }, "Arrays differ" );
    EXPECT_NONFATAL_FAILURE(
        { test_expect_array3d_eq_double_arrays_fail(); }, "Arrays differ" );
}

TEST( NCPAGtestLibraryTest, ExpectComplexEqPassesOnEqual ) {
    test_expect_complex_eq_pass();
}

TEST( NCPAGtestLibraryTest, ExpectComplexEqFailsOnUnequal ) {
    EXPECT_NONFATAL_FAILURES(
        { test_expect_complex_eq_fail(); }, 2 );
}

TEST( NCPAGtestLibraryTest, ExpectComplexArrayEqPassesOnEqual ) {
    test_expect_array_eq_complex_arrays_pass();
}

TEST( NCPAGtestLibraryTest, ExpectComplexArrayEqFailsOnUnequal ) {
    EXPECT_NONFATAL_FAILURE(
        { test_expect_array_eq_complex_arrays_fail(); }, "values differ" );
}

TEST( NCPAGtestLibraryTest, ExpectComplex2DArrayEqPassesOnEqual ) {
    test_expect_array2d_eq_complex_arrays_pass();
}

TEST( NCPAGtestLibraryTest, ExpectComplex2DArrayEqFailsOnUnequal ) {
    EXPECT_NONFATAL_FAILURE(
        { test_expect_array2d_eq_complex_arrays_fail(); }, "values differ" );
}

TEST( NCPAGtestLibraryTest, ExpectComplex3DArrayEqPassesOnEqual ) {
    test_expect_array3d_eq_complex_arrays_pass();
}

TEST( NCPAGtestLibraryTest, ExpectComplex3DArrayEqFailsOnUnequal ) {
    EXPECT_NONFATAL_FAILURE(
        { test_expect_array3d_eq_complex_arrays_fail(); }, "values differ" );
}