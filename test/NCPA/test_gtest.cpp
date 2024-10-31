// #include "NCPA/math.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gtest/gtest-spi.h>
#include "NCPA/gtest.hpp"
#include <vector>
#include <stdexcept>

using namespace std;
using namespace testing;

#define TESTARRAY_SIZE_DIM1 5
#define TESTARRAY_1D { 0, 1, 2, 3, 4 }
#define TESTARRAY_SIZE_DIM2 4
#define TESTARRAY_2D { \
    TESTARRAY_1D, \
    TESTARRAY_1D, \
    TESTARRAY_1D, \
    TESTARRAY_1D, \
}
#define TESTARRAY_SIZE_DIM3 3
#define TESTARRAY_3D { \
    TESTARRAY_2D, \
    TESTARRAY_2D, \
    TESTARRAY_2D, \
}
#define MAKE_TESTARRAY_INT_1D(NAME) int NAME[TESTARRAY_SIZE_DIM1] = TESTARRAY_1D
#define MAKE_TESTARRAY_DOUBLE_1D(NAME) double NAME[TESTARRAY_SIZE_DIM1] = TESTARRAY_1D
#define MAKE_TESTARRAY_INT_2D(NAME) int NAME[TESTARRAY_SIZE_DIM2][TESTARRAY_SIZE_DIM1] = TESTARRAY_2D
#define MAKE_TESTARRAY_DOUBLE_2D(NAME) double NAME[TESTARRAY_SIZE_DIM2][TESTARRAY_SIZE_DIM1] = TESTARRAY_2D
#define MAKE_TESTARRAY_INT_3D(NAME) int NAME[TESTARRAY_SIZE_DIM3][TESTARRAY_SIZE_DIM2][TESTARRAY_SIZE_DIM1] = TESTARRAY_3D
#define MAKE_TESTARRAY_DOUBLE_3D(NAME) double NAME[TESTARRAY_SIZE_DIM3][TESTARRAY_SIZE_DIM2][TESTARRAY_SIZE_DIM1] = TESTARRAY_3D


// functions to perform the tests.  We separate out the functions so that the 
// EXPECT_FATAL_FAILURE and EXPECT_NONFATAL_FAILURE macros work correctly.

void test_assert_array_eq_int_arrays_pass() {
  MAKE_TESTARRAY_INT_1D(a1);
  MAKE_TESTARRAY_INT_1D(a2);
  ASSERT_ARRAY_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_assert_array_eq_double_arrays_pass() {
  MAKE_TESTARRAY_DOUBLE_1D(a1);
  MAKE_TESTARRAY_DOUBLE_1D(a2);
  ASSERT_ARRAY_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_assert_array_eq_int_arrays_fail() {
  MAKE_TESTARRAY_INT_1D(a1);
  MAKE_TESTARRAY_INT_1D(a2);
  a2[TESTARRAY_SIZE_DIM1-1] += 3;
  ASSERT_ARRAY_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_assert_array_eq_double_arrays_fail() {
  MAKE_TESTARRAY_DOUBLE_1D(a1);
  MAKE_TESTARRAY_DOUBLE_1D(a2);
  a2[TESTARRAY_SIZE_DIM1-1] += 3;
  ASSERT_ARRAY_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array2d_eq_int_arrays_pass() {
  MAKE_TESTARRAY_INT_2D(a1);
  MAKE_TESTARRAY_INT_2D(a2);
  ASSERT_ARRAY2D_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_assert_array2d_eq_double_arrays_pass() {
  MAKE_TESTARRAY_DOUBLE_2D(a1);
  MAKE_TESTARRAY_DOUBLE_2D(a2);
  ASSERT_ARRAY2D_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array2d_eq_int_arrays_fail() {
  MAKE_TESTARRAY_INT_2D(a1);
  MAKE_TESTARRAY_INT_2D(a2);
  a2[TESTARRAY_SIZE_DIM2-1][TESTARRAY_SIZE_DIM1-1] += 3;
  ASSERT_ARRAY2D_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_assert_array2d_eq_double_arrays_fail() {
  MAKE_TESTARRAY_DOUBLE_2D(a1);
  MAKE_TESTARRAY_DOUBLE_2D(a2);
  a2[TESTARRAY_SIZE_DIM2-1][TESTARRAY_SIZE_DIM1-1] += 3;
  ASSERT_ARRAY2D_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array3d_eq_int_arrays_pass() {
  MAKE_TESTARRAY_INT_3D(a1);
  MAKE_TESTARRAY_INT_3D(a2);
  ASSERT_ARRAY3D_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_assert_array3d_eq_double_arrays_pass() {
  MAKE_TESTARRAY_DOUBLE_3D(a1);
  MAKE_TESTARRAY_DOUBLE_3D(a2);
  ASSERT_ARRAY3D_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_assert_array3d_eq_int_arrays_fail() {
  MAKE_TESTARRAY_INT_3D(a1);
  MAKE_TESTARRAY_INT_3D(a2);
  a2[TESTARRAY_SIZE_DIM3-1][TESTARRAY_SIZE_DIM2-1][TESTARRAY_SIZE_DIM1-1] += 3;
  ASSERT_ARRAY3D_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_assert_array3d_eq_double_arrays_fail() {
  MAKE_TESTARRAY_DOUBLE_3D(a1);
  MAKE_TESTARRAY_DOUBLE_3D(a2);
  a2[TESTARRAY_SIZE_DIM3-1][TESTARRAY_SIZE_DIM2-1][TESTARRAY_SIZE_DIM1-1] += 3;
  ASSERT_ARRAY3D_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array_eq_int_arrays_pass() {
  MAKE_TESTARRAY_INT_1D(a1);
  MAKE_TESTARRAY_INT_1D(a2);
  EXPECT_ARRAY_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_expect_array_eq_double_arrays_pass() {
  MAKE_TESTARRAY_DOUBLE_1D(a1);
  MAKE_TESTARRAY_DOUBLE_1D(a2);
  EXPECT_ARRAY_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_expect_array_eq_int_arrays_fail() {
  MAKE_TESTARRAY_INT_1D(a1);
  MAKE_TESTARRAY_INT_1D(a2);
  a2[TESTARRAY_SIZE_DIM1-1] += 3;
  EXPECT_ARRAY_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_expect_array_eq_double_arrays_fail() {
  MAKE_TESTARRAY_DOUBLE_1D(a1);
  MAKE_TESTARRAY_DOUBLE_1D(a2);
  a2[TESTARRAY_SIZE_DIM1-1] += 3;
  EXPECT_ARRAY_DOUBLE_EQ( TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array2d_eq_int_arrays_pass() {
  MAKE_TESTARRAY_INT_2D(a1);
  MAKE_TESTARRAY_INT_2D(a2);
  EXPECT_ARRAY2D_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_expect_array2d_eq_double_arrays_pass() {
  MAKE_TESTARRAY_DOUBLE_2D(a1);
  MAKE_TESTARRAY_DOUBLE_2D(a2);
  EXPECT_ARRAY2D_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array2d_eq_int_arrays_fail() {
  MAKE_TESTARRAY_INT_2D(a1);
  MAKE_TESTARRAY_INT_2D(a2);
  a2[TESTARRAY_SIZE_DIM2-1][TESTARRAY_SIZE_DIM1-1] += 3;
  EXPECT_ARRAY2D_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_expect_array2d_eq_double_arrays_fail() {
  MAKE_TESTARRAY_DOUBLE_2D(a1);
  MAKE_TESTARRAY_DOUBLE_2D(a2);
  a2[TESTARRAY_SIZE_DIM2-1][TESTARRAY_SIZE_DIM1-1] += 3;
  EXPECT_ARRAY2D_DOUBLE_EQ( TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array3d_eq_int_arrays_pass() {
  MAKE_TESTARRAY_INT_3D(a1);
  MAKE_TESTARRAY_INT_3D(a2);
  EXPECT_ARRAY3D_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_expect_array3d_eq_double_arrays_pass() {
  MAKE_TESTARRAY_DOUBLE_3D(a1);
  MAKE_TESTARRAY_DOUBLE_3D(a2);
  EXPECT_ARRAY3D_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}

void test_expect_array3d_eq_int_arrays_fail() {
  MAKE_TESTARRAY_INT_3D(a1);
  MAKE_TESTARRAY_INT_3D(a2);
  a2[TESTARRAY_SIZE_DIM3-1][TESTARRAY_SIZE_DIM2-1][TESTARRAY_SIZE_DIM1-1] += 3;
  EXPECT_ARRAY3D_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}
void test_expect_array3d_eq_double_arrays_fail() {
  MAKE_TESTARRAY_DOUBLE_3D(a1);
  MAKE_TESTARRAY_DOUBLE_3D(a2);
  a2[TESTARRAY_SIZE_DIM3-1][TESTARRAY_SIZE_DIM2-1][TESTARRAY_SIZE_DIM1-1] += 3;
  EXPECT_ARRAY3D_DOUBLE_EQ( TESTARRAY_SIZE_DIM3, TESTARRAY_SIZE_DIM2, TESTARRAY_SIZE_DIM1, a1, a2 );
}


// Tests begin here

TEST(NCPAgtestTest,AssertArrayEqPassesOnEqual) {
    test_assert_array_eq_int_arrays_pass();
    test_assert_array_eq_double_arrays_pass();
}

TEST(NCPAgtestTest,AssertArrayEqFailsOnUnequal) {
  EXPECT_FATAL_FAILURE({test_assert_array_eq_int_arrays_fail();}, "Arrays differ" );
  EXPECT_FATAL_FAILURE({test_assert_array_eq_double_arrays_fail();}, "Arrays differ" );
}

TEST(NCPAgtestTest,AssertArray2DEqPassesOnEqual) {
    test_assert_array2d_eq_int_arrays_pass();
    test_assert_array2d_eq_double_arrays_pass();
}

TEST(NCPAgtestTest,AssertArray2DEqFailsOnUnequal) {
  EXPECT_FATAL_FAILURE({test_assert_array2d_eq_int_arrays_fail();}, "Arrays differ" );
  EXPECT_FATAL_FAILURE({test_assert_array2d_eq_double_arrays_fail();}, "Arrays differ" );
}


TEST(NCPAgtestTest,AssertArray3DEqPassesOnEqual) {
    test_assert_array3d_eq_int_arrays_pass();
    test_assert_array3d_eq_double_arrays_pass();
}

TEST(NCPAgtestTest,AssertArray3DEqFailsOnUnequal) {
  EXPECT_FATAL_FAILURE({test_assert_array3d_eq_int_arrays_fail();}, "Arrays differ" );
  EXPECT_FATAL_FAILURE({test_assert_array3d_eq_double_arrays_fail();}, "Arrays differ" );
}

TEST(NCPAgtestTest,ExpectArrayEqPassesOnEqual) {
    test_expect_array_eq_int_arrays_pass();
    test_expect_array_eq_double_arrays_pass();
}

TEST(NCPAgtestTest,ExpectArrayEqFailsOnUnequal) {
  EXPECT_NONFATAL_FAILURE({test_expect_array_eq_int_arrays_fail();}, "Arrays differ" );
  EXPECT_NONFATAL_FAILURE({test_expect_array_eq_double_arrays_fail();}, "Arrays differ" );
}

TEST(NCPAgtestTest,ExpectArray2DEqPassesOnEqual) {
    test_expect_array2d_eq_int_arrays_pass();
    test_expect_array2d_eq_double_arrays_pass();
}

TEST(NCPAgtestTest,ExpectArray2DEqFailsOnUnequal) {
  EXPECT_NONFATAL_FAILURE({test_expect_array2d_eq_int_arrays_fail();}, "Arrays differ" );
  EXPECT_NONFATAL_FAILURE({test_expect_array2d_eq_double_arrays_fail();}, "Arrays differ" );
}


TEST(NCPAgtestTest,ExpectArray3DEqPassesOnEqual) {
    test_expect_array3d_eq_int_arrays_pass();
    test_expect_array3d_eq_double_arrays_pass();
}

TEST(NCPAgtestTest,ExpectArray3DEqFailsOnUnequal) {
  EXPECT_NONFATAL_FAILURE({test_expect_array3d_eq_int_arrays_fail();}, "Arrays differ" );
  EXPECT_NONFATAL_FAILURE({test_expect_array3d_eq_double_arrays_fail();}, "Arrays differ" );
}

