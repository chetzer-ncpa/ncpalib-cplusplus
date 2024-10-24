#include "NCPA/math.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gtest/gtest-spi.h>
#include "NCPA/gtest.hpp"

void test_assert_array_eq_int_arrays_pass() {
  int a1[ 5 ] = { 0, 1, 2, 3, 4 };
  int a2[ 5 ] = { 0, 1, 2, 3, 4 };
  ASSERT_ARRAY_EQ( 5, a1, a2 );
}
void test_assert_array_eq_double_arrays_pass() {
  double a1[ 5 ] = { 0, 1, 2, 3, 4 };
  double a2[ 5 ] = { 0, 1, 2, 3, 4 };
  ASSERT_DOUBLE_ARRAY_EQ( 5, a1, a2 );
}
void test_assert_array_eq_int_arrays_fail() {
  int a1[ 5 ] = { 0, 1, 2, 3, 4 };
  int a2[ 5 ] = { 0, 1, 5, 3, 4 };
  ASSERT_ARRAY_EQ( 5, a1, a2 );
}
void test_assert_array_eq_double_arrays_fail() {
  double a1[ 5 ] = { 0, 1, 2, 3, 4 };
  double a2[ 5 ] = { 0, 1, 5, 3, 4 };
  ASSERT_DOUBLE_ARRAY_EQ( 5, a1, a2 );
}

TEST(NCPAgtestTest,ExpectArrayEqPassesOnEqual) {
    test_assert_array_eq_int_arrays_pass();
    test_assert_array_eq_double_arrays_pass();
}

TEST(NCPAgtestTest,ExpectArrayEqFailsOnUnequal) {
  EXPECT_FATAL_FAILURE({test_assert_array_eq_int_arrays_fail();}, "Arrays differ" );
  EXPECT_FATAL_FAILURE({test_assert_array_eq_double_arrays_fail();}, "Arrays differ" );
}

class ThrowListener : public testing::EmptyTestEventListener {
  void OnTestPartResult(const testing::TestPartResult& result) override {
    if (result.type() == testing::TestPartResult::kFatalFailure) {
      throw testing::AssertionException(result);
    }
  }
};
int main(int argc, char** argv) {
  testing::InitGoogleTest();
  testing::UnitTest::GetInstance()->listeners().Append(new ThrowListener);
  return RUN_ALL_TESTS();
}

