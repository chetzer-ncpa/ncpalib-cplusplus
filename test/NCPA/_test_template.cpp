#include "NCPA/math.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;

typedef double test_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    TestSuiteName

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            
        }  // void TearDown() override {}

        // declare stuff here
        const test_t zero = NCPA::math::zero<test_t>(),
                     one  = NCPA::math::one<test_t>();
};