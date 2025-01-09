#include "NCPA/parameters.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::params;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    ParametersTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            intParamDefault = Parameter<int>( 1 );
            doubleParamDefault = Parameter<double>( 1.0 );
            stringParamDefault = Parameter<string>( "1" );
        }  // void TearDown() override {}

        // declare stuff here
        Parameter<int> intParam;
        Parameter<double> doubleParam;
        Parameter<string> stringParam;
        Parameter<int> intParamDefault;
        Parameter<double> doubleParamDefault;
        Parameter<string> stringParamDefault;
};

TEST_F( _TEST_TITLE_, DefaultConstructorsSetInternals ) {
    EXPECT_TRUE( intParam.required() );
    EXPECT_FALSE( intParamDefault.required() );
    EXPECT_TRUE( doubleParam.required() );
    EXPECT_FALSE( doubleParamDefault.required() );
    EXPECT_TRUE( stringParam.required() );
    EXPECT_FALSE( stringParamDefault.required() );
}