#include "NCPA/configuration.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <vector>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::config;
using namespace NCPA::config::validation;

typedef double test_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    TestSuiteName

class ConfigurableVector
    : public std::vector<double>,
      public Configurable<ConfigurableVector, std::string> {};

enum class testenum { NO, YES, FIRST_ONLY };

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            sparam.set( "3.5" );
            iparam.set( 4 );
            dparam.set( 3.5 );
            bparam.set( false );
        }  // void TearDown() override {}

        // declare stuff here
        StringParameter sparam;
        IntegerParameter iparam;
        DoubleParameter dparam;
        BooleanParameter bparam;

        ConfigurableVector cvec1, cvec2;
};

TEST_F( _TEST_TITLE_, StringParameterReturnsExpectedOutputs ) {
    EXPECT_TRUE( sparam.as_string() == "3.5" );
    EXPECT_EQ( sparam.as_int(), 3 );
    EXPECT_TRUE( sparam.as_bool() );
    EXPECT_NEAR( sparam.as_double(), 3.5, 1e-12 );
}

TEST_F( _TEST_TITLE_, IntegerParameterReturnsExpectedOutputs ) {
    EXPECT_TRUE( iparam.as_string() == "4" );
    EXPECT_EQ( iparam.as_int(), 4 );
    EXPECT_TRUE( iparam.as_bool() );
    EXPECT_NEAR( iparam.as_double(), 4.0, 1e-12 );
}

TEST_F( _TEST_TITLE_, DoubleParameterReturnsExpectedOutputs ) {
    EXPECT_TRUE( dparam.as_string() == "3.5" );
    EXPECT_EQ( dparam.as_int(), 3 );
    EXPECT_TRUE( dparam.as_bool() );
    EXPECT_NEAR( dparam.as_double(), 3.5, 1e-12 );
}

TEST_F( _TEST_TITLE_, BooleanParameterReturnsExpectedOutputs ) {
    EXPECT_TRUE( bparam.as_string() == "false" );
    EXPECT_EQ( bparam.as_int(), 0 );
    EXPECT_FALSE( bparam.as_bool() );
    EXPECT_NEAR( bparam.as_double(), 0.0, 1e-12 );
}

TEST_F( _TEST_TITLE_, IntegerParameterAllowsTests ) {
    EXPECT_TRUE( iparam.passed() );
    EXPECT_FALSE( iparam.failed() );
    iparam.append_test( IsNotZero<int>() );
    EXPECT_TRUE( iparam.pending() );
    EXPECT_FALSE( iparam.validate().pending() );
    EXPECT_TRUE( iparam.passed() );
    EXPECT_FALSE( iparam.failed() );
    iparam.append_test( IsZero<int>() );
    EXPECT_FALSE( iparam.validate().pending() );
    EXPECT_TRUE( iparam.failed() );
    EXPECT_FALSE( iparam.passed() );
}

TEST_F( _TEST_TITLE_, CanAddTestsAtCreateTime ) {
    DoubleParameter d2( { IsNotZero<double>(), IsLessThan<double>( 5.0 ) } );
    d2.set( 3.0 );
    EXPECT_TRUE( d2.validate().passed() );
    d2.set( 0.0 );
    EXPECT_FALSE( d2.validate().passed() );
}

TEST_F( _TEST_TITLE_, ConfigurableClassWorksAsExpected ) {
    cvec1.push_back( 1.0 );
    cvec1.push_back( 2.0 );
    cvec1.add_parameter( "par1",
                         DoubleParameter( 1.0, { IsNotZero<double>() } ) );
    EXPECT_TRUE( cvec1.validate_parameters().passed() );
    cvec1.add_parameter( "par2",
                         DoubleParameter( 0.0, { IsNotZero<double>() } ) );
    EXPECT_FALSE( cvec1.validate_parameters().passed() );
}

TEST_F( _TEST_TITLE_, ConfigurableClassWorksWithEnum ) {
    cvec1.add_parameter( "write_atmosphere",
                         Parameter<testenum>( testenum::FIRST_ONLY ) );
    cvec1.add_parameter( "par1",
                         DoubleParameter( 1.0, { IsNotZero<double>() } ) );
    EXPECT_TRUE( cvec1.validate_parameters().passed() );
    EXPECT_TRUE( cvec1.get<testenum>( "write_atmosphere" )
                 == testenum::FIRST_ONLY );
}

TEST_F( _TEST_TITLE_, ConfigurableClassCopiesParameters ) {
    cvec1.add_parameter( "write_atmosphere",
                         Parameter<testenum>( testenum::FIRST_ONLY ) );
    cvec2.copy_parameters( cvec1.parameters() );
    EXPECT_TRUE( cvec2.has_parameter( "write_atmosphere" ) );
    EXPECT_TRUE( cvec2.get<testenum>( "write_atmosphere" ) == testenum::FIRST_ONLY );
}