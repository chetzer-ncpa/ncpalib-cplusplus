#include "NCPA/configuration.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <vector>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::config;

typedef double test_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    ConfigurationTest

class ConfigurableVector : public std::vector<double>,
                           public Configurable<std::string> {};

enum class testenum { NO, YES, FIRST_ONLY };

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            sparam = param_ptr_t( new StringParameter( "3.5" ) );
            iparam = param_ptr_t( new IntegerParameter( 4 ) );
            dparam = param_ptr_t( new DoubleParameter( 3.5 ) );
            bparam = param_ptr_t( new BooleanParameter( false ) );
        }  // void TearDown() override {}

        // declare stuff here
        param_ptr_t sparam, iparam, dparam, bparam;

        ConfigurableVector cvec1, cvec2;
};

TEST_F( _TEST_TITLE_, StringParameterReturnsExpectedOutputs ) {
    EXPECT_TRUE( sparam->as_string() == "3.5" );
    EXPECT_EQ( sparam->as_int(), 3 );
    EXPECT_FALSE( sparam->as_bool() );
    EXPECT_NEAR( sparam->as_double(), 3.5, 1e-12 );
}

TEST_F( _TEST_TITLE_, IntegerParameterReturnsExpectedOutputs ) {
    EXPECT_TRUE( iparam->as_string() == "4" );
    EXPECT_EQ( iparam->as_int(), 4 );
    EXPECT_TRUE( iparam->as_bool() );
    EXPECT_NEAR( iparam->as_double(), 4.0, 1e-12 );
}

TEST_F( _TEST_TITLE_, DoubleParameterReturnsExpectedOutputs ) {
    EXPECT_TRUE( dparam->as_string() == "3.500000" );
    EXPECT_EQ( dparam->as_int(), 3 );
    EXPECT_TRUE( dparam->as_bool() );
    EXPECT_NEAR( dparam->as_double(), 3.5, 1e-12 );
}

TEST_F( _TEST_TITLE_, BooleanParameterReturnsExpectedOutputs ) {
    EXPECT_TRUE( bparam->as_string() == "false" );
    EXPECT_EQ( bparam->as_int(), 0 );
    EXPECT_FALSE( bparam->as_bool() );
    EXPECT_NEAR( bparam->as_double(), 0.0, 1e-12 );
}

TEST_F( _TEST_TITLE_, IntegerParameterAllowsTests ) {
    EXPECT_TRUE( iparam->passed() );
    EXPECT_FALSE( iparam->failed() );
    iparam->append_test( IsNotZero<int>() );
    EXPECT_TRUE( iparam->pending() );
    EXPECT_FALSE( iparam->validate().pending() );
    EXPECT_TRUE( iparam->passed() );
    EXPECT_FALSE( iparam->failed() );
    iparam->append_test( IsZero<int>() );
    EXPECT_FALSE( iparam->validate().pending() );
    EXPECT_TRUE( iparam->failed() );
    EXPECT_FALSE( iparam->passed() );
}

TEST_F( _TEST_TITLE_, CanAddTestsAtCreateTime ) {
    DoubleParameter d2( { IsNotZero<double>(), IsLessThan<double>( 5.0 ) } );
    d2.from_double( 3.0 );
    EXPECT_TRUE( d2.validate().passed() );
    d2.from_double( 0.0 );
    EXPECT_FALSE( d2.validate().passed() );
}



TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor1Works ) {
    param_ptr_t suparam( new ScalarParameterWithUnits<double>() );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 0.0 );
    EXPECT_TRUE( suparam->get_units() == nullptr );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor2Works ) {
    param_ptr_t suparam( new ScalarParameterWithUnits<double>( 2.5 ) );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_TRUE( suparam->get_units() == nullptr );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor3Works ) {
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( 2.5, &KILOMETERS ) );
    EXPECT_FALSE( suparam->get_units() == nullptr );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor4Works ) {
    ScalarWithUnits<double> su( 2.5, &KILOMETERS );
    param_ptr_t suparam( new ScalarParameterWithUnits<double>( su ) );
    EXPECT_FALSE( suparam->get_units() == nullptr );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor5Works ) {
    VectorWithUnits<double> su( std::vector<double> { 2.5 }, &KILOMETERS );
    param_ptr_t suparam( new ScalarParameterWithUnits<double>( su ) );
    EXPECT_FALSE( suparam->get_units() == nullptr );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsCloneWorks ) {
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( 2.5, &KILOMETERS ) );
    param_ptr_t suparam2( suparam->clone() );
    EXPECT_FALSE( suparam2->get_units() == nullptr );
    EXPECT_DOUBLE_EQ( suparam2->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam2->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsHasUnitsWorks ) {
    param_ptr_t suparam( new ScalarParameterWithUnits<double>( 2.5 ) );
    EXPECT_FALSE( suparam->has_units() );
    param_ptr_t suparam2(
        new ScalarParameterWithUnits<double>( 2.5, &KILOMETERS ) );
    EXPECT_TRUE( suparam2->has_units() );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsSetUnitsWorks ) {
    param_ptr_t suparam( new ScalarParameterWithUnits<double>( 2.5 ) );
    EXPECT_FALSE( suparam->has_units() );
    suparam->set_units( &KILOMETERS );
    EXPECT_TRUE( suparam->has_units() );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConvertUnitsWorks ) {
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( 2.5, &KILOMETERS ) );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( "mm" ).as_double(), 2500000.0 );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConstructor1Works ) {
    param_ptr_t suparam( new VectorParameterWithUnits<double>() );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 0.0 );
    EXPECT_TRUE( suparam->get_units() == nullptr );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConstructor2Works ) {
    param_ptr_t suparam( new VectorParameterWithUnits<double>( 2.5 ) );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_TRUE( suparam->get_units() == nullptr );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConstructor3Works ) {
    param_ptr_t suparam(
        new VectorParameterWithUnits<double>( 2.5, &KILOMETERS ) );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConstructor4Works ) {
    VectorWithUnits<double> su( std::vector<double> { 2.5, 3.5, 4.5 },
                                &KILOMETERS );
    param_ptr_t suparam( new VectorParameterWithUnits<double>( su ) );
    ASSERT_EQ( suparam->size(), 3 );
    EXPECT_DOUBLE_EQ( suparam->as_double( 0 ), 2.5 );
    EXPECT_DOUBLE_EQ( suparam->as_double( 1 ), 3.5 );
    EXPECT_DOUBLE_EQ( suparam->as_double( 2 ), 4.5 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double( 0 ),
                      2500.0 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double( 1 ),
                      3500.0 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double( 2 ),
                      4500.0 );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsCloneWorks ) {
    param_ptr_t suparam(
        new VectorParameterWithUnits<double>( 2.5, &KILOMETERS ) );
    param_ptr_t suparam2( suparam->clone() );
    EXPECT_DOUBLE_EQ( suparam2->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam2->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsHasUnitsWorks ) {
    param_ptr_t suparam( new VectorParameterWithUnits<double>( 2.5 ) );
    EXPECT_FALSE( suparam->has_units() );
    param_ptr_t suparam2(
        new ScalarParameterWithUnits<double>( 2.5, &KILOMETERS ) );
    EXPECT_TRUE( suparam2->has_units() );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsSetUnitsWorks ) {
    param_ptr_t suparam( new VectorParameterWithUnits<double>( 2.5 ) );
    EXPECT_FALSE( suparam->has_units() );
    suparam->set_units( &KILOMETERS );
    EXPECT_TRUE( suparam->has_units() );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConvertUnitsWorks ) {
    param_ptr_t suparam(
        new VectorParameterWithUnits<double>( 2.5, &KILOMETERS ) );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( "mm" ).as_double(), 2500000.0 );
}

TEST_F( _TEST_TITLE_, ParameterScalarMethodWorksForDouble ) {
    param_ptr_t ptr = Parameter::scalar<double>();
    EXPECT_EQ( ptr->size(), 1 );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 0.0 );
    EXPECT_FALSE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterScalarMethodWorksForDoubleWithUnits ) {
    param_ptr_t ptr = Parameter::scalar<double>( &KILOMETERS );
    EXPECT_EQ( ptr->size(), 1 );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 0.0 );
    EXPECT_TRUE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterScalarMethodWorksForDoubleWithUnitsAndValue ) {
    param_ptr_t ptr = Parameter::scalar<double>( 2.5, &KILOMETERS );
    EXPECT_EQ( ptr->size(), 1 );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 2.5 );
    EXPECT_TRUE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterVectorMethodWorksForDouble ) {
    param_ptr_t ptr = Parameter::vector<double>();
    EXPECT_DOUBLE_EQ( ptr->as_double(), 0.0 );
    EXPECT_EQ( ptr->size(), 0 );
    EXPECT_FALSE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterVectorMethodWorksForDoubleWithUnits ) {
    param_ptr_t ptr = Parameter::vector<double>( &KILOMETERS );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 0.0 );
    EXPECT_EQ( ptr->size(), 0 );
    EXPECT_TRUE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterVectorMethodWorksForDoubleWithUnitsAndValue ) {
    param_ptr_t ptr = Parameter::vector<double>( 2.5, &KILOMETERS );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 2.5 );
    EXPECT_EQ( ptr->size(), 1 );
    EXPECT_TRUE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterVectorMethodWorksForDoubleWithUnitsAndValues ) {
    param_ptr_t ptr
        = Parameter::vector<double>( { 2.5, 3.5, 4.5 }, &KILOMETERS );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( ptr->as_double( 0 ), 2.5 );
    EXPECT_DOUBLE_EQ( ptr->as_double( 1 ), 3.5 );
    EXPECT_DOUBLE_EQ( ptr->as_double( 2 ), 4.5 );
    EXPECT_EQ( ptr->size(), 3 );
    EXPECT_TRUE( ptr->has_units() );
    EXPECT_DOUBLE_EQ( ptr->convert_units(&METERS).as_double( 0 ), 2.5e3 );
    EXPECT_DOUBLE_EQ( ptr->convert_units(&METERS).as_double( 1 ), 3.5e3 );
    EXPECT_DOUBLE_EQ( ptr->convert_units(&METERS).as_double( 2 ), 4.5e3 );
}

TEST_F( _TEST_TITLE_, ConfigurableClassWorksAsExpected ) {
    cvec1.push_back( 1.0 );
    cvec1.push_back( 2.0 );
    cvec1.add_parameter( "par1",
                         DoubleParameter( 1.0, { IsNotZero<double>() } ) );
    EXPECT_EQ( cvec1.get<int>( "par1" ), 1 );
    DoubleParameter par2( 0.0, { IsNotZero<double>() } );
    cvec1.add_parameter( "par2", par2 );
    EXPECT_EQ( cvec1.get<int>( "par2" ), 0 );
}

TEST_F( _TEST_TITLE_, ConfigurableClassWorksWithEnum ) {
    cvec1.add_parameter( "write_atmosphere",
                         ScalarParameter<testenum>( testenum::FIRST_ONLY ) );
    cvec1.add_parameter( "par1",
                         DoubleParameter( 1.0, { IsNotZero<double>() } ) );
    cvec1.validate_parameters();
    EXPECT_TRUE( cvec1.passed() );
    Configurable<std::string> c;
    c.add_parameter( "write_atmosphere",
                     ScalarParameter<testenum>( testenum::FIRST_ONLY ) );
    std::string key = "write_atmosphere";
    testenum te     = c.get<testenum>( key );
    EXPECT_TRUE( cvec1.get<testenum>( "write_atmosphere" )
                 == testenum::FIRST_ONLY );
}

TEST_F( _TEST_TITLE_, ConfigurableClassCopiesParameters ) {
    cvec1.add_parameter( "write_atmosphere",
                         ScalarParameter<testenum>( testenum::FIRST_ONLY ) );
    cvec2.copy_parameters_from( cvec1 );
    EXPECT_TRUE( cvec2.has_parameter( "write_atmosphere" ) );
    EXPECT_TRUE( cvec2.get<testenum>( "write_atmosphere" )
                 == testenum::FIRST_ONLY );
}

TEST_F( _TEST_TITLE_, ConfigurableClassHandlesScalarUnits) {
    cvec1.add_parameter( "speed", Parameter::scalar<double>( 360.0, &METERS_PER_SECOND ) );
    EXPECT_TRUE( cvec1.has_parameter( "speed" ) );
    EXPECT_EQ( cvec1.parameter("speed").size(), 1 );
    EXPECT_DOUBLE_EQ( cvec1.get<double>( "speed" ), 360.0 );
    EXPECT_DOUBLE_EQ( cvec1.convert_parameter( "speed", &KILOMETERS_PER_SECOND ).as_double(), 0.36 );
}

TEST_F( _TEST_TITLE_, ConfigurableClassHandlesVectorUnits) {
    cvec1.add_parameter( "length", Parameter::vector<double>( {-5.0, 2.0, 8.3 }, &KILOMETERS ) );
    EXPECT_TRUE( cvec1.has_parameter( "length" ) );
    EXPECT_EQ( cvec1.parameter("length").size(), 3 );
    EXPECT_DOUBLE_EQ( cvec1.get<double>( "length" ), -5.0 );
    EXPECT_DOUBLE_EQ( cvec1.convert_parameter( "length", &METERS ).as_double(), -5000.0 );
}

TEST_F( _TEST_TITLE_, MappingWorksCorrectly ) {
    cvec1.add_parameter( "par1",
                         DoubleParameter( 1.0, { IsNotZero<double>() } ) );
    Mapping<double>* mapping = new ConfigurationMapping<double,double>( "par1", 
        [](double d){ return d; }, cvec1 );
    mapping->apply( 4.2 );
    EXPECT_DOUBLE_EQ( cvec1.get<double>( "par1" ), 4.2 );
}