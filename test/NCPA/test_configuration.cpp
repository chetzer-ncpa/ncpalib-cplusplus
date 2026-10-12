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
            sparam = param_ptr_t( new StringParameter( "sparam", "3.5" ) );
            iparam = param_ptr_t( new IntegerParameter( "iparam", 4 ) );
            dparam = param_ptr_t( new DoubleParameter( "dparam", 3.5 ) );
            bparam = param_ptr_t( new BooleanParameter( "bparam", false ) );
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
    iparam->add_test( is_nonzero<int>() );
    EXPECT_TRUE( iparam->passed() );
    EXPECT_FALSE( iparam->failed() );
    iparam->add_test( is_positive<int>() );
    EXPECT_TRUE( iparam->passed() );
    EXPECT_FALSE( iparam->failed() );
    iparam->add_test( is_greater_than<int>( 5 ) );
    EXPECT_TRUE( iparam->failed() );
    EXPECT_FALSE( iparam->passed() );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor1Works ) {
    param_ptr_t suparam( new ScalarParameterWithUnits<double>() );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 0.0 );
    EXPECT_TRUE( suparam->get_units() == nullptr );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor2Works ) {
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( "suparam", 2.5 ) );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_TRUE( suparam->get_units() == nullptr );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor3Works ) {
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( "suparam", 2.5, &KILOMETERS ) );
    EXPECT_FALSE( suparam->get_units() == nullptr );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor4Works ) {
    ScalarWithUnits<double> su( 2.5, &KILOMETERS );
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( "suparam", su ) );
    EXPECT_FALSE( suparam->get_units() == nullptr );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConstructor5Works ) {
    VectorWithUnits<double> su( std::vector<double> { 2.5 }, &KILOMETERS );
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( "suparam", su ) );
    EXPECT_FALSE( suparam->get_units() == nullptr );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsCloneWorks ) {
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( "suparam", 2.5, &KILOMETERS ) );
    param_ptr_t suparam2( suparam->clone() );
    EXPECT_FALSE( suparam2->get_units() == nullptr );
    EXPECT_DOUBLE_EQ( suparam2->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam2->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsHasUnitsWorks ) {
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( "suparam", 2.5 ) );
    EXPECT_FALSE( suparam->has_units() );
    param_ptr_t suparam2(
        new ScalarParameterWithUnits<double>( "suparam", 2.5, &KILOMETERS ) );
    EXPECT_TRUE( suparam2->has_units() );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsSetUnitsWorks ) {
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( "suparam", 2.5 ) );
    EXPECT_FALSE( suparam->has_units() );
    suparam->set_units( &KILOMETERS );
    EXPECT_TRUE( suparam->has_units() );
}

TEST_F( _TEST_TITLE_, ScalarParameterWithUnitsConvertUnitsWorks ) {
    param_ptr_t suparam(
        new ScalarParameterWithUnits<double>( "suparam", 2.5, &KILOMETERS ) );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( "mm" ).as_double(), 2500000.0 );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConstructor1Works ) {
    param_ptr_t suparam( new VectorParameterWithUnits<double>( "suparam" ) );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 0.0 );
    EXPECT_TRUE( suparam->get_units() == nullptr );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConstructor2Works ) {
    param_ptr_t suparam(
        new VectorParameterWithUnits<double>( "suparam", 2.5 ) );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_TRUE( suparam->get_units() == nullptr );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConstructor3Works ) {
    param_ptr_t suparam(
        new VectorParameterWithUnits<double>( "suparam", 2.5, &KILOMETERS ) );
    EXPECT_DOUBLE_EQ( suparam->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConstructor4Works ) {
    VectorWithUnits<double> su( std::vector<double> { 2.5, 3.5, 4.5 },
                                &KILOMETERS );
    param_ptr_t suparam(
        new VectorParameterWithUnits<double>( "suparam", su ) );
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
        new VectorParameterWithUnits<double>( "suparam", 2.5, &KILOMETERS ) );
    param_ptr_t suparam2( suparam->clone() );
    EXPECT_DOUBLE_EQ( suparam2->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( suparam2->convert_units( &METERS ).as_double(), 2500.0 );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsHasUnitsWorks ) {
    param_ptr_t suparam(
        new VectorParameterWithUnits<double>( "suparam", 2.5 ) );
    EXPECT_FALSE( suparam->has_units() );
    param_ptr_t suparam2(
        new ScalarParameterWithUnits<double>( "suparam", 2.5, &KILOMETERS ) );
    EXPECT_TRUE( suparam2->has_units() );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsSetUnitsWorks ) {
    param_ptr_t suparam(
        new VectorParameterWithUnits<double>( "suparam", 2.5 ) );
    EXPECT_FALSE( suparam->has_units() );
    suparam->set_units( &KILOMETERS );
    EXPECT_TRUE( suparam->has_units() );
}

TEST_F( _TEST_TITLE_, VectorParameterWithUnitsConvertUnitsWorks ) {
    param_ptr_t suparam(
        new VectorParameterWithUnits<double>( "suparam", 2.5, &KILOMETERS ) );
    EXPECT_DOUBLE_EQ( suparam->convert_units( &METERS ).as_double(), 2500.0 );
    EXPECT_DOUBLE_EQ( suparam->convert_units( "mm" ).as_double(), 2500000.0 );
}

TEST_F( _TEST_TITLE_, ParameterScalarMethodWorksForDouble ) {
    param_ptr_t ptr = Parameter::scalar<double>( "suparam" );
    EXPECT_EQ( ptr->size(), 1 );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 0.0 );
    EXPECT_FALSE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterScalarMethodWorksForDoubleWithUnits ) {
    param_ptr_t ptr = Parameter::scalar<double>( "suparam", &KILOMETERS );
    EXPECT_EQ( ptr->size(), 1 );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 0.0 );
    EXPECT_TRUE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterScalarMethodWorksForDoubleWithUnitsAndValue ) {
    param_ptr_t ptr = Parameter::scalar<double>( "suparam", 2.5, &KILOMETERS );
    EXPECT_EQ( ptr->size(), 1 );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 2.5 );
    EXPECT_TRUE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterVectorMethodWorksForDouble ) {
    param_ptr_t ptr = Parameter::vector<double>( "suparam" );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 0.0 );
    EXPECT_EQ( ptr->size(), 0 );
    EXPECT_FALSE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterVectorMethodWorksForDoubleWithUnits ) {
    param_ptr_t ptr = Parameter::vector<double>( "suparam", &KILOMETERS );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 0.0 );
    EXPECT_EQ( ptr->size(), 0 );
    EXPECT_TRUE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterVectorMethodWorksForDoubleWithUnitsAndValue ) {
    param_ptr_t ptr = Parameter::vector<double>( "suparam", 2.5, &KILOMETERS );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 2.5 );
    EXPECT_EQ( ptr->size(), 1 );
    EXPECT_TRUE( ptr->has_units() );
}

TEST_F( _TEST_TITLE_, ParameterVectorMethodWorksForDoubleWithUnitsAndValues ) {
    param_ptr_t ptr = Parameter::vector<double>( "suparam", { 2.5, 3.5, 4.5 },
                                                 &KILOMETERS );
    EXPECT_DOUBLE_EQ( ptr->as_double(), 2.5 );
    EXPECT_DOUBLE_EQ( ptr->as_double( 0 ), 2.5 );
    EXPECT_DOUBLE_EQ( ptr->as_double( 1 ), 3.5 );
    EXPECT_DOUBLE_EQ( ptr->as_double( 2 ), 4.5 );
    EXPECT_EQ( ptr->size(), 3 );
    EXPECT_TRUE( ptr->has_units() );
    EXPECT_DOUBLE_EQ( ptr->convert_units( &METERS ).as_double( 0 ), 2.5e3 );
    EXPECT_DOUBLE_EQ( ptr->convert_units( &METERS ).as_double( 1 ), 3.5e3 );
    EXPECT_DOUBLE_EQ( ptr->convert_units( &METERS ).as_double( 2 ), 4.5e3 );
}

TEST_F( _TEST_TITLE_, ConfigurableClassWorksAsExpected ) {
    cvec1.push_back( 1.0 );
    cvec1.push_back( 2.0 );
    cvec1.add_parameter( DoubleParameter( "par1", 1.0 ) );
    EXPECT_EQ( cvec1.get<int>( "par1" ), 1 );
    DoubleParameter par2( "par2", 0.0 );
    cvec1.add_parameter( par2 );
    EXPECT_EQ( cvec1.get<int>( "par2" ), 0 );
}

TEST_F( _TEST_TITLE_, ConfigurableClassWorksWithEnum ) {
    cvec1.add_parameter( ScalarParameter<testenum>( "write_atmosphere",
                                                    testenum::FIRST_ONLY ) );
    cvec1.add_parameter( DoubleParameter( "par1", 1.0 ) );
    cvec1.validate_parameters();
    EXPECT_TRUE( cvec1.passed() );
    Configurable<std::string> c;
    c.add_parameter( ScalarParameter<testenum>( "write_atmosphere",
                                                testenum::FIRST_ONLY ) );
    std::string key = "write_atmosphere";
    testenum te     = c.get<testenum>( key );
    EXPECT_TRUE( cvec1.get<testenum>( "write_atmosphere" )
                 == testenum::FIRST_ONLY );
}

TEST_F( _TEST_TITLE_, ConfigurableClassCopiesParameters ) {
    cvec1.add_parameter( ScalarParameter<testenum>( "write_atmosphere",
                                                    testenum::FIRST_ONLY ) );
    cvec2.copy_parameters_from( cvec1 );
    EXPECT_TRUE( cvec2.has_parameter( "write_atmosphere" ) );
    EXPECT_TRUE( cvec2.get<testenum>( "write_atmosphere" )
                 == testenum::FIRST_ONLY );
}

TEST_F( _TEST_TITLE_, ConfigurableClassHandlesScalarUnits ) {
    cvec1.add_parameter(
        Parameter::scalar<double>( "speed", 360.0, &METERS_PER_SECOND ) );
    EXPECT_TRUE( cvec1.has_parameter( "speed" ) );
    EXPECT_EQ( cvec1.parameter( "speed" ).size(), 1 );
    EXPECT_DOUBLE_EQ( cvec1.get<double>( "speed" ), 360.0 );
    EXPECT_DOUBLE_EQ(
        cvec1.convert_parameter( "speed", &KILOMETERS_PER_SECOND ).as_double(),
        0.36 );
}

TEST_F( _TEST_TITLE_, ConfigurableClassHandlesVectorUnits ) {
    cvec1.add_parameter( Parameter::vector<double>(
        "length", { -5.0, 2.0, 8.3 }, &KILOMETERS ) );
    EXPECT_TRUE( cvec1.has_parameter( "length" ) );
    EXPECT_EQ( cvec1.parameter( "length" ).size(), 3 );
    EXPECT_DOUBLE_EQ( cvec1.get<double>( "length" ), -5.0 );
    EXPECT_DOUBLE_EQ( cvec1.convert_parameter( "length", &METERS ).as_double(),
                      -5000.0 );
}

// TEST_F( _TEST_TITLE_, MappingWorksCorrectly ) {
//     cvec1.add_parameter( "par1", DoubleParameter( 1.0 ) );
//     Mapping<double> *mapping = new ConfigurationMapping<double, double>(
//         "par1", []( double d ) { return d; }, cvec1 );
//     mapping->apply( 4.2 );
//     EXPECT_DOUBLE_EQ( cvec1.get<double>( "par1" ), 4.2 );
// }

// TEST_F( _TEST_TITLE_, MappingCloneWorksCorrectly ) {
//     cvec1.add_parameter( "par1", DoubleParameter( 1.0 ) );
//     Mapping<double> *mapping = new ConfigurationMapping<double, double>(
//         "par1", []( double d ) { return d * 2.0; }, cvec1 );
//     mapping_ptr_t<double> ptr( mapping->clone() );
//     mapping->apply( 4.2 );
//     EXPECT_DOUBLE_EQ( cvec1.get<double>( "par1" ), 8.4 );
// }

TEST_F( _TEST_TITLE_, TrueValidationValidateReturnsCorrect ) {
    dparam->add_test( is_greater_than<double>( 2.0 ) );
    EXPECT_EQ( ( dparam->validate().result ), test_result_t::PASSED );
}

TEST_F( _TEST_TITLE_, TrueValidationPassedReturnsTrue ) {
    dparam->add_test( is_greater_than<double>( 2.0 ) );
    EXPECT_TRUE( dparam->passed() );
}

TEST_F( _TEST_TITLE_, FalseValidationValidateReturnsCorrect ) {
    dparam->add_test( is_greater_than<double>( 6.0 ) );
    EXPECT_EQ( ( dparam->validate().result ), test_result_t::FAILED );
}

TEST_F( _TEST_TITLE_, FalseValidationPassedReturnsFalse ) {
    dparam->add_test( is_greater_than<double>( 6.0 ) );
    EXPECT_TRUE( dparam->failed() );
}

TEST_F( _TEST_TITLE_, MultipleTrueValidationsReturnCorrect ) {
    DoubleParameter dparam2(
        "dparam2", 4.5,
        vector<TypedValidation<double>> { is_greater_than<double>( 2.0 ),
                                          is_less_than<double>( 6.0 ) } );
    EXPECT_TRUE( dparam2.passed() );
}

TEST_F( _TEST_TITLE_, FlagArgumentSetsProperly ) {
    FlagArgument flag( "testflag", "Test flag for unit tests" );
    EXPECT_FALSE( flag.was_set() );
    flag.set_value( true );
    EXPECT_TRUE( flag.was_set() );
}

TEST_F( _TEST_TITLE_, FlagArgumentThrowsIfSettingVector ) {
    FlagArgument flag( "testflag", "Test flag for unit tests" );
    EXPECT_THROW(
        { flag.set_value( std::vector<bool> { true, true, true } ); },
        TooManyValues );
}

TEST_F( _TEST_TITLE_, SingleArgumentSetsProperly ) {
    SingleValueArgument<double> arg( "testarg", "Test flag for unit tests" );
    EXPECT_FALSE( arg.was_set() );
    arg.set_value( 3.0 );
    EXPECT_TRUE( arg.was_set() );
    EXPECT_DOUBLE_EQ( arg.as_double(), 3.0 );
    std::vector<double> tv = arg.as_doubles();
    EXPECT_EQ( tv.size(), 1 );
    EXPECT_DOUBLE_EQ( tv[ 0 ], 3.0 );
}

TEST_F( _TEST_TITLE_, SingleArgumentThrowsIfSettingVector ) {
    SingleValueArgument<int> arg( "testarg", "Test flag for unit tests" );
    EXPECT_THROW(
        { arg.set_value( std::vector<int> { 1, 2, 3 } ); }, TooManyValues );
}

TEST_F( _TEST_TITLE_, MultiArgumentSetsProperly ) {
    MultiValueArgument<double> arg( "testarg", "Test flag for unit tests" );
    EXPECT_FALSE( arg.was_set() );
    arg.set_value( 3.0 );
    EXPECT_TRUE( arg.was_set() );
    EXPECT_DOUBLE_EQ( arg.as_double(), 3.0 );
    std::vector<double> tv = arg.as_doubles();
    EXPECT_EQ( tv.size(), 1 );
    EXPECT_DOUBLE_EQ( tv[ 0 ], 3.0 );
}

TEST_F( _TEST_TITLE_, MultiArgumentDoesNotThrowIfSettingVector ) {
    MultiValueArgument<int> arg( "testarg", "Test flag for unit tests" );
    EXPECT_NO_THROW( { arg.set_value( std::vector<int> { 1, 2, 3 } ); } );
}

TEST_F( _TEST_TITLE_, ArgumentSetKeepsArguments ) {
    ArgumentSet argset;
    ASSERT_EQ( argset.count(), 0 );
    argset.add_argument(
        SingleValueArgument<int>( "testarg", "Test argument", -4 ) );
    EXPECT_EQ( argset.count(), 1 );
    EXPECT_EQ( argset.argument( "testarg" ).as_int(), -4 );
    argset.add_argument(
        SingleValueArgument<string>( "teststr", "Test string", "5" ) );
    EXPECT_EQ( argset.count(), 2 );
    EXPECT_EQ( argset.argument( "teststr" ).as_string(), "5" );
    EXPECT_EQ( argset.argument( "teststr" ).as_int(), 5 );
}

TEST_F( _TEST_TITLE_, ParserParsesProperly ) {
    vector<string> tokens { "--testarg", "-4", "--teststr", "5" };
    ArgumentSet argset;
    argset.add_argument(
        SingleValueArgument<int>( "testarg", "Test argument" ) );
    argset.add_argument(
        SingleValueArgument<string>( "teststr", "Test string" ) );
    CommandLineArgumentParser parser;
    parser.load( tokens ).parse( argset );
    EXPECT_EQ( argset.count(), 2 );
    EXPECT_EQ( argset.argument( "testarg" ).as_int(), -4 );
    EXPECT_EQ( argset.argument( "teststr" ).as_string(), "5" );
}

TEST_F( _TEST_TITLE_, TestHelpText ) {
    HelpTextOrganizerSection base( "Help Text Tester" );
    // base.add_subsection( HelpTextParagraphSection(
    //     "Paragraph Section 1",
    //     "This is a long string of text.  This is a long string of text.  This "
    //     "is a long string of text.  This is a long string of text.  This is a "
    //     "long string of text.  This is a long string of text.  This is a long "
    //     "string of text.  This is a long string of text.  This is a long "
    //     "string of text.  This is a long string of text.<br>This is after a "
    //     "newline." ) );
    // base.add_subsection( HelpTextParagraphSection(
    //     "Paragraph Section 2",
    //     "This is a long string of text.  This is a long string of text.  This "
    //     "is a long string of text.  This is a long string of text.<br>This is "
    //     "after a newline.  This is a "
    //     "long string of text.  This is a long string of text.  This is a long "
    //     "string of text.  This is a long string of text.  This is a long "
    //     "string of text.  This is a long string of text." ) );
    
    ArgumentSet argset;
    argset.add_required_argument<int>( "arg1", "This is an integer argument that should receive one value." );
    argset.add_multivalue_argument<double>( "arg2", "This is a floating-point argument that can receive multiple values." );
    argset.add_flag( "flag1", "This is a flag argument" );
    base.add_subsection( HelpTextArgumentSection( "Arguments", argset ) );

    HelpTextFormatter formatter;
    formatter.format( base );
}
