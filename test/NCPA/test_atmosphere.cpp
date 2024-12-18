#include "NCPA/atmosphere.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/units.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::atmos;

typedef double test_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPAAtmosphereLibraryTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            test_atmos    = factory.build( filename );
            first5[ "Z" ] = vector<double> { 0.0, 0.1, 0.2, 0.3, 0.4 };
            first5[ "T" ]
                = vector<double> { 2.881495e+02, 2.875774e+02, 2.870075e+02,
                                   2.864399e+02, 2.858728e+02 };
            first5[ "U" ]
                = vector<double> { 4.709510e-04, 4.897540e-04, 5.092760e-04,
                                   5.295410e-04, 5.505760e-04 };
            first5[ "V" ] = vector<double> { 0, 0, 0, 0, 0 };
        }  // void TearDown() override {}

        // declare stuff here
        Atmosphere1D test_atmos, test_atmos2;
        const test_t zero    = NCPA::math::zero<test_t>(),
                     one     = NCPA::math::one<test_t>();
        std::string filename = "testdata/test_atmosphere_1d.ncpaprop";
        AtmosphereFactory factory;
        map<string, vector<double>> first5;
};

TEST_F( _TEST_TITLE_, DefaultConstructorWorks ) {
    EXPECT_FALSE( test_atmos2 );
}

TEST_F( _TEST_TITLE_, BuilderReturnsValidAtmosphere ) {
    EXPECT_TRUE( test_atmos );
}

TEST_F( _TEST_TITLE_, AtmosphereHasCorrectScalarProperties ) {
    EXPECT_DOUBLE_EQ( test_atmos.get( "Z0" ), 0.0 );
    EXPECT_TRUE(
        test_atmos.get_property_units( "Z0" )->equals( NCPA::units::METERS ) );
}

TEST_F( _TEST_TITLE_, AtmosphereHasCorrectNumberOfPoints ) {
    EXPECT_EQ( test_atmos.get_property( "U" ).size(), 2001 );
}

TEST_F( _TEST_TITLE_, AtmosphereThrowsLogicErrorIfNotInitialized ) {
    EXPECT_THROW( { test_atmos2.get_keys(); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, DependentValuesAreCorrect ) {
    EXPECT_TRUE(
        test_atmos.get_property_units( "T" )->equals( NCPA::units::KELVIN ) );
    EXPECT_TRUE( test_atmos.get_property_units( "U" )->equals(
        NCPA::units::METERS_PER_SECOND ) );
    EXPECT_TRUE( test_atmos.get_property_units( "V" )->equals(
        NCPA::units::METERS_PER_SECOND ) );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( test_atmos.get( "T", first5[ "Z" ][ i ] ),
                          first5[ "T" ][ i ] );
        EXPECT_DOUBLE_EQ( test_atmos.get( "U", first5[ "Z" ][ i ] ),
                          first5[ "U" ][ i ] );
        EXPECT_DOUBLE_EQ( test_atmos.get( "V", first5[ "Z" ][ i ] ),
                          first5[ "V" ][ i ] );
    }
}

TEST_F( _TEST_TITLE_, UnitsCanBeConverted ) {
    test_atmos.convert_property_units( "U",
                                       &NCPA::units::KILOMETERS_PER_SECOND );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( test_atmos.get( "U", first5[ "Z" ][ i ] ),
                          first5[ "U" ][ i ] * 0.001 );
    }
    EXPECT_TRUE( test_atmos.get_property_units( "U" )->equals(
        NCPA::units::KILOMETERS_PER_SECOND ) );
}

TEST_F( _TEST_TITLE_, ResampleWorksWithDZ ) {
    test_atmos.resample( 0.0499999 );
    EXPECT_EQ( test_atmos.get_property( "U" ).size(), 4001 );
}

TEST_F( _TEST_TITLE_, ResampleWorksWithZ ) {
    vector_t new_z( vector<double>( 4000 ), "km" );
    for ( size_t i = 0; i < 4000; i++ ) {
        double di  = (double)i;
        new_z[ i ] = di * 0.05;
    }
    test_atmos.resample( new_z );
    EXPECT_EQ( test_atmos.get_property( "U" ).size(), 4000 );
}

TEST_F( _TEST_TITLE_, AssignmentOperatorWorks ) {
    test_atmos2 = test_atmos;
    EXPECT_TRUE( test_atmos2.get_property_units( "T" )->equals(
        *test_atmos2.get_property_units( "T" ) ) );
    EXPECT_TRUE( test_atmos2.get_property_units( "U" )->equals(
        *test_atmos2.get_property_units( "U" ) ) );
    EXPECT_TRUE( test_atmos2.get_property_units( "V" )->equals(
        *test_atmos2.get_property_units( "V" ) ) );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( test_atmos2.get( "T", first5[ "Z" ][ i ] ),
                          test_atmos.get( "T", first5[ "Z" ][ i ] ) );
        EXPECT_DOUBLE_EQ( test_atmos2.get( "U", first5[ "Z" ][ i ] ),
                          test_atmos.get( "U", first5[ "Z" ][ i ] ) );
        EXPECT_DOUBLE_EQ( test_atmos2.get( "V", first5[ "Z" ][ i ] ),
                          test_atmos.get( "V", first5[ "Z" ][ i ] ) );
    }
}

TEST_F( _TEST_TITLE_, DefaultInterpolatorIsNearestNeighbor ) {
    double dz      = 0.1;
    double test_dz = dz * 0.1;
    for ( double z = dz; z < 100.0; z += dz ) {
        EXPECT_DOUBLE_EQ( test_atmos.get( "T", z - test_dz ),
                          test_atmos.get( "T", z + test_dz ) );
    }
}

TEST_F( _TEST_TITLE_, NewInterpolatorIsUsedCorrectly ) {
    double dz   = 0.1;
    double half = dz * 0.5;
    test_atmos.set_interpolator(
        NCPA::interpolation::interpolator_type_t::LANL_LINEAR );
    for ( double z = dz; z < 100.0; z += dz ) {
        EXPECT_NEAR(
            test_atmos.get( "T", z - half ),
            0.5 * ( test_atmos.get( "T", z ) + test_atmos.get( "T", z - dz ) ),
            0.5 );
    }
}

TEST_F( _TEST_TITLE_, CopyPropertyWorks ) {
    test_atmos.add_property( "NEW_U", test_atmos.get_property( "U" ) );
    double dz = 0.1;
    for ( double z = dz; z < 100.0; z += dz ) {
        EXPECT_NEAR( test_atmos.get( "U", z ), test_atmos.get( "NEW_U", z ), 1e-10 );
    }
}

TEST_F( _TEST_TITLE_, InvalidKeyCausesException ) {
    EXPECT_THROW( { test_atmos.get( "G", 0.5 ); }, std::range_error );
}

TEST_F( _TEST_TITLE_, TemperatureToSoundSpeedCalculationIsCorrect ) {
    vector<double> t { 273, 283, 293, 303, 313 };
    vector<double> cx;
    for ( auto it = t.begin(); it != t.end(); ++it ) {
        cx.push_back( sqrt( 1.4 * 287 * *it ) );
    }
    vector_t tv( t, NCPA::units::Units::from_string( "K" ) );
    vector_t cv = t2c( tv );
    for (auto i = 0; i < cv.size(); i++) {
        EXPECT_DOUBLE_EQ( cv.get( i ), cx[ i ] );
    }
}

TEST_F( _TEST_TITLE_, PressureAndDensityToSoundSpeedIsCorrect ) {
    EXPECT_DOUBLE_EQ( pd2c( 
        scalar_t( 1.0, "atm" ),
        scalar_t( 1.225, "kg/m3" )
    ).get(), std::sqrt( 101325.0 * 1.4 / 1.225 ) );
}

TEST_F( _TEST_TITLE_, WindVectorsToSpeedAndDirectionAreCorrect ) {
    scalar_t u( 0.0, "m/s" ), v( 0.0, "m/s" );
    EXPECT_DOUBLE_EQ( uv2ws( u, v ).get(), 0.0 );
    u = 2.0;
    EXPECT_DOUBLE_EQ( uv2ws( u, v ).get(), 2.0 );
    EXPECT_DOUBLE_EQ( uv2wd( u, v ).get(), 90.0 );
    v = 2.0;
    EXPECT_DOUBLE_EQ( uv2ws( u, v ).get(), sqrt(8.0) );
    EXPECT_DOUBLE_EQ( uv2wd( u, v ).get(), 45.0 );
    u = -2.0;
    EXPECT_DOUBLE_EQ( uv2ws( u, v ).get(), sqrt(8.0) );
    EXPECT_DOUBLE_EQ( uv2wd( u, v ).get(), 315.0 );
}

TEST_F( _TEST_TITLE_, WindSpeedAndDirectionToComponentIsCorrect ) {
    scalar_t ws( sqrt( 8.0 ), "m/s" ), wd( 45.0, "deg" );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 45.0 ).get(), sqrt( 8.0 ) );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 0.0 ).get(), 2.0 );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 90.0 ).get(), 2.0 );
    EXPECT_NEAR( w2wc( ws, wd, 135.0 ).get(), 0.0, 1e-12 );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 180.0 ).get(), -2.0 );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 225.0 ).get(), -sqrt( 8.0 ) );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 270.0 ).get(), -2.0 );
    EXPECT_NEAR( w2wc( ws, wd, 315.0 ).get(), 0.0, 1e-12 );
}