#include "NCPA/atmosphere.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/ndvector.hpp"
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
using namespace NCPA::units;
using namespace NCPA::arrays;

typedef double test_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPAAtmosphereLibraryTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            AtmosphereReader1D reader
                = AtmosphereFactory::build( reader_1d_t::NCPAPROP );
            test_atmos                  = reader.read( filename );
            AtmosphereReader2D reader2d = AtmosphereFactory::build(
                reader_2d_t::NCPAPROP_PIECEWISE_STRATIFIED );
            test_2datmos = reader2d.read( filename2d );


            first5[ "Z" ] = vector<double> { 0.0, 0.1, 0.2, 0.3, 0.4 };
            first5[ "T" ]
                = vector<double> { 2.881495e+02, 2.875774e+02, 2.870075e+02,
                                   2.864399e+02, 2.858728e+02 };
            first5[ "U" ]
                = vector<double> { 4.709510e-04, 4.897540e-04, 5.092760e-04,
                                   5.295410e-04, 5.505760e-04 };
            first5[ "V" ] = vector<double> { 0, 0, 0, 0, 0 };

            vector<double> upto5 { 0, 1, 2, 3, 4 };
            VectorWithUnits<double> km5( upto5, KILOMETERS );
            Vector2DWithUnits<double> kmgrid( 5, 5, KILOMETERS );
            for (auto it = kmgrid.begin(); it != kmgrid.end(); ++it) {
                *it = upto5;
            }

            prop2d
                = AtmosphereFactory::build( atmospheric_property_2d_t::GRID );
            prop2d.set( km5, km5, kmgrid );

            strat.add_property( "U", prop2d );
            grid.add_property( "U", prop2d );

            stratwrap
                = AtmosphereFactory::build( atmosphere_2d_t::STRATIFIED );
            stratwrap.add_property( "U", prop2d );

        }  // void TearDown() override {}

        // declare stuff here
        Atmosphere1D test_atmos, test_atmos2;
        Atmosphere2D test_2datmos;
        Atmosphere3D test_3datmos;
        const test_t zero      = NCPA::math::zero<test_t>(),
                     one       = NCPA::math::one<test_t>();
        std::string filename   = "testdata/test_atmosphere_1d.ncpaprop";
        std::string filename2d = "testdata/test_atmosphere_2d.ncpaprop";
        std::string filename3d = "testdata/test_atmosphere_3d.ncpaprop";
        map<string, vector<double>> first5;

        AtmosphericProperty2D prop2d;
        stratified_atmosphere_2d strat;
        grid_atmosphere_2d grid;
        Atmosphere2D stratwrap, gridwrap;
};

TEST_F( _TEST_TITLE_, DefaultConstructorWorks ) {
    EXPECT_FALSE( test_atmos2 );
}

TEST_F( _TEST_TITLE_, BuilderReturnsValidAtmosphere ) {
    EXPECT_TRUE( test_atmos );
}

TEST_F( _TEST_TITLE_, AtmosphereHasCorrectScalarProperties ) {
    EXPECT_DOUBLE_EQ( test_atmos.get( "Z0" ), 0.0 );
    EXPECT_TRUE( test_atmos.get_units( "Z0" )->equals( NCPA::units::METERS ) );
}

TEST_F( _TEST_TITLE_, AtmosphereHasCorrectNumberOfPoints ) {
    EXPECT_EQ( test_atmos.get_property( "U" ).size(), 2001 );
}

TEST_F( _TEST_TITLE_, AtmosphereThrowsLogicErrorIfNotInitialized ) {
    EXPECT_THROW( { test_atmos2.get_keys(); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, DependentValuesAreCorrect ) {
    EXPECT_TRUE( test_atmos.get_units( "T" )->equals( NCPA::units::KELVIN ) );
    EXPECT_TRUE( test_atmos.get_units( "U" )->equals(
        NCPA::units::METERS_PER_SECOND ) );
    EXPECT_TRUE( test_atmos.get_units( "V" )->equals(
        NCPA::units::METERS_PER_SECOND ) );
    for (size_t i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( test_atmos.get( "T", first5[ "Z" ][ i ] ),
                          first5[ "T" ][ i ] );
        EXPECT_DOUBLE_EQ( test_atmos.get( "U", first5[ "Z" ][ i ] ),
                          first5[ "U" ][ i ] );
        EXPECT_DOUBLE_EQ( test_atmos.get( "V", first5[ "Z" ][ i ] ),
                          first5[ "V" ][ i ] );
    }
}

TEST_F( _TEST_TITLE_, UnitsCanBeConverted ) {
    test_atmos.convert_units( "U", &NCPA::units::KILOMETERS_PER_SECOND );
    for (size_t i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( test_atmos.get( "U", first5[ "Z" ][ i ] ),
                          first5[ "U" ][ i ] * 0.001 );
    }
    EXPECT_TRUE( test_atmos.get_units( "U" )->equals(
        NCPA::units::KILOMETERS_PER_SECOND ) );
}

TEST_F( _TEST_TITLE_, ResampleWorksWithDZ ) {
    test_atmos.resample( 0.0499999 );
    EXPECT_EQ( test_atmos.get_property( "U" ).size(), 4001 );
}

TEST_F( _TEST_TITLE_, ResampleWorksWithZ ) {
    vector_u_t new_z( vector<double>( 4000 ), "km" );
    for (size_t i = 0; i < 4000; i++) {
        double di  = (double)i;
        new_z[ i ] = di * 0.05;
    }
    test_atmos.resample( new_z );
    EXPECT_EQ( test_atmos.get_property( "U" ).size(), 4000 );
}

TEST_F( _TEST_TITLE_, AssignmentOperatorWorks ) {
    test_atmos2 = test_atmos;
    EXPECT_TRUE( test_atmos2.get_units( "T" )->equals(
        *test_atmos2.get_units( "T" ) ) );
    EXPECT_TRUE( test_atmos2.get_units( "U" )->equals(
        *test_atmos2.get_units( "U" ) ) );
    EXPECT_TRUE( test_atmos2.get_units( "V" )->equals(
        *test_atmos2.get_units( "V" ) ) );
    for (size_t i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( test_atmos2.get( "T", first5[ "Z" ][ i ] ),
                          test_atmos.get( "T", first5[ "Z" ][ i ] ) );
        EXPECT_DOUBLE_EQ( test_atmos2.get( "U", first5[ "Z" ][ i ] ),
                          test_atmos.get( "U", first5[ "Z" ][ i ] ) );
        EXPECT_DOUBLE_EQ( test_atmos2.get( "V", first5[ "Z" ][ i ] ),
                          test_atmos.get( "V", first5[ "Z" ][ i ] ) );
    }
}

TEST_F( _TEST_TITLE_, NewInterpolatorIsUsedCorrectly ) {
    double dz   = 0.1;
    double half = dz * 0.5;
    test_atmos.set_interpolator(
        NCPA::interpolation::interpolator_1d_type_t::LANL_LINEAR );
    for (double z = dz; z < 100.0; z += dz) {
        EXPECT_NEAR(
            test_atmos.get( "T", z - half ),
            0.5 * ( test_atmos.get( "T", z ) + test_atmos.get( "T", z - dz ) ),
            0.5 );
    }
}

TEST_F( _TEST_TITLE_, CopyPropertyWorks ) {
    test_atmos.add_property( "NEW_U", test_atmos.get_property( "U" ) );
    double dz = 0.1;
    for (double z = dz; z < 100.0; z += dz) {
        EXPECT_NEAR( test_atmos.get( "U", z ), test_atmos.get( "NEW_U", z ),
                     1e-10 );
    }
}

TEST_F( _TEST_TITLE_, InvalidKeyCausesException ) {
    EXPECT_THROW( { test_atmos.get( "G", 0.5 ); }, std::range_error );
}

TEST_F( _TEST_TITLE_, TemperatureToSoundSpeedCalculationIsCorrect ) {
    vector<double> t { 273, 283, 293, 303, 313 };
    vector<double> cx;
    for (auto it = t.begin(); it != t.end(); ++it) {
        cx.push_back( sqrt( 1.4 * 287 * *it ) );
    }
    vector_u_t tv( t, NCPA::units::Units::from_string( "K" ) );
    vector_u_t cv = t2c( tv );
    for (auto i = 0; i < cv.size(); i++) {
        EXPECT_DOUBLE_EQ( cv.get( i ), cx[ i ] );
    }
}

TEST_F( _TEST_TITLE_, PressureAndDensityToSoundSpeedIsCorrect ) {
    EXPECT_DOUBLE_EQ(
        pd2c( scalar_u_t( 1.0, "atm" ), scalar_u_t( 1.225, "kg/m3" ) ).get(),
        std::sqrt( 101325.0 * 1.4 / 1.225 ) );
}

TEST_F( _TEST_TITLE_, WindVectorsToSpeedAndDirectionAreCorrect ) {
    scalar_u_t u( 0.0, "m/s" ), v( 0.0, "m/s" );
    EXPECT_DOUBLE_EQ( uv2ws( u, v ).get(), 0.0 );
    u = 2.0;
    EXPECT_DOUBLE_EQ( uv2ws( u, v ).get(), 2.0 );
    EXPECT_DOUBLE_EQ( uv2wd( u, v ).get(), 90.0 );
    v = 2.0;
    EXPECT_DOUBLE_EQ( uv2ws( u, v ).get(), sqrt( 8.0 ) );
    EXPECT_DOUBLE_EQ( uv2wd( u, v ).get(), 45.0 );
    u = -2.0;
    EXPECT_DOUBLE_EQ( uv2ws( u, v ).get(), sqrt( 8.0 ) );
    EXPECT_DOUBLE_EQ( uv2wd( u, v ).get(), 315.0 );
}

TEST_F( _TEST_TITLE_, WindSpeedAndDirectionToComponentIsCorrect ) {
    scalar_u_t ws( sqrt( 8.0 ), "m/s" ), wd( 45.0, "deg" );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 45.0 ).get(), sqrt( 8.0 ) );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 0.0 ).get(), 2.0 );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 90.0 ).get(), 2.0 );
    EXPECT_NEAR( w2wc( ws, wd, 135.0 ).get(), 0.0, 1e-12 );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 180.0 ).get(), -2.0 );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 225.0 ).get(), -sqrt( 8.0 ) );
    EXPECT_DOUBLE_EQ( w2wc( ws, wd, 270.0 ).get(), -2.0 );
    EXPECT_NEAR( w2wc( ws, wd, 315.0 ).get(), 0.0, 1e-12 );
}

TEST_F( _TEST_TITLE_, Property2DReportsCorrectly ) {
    EXPECT_DOUBLE_EQ( prop2d.get( 0.0, 0.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( prop2d.get( 1.0, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( prop2d.get( 1.1, 1.1 ), 1.0 );
    EXPECT_DOUBLE_EQ( prop2d.get( 1.6, 1.6 ), 2.0 );
    prop2d.set_interpolator(
        NCPA::interpolation::interpolator_2d_type_t::LANL_NATURAL );
    EXPECT_DOUBLE_EQ( prop2d.get( 1.6, 1.6 ), 1.6 );
    EXPECT_DOUBLE_EQ( prop2d.get_first_derivative( 1.6, 1.6, 0 ), 0.0 );
    EXPECT_DOUBLE_EQ( prop2d.get_first_derivative( 1.6, 1.6, 1 ), 1.0 );
    EXPECT_DOUBLE_EQ( prop2d.get_second_derivative( 1.6, 1.6, 0, 0 ), 0.0 );
    EXPECT_DOUBLE_EQ( prop2d.get_second_derivative( 1.6, 1.6, 1, 0 ), 0.0 );
    EXPECT_DOUBLE_EQ( prop2d.get_second_derivative( 1.6, 1.6, 1, 1 ), 0.0 );
}

TEST_F( _TEST_TITLE_, Property2DConvertsUnitsCorrectly ) {
    prop2d.convert_units( METERS );
    EXPECT_DOUBLE_EQ( prop2d.get( 1.0, 1.0 ), 1000.0 );
}

TEST_F( _TEST_TITLE_, Property2DConvertsAxisUnitsCorrectly ) {
    prop2d.convert_axis_units( 0, METERS );
    EXPECT_DOUBLE_EQ( prop2d.get( 1.0, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( prop2d.get( 1499.0, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( prop2d.get( 1501.0, 1.0 ), 1.0 );
    prop2d.convert_axis_units( 1, METERS );
    EXPECT_DOUBLE_EQ( prop2d.get( 1.0, 1.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( prop2d.get( 1499.0, 1499.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( prop2d.get( 1501.0, 1501.0 ), 2.0 );
}

TEST_F( _TEST_TITLE_, Property2DResampleWorksCorrectly ) {
    EXPECT_EQ( prop2d.size( 0 ), 5 );
    EXPECT_EQ( prop2d.size( 1 ), 5 );
    vector_u_t new_r( 4000, prop2d.get_axis_units( 0 ) );
    for (size_t i = 0; i < 4000; ++i) {
        new_r[ i ] = (double)i * 0.001;
    }
    vector_u_t new_z = prop2d.axis( 1 );
    prop2d.resample( new_r, new_z );
    EXPECT_EQ( prop2d.size( 0 ), 4000 );
    EXPECT_EQ( prop2d.size( 1 ), 5 );
}

TEST_F( _TEST_TITLE_, Stratified2DWorksAsExpected ) {
    // strat.set_interpolator(
    // NCPA::interpolation::interpolator_1d_type_t::LANL_LINEAR );
    EXPECT_DOUBLE_EQ( strat.get( "U", 0.0, 0.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.0, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.1, 1.1 ), 1.1 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.6, 1.6 ), 1.6 );
    strat.set_interpolator(
        NCPA::interpolation::interpolator_1d_type_t::NEAREST_NEIGHBOR );
    EXPECT_DOUBLE_EQ( strat.get( "U", 0.0, 0.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.0, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.1, 1.1 ), 1.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.6, 1.6 ), 2.0 );
}

TEST_F( _TEST_TITLE_, Stratified2DWrapperWorksAsExpected ) {
    EXPECT_DOUBLE_EQ( stratwrap.get( "U", 0.0, 0.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( stratwrap.get( "U", 1.0, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( stratwrap.get( "U", 1.1, 1.1 ), 1.1 );
    EXPECT_DOUBLE_EQ( stratwrap.get( "U", 1.6, 1.6 ), 1.6 );
    stratwrap.set_interpolator(
        NCPA::interpolation::interpolator_1d_type_t::NEAREST_NEIGHBOR );
    EXPECT_DOUBLE_EQ( stratwrap.get( "U", 0.0, 0.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( stratwrap.get( "U", 1.0, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( stratwrap.get( "U", 1.1, 1.1 ), 1.0 );
    EXPECT_DOUBLE_EQ( stratwrap.get( "U", 1.6, 1.6 ), 2.0 );
}

TEST_F( _TEST_TITLE_, Stratified2DCanBeBuiltFrom1D ) {
    stratwrap.set( test_atmos );
    // cout << "Set OK" << std::endl;
    EXPECT_DOUBLE_EQ( stratwrap.get( "Z0" ), 0.0 );
    // cout << "get OK," << std::endl;
    EXPECT_TRUE( stratwrap.get_units( "Z0" )->equals( NCPA::units::METERS ) );
    // cout << "getz OK" << std::endl;
    EXPECT_TRUE( stratwrap.get_units( "T" )->equals( NCPA::units::KELVIN ) );
    // cout << "gett OK" << std::endl;
    EXPECT_TRUE(
        stratwrap.get_units( "U" )->equals( NCPA::units::METERS_PER_SECOND ) );
    // cout << "getu OK" << std::endl;
    EXPECT_TRUE(
        stratwrap.get_units( "V" )->equals( NCPA::units::METERS_PER_SECOND ) );
    // cout << "getv OK" << std::endl;
    for (size_t i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( stratwrap.get( "T", 0.0, first5[ "Z" ][ i ] ),
                          first5[ "T" ][ i ] );
        EXPECT_DOUBLE_EQ( stratwrap.get( "U", 0.0, first5[ "Z" ][ i ] ),
                          first5[ "U" ][ i ] );
        EXPECT_DOUBLE_EQ( stratwrap.get( "V", 0.0, first5[ "Z" ][ i ] ),
                          first5[ "V" ][ i ] );
    }
}

TEST_F( _TEST_TITLE_, Grid2DWorksAsExpected ) {
    grid.set_interpolator(
        NCPA::interpolation::interpolator_1d_type_t::LANL_LINEAR );
    EXPECT_DOUBLE_EQ( strat.get( "U", 0.0, 0.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.0, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.1, 1.1 ), 1.1 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.6, 1.6 ), 1.6 );
    strat.set_interpolator(
        NCPA::interpolation::interpolator_1d_type_t::NEAREST_NEIGHBOR );
    EXPECT_DOUBLE_EQ( strat.get( "U", 0.0, 0.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.0, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.1, 1.1 ), 1.0 );
    EXPECT_DOUBLE_EQ( strat.get( "U", 1.6, 1.6 ), 2.0 );
}

TEST_F( _TEST_TITLE_, LocallyStratified2DReadProperly ) {
    EXPECT_FALSE( test_2datmos.same( 100.0, 300.0 ) );
    EXPECT_DOUBLE_EQ( test_2datmos.get( "U", 100.0, 2.1 ), 8.808512e-04 );
    EXPECT_DOUBLE_EQ( test_2datmos.get( "U", 400.0, 52.2 ), 40.99138 );
    EXPECT_DOUBLE_EQ( test_2datmos.get_first_derivative( "U", 400.0, 52.2, 0 ),
                      0.0 );
    EXPECT_FALSE( test_2datmos.get_first_derivative( "U", 400.0, 52.2, 1 )
                  == 0.0 );
    EXPECT_DOUBLE_EQ(
        test_2datmos.get_second_derivative( "U", 400.0, 52.2, 0, 1 ), 0.0 );
    EXPECT_DOUBLE_EQ(
        test_2datmos.get_second_derivative( "U", 400.0, 52.2, 0, 0 ), 0.0 );
    EXPECT_DOUBLE_EQ(
        test_2datmos.get_second_derivative( "U", 400.0, 52.2, 1, 0 ), 0.0 );
    EXPECT_FALSE( test_2datmos.get_second_derivative( "U", 400.0, 52.2, 1, 1 )
                  == 0.0 );
}

TEST_F( _TEST_TITLE_, Atmosphere3DReadsProperly ) {
    AtmosphereReader3D reader3d
        = AtmosphereFactory::build( reader_3d_t::NCPAPROP );
    reader3d.set_axis_units( 1, &NCPA::units::RADIANS );
    test_3datmos = reader3d.read( filename3d );
    EXPECT_FALSE( test_3datmos.same( 100.0, -0.1, 300.0, 0.1 ) );

    // test atmosphere is stratified in theta
    EXPECT_DOUBLE_EQ( test_3datmos.get( "U", 100.0, 0.0, 2.1 ), 
        test_3datmos.get( "U", 100.0, 0.1, 2.1 ) );

    EXPECT_DOUBLE_EQ( test_3datmos.get( "U", 0.0, 0.0, 2.1 ), 8.808512e-04 );
    // EXPECT_DOUBLE_EQ( test_3datmos.get( "U", 400.0, 0.2, 52.2 ), 40.99138 );
}
