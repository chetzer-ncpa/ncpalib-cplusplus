#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/math.hpp"
#include "NCPA/units.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstring>
#include <list>
#include <stdexcept>
#include <vector>

#include <gtest/gtest-spi.h>

using namespace std;
using namespace testing;
using namespace NCPA::units;

#define _TEST_TITLE_ NCPAUnitsLibraryTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            M  = Unit( "m", { "meters" } );
            MM = Unit( "mm", { "millimeters" }, &M, 1.0e-3 );
            KM = Unit( "km", { "kilometers" }, &M, 1.0e3 );

            K    = Unit( "K", { "Kelvin" } );
            C    = Unit( "C", { "degrees C", "Celsius", "degC" }, &K, 1.0,
                         273.15 );
            F    = Unit( "F", { "degrees F", "Fahrenheit", "degF" }, &C,
                         5.0 / 9.0, -32.0 );
            FAKE = Unit( "FAKE", {}, &F, 1.0, 0.0, -5.0 );

            testval1 = 4.35;

            Cvec = { 0.0, 100.0, 200.0 };
            Fvec = { 32.0, 212.0, 392.0 };
            Carr = NCPA::arrays::zeros<double>( 3 );
            std::memcpy( Carr, &Cvec[ 0 ], 3 * sizeof( double ) );
            Farr = NCPA::arrays::zeros<double>( 3 );
            std::memcpy( Farr, &Fvec[ 0 ], 3 * sizeof( double ) );

            freezing  = ScalarWithUnits<>( 0.0, CELSIUS );
            boiling   = ScalarWithUnits<>( 100.0, CELSIUS );
            halfway   = ScalarWithUnits<>( 50.0, CELSIUS );
            s1        = ScalarWithUnits<>( boiling );
            s2        = ScalarWithUnits<>( freezing );
            firstLeg  = ScalarWithUnits<>( 10.0, METERS );
            secondLeg = ScalarWithUnits<>( 20.0, METERS );
            thirdLeg  = ScalarWithUnits<>( 0.03, KILOMETERS );

            v2 = VectorWithUnits<>( 10, kms, KILOMETERS );
            v3 = VectorWithUnits<>( 10, kms, "kilometers" );
            ScalarWithUnits<> s( 10.0, CELSIUS );
            v4 = VectorWithUnits<>( 10, s );
            v5 = VectorWithUnits<>( 10, 10.0, CELSIUS );
            v6 = VectorWithUnits<>( v4 );
            // v7 = VectorWithUnits<>( 5, svec );
            v8 = VectorWithUnits<>( hottemps, CELSIUS );
            v9 = VectorWithUnits<>( hottemps, "F" );
            v0 = VectorWithUnits<>();

            size_t nd = 5;
            // double **testnums = NCPA::arrays::zeros<double>( 5, 5 );
            NCPA::arrays::vector2d_t<double> testnums( 5, 5 );
            for ( auto i = 0; i < 5; ++i ) {
                for ( auto j = 0; j < 5; ++j ) {
                    testnums[ i ][ j ] = (double)i * 10.0;
                }
            }
            v2d_2 = Vector2DWithUnits<double>( testnums, CELSIUS );
            v2d_3 = Vector2DWithUnits<double>( 5, 5, 10.0, CELSIUS );
        }  // void TearDown() override {}

        // declare stuff here
        Unit M, MM, KM, K, C, F, FAKE;
        double testval1;
        std::vector<double> Cvec, Fvec;
        double *Carr, *Farr;

        ScalarWithUnits<> freezing, boiling, halfway, s0, s1, s2;
        ScalarWithUnits<> firstLeg, secondLeg, thirdLeg;

        VectorWithUnits<> v1, v2, v3, v4, v5, v6, v8, v9, v0;
        double kms[ 10 ]
            = { 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0 };
        vector<double> temps        = { 10, 20, 30, 40, 50 };
        vector<double> hottemps     = { 50, 60, 70, 80, 90, 100 };
        ScalarWithUnits<> svec[ 5 ] = {
            ScalarWithUnits<>( 10.0, CELSIUS ),
            ScalarWithUnits<>( 20.0, CELSIUS ),
            ScalarWithUnits<>( 30.0, CELSIUS ),
            ScalarWithUnits<>( 40.0, CELSIUS ),
            ScalarWithUnits<>( 50.0, CELSIUS ),
        };
        double all10[ 10 ] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 };

        Vector2DWithUnits<double> v2d_1, v2d_2, v2d_3;
};

TEST_F( _TEST_TITLE_, ReferenceReturnedCorrectly ) {
    EXPECT_TRUE( K.reference()->equals( K ) );
    EXPECT_TRUE( F.reference()->equals( C ) );
    EXPECT_FALSE( F.reference()->equals( K ) );
}

TEST_F( _TEST_TITLE_, IsConvertibleToReturnsTrueForCommonReference ) {
    EXPECT_TRUE( M.is_convertible_to( MM ) );
    EXPECT_TRUE( F.is_convertible_to( K ) );
}

TEST_F( _TEST_TITLE_, IsConvertibleToReturnsFalseForNoCommonReference ) {
    EXPECT_FALSE( M.is_convertible_to( C ) );
    EXPECT_FALSE( K.is_convertible_to( M ) );
}

TEST_F( _TEST_TITLE_, BaseReferencesSelfIdentify ) {
    EXPECT_TRUE( M.is_base_reference() );
    EXPECT_FALSE( MM.is_base_reference() );
}

TEST_F( _TEST_TITLE_, EqualsReturnsTrueForEquality ) {
    EXPECT_TRUE( C.equals( C ) );
    Unit Cdup
        = Unit( "Cdup", { "degrees C", "Celsius", "degC" }, &K, 1.0, 273.15 );
    EXPECT_TRUE( Cdup.equals( C ) );
}

TEST_F( _TEST_TITLE_, EqualsReturnsFalseForInequality ) {
    EXPECT_FALSE( C.equals( F ) );
}

TEST_F( _TEST_TITLE_, ToReferenceWorksCorrectly ) {
    EXPECT_DOUBLE_EQ( M.to_reference( 3.0 ), 3.0 );
    EXPECT_DOUBLE_EQ( MM.to_reference( 10000.0 ), 10.0 );
    EXPECT_DOUBLE_EQ( C.to_reference( 0.0 ), 273.15 );
    EXPECT_DOUBLE_EQ( F.to_reference( 212.0 ), 100.0 );
    EXPECT_DOUBLE_EQ( FAKE.to_reference( 5.0 ), 0.0 );
}

TEST_F( _TEST_TITLE_, FromReferenceWorksCorrectly ) {
    EXPECT_DOUBLE_EQ( M.from_reference( 3.0 ), 3.0 );
    EXPECT_DOUBLE_EQ( MM.from_reference( 10000.0 ), 10000000.0 );
    EXPECT_DOUBLE_EQ( C.from_reference( 0.0 ), -273.15 );
    EXPECT_DOUBLE_EQ( F.from_reference( 0.0 ), 32.0 );
    EXPECT_DOUBLE_EQ( FAKE.from_reference( 5.0 ), 10.0 );
}

TEST_F( _TEST_TITLE_, ToBaseReferenceWorksCorrectly ) {
    EXPECT_DOUBLE_EQ( M.to_base_reference( 3.0 ), 3.0 );
    EXPECT_DOUBLE_EQ( MM.to_base_reference( 10000.0 ), 10.0 );
    EXPECT_DOUBLE_EQ( C.to_base_reference( 0.0 ), 273.15 );
    EXPECT_DOUBLE_EQ( F.to_base_reference( 212.0 ), 373.15 );
    EXPECT_DOUBLE_EQ( FAKE.to_base_reference( 217.0 ), 373.15 );
}

TEST_F( _TEST_TITLE_, FromBaseReferenceWorksCorrectly ) {
    EXPECT_DOUBLE_EQ( M.from_base_reference( 3.0 ), 3.0 );
    EXPECT_DOUBLE_EQ( MM.from_base_reference( 10000.0 ), 10000000.0 );
    EXPECT_DOUBLE_EQ( C.from_base_reference( 0.0 ), -273.15 );
    EXPECT_DOUBLE_EQ( F.from_base_reference( 9.0 ), -443.47 );
    EXPECT_DOUBLE_EQ( FAKE.from_base_reference( 9.0 ), -438.47 );
}

TEST_F( _TEST_TITLE_, ConvertingToSameUnitGivesSameValue ) {
    EXPECT_DOUBLE_EQ( M.convert_to( testval1, M ), testval1 );
    EXPECT_DOUBLE_EQ( MM.convert_to( testval1, MM ), testval1 );
    EXPECT_DOUBLE_EQ( K.convert_to( testval1, K ), testval1 );
    EXPECT_DOUBLE_EQ( F.convert_to( testval1, F ), testval1 );
}

TEST_F( _TEST_TITLE_, ConvertingToInvalidUnitThrowsInvalidArgument ) {
    EXPECT_THROW( { M.convert_to( testval1, K ); }, invalid_conversion<> );
    EXPECT_THROW( { M.convert_to( testval1, F ); }, invalid_conversion<> );
    EXPECT_THROW( { K.convert_to( testval1, M ); }, invalid_conversion<> );
    EXPECT_THROW( { F.convert_to( testval1, M ); }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, TestConversionsGiveCorrectValues ) {
    EXPECT_DOUBLE_EQ( M.convert_to( testval1, KM ), testval1 * 0.001 );
    EXPECT_DOUBLE_EQ( KM.convert_to( testval1, M ), testval1 * 1000.0 );
    EXPECT_DOUBLE_EQ( M.convert_to( testval1, MM ), testval1 * 1000.0 );
    EXPECT_DOUBLE_EQ( MM.convert_to( testval1, M ), testval1 * 0.001 );
    EXPECT_DOUBLE_EQ( KM.convert_to( testval1, MM ), testval1 * 1.0e6 );
    EXPECT_DOUBLE_EQ( MM.convert_to( testval1, KM ), testval1 * 1.0e-6 );

    EXPECT_DOUBLE_EQ( K.convert_to( 0.0, C ), -273.15 );
    EXPECT_DOUBLE_EQ( F.convert_to( 212.0, C ), 100.0 );
    EXPECT_DOUBLE_EQ( F.convert_to( 32.0, K ), 273.15 );
}

TEST_F( _TEST_TITLE_, ConvertingToWorksWithVectors ) {
    std::vector<double> Ftest = C.convert_to( Cvec, F );
    EXPECT_ARRAY_DOUBLE_EQ( 3, Ftest, Fvec );
    std::vector<double> Ctest = F.convert_to( Fvec, C );
    EXPECT_ARRAY_DOUBLE_EQ( 3, Ctest, Cvec );
}

TEST_F( _TEST_TITLE_, ConvertingToWorksWithArrays ) {
    double *Ftest = NCPA::arrays::zeros<double>( 3 );
    double *Ctest = NCPA::arrays::zeros<double>( 3 );
    F.convert_to( 3, Farr, C, Ctest );
    C.convert_to( 3, Carr, F, Ftest );
    EXPECT_ARRAY_DOUBLE_EQ( 3, Ftest, Farr );
    EXPECT_ARRAY_DOUBLE_EQ( 3, Ctest, Carr );
}

/////////////////////////////////////////////////
// Tests for Units static framework

TEST_F( _TEST_TITLE_, UnitsStringConvertWorksForScalar ) {
    EXPECT_DOUBLE_EQ( Units::convert( 10.0, "km", "m" ), 10000.0 );
}

TEST_F( _TEST_TITLE_, UnitsStringConvertWorksForIterable ) {
    EXPECT_ARRAY_DOUBLE_EQ( 3, Units::convert( Cvec, "C", "F" ), Fvec );
}

TEST_F( _TEST_TITLE_, UnitsStringConvertWorksForArray ) {
    Units::convert( 3, Carr, "C", "F", Carr );
    EXPECT_ARRAY_DOUBLE_EQ( 3, Carr, Farr );
}

TEST_F( _TEST_TITLE_, UnitsLookupWithAliasReturnsCorrectUnit ) {
    EXPECT_TRUE( Units::from_string( "C" )->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, UnitsLookupWithAliasReturnsSameUnit ) {
    EXPECT_TRUE( Units::from_string( "C" )->equals(
        *( Units::from_string( "Celsius" ) ) ) );
}

TEST_F( _TEST_TITLE_, CanConvertReturnsCorrectlyWithStrings ) {
    EXPECT_TRUE( Units::can_convert( "C", "F" ) );
    EXPECT_TRUE( Units::can_convert( "m", "km" ) );
    EXPECT_FALSE( Units::can_convert( "C", "m" ) );
}

///////////////////////////////////////////////////
// Tests for ScalarWithUnits class
TEST_F( _TEST_TITLE_, ConstructorTests ) {
    ASSERT_TRUE( s0.get_units() == nullptr );
    ASSERT_DOUBLE_EQ( s0.get(), 0.0 );

    ASSERT_TRUE( boiling.get_units()->equals( CELSIUS ) );
    ASSERT_DOUBLE_EQ( boiling.get(), 100.0 );

    ASSERT_TRUE( s1.get_units()->equals( CELSIUS ) );
    ASSERT_DOUBLE_EQ( s1.get(), 100.0 );

    // copy constructor
    ScalarWithUnits<> s4( s0 );
    ASSERT_TRUE( s4.get_units() == nullptr );
    ASSERT_DOUBLE_EQ( s4.get(), 0.0 );
    s0.set( 50.0, CELSIUS );
    ASSERT_TRUE( s0.get_units()->equals( CELSIUS ) );
    ASSERT_DOUBLE_EQ( s0.get(), 50.0 );
    ASSERT_TRUE( s4.get_units() == nullptr );
    ASSERT_DOUBLE_EQ( s4.get(), 0.0 );

    // move constructor?
    s2 = std::move( boiling );
    ASSERT_TRUE( s2.get_units()->equals( CELSIUS ) );
    ASSERT_DOUBLE_EQ( s2.get(), 100.0 );
}

TEST_F( _TEST_TITLE_, GetTests ) {
    ASSERT_DOUBLE_EQ( freezing.get(), 0.0 );
    ASSERT_TRUE( freezing.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, GetAsReturnsCorrectValue ) {
    ASSERT_DOUBLE_EQ( freezing.get_as( FAHRENHEIT ), 32.0 );
    ASSERT_DOUBLE_EQ( freezing.as( "Kelvin" ).get(), 273.15 );
}

TEST_F( _TEST_TITLE_, GetAsThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW(
        { double d = freezing.get_as( METERS ); }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, SetValueSetsCorrectValue ) {
    freezing.set_value( 212.0 );
    ASSERT_DOUBLE_EQ( freezing.get(), 212.0 );
}

TEST_F( _TEST_TITLE_, SetValueDoesNotChangeUnits ) {
    const Unit *u = freezing.get_units();
    freezing.set_value( 212.0 );
    ASSERT_TRUE( freezing.get_units()->equals( *u ) );
}

TEST_F( _TEST_TITLE_, SetUnitsSetsCorrectUnits ) {
    freezing.set_units( FAHRENHEIT );
    ASSERT_TRUE( freezing.get_units()->equals( FAHRENHEIT ) );
}

TEST_F( _TEST_TITLE_, SetUnitsDoesNotChangeValue ) {
    freezing.set_units( FAHRENHEIT );
    ASSERT_DOUBLE_EQ( freezing.get(), 0.0 );
}

TEST_F( _TEST_TITLE_, SetUnitsWorksWithString ) {
    freezing.set_units( "F" );
    ASSERT_TRUE( freezing.get_units()->equals( FAHRENHEIT ) );
}

TEST_F( _TEST_TITLE_, SetOverwritesBothValueAndUnits ) {
    // Set units to 32 degrees fahrenheit
    firstLeg.set( 32.0, FAHRENHEIT );
    ASSERT_DOUBLE_EQ( firstLeg.get(), 32.0 );
    ASSERT_TRUE( firstLeg.get_units()->equals( FAHRENHEIT ) );
}

TEST_F( _TEST_TITLE_, ConvertUnitsConvertsCorrectly ) {
    freezing.convert_units( FAHRENHEIT );
    ASSERT_DOUBLE_EQ( freezing.get(), 32.0 );
}

TEST_F( _TEST_TITLE_, ConvertUnitsThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW(
        { freezing.convert_units( METERS ); }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, EqualityOperatorIsTrueForEqual ) {
    ASSERT_TRUE( s1 == boiling );
}

TEST_F( _TEST_TITLE_, EqualityOperatorIsFalseForUnequal ) {
    ASSERT_FALSE( s1 == freezing );
}

TEST_F( _TEST_TITLE_,
        EqualityOperatorThrowsInvalidConversionForEqualWithNullUnits ) {
    boiling.set_units( (Unit *)nullptr );
    EXPECT_THROW( { bool tf = s1 == boiling; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, InequalityOperatorIsTrueForUnequal ) {
    ASSERT_TRUE( freezing != boiling );
}

TEST_F( _TEST_TITLE_, InequalityOperatorIsFalseForEqual ) {
    ASSERT_FALSE( boiling != boiling );
}

TEST_F( _TEST_TITLE_,
        InequalityOperatorThrowsInvalidConversionForEqualWithNullUnits ) {
    boiling.set_units( (Unit *)nullptr );
    // ASSERT_FALSE( freezing == boiling );
    EXPECT_THROW( { bool tf = freezing != boiling; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, GreaterThanOperatorTrueForGreater ) {
    ASSERT_TRUE( boiling > freezing );
}

TEST_F( _TEST_TITLE_, GreaterThanOperatorFalseForLess ) {
    ASSERT_FALSE( freezing > boiling );
}

TEST_F( _TEST_TITLE_, GreaterThanOperatorFalseForEqual ) {
    ASSERT_FALSE( boiling > boiling );
}

TEST_F( _TEST_TITLE_,
        GreaterThanOperatorThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW( { bool b = freezing > firstLeg; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, LessThanOperatorFalseForGreater ) {
    ASSERT_FALSE( boiling < freezing );
}

TEST_F( _TEST_TITLE_, LessThanOperatorTrueForLesser ) {
    ASSERT_TRUE( freezing < boiling );
}

TEST_F( _TEST_TITLE_, LessThanOperatorFalseForEqual ) {
    ASSERT_FALSE( boiling < boiling );
}

TEST_F( _TEST_TITLE_, LessThanOperatorThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW( { bool b = freezing < firstLeg; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, GreaterThanOrEqualOperatorTrueForGreater ) {
    ASSERT_TRUE( boiling >= freezing );
}

TEST_F( _TEST_TITLE_, GreaterThanOrEqualOperatorFalseForLess ) {
    ASSERT_FALSE( freezing >= boiling );
}

TEST_F( _TEST_TITLE_, GreaterThanOrEqualOperatorTrueForEqual ) {
    ASSERT_TRUE( boiling >= boiling );
}

TEST_F( _TEST_TITLE_,
        GreaterThanOrEqualOperatorThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW( { bool b = freezing >= firstLeg; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, LessThanOrEqualOperatorFalseForGreater ) {
    ASSERT_FALSE( boiling <= freezing );
}

TEST_F( _TEST_TITLE_, LessThanOrEqualOperatorTrueForLesser ) {
    ASSERT_TRUE( freezing <= boiling );
}

TEST_F( _TEST_TITLE_, LessThanOrEqualOperatorTrueForEqual ) {
    ASSERT_TRUE( boiling <= boiling );
}

TEST_F( _TEST_TITLE_,
        LessThanOrEqualOperatorThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW( { bool b = freezing <= firstLeg; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, NegationOperatorMakesValueNegative ) {
    ASSERT_DOUBLE_EQ( ( -boiling ).get(), -( boiling.get() ) );
}

TEST_F( _TEST_TITLE_, NegationOperatorDoesNotChangeUnits ) {
    ASSERT_TRUE(
        ( -boiling ).get_units()->equals( *( boiling.get_units() ) ) );
}

TEST_F( _TEST_TITLE_, PlusOperatorGivesCorrectSum ) {
    ScalarWithUnits<> firstAndSecond = firstLeg + secondLeg;
    ASSERT_DOUBLE_EQ( firstAndSecond.get(), 30.0 );
}

TEST_F( _TEST_TITLE_, PlusOperatorGivesCorrectSumWithUnitConversion ) {
    ScalarWithUnits<> secondAndThird = secondLeg + thirdLeg;
    ASSERT_DOUBLE_EQ( secondAndThird.get(), 50.0 );
}

TEST_F( _TEST_TITLE_, PlusOperatorKeepsCorrectUnitsWithUnitConversion ) {
    ScalarWithUnits<> secondAndThird = secondLeg + thirdLeg;
    ASSERT_TRUE( secondAndThird.get_units()->equals( METERS ) );
}

TEST_F( _TEST_TITLE_, PlusOperatorChainsCorrectly ) {
    ScalarWithUnits<> fullTrip = firstLeg + secondLeg + thirdLeg;
    ASSERT_DOUBLE_EQ( fullTrip.get(), 60.0 );
}

TEST_F( _TEST_TITLE_, PlusOperatorChainKeepsCorrectUnits ) {
    ScalarWithUnits<> fullTrip = firstLeg + secondLeg + thirdLeg;
    ASSERT_TRUE( fullTrip.get_units()->equals( METERS ) );
}

TEST_F( _TEST_TITLE_, PlusOperatorWorksWithDouble ) {
    ScalarWithUnits<> superHeated = boiling + 50.0;
    ASSERT_DOUBLE_EQ( superHeated.get(), 150.0 );
    ASSERT_TRUE( superHeated.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, PlusOperatorThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW(
        { ScalarWithUnits<> b = freezing + firstLeg; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, MinusOperatorGivesCorrectDifference ) {
    ScalarWithUnits<> firstAndSecond = firstLeg - secondLeg;
    ASSERT_DOUBLE_EQ( firstAndSecond.get(), -10.0 );
}

TEST_F( _TEST_TITLE_, MinusOperatorGivesCorrectDifferenceWithUnitConversion ) {
    ScalarWithUnits<> secondAndThird = secondLeg - thirdLeg;
    ASSERT_DOUBLE_EQ( secondAndThird.get(), -10.0 );
}

TEST_F( _TEST_TITLE_, MinusOperatorKeepsCorrectUnitsWithUnitConversion ) {
    ScalarWithUnits<> secondAndThird = secondLeg - thirdLeg;
    ASSERT_TRUE( secondAndThird.get_units()->equals( METERS ) );
}

TEST_F( _TEST_TITLE_, MinusOperatorChainsCorrectly ) {
    ScalarWithUnits<> fullTrip = firstLeg - secondLeg - thirdLeg;
    ASSERT_DOUBLE_EQ( fullTrip.get(), -40.0 );
}

TEST_F( _TEST_TITLE_, MinusOperatorChainKeepsCorrectUnits ) {
    ScalarWithUnits<> fullTrip = firstLeg - secondLeg - thirdLeg;
    ASSERT_TRUE( fullTrip.get_units()->equals( METERS ) );
}

TEST_F( _TEST_TITLE_, MinusOperatorWorksWithDouble ) {
    ScalarWithUnits<> superCooled = freezing - 50.0;
    ASSERT_DOUBLE_EQ( superCooled.get(), -50.0 );
    ASSERT_TRUE( superCooled.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, MinusOperatorThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW(
        { ScalarWithUnits<> b = freezing - firstLeg; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, PlusEqualsOperatorGivesCorrectSum ) {
    firstLeg += secondLeg;
    ASSERT_DOUBLE_EQ( firstLeg.get(), 30.0 );
}

TEST_F( _TEST_TITLE_, PlusEqualsOperatorGivesCorrectSumWithUnitConversion ) {
    secondLeg += thirdLeg;
    ASSERT_DOUBLE_EQ( secondLeg.get(), 50.0 );
}

TEST_F( _TEST_TITLE_, PlusEqualsOperatorKeepsCorrectUnitsWithUnitConversion ) {
    secondLeg += thirdLeg;
    ASSERT_TRUE( secondLeg.get_units()->equals( METERS ) );
}

TEST_F( _TEST_TITLE_, PlusEqualsOperatorWorksWithDouble ) {
    boiling += 50.0;
    ASSERT_DOUBLE_EQ( boiling.get(), 150.0 );
    ASSERT_TRUE( boiling.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, PlusEqualsOperatorThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW( { freezing += firstLeg; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, MinusEqualsOperatorGivesCorrectDifference ) {
    firstLeg -= secondLeg;
    ASSERT_DOUBLE_EQ( firstLeg.get(), -10.0 );
}

TEST_F( _TEST_TITLE_,
        MinusEqualsOperatorGivesCorrectDifferenceWithUnitConversion ) {
    secondLeg -= thirdLeg;
    ASSERT_DOUBLE_EQ( secondLeg.get(), -10.0 );
}

TEST_F( _TEST_TITLE_, MinusEqualsOperatorWorksWithDouble ) {
    freezing -= 50.0;
    ASSERT_DOUBLE_EQ( freezing.get(), -50.0 );
    ASSERT_TRUE( freezing.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_,
        MinusEqualsOperatorThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW( { freezing -= firstLeg; }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, MultiplyEqualsOperatorScalesCorrectly ) {
    boiling *= 3.0;
    ASSERT_DOUBLE_EQ( boiling.get(), 300.0 );
    ASSERT_TRUE( boiling.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, DivideEqualsOperatorScalesCorrectly ) {
    boiling /= 4.0;
    ASSERT_DOUBLE_EQ( boiling.get(), 25.0 );
    ASSERT_TRUE( boiling.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, MultiplyOperatorScalesCorrectly ) {
    ScalarWithUnits<> temp = boiling * 3.0;
    ASSERT_DOUBLE_EQ( temp.get(), 300.0 );
    ASSERT_TRUE( temp.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, DivideOperatorScalesCorrectly ) {
    ScalarWithUnits<> temp = boiling / 4.0;
    ASSERT_DOUBLE_EQ( temp.get(), 25.0 );
    ASSERT_TRUE( temp.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, OverGivesCorrectRatio ) {
    ASSERT_DOUBLE_EQ( boiling.over( halfway ), 2.0 );
}

TEST_F( _TEST_TITLE_, OverThrowsOutOfRangeOnInvalidConversion ) {
    EXPECT_THROW(
        { double d = boiling.over( firstLeg ); }, invalid_conversion<> );
}

//////////////////////////////////////
// Tests for VectorWithUnits class
TEST_F( _TEST_TITLE_, DefaultConstructorIsEmpty ) {
    ASSERT_EQ( v1.size(), 0 );
}

TEST_F( _TEST_TITLE_, ConstructorWithDoubleArrayAndUnitsObjectWorks ) {
    EXPECT_EQ( v2.size(), 10 );
    EXPECT_TRUE( v2.get_units()->equals( KILOMETERS ) );
    double buffer[ 10 ];
    v2.get_values( buffer );
    EXPECT_ARRAY_DOUBLE_EQ( 10, buffer, kms );
}

TEST_F( _TEST_TITLE_, ConstructorWithDoubleArrayAndUnitsStringWorks ) {
    EXPECT_EQ( v3.size(), 10 );
    EXPECT_TRUE( v3.get_units()->equals( KILOMETERS ) );
    double buffer[ 10 ];
    v3.get_values( buffer );
    EXPECT_ARRAY_DOUBLE_EQ( 10, buffer, kms );
}

TEST_F( _TEST_TITLE_, ConstructorWithConstantScalarWithUnitsWorks ) {
    // EXPECT_THAT( v4, SizeIs( 10 ) );
    EXPECT_EQ( v4.size(), 10 );
    EXPECT_TRUE( v4.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < v4.size(); i++ ) {
        EXPECT_DOUBLE_EQ( v4[ i ], 10.0 );
    }
}

TEST_F( _TEST_TITLE_, ConstructorWithConstantDoubleWorks ) {
    // EXPECT_THAT( v5, SizeIs( 10 ) );
    EXPECT_EQ( v5.size(), 10 );
    EXPECT_TRUE( v5.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < v5.size(); i++ ) {
        EXPECT_DOUBLE_EQ( v5[ i ], 10.0 );
    }
}

TEST_F( _TEST_TITLE_, ConstructorWithVectorWorks ) {
    EXPECT_EQ( v8.size(), 6 );
    // EXPECT_THAT( v8, SizeIs( 6 ) );
    EXPECT_TRUE( v8.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < v8.size(); i++ ) {
        EXPECT_DOUBLE_EQ( v8[ i ], hottemps[ i ] );
    }
}

TEST_F( _TEST_TITLE_, ConstructorWithVectorAndStringWorks ) {
    EXPECT_EQ( v9.size(), 6 );
    // EXPECT_THAT( v9, SizeIs( 6 ) );
    EXPECT_TRUE( v9.get_units()->equals( FAHRENHEIT ) );
    for ( size_t i = 0; i < v9.size(); i++ ) {
        EXPECT_DOUBLE_EQ( v9[ i ], hottemps[ i ] );
    }
}

TEST_F( _TEST_TITLE_, CopyConstructorWorks ) {
    EXPECT_EQ( v6.size(), 10 );
    // EXPECT_THAT( v6, SizeIs( 10 ) );
    EXPECT_TRUE( v6.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < v6.size(); i++ ) {
        EXPECT_DOUBLE_EQ( v6[ i ], 10.0 );
    }
}

TEST_F( _TEST_TITLE_, AssignmentOperatorWorks ) {
    v3 = v6;
    EXPECT_EQ( v3.size(), v6.size() );
    EXPECT_TRUE( v3.get_units()->equals( *v6.get_units() ) );
    for ( auto v3it = v3.cbegin(), v6it = v6.cbegin();
          v3it != v3.cend() && v6it != v6.cend(); ++v3it, ++v6it ) {
        EXPECT_TRUE( *v3it == *v6it );
    }
}

TEST_F( _TEST_TITLE_, SwapWorks ) {
    size_t v6size = v6.size(), v3size = v3.size();
    const Unit *v6u = v6.get_units();
    const Unit *v3u = v3.get_units();
    swap( v6, v3 );

    EXPECT_EQ( v6.size(), v3size );
    EXPECT_EQ( v6.get_units(), v3u );
    double *buffer = new double[ v6.size() ];
    v6.get_values( buffer );
    EXPECT_ARRAY_DOUBLE_EQ( 5, buffer, kms );
    delete[] buffer;
    buffer = nullptr;

    // EXPECT_THAT( v7, SizeIs( v6size ) );
    EXPECT_EQ( v3.size(), v6size );
    EXPECT_EQ( v3.get_units(), v6u );
    buffer = new double[ v3.size() ];
    v3.get_values( buffer );
    EXPECT_ARRAY_DOUBLE_EQ( 10, buffer, all10 );
    delete[] buffer;
}

TEST_F( _TEST_TITLE_, AsArrayReturnsCorrectArray ) {
    ScalarWithUnits<double> *buffer = nullptr;
    v4.as_array( buffer );
    for ( auto i = 0; i < v4.size(); i++ ) {
        EXPECT_DOUBLE_EQ( buffer[ i ].get(), 10.0 );
        EXPECT_TRUE( buffer[ i ].get_units()->equals( CELSIUS ) );
    }
    delete[] buffer;
}

TEST_F( _TEST_TITLE_, AsArrayReturnsCorrectDoubleArray ) {
    double *buffer = nullptr;
    v4.as_array( buffer );
    for ( auto i = 0; i < v4.size(); i++ ) {
        EXPECT_DOUBLE_EQ( buffer[ i ], 10.0 );
    }
    delete[] buffer;
}

TEST_F( _TEST_TITLE_, AsArrayReturnsCorrectDoubleVector ) {
    vector<double> v = v4;
    for ( auto it = v.cbegin(); it != v.end(); ++it ) {
        EXPECT_DOUBLE_EQ( *it, 10.0 );
    }
}

TEST_F( _TEST_TITLE_, ConvertUnitsCreatesCorrectValues ) {
    v2.convert_units( METERS );
    for ( size_t i = 0; i < 10; i++ ) {
        EXPECT_DOUBLE_EQ( v2[ i ], kms[ i ] * 1000.0 );
    }
}

TEST_F( _TEST_TITLE_, ConvertUnitsStoresCorrectUnits ) {
    v2.convert_units( METERS );
    EXPECT_TRUE( v2.get_units()->equals( METERS ) );
}

TEST_F( _TEST_TITLE_, ConvertUnitsWithStringCreatesCorrectValues ) {
    v2.convert_units( "m" );
    for ( size_t i = 0; i < 10; i++ ) {
        EXPECT_DOUBLE_EQ( v2[ i ], kms[ i ] * 1000.0 );
    }
}

TEST_F( _TEST_TITLE_, ConvertUnitsWithStringStoresCorrectUnits ) {
    v2.convert_units( "m" );

    EXPECT_TRUE( v2.get_units()->equals( METERS ) );
}

TEST_F( _TEST_TITLE_, GetUnitsReturnsCorrectUnits ) {
    EXPECT_TRUE( v4.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, ConvertUnitsThrowsInvalidConversionCorrectly ) {
    EXPECT_THROW( { v2.convert_units( CELSIUS ); }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, FillSetsValuesCorrectly ) {
    size_t n = v2.size();
    v2.fill( 12.0, Units::from_string( "kg/m3" ) );
    for ( size_t i = 0; i < n; i++ ) {
        EXPECT_DOUBLE_EQ( v2[ i ], 12.0 );
    }
}

TEST_F( _TEST_TITLE_, FillSetsUnitsCorrectly ) {
    size_t n      = v2.size();
    const Unit *u = Units::from_string( "kg/m3" );
    v2.fill( 12.0, u );

    EXPECT_TRUE( v2.get_units()->equals( u ) );
}

TEST_F( _TEST_TITLE_, FillSetsUnitsCorrectlyFromString ) {
    size_t n      = v2.size();
    const Unit *u = Units::from_string( "kg/m3" );
    v2.fill( 12.0, "kg/m3" );

    EXPECT_TRUE( v2.get_units()->equals( u ) );
}

TEST_F( _TEST_TITLE_, FillSetsValuesCorrectlyFromScalarObject ) {
    size_t n      = v2.size();
    const Unit *u = Units::from_string( "kg/m3" );
    ScalarWithUnits<double> s( 12.0, *u );
    v2.fill( s );
    for ( size_t i = 0; i < n; i++ ) {
        EXPECT_DOUBLE_EQ( v2[ i ], s.get() );
    }
}

TEST_F( _TEST_TITLE_, FillSetsUnitsCorrectlyFromScalarObject ) {
    size_t n      = v2.size();
    const Unit *u = Units::from_string( "kg/m3" );
    ScalarWithUnits<double> s( 12.0, *u );
    v2.fill( s );

    EXPECT_TRUE( v2.get_units()->equals( *s.get_units() ) );
}

TEST_F( _TEST_TITLE_, GetValuesReturnsCorrectArray ) {
    double buffer[ 10 ];
    size_t n = 10;
    v2.get_values( n, buffer );
    EXPECT_ARRAY_DOUBLE_EQ( n, buffer, kms );
}

TEST_F( _TEST_TITLE_, GetValuesWithoutSizeReturnsCorrectArray ) {
    double buffer[ 10 ];
    v2.get_values( buffer );
    EXPECT_ARRAY_DOUBLE_EQ( 10, buffer, kms );
}

TEST_F( _TEST_TITLE_, SetResizesCorrectly ) {
    ASSERT_EQ( v2.size(), 10 );
    v2.set( 5, &temps[ 0 ], &CELSIUS );
    ASSERT_EQ( v2.size(), 5 );
}

TEST_F( _TEST_TITLE_, SetSetsNewVectorValuesCorrectly ) {
    v2.set( 5, &temps[ 0 ], &CELSIUS );
    for ( size_t i = 0; i < 5; i++ ) {
        ASSERT_DOUBLE_EQ( v2[ i ], temps[ i ] );
    }
}

TEST_F( _TEST_TITLE_, SetSetsNewVectorUnitsCorrectly ) {
    v2.set( 5, &temps[ 0 ], &CELSIUS );

    ASSERT_TRUE( v2.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, SetSetsNewVectorUnitsCorrectlyFromString ) {
    v2.set( temps, "C" );

    ASSERT_TRUE( v2.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_,
        SetSetsNewVectorValuesCorrectlyFromArrayOfScalarObjects ) {
    v2.set( 5, svec );
    ASSERT_TRUE( v2.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < 5; i++ ) {
        ASSERT_DOUBLE_EQ( v2[ i ], temps[ i ] );
    }
}

TEST_F( _TEST_TITLE_, SetWorksWithVector ) {
    v2.set( hottemps, &CELSIUS );
    ASSERT_EQ( v2.size(), 6 );
    ASSERT_TRUE( v2.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < 6; i++ ) {
        ASSERT_DOUBLE_EQ( v2[ i ], hottemps[ i ] );
    }
}

TEST_F( _TEST_TITLE_, SetWorksWithVectorAAndString ) {
    v2.set( hottemps, "F" );
    ASSERT_EQ( v2.size(), 6 );
    ASSERT_TRUE( v2.get_units()->equals( FAHRENHEIT ) );
    for ( size_t i = 0; i < 6; i++ ) {
        ASSERT_DOUBLE_EQ( v2[ i ], hottemps[ i ] );
    }
}

TEST_F( _TEST_TITLE_, SetUnitsSetsAllUnitsCorrectly ) {
    v2.set_units( METERS );

    ASSERT_TRUE( v2.get_units()->equals( METERS ) );
}

TEST_F( _TEST_TITLE_, SetUnitsSetsAllUnitsCorrectlyFromString ) {
    v2.set_units( "C" );

    ASSERT_TRUE( v2.get_units()->equals( CELSIUS ) );
}

TEST_F( _TEST_TITLE_, FakeUnitsRegisterProperly ) {
    Unit HELENS( "helens", { "hn" } );
    Unit MILLIHELENS( "mh", { "millihelens" }, &HELENS, 0.001 );
    Units::register_unit( HELENS );
    Units::register_unit( MILLIHELENS );

    ScalarWithUnits<> ships( 1.0, HELENS );
    EXPECT_DOUBLE_EQ( ships.get_as( MILLIHELENS ), 1000.0 );
}

TEST_F( _TEST_TITLE_, RegisteringExistingUnitRaisesInvalidArgument ) {
    EXPECT_THROW(
        { Units::register_unit( CELSIUS ); }, std::invalid_argument );
    Unit HELENS( "helens", { "hn", "m" } );
    EXPECT_THROW( { Units::register_unit( HELENS ); }, std::invalid_argument );
}

// actual unit conversion tests
TEST_F( _TEST_TITLE_, DefinedUnitConversionsAreCorrect ) {
    EXPECT_DOUBLE_EQ( KILOMETERS.convert_to( 1.0, METERS ), 1000.0 );
    EXPECT_DOUBLE_EQ( METERS.convert_to( 1.0, KILOMETERS ), 1.0e-3 );

    EXPECT_DOUBLE_EQ( KILOMETERS.convert_to( 1.0, MILLIMETERS ), 1.0e6 );
    EXPECT_DOUBLE_EQ( MILLIMETERS.convert_to( 1.0, KILOMETERS ), 1.0e-6 );

    EXPECT_DOUBLE_EQ( METERS.convert_to( 1.0, MILLIMETERS ), 1.0e3 );
    EXPECT_DOUBLE_EQ( MILLIMETERS.convert_to( 1.0, METERS ), 1.0e-3 );

    EXPECT_DOUBLE_EQ( KELVIN.convert_to( 1.0, CELSIUS ), -272.15 );
    EXPECT_DOUBLE_EQ( CELSIUS.convert_to( 1.0, KELVIN ), 274.15 );

    EXPECT_DOUBLE_EQ( KELVIN.convert_to( 1.0, FAHRENHEIT ), -457.87 );
    EXPECT_DOUBLE_EQ( FAHRENHEIT.convert_to( 41.0, KELVIN ), 278.15 );

    EXPECT_DOUBLE_EQ( CELSIUS.convert_to( 100.0, FAHRENHEIT ), 212.0 );
    EXPECT_DOUBLE_EQ( FAHRENHEIT.convert_to( 41.0, CELSIUS ), 5.0 );

    EXPECT_DOUBLE_EQ(
        METERS_PER_SECOND.convert_to( 1.0, KILOMETERS_PER_SECOND ), 1.0e-3 );
    EXPECT_DOUBLE_EQ(
        KILOMETERS_PER_SECOND.convert_to( 1.0, METERS_PER_SECOND ), 1000.0 );

    EXPECT_DOUBLE_EQ( PASCALS.convert_to( 1.0, HECTOPASCALS ), 0.01 );
    EXPECT_DOUBLE_EQ( HECTOPASCALS.convert_to( 1.0, PASCALS ), 100.0 );

    EXPECT_DOUBLE_EQ( PASCALS.convert_to( 1.0, MILLIBARS ), 0.01 );
    EXPECT_DOUBLE_EQ( MILLIBARS.convert_to( 1.0, PASCALS ), 100.0 );

    EXPECT_DOUBLE_EQ( HECTOPASCALS.convert_to( 1.0, MILLIBARS ), 1.0 );
    EXPECT_DOUBLE_EQ( MILLIBARS.convert_to( 1.0, HECTOPASCALS ), 1.0 );

    EXPECT_NEAR( PASCALS.convert_to( 1.0, ATMOSPHERES ), 0.00000986923, 1e-9 );
    EXPECT_DOUBLE_EQ( ATMOSPHERES.convert_to( 1.0, PASCALS ), 101325.0 );

    EXPECT_NEAR( HECTOPASCALS.convert_to( 1.0, ATMOSPHERES ), 0.000986923,
                 1e-7 );
    EXPECT_DOUBLE_EQ( ATMOSPHERES.convert_to( 1.0, HECTOPASCALS ), 1013.250 );

    EXPECT_NEAR( MILLIBARS.convert_to( 1.0, ATMOSPHERES ), 0.000986923, 1e-7 );
    EXPECT_DOUBLE_EQ( ATMOSPHERES.convert_to( 1.0, MILLIBARS ), 1013.250 );

    EXPECT_DOUBLE_EQ( KILOGRAMS_PER_CUBIC_METER.convert_to(
                          1.0, GRAMS_PER_CUBIC_CENTIMETER ),
                      0.001 );
    EXPECT_DOUBLE_EQ( GRAMS_PER_CUBIC_CENTIMETER.convert_to(
                          1.0, KILOGRAMS_PER_CUBIC_METER ),
                      1000.0 );

    EXPECT_DOUBLE_EQ(
        DECIBELS_PER_METER.convert_to( 1.0, DECIBELS_PER_KILOMETER ), 1.0e-3 );
    EXPECT_DOUBLE_EQ(
        DECIBELS_PER_KILOMETER.convert_to( 1.0, DECIBELS_PER_METER ), 1.0e3 );

    EXPECT_NEAR( DECIBELS_PER_METER.convert_to( 1.0, NEPERS_PER_METER ),
                 0.115129254650564, 1e-7 );
    EXPECT_NEAR( NEPERS_PER_METER.convert_to( 1.0, DECIBELS_PER_METER ),
                 8.685889638, 1e-7 );

    EXPECT_NEAR( DECIBELS_PER_KILOMETER.convert_to( 1.0, NEPERS_PER_METER ),
                 115.129254650564, 1e-5 );
    EXPECT_NEAR( NEPERS_PER_METER.convert_to( 1.0, DECIBELS_PER_KILOMETER ),
                 0.008685889638, 1e-9 );

    EXPECT_DOUBLE_EQ( KILOGRAMS.convert_to( 1.0, GRAMS ), 1000.0 );
    EXPECT_DOUBLE_EQ( GRAMS.convert_to( 1.0, KILOGRAMS ), 1.0e-3 );

    EXPECT_DOUBLE_EQ( MINUTES.convert_to( 1.0, SECONDS ), 60.0 );
    EXPECT_DOUBLE_EQ( HOURS.convert_to( 1.0, SECONDS ), 3600.0 );
    EXPECT_DOUBLE_EQ( DAYS.convert_to( 1.0, SECONDS ), 86400.0 );

    EXPECT_DOUBLE_EQ( SECONDS.convert_to( 120.0, MINUTES ), 2.0 );
    EXPECT_DOUBLE_EQ( HOURS.convert_to( 1.0, MINUTES ), 60.0 );
    EXPECT_DOUBLE_EQ( DAYS.convert_to( 1.0, MINUTES ), 24.0 * 60.0 );

    EXPECT_DOUBLE_EQ( SECONDS.convert_to( 1.0, HOURS ), 1.0 / 3600.0 );
    EXPECT_DOUBLE_EQ( MINUTES.convert_to( 30.0, HOURS ), 0.5 );
    EXPECT_DOUBLE_EQ( DAYS.convert_to( 1.0, HOURS ), 24.0 );

    EXPECT_DOUBLE_EQ( SECONDS.convert_to( 1.0, DAYS ), 1.0 / 86400.0 );
    EXPECT_DOUBLE_EQ( MINUTES.convert_to( 30.0, DAYS ), 1.0 / 48.0 );
    EXPECT_DOUBLE_EQ( HOURS.convert_to( 72.0, DAYS ), 3.0 );
}

//////////////////////////////////////
// Tests for Vector2DWithUnits class
TEST_F( _TEST_TITLE_, Default2DConstructorIsEmpty ) {
    ASSERT_EQ( v2d_1.dim( 0 ), 0 );
    ASSERT_EQ( v2d_1.dim( 1 ), 0 );
}

TEST_F( _TEST_TITLE_, Constructor2DWithConstantDoubleWorks ) {
    // EXPECT_THAT( v4, SizeIs( 10 ) );
    EXPECT_EQ( v2d_3.dim( 0 ), 5 );
    EXPECT_EQ( v2d_3.dim( 1 ), 5 );
    EXPECT_TRUE( v2d_3.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < v2d_3.dim( 0 ); i++ ) {
        for ( size_t j = 0; j < v2d_3.dim( 1 ); j++ ) {
            EXPECT_DOUBLE_EQ( v2d_3[ i ][ j ], 10.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, ConstructorWithVector2DWorks ) {
    EXPECT_EQ( v2d_2.dim( 0 ), 5 );
    EXPECT_EQ( v2d_2.dim( 1 ), 5 );
    EXPECT_TRUE( v2d_2.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < v2d_2.dim( 0 ); i++ ) {
        for ( size_t j = 0; j < v2d_2.dim( 1 ); j++ ) {
            EXPECT_DOUBLE_EQ( v2d_2[ i ][ j ], (double)i * 10.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, CopyConstructor2DWorks ) {
    Vector2DWithUnits<double> tester(v2d_2);
    EXPECT_EQ( tester.dim( 0 ), 5 );
    EXPECT_EQ( tester.dim( 1 ), 5 );
    EXPECT_TRUE( tester.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < tester.dim( 0 ); i++ ) {
        for ( size_t j = 0; j < tester.dim( 1 ); j++ ) {
            EXPECT_DOUBLE_EQ( tester[ i ][ j ], (double)i * 10.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, AssignmentOperatorWorks2D ) {
    Vector2DWithUnits<double> tester = v2d_2;
    EXPECT_EQ( tester.dim( 0 ), 5 );
    EXPECT_EQ( tester.dim( 1 ), 5 );
    EXPECT_TRUE( tester.get_units()->equals( CELSIUS ) );
    for ( size_t i = 0; i < tester.dim( 0 ); i++ ) {
        for ( size_t j = 0; j < tester.dim( 1 ); j++ ) {
            EXPECT_DOUBLE_EQ( tester[ i ][ j ], (double)i * 10.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, SwapWorks2D ) {
    v2d_1.set_units( KILOGRAMS );
    swap( v2d_1, v2d_2 );

    EXPECT_EQ( v2d_1.dim(0), 5 );
    EXPECT_EQ( v2d_1.dim(1), 5 );
    EXPECT_EQ( v2d_2.dim(0), 0 );
    EXPECT_EQ( v2d_2.dim(1), 0 );
    EXPECT_TRUE( v2d_1.get_units()->equals( CELSIUS ) );
    EXPECT_TRUE( v2d_2.get_units()->equals( KILOGRAMS ));

    for ( size_t i = 0; i < v2d_1.dim( 0 ); i++ ) {
        for ( size_t j = 0; j < v2d_1.dim( 1 ); j++ ) {
            EXPECT_DOUBLE_EQ( v2d_1[ i ][ j ], (double)i * 10.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, AsArrayReturnsCorrectArray2D ) {
    ScalarWithUnits<double> **tester = nullptr;
    v2d_2.as_array( tester );
    for ( size_t i = 0; i < 5; i++ ) {
        for ( size_t j = 0; j < 5; j++ ) {
            EXPECT_DOUBLE_EQ( tester[ i ][ j ].get(), (double)i * 10.0 );
        }
    }
    delete[] tester;
}

TEST_F( _TEST_TITLE_, AsArrayReturnsCorrectDoubleArray2D ) {
    double **tester = nullptr;
    v2d_2.as_array( tester );
    for ( size_t i = 0; i < 5; i++ ) {
        for ( size_t j = 0; j < 5; j++ ) {
            EXPECT_DOUBLE_EQ( tester[ i ][ j ], (double)i * 10.0 );
        }
    }
    delete[] tester;
}

TEST_F( _TEST_TITLE_, AsArrayReturnsCorrectDoubleVector2D ) {
    NCPA::arrays::vector2d_t<double> tester = v2d_2;
    for ( size_t i = 0; i < 5; i++ ) {
        for ( size_t j = 0; j < 5; j++ ) {
            EXPECT_DOUBLE_EQ( tester[ i ][ j ], (double)i * 10.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, ConvertUnitsCreatesCorrectValues2D ) {
    v2d_2.set_units( KILOMETERS );
    v2d_2.convert_units( METERS );
    for ( size_t i = 0; i < 5; i++ ) {
        for ( size_t j = 0; j < 5; j++ ) {
            EXPECT_DOUBLE_EQ( v2d_2[ i ][ j ], (double)i * 10000.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, ConvertUnitsWithStringCreatesCorrectValues2D ) {
    v2d_2.set_units( "km" );
    v2d_2.convert_units( "m" );
    for ( size_t i = 0; i < 5; i++ ) {
        for ( size_t j = 0; j < 5; j++ ) {
            EXPECT_DOUBLE_EQ( v2d_2[ i ][ j ], (double)i * 10000.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, ConvertUnitsThrowsInvalidConversionCorrectly2D ) {
    EXPECT_THROW( { v2d_2.convert_units( KILOGRAMS ); }, invalid_conversion<> );
}

TEST_F( _TEST_TITLE_, FillSetsValuesCorrectly2D ) {
    v2d_2.fill( 12.0, Units::from_string( "kg/m3" ) );
    for ( size_t i = 0; i < 5; i++ ) {
        for ( size_t j = 0; j < 5; j++ ) {
            EXPECT_DOUBLE_EQ( v2d_2[ i ][ j ], 12.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, FillSetsUnitsCorrectly2D ) {
    const Unit *u = Units::from_string( "kg/m3" );
    v2d_2.fill( 12.0, u );
    EXPECT_TRUE( v2d_2.get_units()->equals( u ) );
}

TEST_F( _TEST_TITLE_, SetWorksCorrectly2D ) {
    ASSERT_EQ( v2d_2.dim(0), 5 );
    ASSERT_EQ( v2d_2.dim(1), 5 );
    NCPA::arrays::vector2d_t<double> z( 10, 10, 0.0 );
    v2d_2.set( z, CELSIUS );
    EXPECT_EQ( v2d_2.dim(0), 10 );
    EXPECT_EQ( v2d_2.dim(1), 10 );
    for ( size_t i = 0; i < v2d_2.dim( 0 ); i++ ) {
        for ( size_t j = 0; j < v2d_2.dim( 1 ); j++ ) {
            EXPECT_DOUBLE_EQ( v2d_2[ i ][ j ], 0.0 );
        }
    }
}
