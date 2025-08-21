#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <limits>
#include <numbers>
#include <utility>

using namespace std;
using namespace NCPA::math;
using namespace NCPA::arrays;
using namespace testing;

TEST( NCPAMathLibraryTest, SignIsCorrect ) {
    EXPECT_EQ( sign<int>( 3 ), 1 );
    EXPECT_EQ( sign<int>( -3 ), -1 );
    EXPECT_EQ( sign<int>( 0 ), 0 );
    EXPECT_EQ( sign<double>( 3.0 ), 1 );
    EXPECT_EQ( sign<double>( -3.2 ), -1 );
}

TEST( NCPAMathLibraryTest, FindIntervalInclusiveFindsInterval ) {
    vector<double> testArray = NCPA::arrays::index_vector<double>( 5 );
    size_t below, above;
    for (size_t i = 1; i < 5; i++) {
        double target = 0.5 * ( (double)i + (double)( i - 1 ) );
        EXPECT_TRUE( find_interval_inclusive<double>( as_array( testArray ), 5,
                                                      target, below, above ) );
        EXPECT_EQ( below, i - 1 );
        EXPECT_EQ( above, i );
    }
}

TEST( NCPAMathLibraryTest, FindIntervalInclusiveRecognizesOutOfRangeBelow ) {
    vector<double> testArray = NCPA::arrays::index_vector<double>( 5 );
    size_t below, above;
    double target = testArray[ 0 ] - 1.0;
    EXPECT_FALSE( find_interval_inclusive<double>( as_array( testArray ), 5,
                                                   target, below, above ) );
    EXPECT_EQ( below, 0 );
    EXPECT_EQ( above, 0 );
}

TEST( NCPAMathLibraryTest, FindIntervalInclusiveRecognizesOutOfRangeAbove ) {
    vector<double> testArray = NCPA::arrays::index_vector<double>( 5 );
    size_t below, above;
    double target = testArray[ 4 ] + 1.0;
    EXPECT_FALSE( find_interval_inclusive<double>( as_array( testArray ), 5,
                                                   target, below, above ) );
    EXPECT_EQ( below, 5 );
    EXPECT_EQ( above, 5 );
}

TEST( NCPAMathLibraryTest, FindClosestIndexWorksCorrectlyForArrays ) {
    vector<double> testArray = NCPA::arrays::index_vector<double>( 5 );
    for (size_t i = 0; i < 5; i++) {
        double target = testArray[ i ];
        EXPECT_EQ(
            find_closest_index<double>( as_array( testArray ), 5, target ),
            i );
        target -= 0.25;
        EXPECT_EQ(
            find_closest_index<double>( as_array( testArray ), 5, target ),
            i );
        target += 0.5;
        EXPECT_EQ(
            find_closest_index<double>( as_array( testArray ), 5, target ),
            i );
        if (i == 0) {
            target = -20000.0;
            EXPECT_EQ(
                find_closest_index<double>( as_array( testArray ), 5, target ),
                i );
        } else if (i == 4) {
            target = 20000.0;
            EXPECT_EQ(
                find_closest_index<double>( as_array( testArray ), 5, target ),
                i );
        } else {
            target += 0.5;
            EXPECT_EQ(
                find_closest_index<double>( as_array( testArray ), 5, target ),
                i + 1 );
            target -= 1.5;
            EXPECT_EQ(
                find_closest_index<double>( as_array( testArray ), 5, target ),
                i - 1 );
        }
    }
}

TEST( NCPAMathLibraryTest, FindClosestIndexWorksCorrectlyForVectors ) {
    vector<double> testArray = NCPA::arrays::index_vector<double>( 5 );
    for (size_t i = 0; i < 5; i++) {
        double target = testArray[ i ];
        EXPECT_EQ( find_closest_index<double>( testArray, target ), i );
        target -= 0.25;
        EXPECT_EQ( find_closest_index<double>( testArray, target ), i );
        target += 0.5;
        EXPECT_EQ( find_closest_index<double>( testArray, target ), i );
        if (i == 0) {
            target = -20000.0;
            EXPECT_EQ( find_closest_index<double>( testArray, target ), i );
        } else if (i == 4) {
            target = 20000.0;
            EXPECT_EQ( find_closest_index<double>( testArray, target ), i );
        } else {
            target += 0.5;
            EXPECT_EQ( find_closest_index<double>( testArray, target ),
                       i + 1 );
            target -= 1.5;
            EXPECT_EQ( find_closest_index<double>( testArray, target ),
                       i - 1 );
        }
    }
}

TEST( NCPAMathLibraryTest, FindClosestPointWorksCorrectly ) {
    vector<double> v1 = NCPA::arrays::index_vector<double>( 5 ),
                   v2 = NCPA::arrays::index_vector<double>( 5 );
    std::pair<size_t, size_t> coords;
    coords = find_closest_point( v1, v2, -0.25, -0.25 );
    EXPECT_EQ( coords.first, 0 );
    EXPECT_EQ( coords.second, 0 );

    coords = find_closest_point( v1, v2, -0.25, 2.75 );
    EXPECT_EQ( coords.first, 0 );
    EXPECT_EQ( coords.second, 3 );

    coords = find_closest_point( v1, v2, 2.75, -0.25 );
    EXPECT_EQ( coords.first, 3 );
    EXPECT_EQ( coords.second, 0 );

    coords = find_closest_point( v1, v2, 0.25, 0.25 );
    EXPECT_EQ( coords.first, 0 );
    EXPECT_EQ( coords.second, 0 );

    coords = find_closest_point( v1, v2, 1.25, 1.25 );
    EXPECT_EQ( coords.first, 1 );
    EXPECT_EQ( coords.second, 1 );

    coords = find_closest_point( v1, v2, 1.25, 2.999 );
    EXPECT_EQ( coords.first, 1 );
    EXPECT_EQ( coords.second, 3 );

    coords = find_closest_point( v1, v2, 6.0, 6.0 );
    EXPECT_EQ( coords.first, 4 );
    EXPECT_EQ( coords.second, 4 );
}

TEST( NCPAMathLibraryTest, Cart2PolConvertsCorrectly ) {
    double r, theta;

    cart2pol<double>( 1.0, 0.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 1.0 );
    EXPECT_DOUBLE_EQ( theta, 0.0 );

    cart2pol<double>( 1.0, 1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( theta, 0.25 * NCPA::math::PI );

    cart2pol<double>( 0.0, 1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 1.0 );
    EXPECT_DOUBLE_EQ( theta, 0.5 * NCPA::math::PI );

    cart2pol<double>( -1.0, 1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( theta, 0.75 * NCPA::math::PI );

    cart2pol<double>( -1.0, 0.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 1.0 );
    EXPECT_DOUBLE_EQ( std::fabs( theta ), NCPA::math::PI );

    cart2pol<double>( -1.0, -1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( theta, -0.75 * NCPA::math::PI );

    cart2pol<double>( 0.0, -1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 1.0 );
    EXPECT_DOUBLE_EQ( theta, -0.5 * NCPA::math::PI );

    cart2pol<double>( 1.0, -1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( theta, -0.25 * NCPA::math::PI );

    cart2pol<double>( 0.0, 0.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 0.0 );
    EXPECT_DOUBLE_EQ( theta, 0.0 );
}

TEST( NCPAMathLibraryTest, Pol2CartConvertsCorrectly ) {
    double x, y;

    pol2cart<double>( 1.0, 0.0, x, y );
    EXPECT_DOUBLE_EQ( x, 1.0 );
    EXPECT_NEAR( y, 0.0, 1e-12 );

    pol2cart<double>( 1.0, 0.25 * NCPA::math::PI, x, y );
    EXPECT_DOUBLE_EQ( x, 1.0 / std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( y, 1.0 / std::sqrt( 2.0 ) );

    pol2cart<double>( 1, 0.5 * NCPA::math::PI, x, y );
    EXPECT_NEAR( x, 0.0, 1e-12 );
    EXPECT_DOUBLE_EQ( y, 1.0 );

    pol2cart<double>( 1, 0.75 * NCPA::math::PI, x, y );
    EXPECT_DOUBLE_EQ( x, -1.0 / std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( y, 1.0 / std::sqrt( 2.0 ) );

    pol2cart<double>( 1.0, NCPA::math::PI, x, y );
    EXPECT_DOUBLE_EQ( x, -1.0 );
    EXPECT_NEAR( y, 0.0, 1e-12 );

    pol2cart<double>( 1, -0.75 * NCPA::math::PI, x, y );
    EXPECT_DOUBLE_EQ( x, -1.0 / std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( y, -1.0 / std::sqrt( 2.0 ) );

    pol2cart<double>( 1, -0.5 * NCPA::math::PI, x, y );
    EXPECT_NEAR( x, 0.0, 1e-12 );
    EXPECT_DOUBLE_EQ( y, -1.0 );

    pol2cart<double>( 1, -0.25 * NCPA::math::PI, x, y );
    EXPECT_DOUBLE_EQ( x, 1.0 / std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( y, -1.0 / std::sqrt( 2.0 ) );

    pol2cart<double>( 0.0, 0.0, x, y );
    EXPECT_DOUBLE_EQ( x, 0.0 );
    EXPECT_DOUBLE_EQ( y, 0.0 );
}

TEST( NCPAMathLibraryTest, Deg2RadConvertsCorrectly ) {
    for (int i = -8; i <= 8; i++) {
        double d = (double)i * 45.0;
        double r = (double)i * NCPA::math::PI / 4.0;
        EXPECT_DOUBLE_EQ( deg2rad<double>( d ), r );
    }
}

TEST( NCPAMathLibraryTest, Rad2DegConvertsCorrectly ) {
    for (int i = -8; i <= 8; i++) {
        double d = (double)i * 45.0;
        double r = (double)i * NCPA::math::PI / 4.0;
        EXPECT_DOUBLE_EQ( rad2deg<double>( r ), d );
    }
}

TEST( NCPAMathLibraryTest, Real2ComplexConvertsRealVectorsCorrectly ) {
    vector<double> d;
    for (int i = 0; i < 9; i++) {
        d.push_back( (double)i * NCPA::math::PI / 4.0 - NCPA::math::PI );
    }
    vector<complex<double>> c = real2complex( d );
    for (int i = 0; i < 9; i++) {
        EXPECT_EQ( d[ i ], c[ i ].real() );
        EXPECT_EQ( c[ i ].imag(), 0.0 );
    }
}

TEST( NCPAMathLibraryTest, real2complexConvertsRealArraysCorrectly ) {
    double *d          = NCPA::arrays::zeros<double>( 9 );
    complex<double> *c = NCPA::arrays::zeros<complex<double>>( 9 );

    for (int i = 0; i < 9; i++) {
        d[ i ] = (double)i * NCPA::math::PI / 4.0 - NCPA::math::PI;
    }
    real2complex( 9, d, c );
    for (int i = 0; i < 9; i++) {
        EXPECT_EQ( d[ i ], c[ i ].real() );
        EXPECT_EQ( c[ i ].imag(), 0.0 );
    }
}

TEST( NCPAMathLibraryTest, real2complexConvertsRealAndImagVectorsCorrectly ) {
    vector<double> r, i;
    for (int ii = 0; ii < 9; ii++) {
        r.push_back( (double)ii * NCPA::math::PI / 4.0 - NCPA::math::PI );
        i.push_back( (double)ii * NCPA::math::PI );
    }
    vector<complex<double>> c = real2complex( r, i );
    for (int ii = 0; ii < 9; ii++) {
        EXPECT_EQ( c[ ii ].real(), r[ ii ] );
        EXPECT_EQ( c[ ii ].imag(), i[ ii ] );
    }
}

TEST( NCPAMathLibraryTest, real2complexConvertsRealAndImagArraysCorrectly ) {
    double *r          = NCPA::arrays::zeros<double>( 9 );
    double *i          = NCPA::arrays::zeros<double>( 9 );
    complex<double> *c = NCPA::arrays::zeros<complex<double>>( 9 );

    for (int ii = 0; ii < 9; ii++) {
        r[ ii ] = (double)ii * NCPA::math::PI / 4.0 - NCPA::math::PI;
        i[ ii ] = (double)ii * NCPA::math::PI;
    }
    real2complex( 9, r, i, c );
    for (int ii = 0; ii < 9; ii++) {
        EXPECT_EQ( c[ ii ].real(), r[ ii ] );
        EXPECT_EQ( c[ ii ].imag(), i[ ii ] );
    }
}

TEST( NCPAMathLibraryTest, complex2realConvertsVectorsCorrectly ) {
    vector<complex<double>> c;
    for (int ii = 0; ii < 9; ii++) {
        c.push_back( complex<double>( (double)ii * NCPA::math::PI / 4.0
                                          - NCPA::math::PI,
                                      (double)ii * NCPA::math::PI ) );
    }
    vector<double> r, i;
    complex2real( c, r, i );
    for (int ii = 0; ii < 9; ii++) {
        EXPECT_EQ( c[ ii ].real(), r[ ii ] );
        EXPECT_EQ( c[ ii ].imag(), i[ ii ] );
    }
}

TEST( NCPAMathLibraryTest, complex2realConvertsArraysCorrectly ) {
    double *r          = NCPA::arrays::zeros<double>( 9 );
    double *i          = NCPA::arrays::zeros<double>( 9 );
    complex<double> *c = NCPA::arrays::zeros<complex<double>>( 9 );

    for (int ii = 0; ii < 9; ii++) {
        c[ ii ].real( (double)ii * NCPA::math::PI / 4.0 - NCPA::math::PI );
        c[ ii ].imag( (double)ii * NCPA::math::PI );
    }
    complex2real( 9, c, r, i );
    for (int ii = 0; ii < 9; ii++) {
        EXPECT_EQ( c[ ii ].real(), r[ ii ] );
        EXPECT_EQ( c[ ii ].imag(), i[ ii ] );
    }
}

TEST( NCPAMathLibraryTest, MaxWorksProperlyOnArrays ) {
    vector<int> testArray = NCPA::arrays::index_vector<int>( 10 );
    EXPECT_EQ( NCPA::math::max<int>( as_array( testArray ), 10 ), 9 );
    testArray[ 5 ] *= 3;
    EXPECT_EQ( NCPA::math::max<int>( as_array( testArray ), 10 ), 15 );
}

TEST( NCPAMathLibraryTest, MaxWorksProperlyOnVectors ) {
    vector<int> testVector = NCPA::arrays::index_vector<int>( 10 );
    EXPECT_EQ( NCPA::math::max<int>( testVector ), 9 );
    testVector[ 5 ] *= 3;
    EXPECT_EQ( NCPA::math::max<int>( testVector ), 15 );
}

TEST( NCPAMathLibraryTest, MinWorksProperlyOnArrays ) {
    vector<int> testArray = NCPA::arrays::index_vector<int>( 10 );
    EXPECT_EQ( NCPA::math::min<int>( as_array( testArray ), 10 ), 0 );
    testArray[ 5 ] *= -2;
    EXPECT_EQ( NCPA::math::min<int>( as_array( testArray ), 10 ), -10 );
}

TEST( NCPAMathLibraryTest, MinWorksProperlyOnVectors ) {
    vector<int> testVector = NCPA::arrays::index_vector<int>( 10 );
    EXPECT_EQ( NCPA::math::min<int>( testVector ), 0 );
    testVector[ 5 ] *= -2;
    EXPECT_EQ( NCPA::math::min<int>( testVector ), -10 );
}

TEST( NCPAMathLibraryTest, RandomNumbersAreInExpectedRange ) {
    vector<double> r = random_numbers<double>( 10, 0.0, 1.0 );
    for (auto i = 0; i < 10; i++) {
        EXPECT_LT( r[ i ], 1.0 );
        EXPECT_GE( r[ i ], 0.0 );
    }
    r = random_numbers<double>( 10, 0.0, 10.0 );
    for (auto i = 0; i < 10; i++) {
        EXPECT_LT( r[ i ], 10.0 );
        EXPECT_GE( r[ i ], 0.0 );
    }

    r = random_numbers<double>( 10, -5.0, 5.0 );
    for (auto i = 0; i < 10; i++) {
        EXPECT_LT( r[ i ], 5.0 );
        EXPECT_GE( r[ i ], -5.0 );
    }
}

TEST( NCPAMathLibraryTest, SingleRandomNumbersAreInExpectedRange ) {
    double minrange = -3.0, maxrange = 12.0;
    for (auto i = 0; i < 100; i++) {
        EXPECT_LT( random_number<double>( minrange, maxrange ), maxrange );
        EXPECT_GE( random_number<double>( minrange, maxrange ), minrange );
    }
}

TEST( NCPAMathLibraryTest, ComplexRandomNumbersAreInExpectedRange ) {
    vector<complex<double>> r
        = random_numbers<complex<double>>( 10, 0.0, 1.0 );
    for (auto i = 0; i < 10; i++) {
        EXPECT_LE( r[ i ].real(), 1.0 );
        EXPECT_GE( r[ i ].real(), 0.0 );
        EXPECT_LE( r[ i ].imag(), 1.0 );
        EXPECT_GE( r[ i ].imag(), 0.0 );
    }
    r = random_numbers<complex<double>>( 10, 0.0, 10.0 );
    for (auto i = 0; i < 10; i++) {
        EXPECT_LE( r[ i ].real(), 10.0 );
        EXPECT_GE( r[ i ].real(), 0.0 );
        EXPECT_LE( r[ i ].imag(), 10.0 );
        EXPECT_GE( r[ i ].imag(), 0.0 );
    }
    r = random_numbers<complex<double>>( 10, -5.0, 5.0 );
    for (auto i = 0; i < 10; i++) {
        EXPECT_LE( r[ i ].real(), 5.0 );
        EXPECT_GE( r[ i ].real(), -5.0 );
        EXPECT_LE( r[ i ].imag(), 5.0 );
        EXPECT_GE( r[ i ].imag(), -5.0 );
    }
}

TEST( NCPAMathLibraryTest, EvalpolyReturnsExpectedResult ) {
    size_t nx               = 20;
    size_t n_coeffs         = 4;
    vector<double> coeffs_r = random_numbers<double>( n_coeffs, -2.0, 2.0 );
    vector<double> x        = random_numbers<double>( nx, -10, 10 );
    for (size_t i = 0; i < nx; i++) {
        double result = 0.0;
        for (size_t j = 0; j < coeffs_r.size(); j++) {
            result += coeffs_r[ j ] * std::pow( x[ i ], (double)j );
        }
        EXPECT_DOUBLE_EQ( result, evalpoly<double>( coeffs_r, x[ i ] ) );
        EXPECT_DOUBLE_EQ(
            result, evalpoly<double>( n_coeffs, &coeffs_r[ 0 ], x[ i ] ) );
    }

    vector<complex<double>> coeffs_c
        = random_numbers<complex<double>>( n_coeffs, -2.0, 2.0 );
    vector<complex<double>> xc
        = random_numbers<complex<double>>( nx, -10, 10 );
    for (size_t i = 0; i < nx; i++) {
        complex<double> result = 0.0;
        for (size_t j = 0; j < coeffs_r.size(); j++) {
            result += coeffs_c[ j ] * std::pow( xc[ i ], (double)j );
        }
        complex<double> evalresult
            = evalpoly<complex<double>>( coeffs_c, xc[ i ] );
        EXPECT_DOUBLE_EQ( result.real(), evalresult.real() );
        EXPECT_DOUBLE_EQ( result.imag(), evalresult.imag() );
        evalresult
            = evalpoly<complex<double>>( n_coeffs, &coeffs_c[ 0 ], xc[ i ] );
        EXPECT_DOUBLE_EQ( result.real(), evalresult.real() );
        EXPECT_DOUBLE_EQ( result.imag(), evalresult.imag() );
    }
}

TEST( NCPAMathLibraryTest, LinearInterpWorksCorrectly ) {
    vector<double> coords = random_numbers<double>( 4, -5, 5 );
    double dy             = coords[ 3 ] - coords[ 1 ];
    double dx             = coords[ 2 ] - coords[ 0 ];
    for (size_t i = 0; i < 5; i++) {
        double frac = (double)i / 4.0;
        EXPECT_NEAR( linear_interpolation<double>( coords[ 0 ], coords[ 1 ],
                                                   coords[ 2 ], coords[ 3 ],
                                                   coords[ 0 ] + frac * dx ),
                     coords[ 1 ] + dy * frac, 1.0e-10 );
    }
}

TEST( NCPAMathLibraryTest, LinspaceWorksAsExpectedForDoubles ) {
    vector<double> tester = linspace<double>( 0.0, 5.0, 21 );
    double *shouldBe      = NCPA::arrays::zeros<double>( 21 );
    for (size_t i = 0; i < 21; i++) {
        shouldBe[ i ] = (double)i * 0.25;
    }
    EXPECT_ARRAY_DOUBLE_EQ( 21, tester, shouldBe );
    linspace<double>( 0.0, 5.0, 21, shouldBe );
    EXPECT_ARRAY_DOUBLE_EQ( 21, tester, shouldBe );
}

TEST( NCPAMathLibraryTest, LinspaceWorksAsExpectedForIntegers ) {
    vector<int> tester = linspace<int>( 0, 10, 11 );
    int *shouldBe      = NCPA::arrays::zeros<int>( 11 );
    for (size_t i = 0; i < 11; i++) {
        shouldBe[ i ] = (int)i;
    }
    EXPECT_ARRAY_EQ( 11, tester, shouldBe );
    linspace<int>( 0, 10, 11, shouldBe );
    EXPECT_ARRAY_EQ( 11, tester, shouldBe );
}

TEST( NCPAMathLibraryTest, LogspaceWorksAsExpected ) {
    vector<double> tester = logspace<double>( 1.0, 1000.0, 21 );
    double *shouldBe      = NCPA::arrays::zeros<double>( 21 );
    linspace<double>( 0.0, 3.0, 21, shouldBe );
    for (size_t i = 0; i < 21; i++) {
        shouldBe[ i ] = std::pow( 10.0, shouldBe[ i ] );
    }
    EXPECT_ARRAY_DOUBLE_EQ( 21, tester, shouldBe );
    logspace<double>( 1.0, 1000.0, 21, shouldBe );
    EXPECT_ARRAY_DOUBLE_EQ( 21, tester, shouldBe );
}

TEST( NCPAMathLibraryTest, MeanWorksAsExpected ) {
    size_t N                 = 10;
    double dN                = (double)N;
    vector<double> testArray = random_numbers<double>( N, 0.0, 5.0 );
    double testmean          = 0.0;
    for (size_t i = 0; i < N; i++) {
        testmean += testArray[ i ] / dN;
    }
    EXPECT_DOUBLE_EQ( testmean, mean( testArray ) );
}

TEST( NCPAMathLibraryTest, NextPow2ReturnsCorrectPowers ) {
    EXPECT_DOUBLE_EQ( nextpow2<double>( 1.5 ), 1.0 );
    EXPECT_DOUBLE_EQ( nextpow2<double>( 2.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( nextpow2<double>( 1024.0 ), 10.0 );
    EXPECT_DOUBLE_EQ( nextpow2<double>( 1024.0000000001 ), 11.0 );
    EXPECT_EQ( nextpow2<int>( 1 ), 0 );
    EXPECT_EQ( nextpow2<int>( 8 ), 3 );
}

TEST( NCPAMathLibraryTest, NextPow2ReturnsNaNForNonPositiveArgument ) {
    EXPECT_TRUE( std::isnan( nextpow2<double>( -1.0 ) ) );
}

TEST( NCPAMathLibraryTest, TrapzComputesIntegralCorrectly ) {
    double x[ 5 ]   = { 0, 1, 3, 4, 5 };
    double y[ 5 ]   = { 0, 2, 6, 0, -2 };
    double expected = 1.0 + 8.0 + 3.0 - 1.0;
    EXPECT_DOUBLE_EQ( expected, trapz<double>( 5, x, y ) );
}

TEST( NCPAMathLibraryTest, IsZeroReturnsTrueForZeroValue ) {
    float f        = 0.0;
    double d       = 0.0;
    long double ld = 0.0;
    complex<float> cf;
    complex<double> cd;
    complex<long double> cld;
    EXPECT_TRUE( is_zero( f ) );
    EXPECT_TRUE( is_zero( d ) );
    EXPECT_TRUE( is_zero( ld ) );
    EXPECT_TRUE( is_zero( cf ) );
    EXPECT_TRUE( is_zero( cd ) );
    EXPECT_TRUE( is_zero( cld ) );
}

TEST( NCPAMathLibraryTest, IsZeroReturnsFalseForNonzeroValue ) {
    float f        = std::numeric_limits<float>::epsilon();
    double d       = std::numeric_limits<double>::epsilon();
    long double ld = std::numeric_limits<long double>::epsilon();
    complex<float> cf( f, f );
    complex<double> cd( d, d );
    complex<long double> cld( ld, ld );
    EXPECT_FALSE( is_zero( f ) );
    EXPECT_FALSE( is_zero( d ) );
    EXPECT_FALSE( is_zero( ld ) );
    EXPECT_FALSE( is_zero( cf ) );
    EXPECT_FALSE( is_zero( cd ) );
    EXPECT_FALSE( is_zero( cld ) );
}

TEST( NCPAMathLibraryTest, ComplexRootFinderWorks ) {
    std::vector<double> testroots_real { 1.0, 3.0, 5.0, 7.0, 9.0 };
    std::vector<double> testcoeffs_real { 1.0,    -25.0,  230.0,
                                          -950.0, 1689.0, -945 };
    ComplexRootFinder finder;
    std::vector<std::complex<double>> roots
        = finder.find( testcoeffs_real, true );
    std::vector<double> roots_real;
    for (size_t i = 0; i < roots.size(); ++i) {
        roots_real.push_back( roots[ i ].real() );
    }
    EXPECT_EQ( roots_real.size(), testroots_real.size() );
    std::sort( roots_real.begin(), roots_real.end() );
    for (size_t i = 0; i < roots_real.size(); ++i) {
        EXPECT_NEAR( testroots_real[ i ], roots_real[ i ], 1.0e-8 );
    }

    std::vector<std::complex<double>> testroots_complex {
        {  1.0, -1.0 },
        {  2.0, -2.0 },
        {  3.5, -4.5 },
        {  0.5, -1.3 },
        { -4.6, -2.3 }
    };
    std::vector<std::complex<double>> testcoeffs_complex {
        {    0.01,     0.0 },
        {  -0.024,   0.111 },
        { -0.6194, -0.1582 },
        {  1.0466, -1.8267 },
        {  1.7451,  2.5583 },
        { -1.6284,  0.1288 }
    };
    roots = finder.find( testcoeffs_complex, true );
    EXPECT_EQ( roots.size(), testroots_complex.size() );
    std::vector<bool> found( roots.size(), false );
    for (size_t i = 0; i < roots.size(); ++i) {
        for (size_t j = 0; j < roots.size(); ++j) {
            if (std::abs( testroots_complex[ i ].real() - roots[ j ].real() )
                    <= 1.0e-8
                && std::abs( testroots_complex[ i ].imag()
                             - roots[ j ].imag() )
                       <= 1.0e-8) {
                        found[i] = true;
                       }
        }
    }
    for (size_t i = 0; i < roots.size(); ++i) {
        EXPECT_TRUE( found[ i ] );
    }
}
