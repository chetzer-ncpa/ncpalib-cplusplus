#include "NCPA/gtest.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <complex>
#include <numbers>
#include <algorithm>

using namespace std;
using namespace NCPA::math;
using namespace testing;

TEST( NCPAmathTest, ZerosCreatesArray ) {
    double *testArray     = zeros<double>( 5 );
    double  shouldBe[ 5 ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    ASSERT_ARRAY_EQ( 5, testArray, shouldBe );
    // for (auto i = 0; i < 10; i++) {
    // 	EXPECT_EQ( testArray[i], 0.0 );
    // }
    delete[] testArray;
}

TEST( NCPAmathTest, ZerosCreates2DArray ) {
    double **testArray = zeros<double>( 10, 5 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            EXPECT_EQ( testArray[ i ][ j ], 0.0 );
        }
    }
}

TEST( NCPAmathTest, FreeArray2DFreesArrayWithoutError ) {
    double **testArray = zeros<double>( 10, 5 );
    free_array<double>( testArray, 10, 5 );
    EXPECT_EQ( testArray, nullptr );
}

TEST( NCPAmathTest, ZerosCreates3DArray ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            for ( auto k = 0; k < 3; k++ ) {
                EXPECT_EQ( testArray[ i ][ j ][ k ], 0.0 );
            }
        }
    }
}

TEST( NCPAmathTest, FreeArray3DFreesArrayWithoutError ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    free_array<double>( testArray, 10, 5, 3 );
    ASSERT_EQ( testArray, nullptr );
}

TEST( NCPAmathTest, IndexVectorCreatesExpectedVectors ) {
    int *testArray = index_vector<int>( 5 );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_EQ( i, testArray[ i ] );
    }
    double *testArray2 = index_vector<double>( 5 );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( (double)i, testArray2[ i ] );
    }
}

TEST( NCPAmathTest, IndexVectorCreatesExpectedVectorsWithOffset ) {
    int  offset    = 2;
    int *testArray = index_vector<int>( 5, offset );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_EQ( i + offset, testArray[ i ] );
    }
    double *testArray2 = index_vector<double>( 5, (double)offset );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( (double)( i + offset ), testArray2[ i ] );
    }
}

TEST( NCPAmathTest, SignIsCorrect ) {
    EXPECT_EQ( sign<int>( 3 ), 1 );
    EXPECT_EQ( sign<int>( -3 ), -1 );
    EXPECT_EQ( sign<int>( 0 ), 0 );
    EXPECT_EQ( sign<double>( 3.0 ), 1 );
    EXPECT_EQ( sign<double>( -3.2 ), -1 );
}

TEST( NCPAmathTest, CircShiftShiftsForward ) {
    int *testArray     = index_vector<int>( 5 );
    int  shouldBe[ 5 ] = { 2, 3, 4, 0, 1 };
    circshift( testArray, 5, 2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAmathTest, CircShiftShiftsBackward ) {
    int *testArray     = index_vector<int>( 5 );
    int  shouldBe[ 5 ] = { 3, 4, 0, 1, 2 };
    circshift( testArray, 5, -2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAmathTest, CircShiftZeroDoesNotShift ) {
    int *testArray     = index_vector<int>( 5 );
    int  shouldBe[ 5 ] = { 0, 1, 2, 3, 4 };
    circshift( testArray, 5, 0, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAmathTest, FillArray2DSetsValuesCorrectly ) {
    double **testArray = zeros<double>( 10, 5 );
    fill_array2d<double>( testArray, 10, 5, 4.5 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            EXPECT_DOUBLE_EQ( testArray[ i ][ j ], 4.5 );
        }
    }
    free_array( testArray, 10, 5 );
}

TEST( NCPAmathTest, FillArray3DSetsValuesCorrectly ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    fill_array3d<double>( testArray, 10, 5, 3, 4.5 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            for ( auto k = 0; k < 3; k++ ) {
                EXPECT_DOUBLE_EQ( testArray[ i ][ j ][ k ], 4.5 );
            }
        }
    }
    free_array( testArray, 10, 5, 3 );
}

TEST( NCPAmathTest, FindIntervalInclusiveFindsInterval ) {
    double *testArray = index_vector<double>( 5 );
    size_t  below, above;
    for ( size_t i = 1; i < 5; i++ ) {
        double target = 0.5 * ( (double)i + (double)( i - 1 ) );
        EXPECT_TRUE( find_interval_inclusive<double>( testArray, 5, target,
                                                      below, above ) );
        EXPECT_EQ( below, i - 1 );
        EXPECT_EQ( above, i );
    }
}

TEST( NCPAmathTest, FindIntervalInclusiveRecognizesOutOfRangeBelow ) {
    double *testArray = index_vector<double>( 5 );
    size_t  below, above;
    double  target = testArray[ 0 ] - 1.0;
    EXPECT_FALSE( find_interval_inclusive<double>( testArray, 5, target, below,
                                                   above ) );
    EXPECT_EQ( below, 0 );
    EXPECT_EQ( above, 0 );
}

TEST( NCPAmathTest, FindIntervalInclusiveRecognizesOutOfRangeAbove ) {
    double *testArray = index_vector<double>( 5 );
    size_t  below, above;
    double  target = testArray[ 4 ] + 1.0;
    EXPECT_FALSE( find_interval_inclusive<double>( testArray, 5, target, below,
                                                   above ) );
    EXPECT_EQ( below, 5 );
    EXPECT_EQ( above, 5 );
}

TEST( NCPAmathTest, FindClosestIndexWorksCorrectlyForArrays ) {
    double *testArray = index_vector<double>( 5 );
    for ( size_t i = 0; i < 5; i++ ) {
        double target = testArray[ i ];
        EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
        target -= 0.25;
        EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
        target += 0.5;
        EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
        if ( i == 0 ) {
            target = -20000.0;
            EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
        } else if ( i == 4 ) {
            target = 20000.0;
            EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
        } else {
            target += 0.5;
            EXPECT_EQ( find_closest_index<double>( testArray, 5, target ),
                       i + 1 );
            target -= 1.5;
            EXPECT_EQ( find_closest_index<double>( testArray, 5, target ),
                       i - 1 );
        }
    }
}

TEST( NCPAmathTest, FindClosestIndexWorksCorrectlyForVectors ) {
    double        *tempArray = index_vector<double>( 5 );
    vector<double> testArray( tempArray, tempArray + 5 );
    for ( size_t i = 0; i < 5; i++ ) {
        double target = testArray[ i ];
        EXPECT_EQ( find_closest_index<double>( testArray, target ), i );
        target -= 0.25;
        EXPECT_EQ( find_closest_index<double>( testArray, target ), i );
        target += 0.5;
        EXPECT_EQ( find_closest_index<double>( testArray, target ), i );
        if ( i == 0 ) {
            target = -20000.0;
            EXPECT_EQ( find_closest_index<double>( testArray, target ), i );
        } else if ( i == 4 ) {
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

TEST( NCPAmathTest, Cart2PolConvertsCorrectly ) {
    double r, theta;

    cart2pol<double>( 1.0, 0.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 1.0 );
    EXPECT_DOUBLE_EQ( theta, 0.0 );

    cart2pol<double>( 1.0, 1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( theta, 0.25 * numbers::pi );

    cart2pol<double>( 0.0, 1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 1.0 );
    EXPECT_DOUBLE_EQ( theta, 0.5 * numbers::pi );

    cart2pol<double>( -1.0, 1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( theta, 0.75 * numbers::pi );

    cart2pol<double>( -1.0, 0.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 1.0 );
    EXPECT_DOUBLE_EQ( std::fabs( theta ), numbers::pi );

    cart2pol<double>( -1.0, -1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( theta, -0.75 * numbers::pi );

    cart2pol<double>( 0.0, -1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 1.0 );
    EXPECT_DOUBLE_EQ( theta, -0.5 * numbers::pi );

    cart2pol<double>( 1.0, -1.0, r, theta );
    EXPECT_DOUBLE_EQ( r, std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( theta, -0.25 * numbers::pi );

    cart2pol<double>( 0.0, 0.0, r, theta );
    EXPECT_DOUBLE_EQ( r, 0.0 );
    EXPECT_DOUBLE_EQ( theta, 0.0 );
}

TEST( NCPAmathTest, Pol2CartConvertsCorrectly ) {
    double x, y;

    pol2cart<double>( 1.0, 0.0, x, y );
    EXPECT_DOUBLE_EQ( x, 1.0 );
    EXPECT_NEAR( y, 0.0, 1e-12 );

    pol2cart<double>( 1.0, 0.25 * numbers::pi, x, y );
    EXPECT_DOUBLE_EQ( x, 1.0 / std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( y, 1.0 / std::sqrt( 2.0 ) );

    pol2cart<double>( 1, 0.5 * numbers::pi, x, y );
    EXPECT_NEAR( x, 0.0, 1e-12 );
    EXPECT_DOUBLE_EQ( y, 1.0 );

    pol2cart<double>( 1, 0.75 * numbers::pi, x, y );
    EXPECT_DOUBLE_EQ( x, -1.0 / std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( y, 1.0 / std::sqrt( 2.0 ) );

    pol2cart<double>( 1.0, numbers::pi, x, y );
    EXPECT_DOUBLE_EQ( x, -1.0 );
    EXPECT_NEAR( y, 0.0, 1e-12 );

    pol2cart<double>( 1, -0.75 * numbers::pi, x, y );
    EXPECT_DOUBLE_EQ( x, -1.0 / std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( y, -1.0 / std::sqrt( 2.0 ) );

    pol2cart<double>( 1, -0.5 * numbers::pi, x, y );
    EXPECT_NEAR( x, 0.0, 1e-12 );
    EXPECT_DOUBLE_EQ( y, -1.0 );

    pol2cart<double>( 1, -0.25 * numbers::pi, x, y );
    EXPECT_DOUBLE_EQ( x, 1.0 / std::sqrt( 2.0 ) );
    EXPECT_DOUBLE_EQ( y, -1.0 / std::sqrt( 2.0 ) );

    pol2cart<double>( 0.0, 0.0, x, y );
    EXPECT_DOUBLE_EQ( x, 0.0 );
    EXPECT_DOUBLE_EQ( y, 0.0 );
}

TEST( NCPAmathTest, Deg2RadConvertsCorrectly ) {
    for ( int i = -8; i <= 8; i++ ) {
        double d = (double)i * 45.0;
        double r = (double)i * numbers::pi / 4.0;
        EXPECT_DOUBLE_EQ( deg2rad<double>( d ), r );
    }
}

TEST( NCPAmathTest, Rad2DegConvertsCorrectly ) {
    for ( int i = -8; i <= 8; i++ ) {
        double d = (double)i * 45.0;
        double r = (double)i * numbers::pi / 4.0;
        EXPECT_DOUBLE_EQ( rad2deg<double>( r ), d );
    }
}

TEST( NCPAmathTest, Double2ComplexConvertsRealVectorsCorrectly ) {
    vector<double> d;
    for ( int i = 0; i < 9; i++ ) {
        d.push_back( (double)i * numbers::pi / 4.0 - numbers::pi );
    }
    vector<complex<double>> c = double2complex( d );
    for ( int i = 0; i < 9; i++ ) {
        EXPECT_EQ( d[ i ], c[ i ].real() );
        EXPECT_EQ( c[ i ].imag(), 0.0 );
    }
}

TEST( NCPAmathTest, Double2ComplexConvertsRealArraysCorrectly ) {
    double          *d = zeros<double>( 9 );
    complex<double> *c = zeros<complex<double>>( 9 );

    for ( int i = 0; i < 9; i++ ) {
        d[ i ] = (double)i * numbers::pi / 4.0 - numbers::pi;
    }
    double2complex( 9, d, c );
    for ( int i = 0; i < 9; i++ ) {
        EXPECT_EQ( d[ i ], c[ i ].real() );
        EXPECT_EQ( c[ i ].imag(), 0.0 );
    }
}

TEST( NCPAmathTest, Double2ComplexConvertsRealAndImagVectorsCorrectly ) {
    vector<double> r, i;
    for ( int ii = 0; ii < 9; ii++ ) {
        r.push_back( (double)ii * numbers::pi / 4.0 - numbers::pi );
        i.push_back( (double)ii * numbers::pi );
    }
    vector<complex<double>> c = double2complex( r, i );
    for ( int ii = 0; ii < 9; ii++ ) {
        EXPECT_EQ( c[ ii ].real(), r[ ii ] );
        EXPECT_EQ( c[ ii ].imag(), i[ ii ] );
    }
}

TEST( NCPAmathTest, Double2ComplexConvertsRealAndImagArraysCorrectly ) {
    double          *r = zeros<double>( 9 );
    double          *i = zeros<double>( 9 );
    complex<double> *c = zeros<complex<double>>( 9 );

    for ( int ii = 0; ii < 9; ii++ ) {
        r[ ii ] = (double)ii * numbers::pi / 4.0 - numbers::pi;
        i[ ii ] = (double)ii * numbers::pi;
    }
    double2complex( 9, r, i, c );
    for ( int ii = 0; ii < 9; ii++ ) {
        EXPECT_EQ( c[ ii ].real(), r[ ii ] );
        EXPECT_EQ( c[ ii ].imag(), i[ ii ] );
    }
}

TEST( NCPAmathTest, Complex2DoubleConvertsVectorsCorrectly ) {
    vector<complex<double>> c;
    for ( int ii = 0; ii < 9; ii++ ) {
        c.push_back(
            complex<double>( (double)ii * numbers::pi / 4.0 - numbers::pi,
                             (double)ii * numbers::pi ) );
    }
    vector<double> r, i;
    complex2double( c, r, i );
    for ( int ii = 0; ii < 9; ii++ ) {
        EXPECT_EQ( c[ ii ].real(), r[ ii ] );
        EXPECT_EQ( c[ ii ].imag(), i[ ii ] );
    }
}

TEST( NCPAmathTest, Complex2DoubleConvertsArraysCorrectly ) {
    double          *r = zeros<double>( 9 );
    double          *i = zeros<double>( 9 );
    complex<double> *c = zeros<complex<double>>( 9 );

    for ( int ii = 0; ii < 9; ii++ ) {
        c[ ii ].real( (double)ii * numbers::pi / 4.0 - numbers::pi );
        c[ ii ].imag( (double)ii * numbers::pi );
    }
    complex2double( 9, c, r, i );
    for ( int ii = 0; ii < 9; ii++ ) {
        EXPECT_EQ( c[ ii ].real(), r[ ii ] );
        EXPECT_EQ( c[ ii ].imag(), i[ ii ] );
    }
}

TEST( NCPAmathTest, MaxWorksProperlyOnArrays ) {
    int *testArray = index_vector<int>( 10 );
    EXPECT_EQ( NCPA::math::max<int>( testArray, 10 ), 9 );
    testArray[ 5 ] *= 3;
    EXPECT_EQ( NCPA::math::max<int>( testArray, 10 ), 15 );
}

TEST( NCPAmathTest, MaxWorksProperlyOnVectors ) {
    int        *testArray = index_vector<int>( 10 );
    vector<int> testVector( testArray, testArray + 10 );
    EXPECT_EQ( NCPA::math::max<int>( testVector ), 9 );
    testVector[ 5 ] *= 3;
    EXPECT_EQ( NCPA::math::max<int>( testVector ), 15 );
}

TEST( NCPAmathTest, MinWorksProperlyOnArrays ) {
    int *testArray = index_vector<int>( 10 );
    EXPECT_EQ( NCPA::math::min<int>( testArray, 10 ), 0 );
    testArray[ 5 ] *= -2;
    EXPECT_EQ( NCPA::math::min<int>( testArray, 10 ), -10 );
}

TEST( NCPAmathTest, MinWorksProperlyOnVectors ) {
    int        *testArray = index_vector<int>( 10 );
    vector<int> testVector( testArray, testArray + 10 );
    EXPECT_EQ( NCPA::math::min<int>( testVector ), 0 );
    testVector[ 5 ] *= -2;
    EXPECT_EQ( NCPA::math::min<int>( testVector ), -10 );
}

TEST( NCPAmathTest, RandomNumbersAreInExpectedRange ) {
    vector<double> r = random_numbers<double>( 10, 0.0, 1.0 );
    for ( auto i = 0; i < 10; i++ ) {
        EXPECT_LT( r[ i ], 1.0 );
        EXPECT_GE( r[ i ], 0.0 );
    }
    r = random_numbers<double>( 10, 0.0, 10.0 );
    for ( auto i = 0; i < 10; i++ ) {
        EXPECT_LT( r[ i ], 10.0 );
        EXPECT_GE( r[ i ], 0.0 );
    }

    r = random_numbers<double>( 10, -5.0, 5.0 );
    for ( auto i = 0; i < 10; i++ ) {
        EXPECT_LT( r[ i ], 5.0 );
        EXPECT_GE( r[ i ], -5.0 );
    }
}

TEST( NCPAmathTest, SingleRandomNumbersAreInExpectedRange ) {
    double minrange = -3.0, maxrange = 12.0;
    for ( auto i = 0; i < 100; i++ ) {
        EXPECT_LT( random_number<double>( minrange, maxrange ), maxrange );
        EXPECT_GE( random_number<double>( minrange, maxrange ), minrange );
    }
}

TEST( NCPAmathTest, ComplexRandomNumbersAreInExpectedRange ) {
    vector<complex<double>> r = random_numbers<double>( 10, 0.0, 1.0, 0.0, 1.0 );
    for ( auto i = 0; i < 10; i++ ) {
        EXPECT_LE( r[ i ].real(), 1.0 );
        EXPECT_GE( r[ i ].real(), 0.0 );
        EXPECT_LE( r[ i ].imag(), 1.0 );
        EXPECT_GE( r[ i ].imag(), 0.0 );
    }
    r = random_numbers<double>( 10, 0.0, 10.0, 0.0, 10.0 );
    for ( auto i = 0; i < 10; i++ ) {
        EXPECT_LE( r[ i ].real(), 10.0 );
        EXPECT_GE( r[ i ].real(), 0.0 );
        EXPECT_LE( r[ i ].imag(), 10.0 );
        EXPECT_GE( r[ i ].imag(), 0.0 );
    }
    r = random_numbers<double>( 10, -5.0, 5.0, -5.0, 5.0 );
    for ( auto i = 0; i < 10; i++ ) {
        EXPECT_LE( r[ i ].real(), 5.0 );
        EXPECT_GE( r[ i ].real(), -5.0 );
        EXPECT_LE( r[ i ].imag(), 5.0 );
        EXPECT_GE( r[ i ].imag(), -5.0 );
    }
}

TEST( NCPAmathTest, EvalpolyReturnsExpectedResult ) {
    size_t         nx       = 20;
    size_t         n_coeffs = 4;
    vector<double> coeffs_r = random_numbers<double>( n_coeffs, -2.0, 2.0 );
    vector<double> x        = random_numbers<double>( nx, -10, 10 );
    for ( size_t i = 0; i < nx; i++ ) {
        double result = 0.0;
        for ( size_t j = 0; j < coeffs_r.size(); j++ ) {
            result += coeffs_r[ j ] * std::pow( x[ i ], (double)j );
        }
        EXPECT_DOUBLE_EQ( result, evalpoly<double>( coeffs_r, x[ i ] ) );
        EXPECT_DOUBLE_EQ(
            result, evalpoly<double>( n_coeffs, &coeffs_r[ 0 ], x[ i ] ) );
    }

    vector<complex<double>> coeffs_c
        = random_numbers<double>( n_coeffs, -2.0, 2.0, -2.0, 2.0 );
    vector<complex<double>> xc
        = random_numbers<double>( nx, -10, 10, -10, 10 );
    for ( size_t i = 0; i < nx; i++ ) {
        complex<double> result = 0.0;
        for ( size_t j = 0; j < coeffs_r.size(); j++ ) {
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

TEST( NCPAmathTest, LinearInterpWorksCorrectly ) {
    vector<double> coords = random_numbers<double>( 4, -5, 5 );
    double         dy     = coords[ 3 ] - coords[ 1 ];
    double         dx     = coords[ 2 ] - coords[ 0 ];
    for ( size_t i = 0; i < 5; i++ ) {
        double frac = (double)i / 4.0;
        EXPECT_NEAR( linear_interpolation<double>( coords[ 0 ], coords[ 1 ],
                                                   coords[ 2 ], coords[ 3 ],
                                                   coords[ 0 ] + frac * dx ),
                     coords[ 1 ] + dy * frac, 1.0e-10 );
    }
}

TEST( NCPAmathTest, LinspaceWorksAsExpectedForDoubles ) {
    vector<double> tester   = linspace<double>( 0.0, 5.0, 21 );
    double        *shouldBe = zeros<double>( 21 );
    for ( size_t i = 0; i < 21; i++ ) {
        shouldBe[ i ] = (double)i * 0.25;
    }
    EXPECT_ARRAY_DOUBLE_EQ( 21, tester, shouldBe );
    linspace<double>( 0.0, 5.0, 21, shouldBe );
    EXPECT_ARRAY_DOUBLE_EQ( 21, tester, shouldBe );
}

TEST( NCPAmathTest, LinspaceWorksAsExpectedForIntegers ) {
    vector<int> tester   = linspace<int>( 0, 10, 11 );
    int        *shouldBe = zeros<int>( 11 );
    for ( size_t i = 0; i < 11; i++ ) {
        shouldBe[ i ] = (int)i;
    }
    EXPECT_ARRAY_EQ( 11, tester, shouldBe );
    linspace<int>( 0, 10, 11, shouldBe );
    EXPECT_ARRAY_EQ( 11, tester, shouldBe );
}

TEST( NCPAmathTest, LogspaceWorksAsExpected ) {
    vector<double> tester   = logspace<double>( 1.0, 1000.0, 21 );
    double        *shouldBe = zeros<double>( 21 );
    linspace<double>( 0.0, 3.0, 21, shouldBe );
    for ( size_t i = 0; i < 21; i++ ) {
        shouldBe[ i ] = std::pow( 10.0, shouldBe[ i ] );
    }
    EXPECT_ARRAY_DOUBLE_EQ( 21, tester, shouldBe );
    logspace<double>( 1.0, 1000.0, 21, shouldBe );
    EXPECT_ARRAY_DOUBLE_EQ( 21, tester, shouldBe );
}

TEST( NCPAmathTest, MeanWorksAsExpected ) {
	size_t N = 10;
	double dN = (double)N;
	vector<double> testArray = random_numbers<double>( N, 0.0, 5.0 );
	double testmean = 0.0;
	for (size_t i = 0; i < N; i++) {
		testmean += testArray[i] / dN;
	}
	EXPECT_DOUBLE_EQ( testmean, mean( testArray ) );
}

TEST(NCPAmathTest, NextPow2ReturnsCorrectPowers) {
	EXPECT_DOUBLE_EQ( nextpow2<double>( 1.5 ), 1.0 );
	EXPECT_DOUBLE_EQ( nextpow2<double>( 2.0 ), 1.0 );
	EXPECT_DOUBLE_EQ( nextpow2<double>( 1024.0 ), 10.0 );
	EXPECT_DOUBLE_EQ( nextpow2<double>( 1024.0000000001 ), 11.0 );
	EXPECT_EQ( nextpow2<int>( 1 ), 0 );
	EXPECT_EQ( nextpow2<int>( 8 ), 3 );
}

TEST(NCPAmathTest, NextPow2ReturnsNaNForNonPositiveArgument) {
	EXPECT_TRUE( std::isnan( nextpow2<double>( -1.0 )));
}

TEST(NCPAmathTest, ReverseReversesOrder) {
	double *testArray = index_vector<double>( 11 );
	vector<double> testVector( testArray, testArray + 11 );
	std::reverse( testVector.begin(), testVector.end() );
	NCPA::math::reverse( testArray, 11, testArray );
	EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST(NCPAmathTest, ReverseReversesOrderPartially) {
	double *testArray = index_vector<double>( 11 );
	vector<double> testVector( testArray, testArray + 11 );
	std::reverse( testVector.begin(), testVector.begin() + 5 );
	NCPA::math::reverse( testArray, 5, testArray );
	EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST(NCPAmathTest,TrapzComputesIntegralCorrectly) {
	double x[ 5 ] = {0, 1, 3, 4, 5 };
	double y[ 5 ] = { 0, 2, 6, 0, -2 };
	double expected = 1.0 + 8.0 + 3.0 - 1.0;
	EXPECT_DOUBLE_EQ( expected, trapz<double>( 5, x, y ) );
}

TEST( NCPAmathTest, AddVectorsWorksCorrectlyForVectors ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z = add_vectors( x, y );
    for ( auto i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( x[ i ] + y[ i ], z[ i ] );
    }
    for ( auto i = 5; i < 8; i++ ) {
        EXPECT_DOUBLE_EQ( y[ i ], z[ i ] );
    }
}

TEST( NCPAmathTest, AddVectorsCommutesForVectors ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z1 = add_vectors( x, y );
    vector<double> z2 = add_vectors( y, x );
    EXPECT_ARRAY_DOUBLE_EQ( 8, z1, z2 );
}

TEST( NCPAmathTest, AddVectorsWorksCorrectlyForArrays ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 5, 0.0, 1.0 );
    double *z = zeros<double>( 5 );
    add_vectors( 5, &x[0], &y[0], z );
    for ( auto i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( x[ i ] + y[ i ], z[ i ] );
    }
}

TEST( NCPAmathTest, AddVectorsCommutesForArrays ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 5, 0.0, 1.0 );
    double *z1 = zeros<double>( 5 );
    double *z2 = zeros<double>( 5 );
    add_vectors( 5, &x[0], &y[0], z1 );
    add_vectors( 5, &y[0], &x[0], z2 );
    EXPECT_ARRAY_DOUBLE_EQ( 5, z1, z2 );
}


TEST( NCPAmathTest, DivideVectorsWorksCorrectlyForVectors ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 8, 0.1, 1.0 );
    vector<double> z = divide_vectors( x, y );
    for ( auto i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( x[ i ] / y[ i ], z[ i ] );
    }
    for ( auto i = 5; i < 8; i++ ) {
        EXPECT_DOUBLE_EQ( z[ i ], 0.0 );
    }
}

TEST( NCPAmathTest, DivideVectorsWorksCorrectlyForArrays ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 5, 0.1, 1.0 );
    double *z = zeros<double>( 5 );
    divide_vectors( 5, &x[0], &y[0], z );
    for ( auto i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( x[ i ] / y[ i ], z[ i ] );
    }
}

TEST( NCPAmathTest, MultiplyVectorsWorksCorrectlyForVectors ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 8, -1.0, 1.0 );
    vector<double> z = multiply_vectors( x, y );
    for ( auto i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( x[ i ] * y[ i ], z[ i ] );
    }
    for ( auto i = 5; i < 8; i++ ) {
        EXPECT_DOUBLE_EQ( z[ i ], 0.0 );
    }
}

TEST( NCPAmathTest, MultiplyVectorsWorksCorrectlyForArrays ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 5, -1.0, 1.0 );
    double *z = zeros<double>( 5 );
    multiply_vectors( 5, &x[0], &y[0], z );
    for ( auto i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( x[ i ] * y[ i ], z[ i ] );
    }
}


TEST( NCPAmathTest, MultiplyVectorsCommutesForVectors ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z1 = multiply_vectors( x, y );
    vector<double> z2 = multiply_vectors( y, x );
    EXPECT_ARRAY_DOUBLE_EQ( 8, z1, z2 );
}

TEST( NCPAmathTest, MultiplyVectorsCommutesForArrays ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 5, 0.0, 1.0 );
    double *z1 = zeros<double>( 5 );
    double *z2 = zeros<double>( 5 );
    multiply_vectors<double>( 5, &x[0], &y[0], z1 );
    multiply_vectors<double>( 5, &y[0], &x[0], z2 );
    EXPECT_ARRAY_DOUBLE_EQ( 5, z1, z2 );
}

TEST( NCPAmathTest, ScaleVectorWorksProperlyForVectors) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 );
    double scalar = random_number<double>(-3.0, 3.0);
    vector<double> y = scale_vector( x, scalar );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( y[i], x[i] * scalar );
    }
}

TEST( NCPAmathTest, ScaleVectorWorksProperlyForArrays) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 );
    double scalar = random_number<double>(-3.0, 3.0);
    double *y = zeros<double>( 5 );
    scale_vector(5, &x[0], scalar, y );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( y[i], x[i] * scalar );
    }
}

TEST( NCPAmathTest, ScaleVectorWorksProperlyInPlace) {
    double *x = zeros<double>( 5 );
    vector<double> xv = random_numbers<double>( 5, 0.0, 1.0 );
    std::copy( xv.cbegin(), xv.cend(), x );
    double scalar = random_number<double>(-3.0, 3.0);
    scale_vector(5, x, scalar );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[i], xv[i] * scalar );
    }
}

TEST( NCPAmathTest, SubtractVectorsWorksCorrectlyForVectors ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z = subtract_vectors( x, y );
    for ( auto i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( x[ i ] - y[ i ], z[ i ] );
    }
    for ( auto i = 5; i < 8; i++ ) {
        EXPECT_DOUBLE_EQ( -y[ i ], z[ i ] );
    }
}

TEST( NCPAmathTest, SubtractVectorsWorksCorrectlyForArrays ) {
    vector<double> x = random_numbers<double>( 5, 0.0, 1.0 ),
                   y = random_numbers<double>( 5, 0.0, 1.0 );
    double *z = zeros<double>( 5 );
    subtract_vectors( 5, &x[0], &y[0], z );
    for ( auto i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( x[ i ] - y[ i ], z[ i ] );
    }
}