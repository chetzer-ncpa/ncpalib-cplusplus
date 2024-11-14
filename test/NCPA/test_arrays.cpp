#include "NCPA/gtest.hpp"
#include "NCPA/arrays.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <complex>
#include <numbers>
#include <algorithm>

using namespace std;
using namespace NCPA::arrays;
using namespace testing;

TEST( NCPAArraysLibraryTest, ZerosCreatesArray ) {
    double *testArray     = zeros<double>( 5 );
    double  shouldBe[ 5 ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    ASSERT_ARRAY_EQ( 5, testArray, shouldBe );
    // for (auto i = 0; i < 10; i++) {
    // 	EXPECT_EQ( testArray[i], 0.0 );
    // }
    delete[] testArray;
}

TEST( NCPAArraysLibraryTest, ZerosCreates2DArray ) {
    double **testArray = zeros<double>( 10, 5 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            EXPECT_EQ( testArray[ i ][ j ], 0.0 );
        }
    }
}

TEST( NCPAArraysLibraryTest, FreeArray2DFreesArrayWithoutError ) {
    double **testArray = zeros<double>( 10, 5 );
    free_array( testArray, 10, 5 );
    EXPECT_EQ( testArray, nullptr );
}

TEST( NCPAArraysLibraryTest, ZerosCreates3DArray ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            for ( auto k = 0; k < 3; k++ ) {
                EXPECT_EQ( testArray[ i ][ j ][ k ], 0.0 );
            }
        }
    }
}

TEST( NCPAArraysLibraryTest, FreeArray3DFreesArrayWithoutError ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    free_array( testArray, 10, 5, 3 );
    ASSERT_EQ( testArray, nullptr );
}

TEST( NCPAArraysLibraryTest, IndexVectorCreatesExpectedVectors ) {
    vector<int> testArray = index_vector<int>( 5 );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_EQ( i, testArray[ i ] );
    }
    vector<double> testArray2 = index_vector<double>( 5 );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( (double)i, testArray2[ i ] );
    }
}

TEST( NCPAArraysLibraryTest, IndexVectorCreatesExpectedVectorsWithOffset ) {
    int  offset    = 2;
    vector<int> testArray = index_vector<int>( 5, offset );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_EQ( i + offset, testArray[ i ] );
    }
    vector<double> testArray2 = index_vector<double>( 5, (double)offset );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( (double)( i + offset ), testArray2[ i ] );
    }
}

TEST( NCPAArraysLibraryTest, CircShiftShiftsForward ) {
    vector<int> inds = index_vector<int>( 5 );
    int *testArray     = &inds[0];
    int  shouldBe[ 5 ] = { 2, 3, 4, 0, 1 };
    circshift( testArray, 5, 2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAArraysLibraryTest, CircShiftShiftsBackward ) {
    vector<int> inds = index_vector<int>( 5 );
    int *testArray     = &inds[0];
    int  shouldBe[ 5 ] = { 3, 4, 0, 1, 2 };
    circshift( testArray, 5, -2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAArraysLibraryTest, CircShiftZeroDoesNotShift ) {
    vector<int> inds = index_vector<int>( 5 );
    int *testArray     = &inds[0];
    int  shouldBe[ 5 ] = { 0, 1, 2, 3, 4 };
    circshift( testArray, 5, 0, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAArraysLibraryTest, FillArray2DSetsValuesCorrectly ) {
    double **testArray = zeros<double>( 10, 5 );
    fill( testArray, 10, 5, 4.5 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            EXPECT_DOUBLE_EQ( testArray[ i ][ j ], 4.5 );
        }
    }
    free_array( testArray, 10, 5 );
}

TEST( NCPAArraysLibraryTest, FillArray3DSetsValuesCorrectly ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    fill( testArray, 10, 5, 3, 4.5 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            for ( auto k = 0; k < 3; k++ ) {
                EXPECT_DOUBLE_EQ( testArray[ i ][ j ][ k ], 4.5 );
            }
        }
    }
    free_array( testArray, 10, 5, 3 );
}

TEST(NCPAArraysLibraryTest, ReverseReversesOrder) {
	vector<double> testVector = index_vector<double>( 11 );
    double *testArray = zeros<double>( 11 );
    std::copy( testVector.cbegin(), testVector.cend(), testArray );
	std::reverse( testVector.begin(), testVector.end() );
	NCPA::arrays::reverse( testArray, 11, testArray );
	EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST(NCPAArraysLibraryTest, ReverseReversesOrderPartially) {
	vector<double> testVector = index_vector<double>( 11 );
    double *testArray = zeros<double>( 11 );
    std::copy( testVector.cbegin(), testVector.cend(), testArray );
	std::reverse( testVector.begin(), testVector.begin() + 5 );
	NCPA::arrays::reverse( testArray, 5, testArray );
	EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST(NCPAArraysLibraryTest, AsArrayLinksToVector) {
    std::vector<int> vi = { 0, 1, 2, 3, 4 };
    int *pi = as_array( vi );
    for (auto i = 0; i < 5; i++) {
        ASSERT_EQ( pi[ i ], vi[ i ] );
    }
    vi[ 0 ] = 5;
    for (auto i = 0; i < 5; i++) {
        ASSERT_EQ( pi[ i ], vi[ i ] );
    }
}