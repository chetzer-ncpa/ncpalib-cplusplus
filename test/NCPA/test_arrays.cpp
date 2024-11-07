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

TEST( NCPAarraysTest, ZerosCreatesArray ) {
    double *testArray     = zeros<double>( 5 );
    double  shouldBe[ 5 ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    ASSERT_ARRAY_EQ( 5, testArray, shouldBe );
    // for (auto i = 0; i < 10; i++) {
    // 	EXPECT_EQ( testArray[i], 0.0 );
    // }
    delete[] testArray;
}

TEST( NCPAarraysTest, ZerosCreates2DArray ) {
    double **testArray = zeros<double>( 10, 5 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            EXPECT_EQ( testArray[ i ][ j ], 0.0 );
        }
    }
}

TEST( NCPAarraysTest, FreeArray2DFreesArrayWithoutError ) {
    double **testArray = zeros<double>( 10, 5 );
    free_array( testArray, 10, 5 );
    EXPECT_EQ( testArray, nullptr );
}

TEST( NCPAarraysTest, ZerosCreates3DArray ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            for ( auto k = 0; k < 3; k++ ) {
                EXPECT_EQ( testArray[ i ][ j ][ k ], 0.0 );
            }
        }
    }
}

TEST( NCPAarraysTest, FreeArray3DFreesArrayWithoutError ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    free_array( testArray, 10, 5, 3 );
    ASSERT_EQ( testArray, nullptr );
}

TEST( NCPAarraysTest, IndexVectorCreatesExpectedVectors ) {
    int *testArray = index_vector<int>( 5 );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_EQ( i, testArray[ i ] );
    }
    double *testArray2 = index_vector<double>( 5 );
    for ( int i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( (double)i, testArray2[ i ] );
    }
}

TEST( NCPAarraysTest, IndexVectorCreatesExpectedVectorsWithOffset ) {
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

TEST( NCPAarraysTest, CircShiftShiftsForward ) {
    int *testArray     = index_vector<int>( 5 );
    int  shouldBe[ 5 ] = { 2, 3, 4, 0, 1 };
    circshift( testArray, 5, 2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAarraysTest, CircShiftShiftsBackward ) {
    int *testArray     = index_vector<int>( 5 );
    int  shouldBe[ 5 ] = { 3, 4, 0, 1, 2 };
    circshift( testArray, 5, -2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAarraysTest, CircShiftZeroDoesNotShift ) {
    int *testArray     = index_vector<int>( 5 );
    int  shouldBe[ 5 ] = { 0, 1, 2, 3, 4 };
    circshift( testArray, 5, 0, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAarraysTest, FillArray2DSetsValuesCorrectly ) {
    double **testArray = zeros<double>( 10, 5 );
    fill( testArray, 10, 5, 4.5 );
    for ( auto i = 0; i < 10; i++ ) {
        for ( auto j = 0; j < 5; j++ ) {
            EXPECT_DOUBLE_EQ( testArray[ i ][ j ], 4.5 );
        }
    }
    free_array( testArray, 10, 5 );
}

TEST( NCPAarraysTest, FillArray3DSetsValuesCorrectly ) {
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

TEST(NCPAarraysTest, ReverseReversesOrder) {
	double *testArray = index_vector<double>( 11 );
	vector<double> testVector( testArray, testArray + 11 );
	std::reverse( testVector.begin(), testVector.end() );
	NCPA::arrays::reverse( testArray, 11, testArray );
	EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST(NCPAarraysTest, ReverseReversesOrderPartially) {
	double *testArray = index_vector<double>( 11 );
	vector<double> testVector( testArray, testArray + 11 );
	std::reverse( testVector.begin(), testVector.begin() + 5 );
	NCPA::arrays::reverse( testArray, 5, testArray );
	EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}
