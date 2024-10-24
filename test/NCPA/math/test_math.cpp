#include "NCPA/math.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "NCPA/gtest.hpp"

using namespace std;
using namespace NCPA::math;
using namespace testing;


TEST(NCPAmathTest,ZerosCreatesArray) {
	double *testArray = zeros<double>( 5 );
	double shouldBe[ 5 ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	ASSERT_ARRAY_EQ( 5, testArray, shouldBe );
	// for (auto i = 0; i < 10; i++) {
	// 	EXPECT_EQ( testArray[i], 0.0 );
	// }
	delete [] testArray;
}

TEST(NCPAmathTest,Zeros2DCreates2DArray) {
	double **testArray = zeros2d<double>( 10, 5 );
	for (auto i = 0; i < 10; i++) {
		for (auto j = 0; j < 5; j++) {
			EXPECT_EQ( testArray[i][j], 0.0 );
		}
	}
}

TEST(NCPAmathTest,FreeArray2DFreesArrayWithoutError) {
	double **testArray = zeros2d<double>( 10, 5 );
	free_array2d<double>( testArray, 10, 5 );
	EXPECT_EQ( testArray, nullptr );
}

TEST(NCPAmathTest,Zeros3DCreates3DArray) {
	double ***testArray = zeros3d<double>( 10, 5, 3 );
	for (auto i = 0; i < 10; i++) {
		for (auto j = 0; j < 5; j++) {
			for (auto k = 0; k < 3; k++) {
				EXPECT_EQ( testArray[i][j][k], 0.0 );
			}	
		}
	}
}

TEST(NCPAmathTest,FreeArray3DFreesArrayWithoutError) {
	double ***testArray = zeros3d<double>( 10, 5, 3 );
	free_array3d<double>( testArray, 10, 5, 3 );
	ASSERT_EQ( testArray, nullptr );
}

TEST(NCPAmathTest,IndexVectorCreatesExpectedVectors) {
	int *testArray = index_vector<int>( 5 );
	for (int i = 0; i < 5; i++) {
		EXPECT_EQ( i, testArray[ i ] );
	}
	double *testArray2 = index_vector<double>( 5 );
	for (int i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( (double)i, testArray2[ i ] );
	}
}

TEST(NCPAmathTest,IndexVectorCreatesExpectedVectorsWithOffset) {
	int offset = 2;
	int *testArray = index_vector<int>( 5, offset );
	for (int i = 0; i < 5; i++) {
		EXPECT_EQ( i+offset, testArray[ i ] );
	}
	double *testArray2 = index_vector<double>( 5, (double)offset );
	for (int i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( (double)(i + offset), testArray2[ i ] );
	}
}

TEST(NCPAmathTest,SignIsCorrect) {
	EXPECT_EQ(sign<int>(3),1);
	EXPECT_EQ(sign<int>(-3),-1);
	EXPECT_EQ(sign<int>(0),0);
	EXPECT_EQ(sign<double>(3.0),1);
	EXPECT_EQ(sign<double>(-3.2),-1);
}

TEST(NCPAmathTest,CircShiftShiftsForward) {
	int *testArray = index_vector<int>( 5 );
	int shouldBe[5] = {2,3,4,0,1};
	circshift( testArray, 5, 2, testArray );
	for (int i = 0; i < 5; ++i) {
  		EXPECT_EQ(testArray[i], shouldBe[i] );
	}
}

TEST(NCPAmathTest,CircShiftShiftsBackward) {
	int *testArray = index_vector<int>( 5 );
	int shouldBe[5] = {3,4,0,1,2};
	circshift( testArray, 5, -2, testArray );
	for (int i = 0; i < 5; ++i) {
  		EXPECT_EQ(testArray[i], shouldBe[i] );
	}
}

TEST(NCPAmathTest,CircShiftZeroDoesNotShift) {
	int *testArray = index_vector<int>( 5 );
	int shouldBe[5] = {0,1,2,3,4};
	circshift( testArray, 5, 0, testArray );
	for (int i = 0; i < 5; ++i) {
  		EXPECT_EQ(testArray[i], shouldBe[i] );
	}
}

TEST(NCPAmathTest,FillArray2DSetsValuesCorrectly) {
	double **testArray = zeros2d<double>( 10, 5 );
	fill_array2d<double>( testArray, 10, 5, 4.5 );
	for (auto i = 0; i < 10; i++) {
		for (auto j = 0; j < 5; j++) {
			EXPECT_DOUBLE_EQ( testArray[i][j], 4.5 );
		}
	}
	free_array2d( testArray, 10, 5 );
}

TEST(NCPAmathTest,FillArray3DSetsValuesCorrectly) {
	double ***testArray = zeros3d<double>( 10, 5, 3 );
	fill_array3d<double>( testArray, 10, 5, 3, 4.5 );
	for (auto i = 0; i < 10; i++) {
		for (auto j = 0; j < 5; j++) {
			for (auto k = 0; k < 3; k++) {
				EXPECT_DOUBLE_EQ( testArray[i][j][k], 4.5 );
			}
		}
	}
	free_array2d( testArray, 10, 5 );
}