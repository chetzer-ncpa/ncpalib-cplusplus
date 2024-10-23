#include "NCPA/math.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>

using namespace std;
using namespace NCPA::math;
using namespace testing;


TEST(NCPAmathTest,ZerosCreatesArray) {
	double *testArray = zeros<double>( 10 );
	for (auto i = 0; i < 10; i++) {
		ASSERT_EQ( testArray[i], 0.0 );
	}
	delete [] testArray;
}

TEST(NCPAmathTest,Zeros2DCreates2DArray) {
	double **testArray = zeros2d<double>( 10, 5 );
	for (auto i = 0; i < 10; i++) {
		for (auto j = 0; j < 5; j++) {
			ASSERT_EQ( testArray[i][j], 0.0 );
		}
	}
	free_array2d<double>( testArray, 10, 5 );
}
