#include "NCPA/math.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "NCPA/gtest.hpp"

#include <numbers>
#include <complex>

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
	EXPECT_ARRAY_EQ(5,testArray,shouldBe);
	// for (int i = 0; i < 5; ++i) {
  	// 	EXPECT_EQ(testArray[i], shouldBe[i] );
	// }
}

TEST(NCPAmathTest,CircShiftShiftsBackward) {
	int *testArray = index_vector<int>( 5 );
	int shouldBe[5] = {3,4,0,1,2};
	circshift( testArray, 5, -2, testArray );
	EXPECT_ARRAY_EQ(5,testArray,shouldBe);
	// for (int i = 0; i < 5; ++i) {
  	// 	EXPECT_EQ(testArray[i], shouldBe[i] );
	// }
}

TEST(NCPAmathTest,CircShiftZeroDoesNotShift) {
	int *testArray = index_vector<int>( 5 );
	int shouldBe[5] = {0,1,2,3,4};
	circshift( testArray, 5, 0, testArray );
	EXPECT_ARRAY_EQ(5,testArray,shouldBe);
	// for (int i = 0; i < 5; ++i) {
  	// 	EXPECT_EQ(testArray[i], shouldBe[i] );
	// }
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
	free_array3d( testArray, 10, 5, 3 );
}

TEST(NCPAmathTest,FindIntervalInclusiveFindsInterval) {
	double *testArray = index_vector<double>( 5 );
	size_t below, above;
	for (size_t i = 1; i < 5; i++) {
		double target = 0.5*((double)i + (double)(i-1));
		EXPECT_TRUE( find_interval_inclusive<double>( testArray, 5, target, below, above ) );
		EXPECT_EQ( below, i-1 );
		EXPECT_EQ( above, i );
	}
}

TEST(NCPAmathTest,FindIntervalInclusiveRecognizesOutOfRangeBelow) {
	double *testArray = index_vector<double>( 5 );
	size_t below, above;
	double target = testArray[0] - 1.0;
	EXPECT_FALSE( find_interval_inclusive<double>( testArray, 5, target, below, above ) );
	EXPECT_EQ( below, 0 );
	EXPECT_EQ( above, 0 );
}

TEST(NCPAmathTest,FindIntervalInclusiveRecognizesOutOfRangeAbove) {
	double *testArray = index_vector<double>( 5 );
	size_t below, above;
	double target = testArray[4] + 1.0;
	EXPECT_FALSE( find_interval_inclusive<double>( testArray, 5, target, below, above ) );
	EXPECT_EQ( below, 5 );
	EXPECT_EQ( above, 5 );
}

TEST(NCPAmathTest,FindClosestIndexWorksCorrectlyForArrays) {
	double *testArray = index_vector<double>( 5 );
	for (size_t i = 0; i < 5; i++) {
		double target = testArray[i];
		EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
		target -= 0.25;
		EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
		target += 0.5;
		EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
		if (i == 0) {
			target = -20000.0;
			EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
		} else if (i == 4) {
			target = 20000.0;
			EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i );
		} else {
			target += 0.5;
			EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i+1 );
			target -= 1.5;
			EXPECT_EQ( find_closest_index<double>( testArray, 5, target ), i-1 );
		}
	}
}

TEST(NCPAmathTest,FindClosestIndexWorksCorrectlyForVectors) {
	double *tempArray = index_vector<double>( 5 );
	vector<double> testArray( tempArray, tempArray + 5 );
	for (size_t i = 0; i < 5; i++) {
		double target = testArray[i];
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
			EXPECT_EQ( find_closest_index<double>( testArray, target ), i+1 );
			target -= 1.5;
			EXPECT_EQ( find_closest_index<double>( testArray, target ), i-1 );
		}
	}
}

TEST(NCPAmathTest,Cart2PolConvertsCorrectly) {
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
	EXPECT_DOUBLE_EQ( std::fabs(theta), numbers::pi );

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

TEST(NCPAmathTest,Pol2CartConvertsCorrectly) {
	double x, y;

	pol2cart<double>( 1.0, 0.0, x, y );
	EXPECT_DOUBLE_EQ( x, 1.0 );
	EXPECT_NEAR( y, 0.0, 1e-12 );

	pol2cart<double>( 1.0, 0.25 * numbers::pi, x, y );
	EXPECT_DOUBLE_EQ( x, 1.0 / std::sqrt(2.0) );
	EXPECT_DOUBLE_EQ( y, 1.0 / std::sqrt(2.0) );

	pol2cart<double>( 1, 0.5 * numbers::pi, x, y );
	EXPECT_NEAR( x, 0.0, 1e-12 );
	EXPECT_DOUBLE_EQ( y, 1.0 );
	
	pol2cart<double>( 1, 0.75 * numbers::pi, x, y );
	EXPECT_DOUBLE_EQ( x, -1.0 / std::sqrt(2.0) );
	EXPECT_DOUBLE_EQ( y, 1.0 / std::sqrt(2.0) );

	pol2cart<double>( 1.0, numbers::pi, x, y );
	EXPECT_DOUBLE_EQ( x, -1.0 );
	EXPECT_NEAR( y, 0.0, 1e-12 );

	pol2cart<double>( 1, -0.75 * numbers::pi, x, y );
	EXPECT_DOUBLE_EQ( x, -1.0 / std::sqrt(2.0) );
	EXPECT_DOUBLE_EQ( y, -1.0 / std::sqrt(2.0) );

	pol2cart<double>( 1, -0.5 * numbers::pi, x, y );
	EXPECT_NEAR( x, 0.0, 1e-12 );
	EXPECT_DOUBLE_EQ( y, -1.0 );

	pol2cart<double>( 1, -0.25 * numbers::pi, x, y );
	EXPECT_DOUBLE_EQ( x, 1.0 / std::sqrt(2.0) );
	EXPECT_DOUBLE_EQ( y, -1.0 / std::sqrt(2.0) );

	pol2cart<double>( 0.0, 0.0, x, y );
	EXPECT_DOUBLE_EQ( x, 0.0 );
	EXPECT_DOUBLE_EQ( y, 0.0 );
}

TEST(NCPAmathTest,Deg2RadConvertsCorrectly) {
	for (int i = -8; i <= 8; i++) {
		double d = (double)i * 45.0;
		double r = (double)i * numbers::pi / 4.0;
		EXPECT_DOUBLE_EQ( deg2rad<double>( d ), r );
	}
}

TEST(NCPAmathTest,Rad2DegConvertsCorrectly) {
	for (int i = -8; i <= 8; i++) {
		double d = (double)i * 45.0;
		double r = (double)i * numbers::pi / 4.0;
		EXPECT_DOUBLE_EQ( rad2deg<double>( r ), d );
	}
}

TEST(NCPAmathTest,Double2ComplexConvertsRealCorrectly) {
	double *d = zeros<double>( 9 );
	complex<double> *c = zeros<complex<double>>( 9 );

	for (int i = 0; i < 9; i++) {
		d[ i ] = (double)i * numbers::pi / 4.0 - numbers::pi;
	}
	double2complex( 9, d, c );
	for (int i = 0; i < 9; i++) {
		EXPECT_EQ( d[i], c[i].real() );
		EXPECT_EQ( c[i].imag(), 0.0 );
	}
}

TEST(NCPAmathTest,Double2ComplexConvertsRealAndImagCorrectly) {
	double *r = zeros<double>( 9 );
	double *i = zeros<double>( 9 );
	complex<double> *c = zeros<complex<double>>( 9 );

	for (int ii = 0; ii < 9; ii++) {
		r[ ii ] = (double)ii * numbers::pi / 4.0 - numbers::pi;
		i[ ii ] = (double)ii * numbers::pi;
	}
	double2complex( 9, r, i, c );
	for (int ii = 0; ii < 9; ii++) {
		EXPECT_EQ( c[ii].real(), r[ii] );
		EXPECT_EQ( c[ii].imag(), i[ii] );
	}
}

TEST(NCPAmathTest,Complex2DoubleConvertsCorrectly) {
	double *r = zeros<double>( 9 );
	double *i = zeros<double>( 9 );
	complex<double> *c = zeros<complex<double>>( 9 );

	for (int ii = 0; ii < 9; ii++) {
		c[ ii ].real( (double)ii * numbers::pi / 4.0 - numbers::pi );
		c[ ii ].imag( (double)ii * numbers::pi );
	}
	complex2double( 9, c, r, i );
	for (int ii = 0; ii < 9; ii++) {
		EXPECT_EQ( c[ii].real(), r[ii] );
		EXPECT_EQ( c[ii].imag(), i[ii] );
	}
}

TEST(NCPAmathTest,MaxWorksProperlyOnArrays) {
	int *testArray = index_vector<int>( 10 );
	EXPECT_EQ( NCPA::math::max<int>( testArray, 10 ), 9 );
	testArray[5] *= 3;
	EXPECT_EQ( NCPA::math::max<int>( testArray, 10 ), 15 );
}

TEST(NCPAmathTest,MaxWorksProperlyOnVectors) {
	int *testArray = index_vector<int>( 10 );
	vector<int> testVector( testArray, testArray+10 );
	EXPECT_EQ( NCPA::math::max<int>( testVector ), 9 );
	testVector[5] *= 3;
	EXPECT_EQ( NCPA::math::max<int>( testVector ), 15 );
}

TEST(NCPAmathTest,MinWorksProperlyOnArrays) {
	int *testArray = index_vector<int>( 10 );
	EXPECT_EQ( NCPA::math::min<int>( testArray, 10 ), 0 );
	testArray[5] *= -2;
	EXPECT_EQ( NCPA::math::min<int>( testArray, 10 ), -10 );
}

TEST(NCPAmathTest,MinWorksProperlyOnVectors) {
	int *testArray = index_vector<int>( 10 );
	vector<int> testVector( testArray, testArray+10 );
	EXPECT_EQ( NCPA::math::min<int>( testVector ), 0 );
	testVector[5] *= -2;
	EXPECT_EQ( NCPA::math::min<int>( testVector ), -10 );
}