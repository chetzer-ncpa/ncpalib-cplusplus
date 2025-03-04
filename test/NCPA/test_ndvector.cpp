#include "NCPA/gtest.hpp"
#include "NCPA/ndvector.hpp"
// #include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <complex>
#include <iostream>
#include <numbers>
#include <vector>

using namespace std;
using namespace NCPA::arrays;
using namespace testing;

#define _TEST_TITLE_ NDVectorTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {
            // cout << "Setting dims" << endl;
            dims = { 5, 4, 3 };
            // cout << "Setting d1" << endl;
            d1   = { 4.2, 3.1, 2.0 };
            // cout << "Setting d2" << endl;
            d2   = { d1, d1 };
            // cout << "Setting d3" << endl;
            d3   = { d2, d2, d2 };
            // cout << "Done" << endl;
        }

        ndvector<1, double> d1;
        ndvector<2, double> d2;
        ndvector<3, double> d3;
        _dimarray<3> dims;
};

TEST_F( _TEST_TITLE_, _dimarrayAssigns ) {
    EXPECT_TRUE( dims == ( std::vector<size_t> { 5, 4, 3 } ) );
}

TEST_F( _TEST_TITLE_, _dimarraySubdimsIsCorrect ) {
    auto sub = dims.subdims();
    EXPECT_EQ( sub.size(), 2 );
    for ( size_t i = 0; i < 2; i++ ) {
        EXPECT_EQ( sub[ i ], dims[ i + 1 ] );
    }
    auto subsub = dims.subdims().subdims();
    EXPECT_EQ( subsub.size(), 1 );
    EXPECT_EQ( subsub[ 0 ], dims[ 2 ] );
}

TEST_F( _TEST_TITLE_, ZeroLengthVectorIsDouble ) {
    double d = 4.2;
    ndvector<0, double> d0( d );
    EXPECT_DOUBLE_EQ( d, d0 );
}

TEST_F( _TEST_TITLE_, Length1ArrayIsStdVector ) {
    std::vector<double> v { 4.2, 3.1, 2.0 };
    EXPECT_TRUE( v == d1 );
    EXPECT_TRUE( d1 == v );
}

TEST_F( _TEST_TITLE_, AssignmentOperatorWorksOnSize1 ) {
    std::vector<double> v  = d1;
    v[ 0 ]                *= 2.0;
    EXPECT_FALSE( v == d1 );
    EXPECT_FALSE( d1 == v );
    d1 = v;
    EXPECT_TRUE( v == d1 );
    EXPECT_TRUE( d1 == v );
}

TEST_F( _TEST_TITLE_, ShapeFunctionWorksOnSize1 ) {
    EXPECT_TRUE( d1.shape() == ( std::vector<size_t> { 3 } ) );
}

TEST_F( _TEST_TITLE_, SquareBracketsWorkOnSize1 ) {
    EXPECT_TRUE( d1[ 0 ] == 4.2 );
    EXPECT_TRUE( d1[ 1 ] == 3.1 );
    EXPECT_TRUE(( d1[ _dimarray<1>( { 2 } ) ] == 2.0 ));
    d1[ 0 ] = 10.0;
    EXPECT_TRUE( d1[ 0 ] == 10.0 );
    d1[ 4 ] = 10.0;
    EXPECT_TRUE( d1[ 4 ] == 10.0 );
    EXPECT_EQ( d1.size(), 4 );
}

TEST_F( _TEST_TITLE_, ResizeWorksOnSize1 ) {
    d1.reshape( { 4 } );
    EXPECT_TRUE( d1[ 0 ] == 4.2 );
    EXPECT_TRUE( d1[ 1 ] == 3.1 );
    EXPECT_TRUE( d1[ 2 ] == 2.0 );
    EXPECT_TRUE( d1[ 3 ] == 0.0 );
    d1.reshape( { 2 } );
    EXPECT_EQ( d1.size(), 2 );
    EXPECT_TRUE( d1[ 0 ] == 4.2 );
    EXPECT_TRUE( d1[ 1 ] == 3.1 );
}

TEST_F( _TEST_TITLE_, ShapeFunctionWorksOnSize2 ) {
    EXPECT_TRUE( d2.shape() == ( std::vector<size_t> { 2, 3 } ) );
}

TEST_F( _TEST_TITLE_, Length2ArrayIsVectorOfVectors ) {
    std::vector<double> v { 4.2, 3.1, 2.0 };
    EXPECT_TRUE( v == d2[ 0 ] );
    EXPECT_TRUE( d2[ 1 ] == v );
}

TEST_F( _TEST_TITLE_, SquareBracketsWorkOnSize2 ) {
    EXPECT_TRUE( d2[ 0 ][ 0 ] == 4.2 );
    EXPECT_TRUE( d2[ 1 ][ 1 ] == 3.1 );
    EXPECT_TRUE( d2[ 0 ][ 2 ] == 2.0 );
    d2[ 0 ][ 0 ] = 10.0;
    EXPECT_TRUE( d2[ 0 ][ 0 ] == 10.0 );
}

TEST_F( _TEST_TITLE_, SquareBracketsWithArrayWorkOnSize2 ) {
    EXPECT_TRUE(( d2[ { 0, 0 } ] == 4.2 ));
    EXPECT_TRUE(( d2[ { 1, 1 } ] == 3.1 ));
    EXPECT_TRUE(( d2[ { 0, 2 } ] == 2.0 ));
    d2[ { 0, 0 } ] = 10.0;
    EXPECT_TRUE(( d2[ { 0, 0 } ] == 10.0 ));
}

TEST_F( _TEST_TITLE_, ReshapeWorksOnSize2 ) {
    d2.reshape( { 3, 3 } );
    EXPECT_TRUE( d2.shape() == ( std::vector<size_t> { 3, 3 } ) );
    EXPECT_TRUE(( d2[ { 0, 0 } ] == 4.2 ));
    EXPECT_TRUE(( d2[ { 1, 1 } ] == 3.1 ));
    EXPECT_TRUE(( d2[ { 2, 2 } ] == 0.0 ));
}

TEST_F( _TEST_TITLE_, ShapeFunctionWorksOnSize3 ) {
    EXPECT_TRUE( d3.shape() == ( std::vector<size_t> { 3, 2, 3 } ) );
}

TEST_F( _TEST_TITLE_, Length3ArrayIsVectorOfVectorOfVectorss ) {
    std::vector<double> v { 4.2, 3.1, 2.0 };
    EXPECT_TRUE( v == d3[ 0 ][ 0 ] );
    EXPECT_TRUE( d3[ 1 ][ 1 ] == v );
}

TEST_F( _TEST_TITLE_, SquareBracketsWorkOnSize3 ) {
    EXPECT_TRUE( d3[ 0 ][ 0 ][ 0 ] == 4.2 );
    EXPECT_TRUE( d3[ 1 ][ 1 ][ 1 ] == 3.1 );
    EXPECT_TRUE( d3[ 2 ][ 0 ][ 2 ] == 2.0 );
    d3[ 0 ][ 0 ][ 0 ] = 10.0;
    EXPECT_TRUE( d3[ 0 ][ 0 ][ 0 ] == 10.0 );
}

TEST_F( _TEST_TITLE_, SquareBracketsWithArrayWorksOnSize3 ) {
    EXPECT_TRUE( (d3[ {0 , 0 , 0} ] == 4.2) );
    EXPECT_TRUE( (d3[ {1 , 1 , 1} ] == 3.1) );
    EXPECT_TRUE( (d3[ {2 , 0 , 2} ] == 2.0) );
    d3[ {0 , 0 , 0} ] = 10.0;
    EXPECT_TRUE( (d3[ {0 , 0 , 0} ] == 10.0) );
}

TEST_F( _TEST_TITLE_, ReshapeWorksOnSize3 ) {
    d3.reshape( { 3, 3, 3 } );
    EXPECT_TRUE( d3.shape() == ( std::vector<size_t> { 3, 3, 3 } ) );
    EXPECT_TRUE(( d3[ { 0, 0, 0 } ] == 4.2 ));
    EXPECT_TRUE(( d3[ { 1, 1, 1 } ] == 3.1 ));
    EXPECT_TRUE(( d3[ { 2, 2, 2 } ] == 0.0 ));
    d3.reshape( { 1, 3, 3 } );
    EXPECT_TRUE(( d3[ { 0, 0, 0 } ] == 4.2 ));
}