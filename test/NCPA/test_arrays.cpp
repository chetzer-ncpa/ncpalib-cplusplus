#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <numbers>

using namespace std;
using namespace NCPA::arrays;
using namespace testing;

#define _TEST_TITLE_ NCPAArraysLibraryTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {
            // cout << "Setting dims" << endl;
            dims = { 5, 4, 3 };
            // cout << "Setting d1" << endl;
            d1   = vector<double> { 4.2, 3.1, 2.0 };
            // cout << "Setting d2" << endl;
            d2   = { d1, d1 };
            // cout << "Setting d3" << endl;
            d3   = { d2, d2, d2 };
            // cout << "Done" << endl;
            sv.resize( dims );
            sa.resize( dims );

            for (auto i = 0; i < dims[ 0 ]; ++i) {
                for (auto j = 0; j < dims[ 1 ]; ++j) {
                    for (auto k = 0; k < dims[ 2 ]; ++k) {
                        sa( i, j, k ) = i * 100 + j * 10 + k;
                        sv( i, j, k ) = i * 100 + j * 10 + k;
                    }
                }
            }
        }

        ndvector<1, double> d1;
        ndvector<2, double> d2;
        ndvector<3, double> d3;
        details::_dimarray dims;

        SmartVector<int> sv;
        SmartArray<3, int> sa;
};

TEST_F( _TEST_TITLE_, ZerosCreatesArray ) {
    double *testArray    = zeros<double>( 5 );
    double shouldBe[ 5 ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    ASSERT_ARRAY_EQ( 5, testArray, shouldBe );
    // for (auto i = 0; i < 10; i++) {
    // 	EXPECT_EQ( testArray[i], 0.0 );
    // }
    delete[] testArray;
}

TEST_F( _TEST_TITLE_, ZerosCreates2DArray ) {
    double **testArray = zeros<double>( 10, 5 );
    for (auto i = 0; i < 10; i++) {
        for (auto j = 0; j < 5; j++) {
            EXPECT_EQ( testArray[ i ][ j ], 0.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, FreeArray2DFreesArrayWithoutError ) {
    double **testArray = zeros<double>( 10, 5 );
    free_array( testArray, 10, 5 );
    EXPECT_EQ( testArray, nullptr );
}

TEST_F( _TEST_TITLE_, ZerosCreates3DArray ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    for (auto i = 0; i < 10; i++) {
        for (auto j = 0; j < 5; j++) {
            for (auto k = 0; k < 3; k++) {
                EXPECT_EQ( testArray[ i ][ j ][ k ], 0.0 );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, FreeArray3DFreesArrayWithoutError ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    free_array( testArray, 10, 5, 3 );
    ASSERT_EQ( testArray, nullptr );
}

TEST_F( _TEST_TITLE_, IndexVectorCreatesExpectedVectors ) {
    vector<int> testArray = index_vector<int>( 5 );
    for (int i = 0; i < 5; i++) {
        EXPECT_EQ( i, testArray[ i ] );
    }
    vector<double> testArray2 = index_vector<double>( 5 );
    for (int i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( (double)i, testArray2[ i ] );
    }
}

TEST_F( _TEST_TITLE_, IndexVectorCreatesExpectedVectorsWithOffset ) {
    int offset            = 2;
    vector<int> testArray = index_vector<int>( 5, offset );
    for (int i = 0; i < 5; i++) {
        EXPECT_EQ( i + offset, testArray[ i ] );
    }
    vector<double> testArray2 = index_vector<double>( 5, (double)offset );
    for (int i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( (double)( i + offset ), testArray2[ i ] );
    }
}

TEST_F( _TEST_TITLE_, CircShiftShiftsForward ) {
    vector<int> inds  = index_vector<int>( 5 );
    int *testArray    = &inds[ 0 ];
    int shouldBe[ 5 ] = { 2, 3, 4, 0, 1 };
    circshift( testArray, 5, 2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST_F( _TEST_TITLE_, CircShiftShiftsBackward ) {
    vector<int> inds  = index_vector<int>( 5 );
    int *testArray    = &inds[ 0 ];
    int shouldBe[ 5 ] = { 3, 4, 0, 1, 2 };
    circshift( testArray, 5, -2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST_F( _TEST_TITLE_, CircShiftZeroDoesNotShift ) {
    vector<int> inds  = index_vector<int>( 5 );
    int *testArray    = &inds[ 0 ];
    int shouldBe[ 5 ] = { 0, 1, 2, 3, 4 };
    circshift( testArray, 5, 0, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST_F( _TEST_TITLE_, FillArray2DSetsValuesCorrectly ) {
    double **testArray = zeros<double>( 10, 5 );
    fill( testArray, 10, 5, 4.5 );
    for (auto i = 0; i < 10; i++) {
        for (auto j = 0; j < 5; j++) {
            EXPECT_DOUBLE_EQ( testArray[ i ][ j ], 4.5 );
        }
    }
    free_array( testArray, 10, 5 );
}

TEST_F( _TEST_TITLE_, FillArray3DSetsValuesCorrectly ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    fill( testArray, 10, 5, 3, 4.5 );
    for (auto i = 0; i < 10; i++) {
        for (auto j = 0; j < 5; j++) {
            for (auto k = 0; k < 3; k++) {
                EXPECT_DOUBLE_EQ( testArray[ i ][ j ][ k ], 4.5 );
            }
        }
    }
    free_array( testArray, 10, 5, 3 );
}

TEST_F( _TEST_TITLE_, ReverseReversesOrder ) {
    vector<double> testVector = index_vector<double>( 11 );
    double *testArray         = zeros<double>( 11 );
    std::copy( testVector.cbegin(), testVector.cend(), testArray );
    std::reverse( testVector.begin(), testVector.end() );
    NCPA::arrays::reverse( testArray, 11, testArray );
    EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST_F( _TEST_TITLE_, ReverseReversesOrderPartially ) {
    vector<double> testVector = index_vector<double>( 11 );
    double *testArray         = zeros<double>( 11 );
    std::copy( testVector.cbegin(), testVector.cend(), testArray );
    std::reverse( testVector.begin(), testVector.begin() + 5 );
    NCPA::arrays::reverse( testArray, 5, testArray );
    EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST_F( _TEST_TITLE_, AsArrayLinksToVector ) {
    std::vector<int> vi = { 0, 1, 2, 3, 4 };
    int *pi             = as_array( vi );
    for (auto i = 0; i < 5; i++) {
        ASSERT_EQ( pi[ i ], vi[ i ] );
    }
    vi[ 0 ] = 5;
    for (auto i = 0; i < 5; i++) {
        ASSERT_EQ( pi[ i ], vi[ i ] );
    }
}

TEST_F( _TEST_TITLE_, SortPermutationCreatedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<size_t> perm_expected { 5, 6, 3, 0, 2, 1, 4 };
    auto perm = sort_permutation(
        unsorted, []( int const& a, int const& b ) { return a < b; } );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST_F( _TEST_TITLE_, SortPermutationIncreasingCreatedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<size_t> perm_expected { 5, 6, 3, 0, 2, 1, 4 };
    auto perm = sort_permutation_increasing( unsorted );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST_F( _TEST_TITLE_, SortPermutationDecreasingCreatedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<size_t> perm_expected { 4, 1, 2, 0, 3, 5, 6 };
    auto perm = sort_permutation_decreasing( unsorted );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST_F( _TEST_TITLE_, SortPermutationAppliedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<double> orig { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    std::vector<double> sorted_expected { 5.0, 6.0, 3.0, 0.0, 2.0, 1.0, 4.0 };
    auto perm = sort_permutation(
        unsorted, []( int const& a, int const& b ) { return a < b; } );
    auto sorted = apply_permutation( orig, perm );
    EXPECT_ARRAY_DOUBLE_EQ( perm.size(), sorted, sorted_expected );
}

TEST_F( _TEST_TITLE_, SortPermutationAppliedProperlyInPlace ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<double> orig { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    std::vector<double> sorted_expected { 5.0, 6.0, 3.0, 0.0, 2.0, 1.0, 4.0 };
    auto perm = sort_permutation(
        unsorted, []( int const& a, int const& b ) { return a < b; } );
    apply_permutation_in_place( orig, perm );
    EXPECT_ARRAY_DOUBLE_EQ( perm.size(), orig, sorted_expected );
}

TEST_F( _TEST_TITLE_, DoubleSortPermutationCreatedProperly ) {
    std::vector<int> unsorted1 { 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5 };
    std::vector<int> unsorted2 { 2, 1, 3, 4, 0, 0, 1, 2, 4, 3, 2, 2, 1 };
    std::vector<size_t> perm_expected { 1, 0, 2, 4,  5,  3, 6,
                                        7, 9, 8, 10, 11, 12 };
    auto perm = sort_permutation_increasing( unsorted1, unsorted2 );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST_F( _TEST_TITLE_, DoubleSortPermutationDecreasingCreatedProperly ) {
    std::vector<int> unsorted1 { 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5 };
    std::vector<int> unsorted2 { 2, 1, 3, 4, 0, 0, 1, 2, 4, 3, 2, 2, 1 };
    std::vector<size_t> perm_expected { 12, 10, 11, 8, 9, 7, 6,
                                        3,  4,  5,  2, 0, 1 };
    auto perm = sort_permutation_decreasing( unsorted1, unsorted2 );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST_F( _TEST_TITLE_, AddVectorsWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z = add_vectors( x, y );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] + y[ i ], z[ i ] );
    }
    for (auto i = 5; i < 8; i++) {
        EXPECT_DOUBLE_EQ( y[ i ], z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, AddVectorsCommutesForVectors ) {
    vector<double> x  = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y  = NCPA::math::random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z1 = add_vectors( x, y );
    vector<double> z2 = add_vectors( y, x );
    EXPECT_ARRAY_DOUBLE_EQ( 8, z1, z2 );
}

TEST_F( _TEST_TITLE_, AddArraysWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double *z        = NCPA::arrays::zeros<double>( 5 );
    add_arrays( 5, &x[ 0 ], &y[ 0 ], z );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] + y[ i ], z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, AddArraysCommutes ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double *z1       = NCPA::arrays::zeros<double>( 5 );
    double *z2       = NCPA::arrays::zeros<double>( 5 );
    add_arrays( 5, &x[ 0 ], &y[ 0 ], z1 );
    add_arrays( 5, &y[ 0 ], &x[ 0 ], z2 );
    EXPECT_ARRAY_DOUBLE_EQ( 5, z1, z2 );
}

TEST_F( _TEST_TITLE_, DivideVectorsWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 8, 0.1, 1.0 );
    vector<double> z = divide_vectors( x, y );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] / y[ i ], z[ i ] );
    }
    for (auto i = 5; i < 8; i++) {
        EXPECT_DOUBLE_EQ( z[ i ], 0.0 );
    }
}

TEST_F( _TEST_TITLE_, DivideArraysWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.1, 1.0 );
    double *z        = NCPA::arrays::zeros<double>( 5 );
    divide_arrays( 5, &x[ 0 ], &y[ 0 ], z );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] / y[ i ], z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, MultiplyVectorsWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 8, -1.0, 1.0 );
    vector<double> z = multiply_vectors( x, y );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] * y[ i ], z[ i ] );
    }
    for (auto i = 5; i < 8; i++) {
        EXPECT_DOUBLE_EQ( z[ i ], 0.0 );
    }
}

TEST_F( _TEST_TITLE_, MultiplyArraysWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, -1.0, 1.0 );
    double *z        = NCPA::arrays::zeros<double>( 5 );
    multiply_arrays( 5, &x[ 0 ], &y[ 0 ], z );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] * y[ i ], z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, MultiplyVectorsCommutes ) {
    vector<double> x  = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y  = NCPA::math::random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z1 = multiply_vectors( x, y );
    vector<double> z2 = multiply_vectors( y, x );
    EXPECT_ARRAY_DOUBLE_EQ( 8, z1, z2 );
}

TEST_F( _TEST_TITLE_, MultiplyArraysCommutes ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double *z1       = NCPA::arrays::zeros<double>( 5 );
    double *z2       = NCPA::arrays::zeros<double>( 5 );
    multiply_arrays<double>( 5, &x[ 0 ], &y[ 0 ], z1 );
    multiply_arrays<double>( 5, &y[ 0 ], &x[ 0 ], z2 );
    EXPECT_ARRAY_DOUBLE_EQ( 5, z1, z2 );
}

TEST_F( _TEST_TITLE_, ScaleVectorWorksProperly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double scalar    = NCPA::math::random_number<double>( -3.0, 3.0 );
    vector<double> y = scale_vector( x, scalar );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( y[ i ], x[ i ] * scalar );
    }
}

TEST_F( _TEST_TITLE_, ScaleArrayWorksProperly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double scalar    = NCPA::math::random_number<double>( -3.0, 3.0 );
    double *y        = NCPA::arrays::zeros<double>( 5 );
    scale_array( 5, &x[ 0 ], scalar, y );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( y[ i ], x[ i ] * scalar );
    }
}

TEST_F( _TEST_TITLE_, ScaleArrayWorksProperlyInPlace ) {
    double *x         = NCPA::arrays::zeros<double>( 5 );
    vector<double> xv = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    std::copy( xv.cbegin(), xv.cend(), x );
    double scalar = NCPA::math::random_number<double>( -3.0, 3.0 );
    scale_array( 5, x, scalar );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ], xv[ i ] * scalar );
    }
}

TEST_F( _TEST_TITLE_, SubtractVectorsWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z = subtract_vectors( x, y );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] - y[ i ], z[ i ] );
    }
    for (auto i = 5; i < 8; i++) {
        EXPECT_DOUBLE_EQ( -y[ i ], z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, SubtractArraysWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double *z        = NCPA::arrays::zeros<double>( 5 );
    subtract_arrays( 5, &x[ 0 ], &y[ 0 ], z );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] - y[ i ], z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, PlusOperatorWorksForVectors ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z = x + y;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] + y[ i ], z[ i ] );
    }
    for (auto i = 5; i < 8; i++) {
        EXPECT_DOUBLE_EQ( y[ i ], z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, PlusOperatorWorksForVectorScalar ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double y         = 3.5;
    vector<double> z = x + y;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] + y, z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, MinusOperatorWorksForVectors ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z = x - y;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] - y[ i ], z[ i ] );
    }
    for (auto i = 5; i < 8; i++) {
        EXPECT_DOUBLE_EQ( -y[ i ], z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, MinusOperatorWorksForVectorScalar ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double y         = 3.5;
    vector<double> z = x - y;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] - y, z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, NegationOperatorWorksForVector ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    vector<double> z = -x;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( -x[ i ], z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, MultiplyOperatorWorksForVectors ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z = x * y;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] * y[ i ], z[ i ] );
    }
    for (auto i = 5; i < 8; i++) {
        EXPECT_DOUBLE_EQ( 0.0, z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, MultiplyOperatorWorksForVectorScalar ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double y         = 3.5;
    vector<double> z = x / y;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] / y, z[ i ] );
    }
}

TEST_F( _TEST_TITLE_, _dimarrayAssigns ) {
    std::vector<size_t> testdims { 5, 4, 3 };
    for (auto i = 0; i < dims.size(); ++i) {
        EXPECT_EQ( dims.at( i ), testdims.at( i ) );
    }
    // EXPECT_TRUE( dims == ( std::vector<size_t> { 5, 4, 3 } ) );
}

TEST_F( _TEST_TITLE_, _dimarraySubdimsIsCorrect ) {
    auto sub = details::subdims( dims );
    EXPECT_EQ( sub.size(), 2 );
    for (size_t i = 0; i < 2; i++) {
        EXPECT_EQ( sub[ i ], dims[ i + 1 ] );
    }
    auto subsub = details::subdims( details::subdims( dims ) );
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
    EXPECT_TRUE(( d1[ 2 ] == 2.0 ));
    d1[ 0 ] = 10.0;
    EXPECT_TRUE( d1[ 0 ] == 10.0 );
    d1[ 4 ] = 10.0;
    EXPECT_TRUE( d1[ 4 ] == 10.0 );
    EXPECT_EQ( d1.size(), 5 );
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
    EXPECT_TRUE(( d3[ { 0, 0, 0 } ] == 4.2 ));
    EXPECT_TRUE(( d3[ { 1, 1, 1 } ] == 3.1 ));
    EXPECT_TRUE(( d3[ { 2, 0, 2 } ] == 2.0 ));
    d3[ { 0, 0, 0 } ] = 10.0;
    EXPECT_TRUE(( d3[ { 0, 0, 0 } ] == 10.0 ));
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

TEST_F( _TEST_TITLE_, ConstructorWorksWithDimarray ) {
    d3 = ndvector<3, double>( dims );
    EXPECT_TRUE( d3.shape() == dims );
}

TEST_F( _TEST_TITLE_, ImplicitConversionCopyConstructorWorks ) {
    ndvector<1, float> fvec1( d1 );
    for (size_t i = 0; i < d1.size(); ++i) {
        EXPECT_FLOAT_EQ( fvec1[ i ], d1[ i ] );
    }

    ndvector<2, float> fvec2( d2 );
    for (size_t i = 0; i < d2.shape()[ 0 ]; ++i) {
        for (size_t j = 0; j < d2.shape()[ 1 ]; ++j) {
            EXPECT_FLOAT_EQ( fvec2[ i ][ j ], d2[ i ][ j ] );
        }
    }

    ndvector<3, float> fvec3( d3 );
    for (size_t i = 0; i < d3.shape()[ 0 ]; ++i) {
        for (size_t j = 0; j < d3.shape()[ 1 ]; ++j) {
            for (size_t k = 0; k < d3.shape()[ 2 ]; ++k) {
                EXPECT_FLOAT_EQ( fvec3[ i ][ j ][ k ], d3[ i ][ j ][ k ] );
            }
        }
    }

    ndvector<3, int> ivec3( d3 );
    for (size_t i = 0; i < d3.shape()[ 0 ]; ++i) {
        for (size_t j = 0; j < d3.shape()[ 1 ]; ++j) {
            for (size_t k = 0; k < d3.shape()[ 2 ]; ++k) {
                EXPECT_EQ( ivec3[ i ][ j ][ k ], (int)d3[ i ][ j ][ k ] );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, CopyMethodWorks ) {
    ndvector<1, float> fvec1;
    ndvector<2, float> fvec2;
    ndvector<3, float> fvec3;
    ndvector<3, int> ivec3;

    fvec1.copy( d1 );
    for (size_t i = 0; i < d1.size(); ++i) {
        EXPECT_FLOAT_EQ( fvec1[ i ], d1[ i ] );
    }

    fvec2.copy( d2 );
    for (size_t i = 0; i < d2.shape()[ 0 ]; ++i) {
        for (size_t j = 0; j < d2.shape()[ 1 ]; ++j) {
            EXPECT_FLOAT_EQ( fvec2[ i ][ j ], d2[ i ][ j ] );
        }
    }

    fvec3.copy( d3 );
    for (size_t i = 0; i < d3.shape()[ 0 ]; ++i) {
        for (size_t j = 0; j < d3.shape()[ 1 ]; ++j) {
            for (size_t k = 0; k < d3.shape()[ 2 ]; ++k) {
                EXPECT_FLOAT_EQ( fvec3[ i ][ j ][ k ], d3[ i ][ j ][ k ] );
            }
        }
    }

    ivec3.copy( d3 );
    for (size_t i = 0; i < d3.shape()[ 0 ]; ++i) {
        for (size_t j = 0; j < d3.shape()[ 1 ]; ++j) {
            for (size_t k = 0; k < d3.shape()[ 2 ]; ++k) {
                EXPECT_EQ( ivec3[ i ][ j ][ k ], (int)d3[ i ][ j ][ k ] );
            }
        }
    }
}

template<typename T>
void check_values( T& container, details::_dimarray dims ) {
    for (auto i = 0; i < dims[ 0 ]; ++i) {
        for (auto j = 0; j < dims[ 1 ]; ++j) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                EXPECT_EQ( container( i, j, k ), i * 100 + j * 10 + k );
            }
        }
    }
}

template<typename T>
void assert_values( T& container, details::_dimarray dims ) {
    for (auto i = 0; i < dims[ 0 ]; ++i) {
        for (auto j = 0; j < dims[ 1 ]; ++j) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                ASSERT_EQ( container( i, j, k ), i * 100 + j * 10 + k );
            }
        }
    }
}

template<typename T>
void check_slices( T& container, details::_dimarray dims ) {
    for (auto i = 0; i < dims[ 0 ]; ++i) {
        auto sav = container.slice( 0, i );
        for (auto j = 0; j < dims[ 1 ]; ++j) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                EXPECT_EQ( sav( j, k ), i * 100 + j * 10 + k );
            }
        }
    }

    for (auto j = 0; j < dims[ 1 ]; ++j) {
        auto sav = container.slice( 1, j );
        for (auto i = 0; i < dims[ 0 ]; ++i) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                EXPECT_EQ( sav( i, k ), i * 100 + j * 10 + k );
            }
        }
    }

    for (auto k = 0; k < dims[ 2 ]; ++k) {
        auto sav = container.slice( 2, k );
        for (auto i = 0; i < dims[ 0 ]; ++i) {
            for (auto j = 0; j < dims[ 1 ]; ++j) {
                EXPECT_EQ( sav( i, j ), i * 100 + j * 10 + k );
            }
        }
    }
}

template<typename T>
void check_iterators( T& container, details::_dimarray dims ) {
    auto it = container.begin();
    for (auto i = 0; i < dims[ 0 ]; ++i) {
        for (auto j = 0; j < dims[ 1 ]; ++j) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                EXPECT_EQ( *( it++ ), i * 100 + j * 10 + k );
            }
        }
    }
}

template<typename T>
void check_iterators_negative( T& container, details::_dimarray dims ) {
    auto it = container.begin();
    for (auto i = 0; i < dims[ 0 ]; ++i) {
        for (auto j = 0; j < dims[ 1 ]; ++j) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                EXPECT_EQ( *( it++ ), -( i * 100 + j * 10 + k ) );
            }
        }
    }
}

template<typename T>
void check_iterators_writeable( T& container, details::_dimarray dims ) {
    auto it = container.begin();
    for (auto i = 0; i < dims[ 0 ]; ++i) {
        for (auto j = 0; j < dims[ 1 ]; ++j) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                *it *= -1;
                ++it;
            }
        }
    }

    check_iterators_negative( container, dims );
}

template<typename T>
void check_slices_writeable( T& container0, details::_dimarray dims ) {
    check_slices( container0, dims );

    for (auto i = 0; i < dims[ 0 ]; ++i) {
        auto container( container0 );
        assert_values( container, dims );
        auto sav = container.slice( 0, i );
        // cout << "Slicing at i=" << i << endl;
        for (auto j = 0; j < dims[ 1 ]; ++j) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                sav( j, k ) *= -1;
                // cout << "container(" << i << "," << j << "," << k
                //      << ") = " << container( i, j, k ) << endl;
            }
        }

        for (auto j = 0; j < dims[ 1 ]; ++j) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                // cout << "checking: container(" << i << "," << j << "," << k
                //      << ") = " << container( i, j, k ) << endl;
                EXPECT_EQ( container( i, j, k ), -( i * 100 + j * 10 + k ) );
            }
        }
        // cout << endl;
    }

    for (auto j = 0; j < dims[ 1 ]; ++j) {
        auto container( container0 );
        assert_values( container, dims );
        auto sav = container.slice( 1, j );
        for (auto i = 0; i < dims[ 0 ]; ++i) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                sav( i, k ) *= -1;
            }
        }
        for (auto i = 0; i < dims[ 0 ]; ++i) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                EXPECT_EQ( container( i, j, k ), -( i * 100 + j * 10 + k ) );
            }
        }
    }

    for (auto k = 0; k < dims[ 2 ]; ++k) {
        auto container = container0;
        assert_values( container, dims );
        auto sav = container.slice( 2, k );
        for (auto i = 0; i < dims[ 0 ]; ++i) {
            for (auto j = 0; j < dims[ 1 ]; ++j) {
                sav( i, j ) *= -1;
            }
        }

        for (auto i = 0; i < dims[ 0 ]; ++i) {
            for (auto j = 0; j < dims[ 1 ]; ++j) {
                EXPECT_EQ( container( i, j, k ), -( i * 100 + j * 10 + k ) );
            }
        }
    }
}

template<typename T>
void check_indexing( T& container, details::_dimarray dims ) {
    for (auto i = 0; i < dims[ 0 ]; ++i) {
        auto sav = container[ i ];
        for (auto j = 0; j < dims[ 1 ]; ++j) {
            for (auto k = 0; k < dims[ 2 ]; ++k) {
                EXPECT_EQ( sav( j, k ), i * 100 + j * 10 + k );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, SmartArraySetCorrectly ) {
    check_values( sa, dims );
}

TEST_F( _TEST_TITLE_, SmartArraySlicesCorrectly ) {
    check_slices( sa, dims );
}

TEST_F( _TEST_TITLE_, SmartArraySlicesWriteable ) {
    check_slices_writeable( sa, dims );
}

TEST_F( _TEST_TITLE_, SmartArrayIndexesCorrectly ) {
    check_indexing( sa, dims );
}

TEST_F( _TEST_TITLE_, SmartArrayIteratorsReadCorrectly ) {
    check_iterators( sa, dims );
}

TEST_F( _TEST_TITLE_, SmartArrayIteratorsWriteCorrectly ) {
    check_iterators_writeable( sa, dims );
}

TEST_F( _TEST_TITLE_, SmartArrayResizesCorrectly ) {
    assert_values( sa, dims );
    sa.resize( 2, 3, 4 );
    EXPECT_EQ( sa.size(), 2*3*4 );
    for (auto val : sa) {
        EXPECT_EQ( val, 0 );
    }
}

TEST_F( _TEST_TITLE_, SmartVectorSetCorrectly ) {
    check_values( sv, dims );
}

TEST_F( _TEST_TITLE_, SmartVectorSlicesCorrectly ) {
    check_slices( sv, dims );
}

TEST_F( _TEST_TITLE_, SmartVectorIndexesCorrectly ) {
    check_indexing( sv, dims );
}

TEST_F( _TEST_TITLE_, SmartVectorIteratorsReadCorrectly ) {
    check_iterators( sv, dims );
}

TEST_F( _TEST_TITLE_, SmartVectorIteratorsWriteCorrectly ) {
    check_iterators_writeable( sv, dims );
}

TEST_F( _TEST_TITLE_, SmartVectorSlicesWriteable ) {
    check_slices_writeable( sv, dims );
}

TEST_F( _TEST_TITLE_, SmartVectorResizesCorrectly ) {
    assert_values( sv, dims );
    sv.resize( 2, 3, 4 );
    EXPECT_EQ( sv.size(), 2*3*4 );
    for (auto val : sv) {
        EXPECT_EQ( val, 0 );
    }
}

TEST_F( _TEST_TITLE_, SmartVectorCanRerank ) {
    assert_values( sa, dims );
    EXPECT_NO_THROW( { sv.resize(2,3); } );
    EXPECT_EQ( sv.size(), 2*3 );
    for (auto val : sv) {
        EXPECT_EQ( val, 0 );
    }
}