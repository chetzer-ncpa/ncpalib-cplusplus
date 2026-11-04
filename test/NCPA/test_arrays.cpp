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

TEST( _TEST_TITLE_, ZerosCreatesArray ) {
    double *testArray    = zeros<double>( 5 );
    double shouldBe[ 5 ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    ASSERT_ARRAY_EQ( 5, testArray, shouldBe );
    // for (auto i = 0; i < 10; i++) {
    // 	EXPECT_EQ( testArray[i], 0.0 );
    // }
    delete[] testArray;
}

TEST( _TEST_TITLE_, ZerosCreates2DArray ) {
    double **testArray = zeros<double>( 10, 5 );
    for (auto i = 0; i < 10; i++) {
        for (auto j = 0; j < 5; j++) {
            EXPECT_EQ( testArray[ i ][ j ], 0.0 );
        }
    }
}

TEST( _TEST_TITLE_, FreeArray2DFreesArrayWithoutError ) {
    double **testArray = zeros<double>( 10, 5 );
    free_array( testArray, 10, 5 );
    EXPECT_EQ( testArray, nullptr );
}

TEST( _TEST_TITLE_, ZerosCreates3DArray ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    for (auto i = 0; i < 10; i++) {
        for (auto j = 0; j < 5; j++) {
            for (auto k = 0; k < 3; k++) {
                EXPECT_EQ( testArray[ i ][ j ][ k ], 0.0 );
            }
        }
    }
}

TEST( _TEST_TITLE_, FreeArray3DFreesArrayWithoutError ) {
    double ***testArray = zeros<double>( 10, 5, 3 );
    free_array( testArray, 10, 5, 3 );
    ASSERT_EQ( testArray, nullptr );
}

TEST( _TEST_TITLE_, IndexVectorCreatesExpectedVectors ) {
    vector<int> testArray = index_vector<int>( 5 );
    for (int i = 0; i < 5; i++) {
        EXPECT_EQ( i, testArray[ i ] );
    }
    vector<double> testArray2 = index_vector<double>( 5 );
    for (int i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( (double)i, testArray2[ i ] );
    }
}

TEST( _TEST_TITLE_, IndexVectorCreatesExpectedVectorsWithOffset ) {
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

TEST( _TEST_TITLE_, CircShiftShiftsForward ) {
    vector<int> inds  = index_vector<int>( 5 );
    int *testArray    = &inds[ 0 ];
    int shouldBe[ 5 ] = { 2, 3, 4, 0, 1 };
    circshift( testArray, 5, 2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( _TEST_TITLE_, CircShiftShiftsBackward ) {
    vector<int> inds  = index_vector<int>( 5 );
    int *testArray    = &inds[ 0 ];
    int shouldBe[ 5 ] = { 3, 4, 0, 1, 2 };
    circshift( testArray, 5, -2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( _TEST_TITLE_, CircShiftZeroDoesNotShift ) {
    vector<int> inds  = index_vector<int>( 5 );
    int *testArray    = &inds[ 0 ];
    int shouldBe[ 5 ] = { 0, 1, 2, 3, 4 };
    circshift( testArray, 5, 0, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( _TEST_TITLE_, FillArray2DSetsValuesCorrectly ) {
    double **testArray = zeros<double>( 10, 5 );
    fill( testArray, 10, 5, 4.5 );
    for (auto i = 0; i < 10; i++) {
        for (auto j = 0; j < 5; j++) {
            EXPECT_DOUBLE_EQ( testArray[ i ][ j ], 4.5 );
        }
    }
    free_array( testArray, 10, 5 );
}

TEST( _TEST_TITLE_, FillArray3DSetsValuesCorrectly ) {
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

TEST( _TEST_TITLE_, ReverseReversesOrder ) {
    vector<double> testVector = index_vector<double>( 11 );
    double *testArray         = zeros<double>( 11 );
    std::copy( testVector.cbegin(), testVector.cend(), testArray );
    std::reverse( testVector.begin(), testVector.end() );
    NCPA::arrays::reverse( testArray, 11, testArray );
    EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST( _TEST_TITLE_, ReverseReversesOrderPartially ) {
    vector<double> testVector = index_vector<double>( 11 );
    double *testArray         = zeros<double>( 11 );
    std::copy( testVector.cbegin(), testVector.cend(), testArray );
    std::reverse( testVector.begin(), testVector.begin() + 5 );
    NCPA::arrays::reverse( testArray, 5, testArray );
    EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST( _TEST_TITLE_, AsArrayLinksToVector ) {
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

TEST( _TEST_TITLE_, SortPermutationCreatedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<size_t> perm_expected { 5, 6, 3, 0, 2, 1, 4 };
    auto perm = sort_permutation(
        unsorted, []( int const& a, int const& b ) { return a < b; } );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST( _TEST_TITLE_, SortPermutationIncreasingCreatedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<size_t> perm_expected { 5, 6, 3, 0, 2, 1, 4 };
    auto perm = sort_permutation_increasing( unsorted );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST( _TEST_TITLE_, SortPermutationDecreasingCreatedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<size_t> perm_expected { 4, 1, 2, 0, 3, 5, 6 };
    auto perm = sort_permutation_decreasing( unsorted );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST( _TEST_TITLE_, SortPermutationAppliedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<double> orig { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    std::vector<double> sorted_expected { 5.0, 6.0, 3.0, 0.0, 2.0, 1.0, 4.0 };
    auto perm = sort_permutation(
        unsorted, []( int const& a, int const& b ) { return a < b; } );
    auto sorted = apply_permutation( orig, perm );
    EXPECT_ARRAY_DOUBLE_EQ( perm.size(), sorted, sorted_expected );
}

TEST( _TEST_TITLE_, SortPermutationAppliedProperlyInPlace ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<double> orig { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    std::vector<double> sorted_expected { 5.0, 6.0, 3.0, 0.0, 2.0, 1.0, 4.0 };
    auto perm = sort_permutation(
        unsorted, []( int const& a, int const& b ) { return a < b; } );
    apply_permutation_in_place( orig, perm );
    EXPECT_ARRAY_DOUBLE_EQ( perm.size(), orig, sorted_expected );
}

TEST( _TEST_TITLE_, DoubleSortPermutationCreatedProperly ) {
    std::vector<int> unsorted1 { 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5 };
    std::vector<int> unsorted2 { 2, 1, 3, 4, 0, 0, 1, 2, 4, 3, 2, 2, 1 };
    std::vector<size_t> perm_expected { 1, 0, 2, 4,  5,  3, 6,
                                        7, 9, 8, 10, 11, 12 };
    auto perm = sort_permutation_increasing( unsorted1, unsorted2 );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST( _TEST_TITLE_, DoubleSortPermutationDecreasingCreatedProperly ) {
    std::vector<int> unsorted1 { 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5 };
    std::vector<int> unsorted2 { 2, 1, 3, 4, 0, 0, 1, 2, 4, 3, 2, 2, 1 };
    std::vector<size_t> perm_expected { 12, 10, 11, 8, 9, 7, 6,
                                        3,  4,  5,  2, 0, 1 };
    auto perm = sort_permutation_decreasing( unsorted1, unsorted2 );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST( _TEST_TITLE_, AddVectorsWorksCorrectly ) {
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

TEST( _TEST_TITLE_, AddVectorsCommutesForVectors ) {
    vector<double> x  = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y  = NCPA::math::random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z1 = add_vectors( x, y );
    vector<double> z2 = add_vectors( y, x );
    EXPECT_ARRAY_DOUBLE_EQ( 8, z1, z2 );
}

TEST( _TEST_TITLE_, AddArraysWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double *z        = NCPA::arrays::zeros<double>( 5 );
    add_arrays( 5, &x[ 0 ], &y[ 0 ], z );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] + y[ i ], z[ i ] );
    }
}

TEST( _TEST_TITLE_, AddArraysCommutes ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double *z1       = NCPA::arrays::zeros<double>( 5 );
    double *z2       = NCPA::arrays::zeros<double>( 5 );
    add_arrays( 5, &x[ 0 ], &y[ 0 ], z1 );
    add_arrays( 5, &y[ 0 ], &x[ 0 ], z2 );
    EXPECT_ARRAY_DOUBLE_EQ( 5, z1, z2 );
}

TEST( _TEST_TITLE_, DivideVectorsWorksCorrectly ) {
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

TEST( _TEST_TITLE_, DivideArraysWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.1, 1.0 );
    double *z        = NCPA::arrays::zeros<double>( 5 );
    divide_arrays( 5, &x[ 0 ], &y[ 0 ], z );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] / y[ i ], z[ i ] );
    }
}

TEST( _TEST_TITLE_, MultiplyVectorsWorksCorrectly ) {
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

TEST( _TEST_TITLE_, MultiplyArraysWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, -1.0, 1.0 );
    double *z        = NCPA::arrays::zeros<double>( 5 );
    multiply_arrays( 5, &x[ 0 ], &y[ 0 ], z );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] * y[ i ], z[ i ] );
    }
}

TEST( _TEST_TITLE_, MultiplyVectorsCommutes ) {
    vector<double> x  = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y  = NCPA::math::random_numbers<double>( 8, 0.0, 1.0 );
    vector<double> z1 = multiply_vectors( x, y );
    vector<double> z2 = multiply_vectors( y, x );
    EXPECT_ARRAY_DOUBLE_EQ( 8, z1, z2 );
}

TEST( _TEST_TITLE_, MultiplyArraysCommutes ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double *z1       = NCPA::arrays::zeros<double>( 5 );
    double *z2       = NCPA::arrays::zeros<double>( 5 );
    multiply_arrays<double>( 5, &x[ 0 ], &y[ 0 ], z1 );
    multiply_arrays<double>( 5, &y[ 0 ], &x[ 0 ], z2 );
    EXPECT_ARRAY_DOUBLE_EQ( 5, z1, z2 );
}

TEST( _TEST_TITLE_, ScaleVectorWorksProperly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double scalar    = NCPA::math::random_number<double>( -3.0, 3.0 );
    vector<double> y = scale_vector( x, scalar );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( y[ i ], x[ i ] * scalar );
    }
}

TEST( _TEST_TITLE_, ScaleArrayWorksProperly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double scalar    = NCPA::math::random_number<double>( -3.0, 3.0 );
    double *y        = NCPA::arrays::zeros<double>( 5 );
    scale_array( 5, &x[ 0 ], scalar, y );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( y[ i ], x[ i ] * scalar );
    }
}

TEST( _TEST_TITLE_, ScaleArrayWorksProperlyInPlace ) {
    double *x         = NCPA::arrays::zeros<double>( 5 );
    vector<double> xv = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    std::copy( xv.cbegin(), xv.cend(), x );
    double scalar = NCPA::math::random_number<double>( -3.0, 3.0 );
    scale_array( 5, x, scalar );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ], xv[ i ] * scalar );
    }
}

TEST( _TEST_TITLE_, SubtractVectorsWorksCorrectly ) {
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

TEST( _TEST_TITLE_, SubtractArraysWorksCorrectly ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 ),
                   y = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double *z        = NCPA::arrays::zeros<double>( 5 );
    subtract_arrays( 5, &x[ 0 ], &y[ 0 ], z );
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] - y[ i ], z[ i ] );
    }
}

TEST( _TEST_TITLE_, PlusOperatorWorksForVectors ) {
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

TEST( _TEST_TITLE_, PlusOperatorWorksForVectorScalar ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double y         = 3.5;
    vector<double> z = x + y;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] + y, z[ i ] );
    }
}

TEST( _TEST_TITLE_, MinusOperatorWorksForVectors ) {
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

TEST( _TEST_TITLE_, MinusOperatorWorksForVectorScalar ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double y         = 3.5;
    vector<double> z = x - y;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] - y, z[ i ] );
    }
}

TEST( _TEST_TITLE_, NegationOperatorWorksForVector ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    vector<double> z = -x;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( -x[ i ], z[ i ] );
    }
}

TEST( _TEST_TITLE_, MultiplyOperatorWorksForVectors ) {
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

TEST( _TEST_TITLE_, MultiplyOperatorWorksForVectorScalar ) {
    vector<double> x = NCPA::math::random_numbers<double>( 5, 0.0, 1.0 );
    double y         = 3.5;
    vector<double> z = x / y;
    for (auto i = 0; i < 5; i++) {
        EXPECT_DOUBLE_EQ( x[ i ] / y, z[ i ] );
    }
}

template<typename Container>
void make_test_array( Container& a ) {
    a.resize( 2, 3, 4 );
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );
}

TEST( _TEST_TITLE_, NDArrayAtWorksCorrectly ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDArray<3, int> a( dims );

    make_test_array( a );

    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( a.at( i, j, k ), int( 100 * i + 10 * j + k ) );
}

TEST( _TEST_TITLE_, NDArrayShapeFunctionIsCorrect ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDArray<3, int> a( dims );
    auto shape = a.shape();
    EXPECT_EQ( shape.size(), dims.size() );
    for (size_t i = 0; i < dims.size(); ++i) {
        EXPECT_EQ( shape.at( i ), dims.at( i ) );
    }
}

TEST( _TEST_TITLE_, NDArraySliceWorksCorrectly ) {
    std::vector<size_t> dims { 2, 3, 4 };
    NDArray<3, int> a( dims );
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );

    for (size_t i = 0; i < 2; ++i) {
        auto plane = a.slice( 0, i );
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( plane.at( j, k ), int( 100 * i + 10 * j + k ) );
    }

    for (size_t j = 0; j < 3; ++j) {
        auto plane = a.slice( 1, j );
        for (size_t i = 0; i < 2; ++i)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( plane.at( i, k ), int( 100 * i + 10 * j + k ) );
    }

    for (size_t k = 0; k < 4; ++k) {
        auto plane = a.slice( 2, k );
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 3; ++j)
                EXPECT_EQ( plane.at( i, j ), int( 100 * i + 10 * j + k ) );
    }
}

TEST( _TEST_TITLE_, NDArrayCopyConstructorsWork ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDArray<3, int> a( dims );

    make_test_array( a );

    NDArray<3, int> b( a );
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( a.at( i, j, k ), b.at( i, j, k ) );

    NDVector<int> v;
    make_test_array( v );

    NDArray<3, int> c( v );
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( a.at( i, j, k ), c.at( i, j, k ) );

    auto plane = a.slice( 1, 1 );
    NDArray<2, int> d( plane );
    for (size_t i = 0; i < 2; ++i)
        for (size_t k = 0; k < 4; ++k)
            EXPECT_EQ( a.at( i, 1, k ), d.at( i, k ) );

    auto vplane = v.slice( 1, 1 );
    NDArray<2, int> e( vplane );
    for (size_t i = 0; i < 2; ++i)
        for (size_t k = 0; k < 4; ++k)
            EXPECT_EQ( v.at( i, 1, k ), e.at( i, k ) );
}

TEST( _TEST_TITLE_, NDArrayBracketOperatorWorksCorrectly ) {
    NDArray<3, int> a( 2, 3, 4 );

    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );

    for (size_t i = 0; i < 2; ++i) {
        auto plane = a[ i ];
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( plane.at( j, k ), int( 100 * i + 10 * j + k ) );
    }

    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            auto row = a[ i ][ j ];
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( row.at( k ), int( 100 * i + 10 * j + k ) );
        }
    }
}

TEST( _TEST_TITLE_, NDArrayIteratorWorks ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDArray<3, int> a( dims );
    std::vector<int> flat;
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k) {
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );
                flat.push_back( a.at( i, j, k ) );
            }

    auto flat_it = flat.begin();
    auto a_it    = a.begin();
    while (flat_it != flat.end() && a_it != a.end()) {
        EXPECT_EQ( *flat_it, *a_it );
        ++flat_it;
        ++a_it;
    }
}

TEST( _TEST_TITLE_, NDArrayReverseIteratorWorks ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDArray<3, int> a( dims );
    std::vector<int> flat;
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k) {
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );
                flat.push_back( a.at( i, j, k ) );
            }

    auto flat_it = flat.rbegin();
    auto a_it    = a.rbegin();
    while (flat_it != flat.rend() && a_it != a.rend()) {
        EXPECT_EQ( *flat_it, *a_it );
        ++flat_it;
        ++a_it;
    }
}

TEST( _TEST_TITLE_, NDArrayClearMethodWorks ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDArray<3, int> a( dims );
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );

    a.clear();
    for (auto x : a) {
        EXPECT_EQ( x, 0 );
    }
}

TEST( _TEST_TITLE_, NDVectorAtWorksCorrectly ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDVector<int> a( dims );

    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );

    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( a.at( i, j, k ), int( 100 * i + 10 * j + k ) );
}

TEST( _TEST_TITLE_, NDVectorShapeFunctionIsCorrect ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDVector<int> a( dims );
    auto shape = a.shape();
    EXPECT_EQ( shape.size(), dims.size() );
    for (size_t i = 0; i < dims.size(); ++i) {
        EXPECT_EQ( shape.at( i ), dims.at( i ) );
    }
}

TEST( _TEST_TITLE_, NDVectorSliceWorksCorrectly ) {
    std::vector<size_t> dims { 2, 3, 4 };
    NDVector<int> a( dims );
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );

    for (size_t i = 0; i < 2; ++i) {
        auto plane = a.slice( 0, i );
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( plane.at( j, k ), int( 100 * i + 10 * j + k ) );
    }

    for (size_t j = 0; j < 3; ++j) {
        auto plane = a.slice( 1, j );
        for (size_t i = 0; i < 2; ++i)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( plane.at( i, k ), int( 100 * i + 10 * j + k ) );
    }

    for (size_t k = 0; k < 4; ++k) {
        auto plane = a.slice( 2, k );
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 3; ++j)
                EXPECT_EQ( plane.at( i, j ), int( 100 * i + 10 * j + k ) );
    }
}


TEST( _TEST_TITLE_, NDVectorCopyConstructorsWork ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDVector<int> a( dims );

    make_test_array( a );

    NDVector<int> b( a );
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( a.at( i, j, k ), b.at( i, j, k ) );

    NDArray<3,int> v( dims );
    make_test_array( v );

    NDVector<int> c( v );
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( a.at( i, j, k ), c.at( i, j, k ) );

    auto plane = a.slice( 1, 1 );
    NDVector<int> d( plane );
    for (size_t i = 0; i < 2; ++i)
        for (size_t k = 0; k < 4; ++k)
            EXPECT_EQ( a.at( i, 1, k ), d.at( i, k ) );

    auto vplane = v.slice( 1, 1 );
    NDVector<int> e( vplane );
    for (size_t i = 0; i < 2; ++i)
        for (size_t k = 0; k < 4; ++k)
            EXPECT_EQ( v.at( i, 1, k ), e.at( i, k ) );
}

TEST( _TEST_TITLE_, NDVectorBracketOperatorWorksCorrectly ) {
    NDVector<int> a( 2, 3, 4 );

    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );

    for (size_t i = 0; i < 2; ++i) {
        auto plane = a[ i ];
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( plane.at( j, k ), int( 100 * i + 10 * j + k ) );
    }

    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            auto row = a[ i ][ j ];
            for (size_t k = 0; k < 4; ++k)
                EXPECT_EQ( row.at( k ), int( 100 * i + 10 * j + k ) );
        }
    }
}

TEST( _TEST_TITLE_, NDVectorIteratorWorks ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDVector<int> a( dims );
    std::vector<int> flat;
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k) {
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );
                flat.push_back( a.at( i, j, k ) );
            }

    auto flat_it = flat.begin();
    auto a_it    = a.begin();
    while (flat_it != flat.end() && a_it != a.end()) {
        EXPECT_EQ( *flat_it, *a_it );
        ++flat_it;
        ++a_it;
    }
}

TEST( _TEST_TITLE_, NDVectorReverseIteratorWorks ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDVector<int> a( dims );
    std::vector<int> flat;
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k) {
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );
                flat.push_back( a.at( i, j, k ) );
            }

    auto flat_it = flat.rbegin();
    auto a_it    = a.rbegin();
    while (flat_it != flat.rend() && a_it != a.rend()) {
        EXPECT_EQ( *flat_it, *a_it );
        ++flat_it;
        ++a_it;
    }
}

TEST( _TEST_TITLE_, NDVectorClearMethodWorks ) {
    std::vector<size_t> dims { 2, 3, 4 };
    std::vector<int> inds( 3 );
    NDVector<int> a( dims );
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 4; ++k)
                a.at( i, j, k ) = int( 100 * i + 10 * j + k );

    a.clear();
    for (auto x : a) {
        EXPECT_EQ( x, 0 );
    }
}
