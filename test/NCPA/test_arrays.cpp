#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <numbers>

using namespace std;
using namespace NCPA::arrays;
using namespace testing;

TEST( NCPAArraysLibraryTest, ZerosCreatesArray ) {
    double *testArray    = zeros<double>( 5 );
    double shouldBe[ 5 ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
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
    int offset            = 2;
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
    vector<int> inds  = index_vector<int>( 5 );
    int *testArray    = &inds[ 0 ];
    int shouldBe[ 5 ] = { 2, 3, 4, 0, 1 };
    circshift( testArray, 5, 2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAArraysLibraryTest, CircShiftShiftsBackward ) {
    vector<int> inds  = index_vector<int>( 5 );
    int *testArray    = &inds[ 0 ];
    int shouldBe[ 5 ] = { 3, 4, 0, 1, 2 };
    circshift( testArray, 5, -2, testArray );
    EXPECT_ARRAY_EQ( 5, testArray, shouldBe );
}

TEST( NCPAArraysLibraryTest, CircShiftZeroDoesNotShift ) {
    vector<int> inds  = index_vector<int>( 5 );
    int *testArray    = &inds[ 0 ];
    int shouldBe[ 5 ] = { 0, 1, 2, 3, 4 };
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

TEST( NCPAArraysLibraryTest, ReverseReversesOrder ) {
    vector<double> testVector = index_vector<double>( 11 );
    double *testArray         = zeros<double>( 11 );
    std::copy( testVector.cbegin(), testVector.cend(), testArray );
    std::reverse( testVector.begin(), testVector.end() );
    NCPA::arrays::reverse( testArray, 11, testArray );
    EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST( NCPAArraysLibraryTest, ReverseReversesOrderPartially ) {
    vector<double> testVector = index_vector<double>( 11 );
    double *testArray         = zeros<double>( 11 );
    std::copy( testVector.cbegin(), testVector.cend(), testArray );
    std::reverse( testVector.begin(), testVector.begin() + 5 );
    NCPA::arrays::reverse( testArray, 5, testArray );
    EXPECT_ARRAY_DOUBLE_EQ( 11, testArray, testVector );
}

TEST( NCPAArraysLibraryTest, AsArrayLinksToVector ) {
    std::vector<int> vi = { 0, 1, 2, 3, 4 };
    int *pi             = as_array( vi );
    for ( auto i = 0; i < 5; i++ ) {
        ASSERT_EQ( pi[ i ], vi[ i ] );
    }
    vi[ 0 ] = 5;
    for ( auto i = 0; i < 5; i++ ) {
        ASSERT_EQ( pi[ i ], vi[ i ] );
    }
}

TEST( NCPAArraysLibraryTest, SortPermutationCreatedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<size_t> perm_expected { 5, 6, 3, 0, 2, 1, 4 };
    auto perm = sort_permutation(
        unsorted, []( int const& a, int const& b ) { return a < b; } );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST( NCPAArraysLibraryTest, SortPermutationIncreasingCreatedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<size_t> perm_expected { 5, 6, 3, 0, 2, 1, 4 };
    auto perm = sort_permutation_increasing( unsorted );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST( NCPAArraysLibraryTest, SortPermutationDecreasingCreatedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<size_t> perm_expected { 4, 1, 2, 0, 3, 5, 6 };
    auto perm = sort_permutation_decreasing( unsorted );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST( NCPAArraysLibraryTest, SortPermutationAppliedProperly ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<double> orig { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    std::vector<double> sorted_expected { 5.0, 6.0, 3.0, 0.0, 2.0, 1.0, 4.0 };
    auto perm = sort_permutation(
        unsorted, []( int const& a, int const& b ) { return a < b; } );
    auto sorted = apply_permutation( orig, perm );
    EXPECT_ARRAY_DOUBLE_EQ( perm.size(), sorted, sorted_expected );
}

TEST( NCPAArraysLibraryTest, SortPermutationAppliedProperlyInPlace ) {
    std::vector<int> unsorted { 1, 4, 2, -1, 5, -3, -3 };
    std::vector<double> orig { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    std::vector<double> sorted_expected { 5.0, 6.0, 3.0, 0.0, 2.0, 1.0, 4.0 };
    auto perm = sort_permutation(
        unsorted, []( int const& a, int const& b ) { return a < b; } );
    apply_permutation_in_place( orig, perm );
    EXPECT_ARRAY_DOUBLE_EQ( perm.size(), orig, sorted_expected );
}

TEST( NCPAArraysLibraryTest, DoubleSortPermutationCreatedProperly ) {
    std::vector<int> unsorted1 { 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5 };
    std::vector<int> unsorted2 { 2, 1, 3, 4, 0, 0, 1, 2, 4, 3, 2, 2, 1 };
    std::vector<size_t> perm_expected { 1, 0, 2, 4, 5, 3, 6, 7, 9, 8, 10, 11, 12 };
    auto perm = sort_permutation_increasing( unsorted1, unsorted2 );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}

TEST( NCPAArraysLibraryTest, DoubleSortPermutationDecreasingCreatedProperly ) {
    std::vector<int> unsorted1 { 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5 };
    std::vector<int> unsorted2 { 2, 1, 3, 4, 0, 0, 1, 2, 4, 3, 2, 2, 1 };
    std::vector<size_t> perm_expected { 12, 10, 11, 8, 9, 7, 6, 3, 4, 5, 2, 0, 1 };
    auto perm = sort_permutation_decreasing( unsorted1, unsorted2 );
    EXPECT_ARRAY_EQ( perm.size(), perm, perm_expected );
}