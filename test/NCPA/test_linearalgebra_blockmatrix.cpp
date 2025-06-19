#define NCPA_DEBUG_ON
#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/linearalgebra.hpp"
#include "NCPA/logging.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::linear;

#define _TEST_EQ_       EXPECT_COMPLEX_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_COMPLEX_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraLibraryBandDiagonalComplexMatrixTest

typedef complex<double> test_t;
typedef band_diagonal_matrix<test_t> mat_t;
typedef BlockMatrix<test_t> bmat_t;
typedef dense_vector<test_t> vec_t;

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here

            blockrows = 3;
            blockcols = 3;
            tilerows = 5;
            tilecols = 5;

            blockmat.block_type( matrix_t::BAND_DIAGONAL ).resize( blockrows, blockcols, tilerows, tilecols ).identity();
            empty.block_type( matrix_t::BAND_DIAGONAL );


        }  // void TearDown() override {}

        // declare stuff here
        bmat_t blockmat, empty;
        mat_t itile, ztile;
        size_t blockrows, blockcols, tilerows, tilecols;
};

TEST_F( _TEST_TITLE_, DefaultConstructorIsCorrect ) {
    EXPECT_EQ( empty.rows(), 0 );
    EXPECT_EQ( empty.columns(), 0 );
}

TEST_F( _TEST_TITLE_, SizedConstructorIsCorrect ) {
    EXPECT_EQ( blockmat.rows(), blockrows * tilerows );
    EXPECT_EQ( blockmat.columns(), blockcols * tilecols );
    EXPECT_EQ( blockmat.block_rows(), blockrows );
    EXPECT_EQ( blockmat.block_columns(), blockcols );
    EXPECT_EQ( blockmat.rows_per_block(), tilerows );
    EXPECT_EQ( blockmat.columns_per_block(), tilecols );
}

TEST_F( _TEST_TITLE_, EqualsReturnsTrueForEqual ) {
    EXPECT_TRUE( blockmat.equals( blockmat ) );
}

TEST_F( _TEST_TITLE_, EqualsReturnsFalseForUnequal ) {
    EXPECT_FALSE( blockmat.equals( empty ) );
}

TEST_F( _TEST_TITLE_, CopyConstructorWorks ) {
    empty = bmat_t( blockmat );
    EXPECT_TRUE( empty.equals( blockmat ) );
    blockmat.clear();
    EXPECT_FALSE( empty.equals( blockmat ) );
}

TEST_F( _TEST_TITLE_, AssignmentOperatorWorks ) {
    empty = blockmat;
    EXPECT_TRUE( empty.equals( blockmat ) );
    blockmat.clear();
    EXPECT_FALSE( empty.equals( blockmat ) );
}
