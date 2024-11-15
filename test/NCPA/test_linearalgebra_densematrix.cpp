#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/linearalgebra.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::linear;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ


typedef double test_t;
typedef details::dense_matrix<test_t> mat_t;

class NCPALinearAlgebraDenseMatrixTest : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here

            square    = mat_t( dim1, dim1 );
            more_rows = mat_t( dim2, dim1 );
            more_cols = mat_t( dim1, dim2 );
            product   = mat_t( dim2, dim2 );
            product2  = mat_t( dim1, dim1 );
            symmetric = mat_t( dim1, dim1 );
            identity  = mat_t( dim1, dim1 );
            for ( size_t i = 0; i < dim1; i++ ) {
                double di = (double)( i + 1 );
                square.set_row( i, { di, di, di } );
                more_rows.set_column( i, { di, di, di, di, di } );
                more_cols.set_row( i, { di, di, di, di, di } );
                symmetric.set_row( i, { 0, 0, 0 } );
                identity.set( i, i, 1.0 );
                for ( size_t j = 0; j < dim1; j++ ) {
                    product2.set( i, j,
                                  (double)( 5 * ( i + 1 ) * ( j + 1 ) ) );
                }
            }
            symmetric.set( 1, 1, testval );
            product.set( 14 );

        }  // void TearDown() override {}

        // declare stuff here
        details::dense_matrix<test_t> empty, square, more_rows, more_cols,
            product, identity, symmetric, product2;
        const size_t dim1 = 3, dim2 = 5;
        test_t testval = -4.2;
};

TEST_F( NCPALinearAlgebraDenseMatrixTest, DefaultConstructorIsCorrect ) {
    EXPECT_EQ( empty.rows(), 0 );
    EXPECT_EQ( empty.columns(), 0 );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, SizedConstructorIsCorrect ) {
    EXPECT_EQ( square.rows(), dim1 );
    EXPECT_EQ( square.columns(), dim1 );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, EqualsReturnsTrueForEqual ) {
    EXPECT_TRUE( square.equals( square ) );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, EqualsReturnsFalseForUnequal ) {
    EXPECT_FALSE( square.equals( more_rows ) );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, CopyConstructorWorks ) {
    empty = mat_t( square );
    EXPECT_TRUE( empty.equals( square ) );
    square.clear();
    EXPECT_FALSE( empty.equals( square ) );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, AssignmentOperatorWorks ) {
    empty = square;
    EXPECT_TRUE( empty.equals( square ) );
    square.clear();
    EXPECT_FALSE( empty.equals( square ) );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, SizeMethodsReturnExpectedValues ) {
    EXPECT_EQ( square.rows(), dim1 );
    EXPECT_EQ( square.columns(), dim1 );
    EXPECT_EQ( more_rows.rows(), dim2 );
    EXPECT_EQ( more_rows.columns(), dim1 );
    EXPECT_EQ( more_cols.rows(), dim1 );
    EXPECT_EQ( more_cols.columns(), dim2 );
    EXPECT_EQ( product.rows(), dim2 );
    EXPECT_EQ( product.columns(), dim2 );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, ResizeWorks ) {
    square.resize( 2, 7 );
    EXPECT_EQ( square.rows(), 2 );
    EXPECT_EQ( square.columns(), 7 );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, ResizePreservesExistingData ) {
    square.resize( dim1 + 1, dim1 + 1 );
    for ( size_t i = 0; i < dim1 + 1; i++ ) {
        for ( size_t j = 0; j < dim1 + 1; j++ ) {
            if ( i < dim1 && j < dim1 ) {
                _TEST_EQ_( more_cols.get( i, j ), square.get( i, j ) );
            } else {
                _TEST_EQ_( square.get( i, j ), 0.0 );
            }
        }
    }

    square.resize( dim1 - 1, dim1 - 1 );
    for ( size_t i = 0; i < dim1 - 1; i++ ) {
        for ( size_t j = 0; j < dim1 - 1; j++ ) {
            _TEST_EQ_( more_cols.get( i, j ), square.get( i, j ) );
        }
    }
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, SwapWorks ) {
    ASSERT_EQ( more_cols.rows(), dim1 );
    ASSERT_EQ( more_cols.columns(), dim2 );
    ASSERT_EQ( more_rows.rows(), dim2 );
    ASSERT_EQ( more_rows.columns(), dim1 );
    swap( more_rows, more_cols );
    EXPECT_EQ( more_rows.rows(), dim1 );
    EXPECT_EQ( more_rows.columns(), dim2 );
    EXPECT_EQ( more_cols.rows(), dim2 );
    EXPECT_EQ( more_cols.columns(), dim1 );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, GetMethodAllowsModification ) {
    square.get( 0, 0 ) = testval;
    _TEST_EQ_( square.get( 0, 0 ), testval );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, GetRowReturnsCopy ) {
    for ( size_t i = 0; i < dim1; i++ ) {
        for ( size_t j = 0; j < dim1; j++ ) {
            _TEST_EQ_( (double)( i + 1 ), square.get_row( i )->get( j ) );
        }
    }
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, GetColumnReturnsCopy ) {
    for ( size_t i = 0; i < dim1; i++ ) {
        for ( size_t j = 0; j < dim1; j++ ) {
            _TEST_EQ_( (double)( j + 1 ), square.get_column( i )->get( j ) );
        }
    }
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, GetDiagonalReturnsExpectedValues ) {
    std::vector<test_t> diag0( { 1.0, 2.0, 3.0 } ), diag1( { 1.0, 2.0 } ),
        diag2( { 1.0 } ), diag_1( { 2.0, 3.0 } ), diag_2( { 3.0 } );
    _TEST_ARRAY_EQ_( 3, square.get_diagonal()->as_std(), diag0 );
    _TEST_ARRAY_EQ_( 2, square.get_diagonal( 1 )->as_std(), diag1 );
    _TEST_ARRAY_EQ_( 1, square.get_diagonal( 2 )->as_std(), diag2 );
    _TEST_ARRAY_EQ_( 2, square.get_diagonal( -1 )->as_std(), diag_1 );
    _TEST_ARRAY_EQ_( 1, square.get_diagonal( -2 )->as_std(), diag_2 );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, ClearWorksProperly ) {
    square.clear();
    EXPECT_EQ( square.rows(), 0 );
    EXPECT_EQ( square.columns(), 0 );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, CloneMethodsWork ) {
    EXPECT_EQ( square.clone()->rows(), dim1 );
    EXPECT_EQ( square.clone()->columns(), dim1 );
    EXPECT_EQ( square.fresh_clone()->rows(), 0 );
    EXPECT_EQ( square.fresh_clone()->columns(), 0 );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, AsArrayWorksWithNullPointer ) {
    test_t **testarr = nullptr;
    size_t nrows = 0, ncols = 0;
    more_rows.as_array( nrows, ncols, testarr );
    EXPECT_EQ( nrows, more_rows.rows() );
    EXPECT_EQ( ncols, more_rows.columns() );
    for ( size_t i = 0; i < more_rows.rows(); i++ ) {
        for ( size_t j = 0; j < more_rows.columns(); j++ ) {
            _TEST_EQ_( testarr[ i ][ j ], more_rows.get( i, j ) );
        }
    }
    NCPA::arrays::free_array( testarr, nrows, ncols );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest,
        AsArrayWorksWithPreallocatedPointer ) {
    size_t nrows = more_rows.rows(), ncols = more_rows.columns();
    test_t **testarr = NCPA::arrays::zeros<test_t>( nrows, ncols );
    more_rows.as_array( nrows, ncols, testarr );
    EXPECT_EQ( nrows, more_rows.rows() );
    EXPECT_EQ( ncols, more_rows.columns() );
    for ( size_t i = 0; i < more_rows.rows(); i++ ) {
        for ( size_t j = 0; j < more_rows.columns(); j++ ) {
            _TEST_EQ_( testarr[ i ][ j ], more_rows.get( i, j ) );
        }
    }
    NCPA::arrays::free_array( testarr, nrows, ncols );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, AsArrayThrowsIfDimensionsWrong ) {
    size_t nrows = more_rows.rows(), ncols = more_rows.columns() + 2;
    test_t **testarr = NCPA::arrays::zeros<test_t>( nrows, ncols );
    EXPECT_THROW(
        { more_rows.as_array( nrows, ncols, testarr ); },
        std::invalid_argument );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, SetMethodsWork ) {
    _TEST_EQ_( square.get( 0, 0 ), 1.0 );
    square.set( 0, 0, testval );
    _TEST_EQ_( square.get( 0, 0 ), testval );
    square.set( testval );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), testval );
        }
    }
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, SetRowMethodsWork ) {
    test_t *row0          = NCPA::arrays::zeros<test_t>( dim1 );
    size_t row0_inds[ 3 ] = { 0, 1, 2 };
    NCPA::arrays::fill( row0, dim1, testval );
    size_t row = 0;
    more_rows.set_row( row, dim1, row0_inds, row0 );
    for ( size_t col = 0; col < dim1; col++ ) {
        _TEST_EQ_( more_rows.get( row, col ), testval );
    }

    row           = 1;
    test_t oldval = more_rows.get( row, 1 );
    std::vector<test_t> row1( 2, testval );
    std::vector<size_t> row1_inds( { 0, 2 } );
    more_rows.set_row( row, row1_inds, row1 );
    _TEST_EQ_( more_rows.get( row, 0 ), testval );
    _TEST_EQ_( more_rows.get( row, 1 ), oldval );
    _TEST_EQ_( more_rows.get( row, 2 ), testval );

    row    = 2;
    oldval = more_rows.get( row, 1 );
    more_rows.set_row( row, { 0, 2 }, { testval, testval } );
    _TEST_EQ_( more_rows.get( row, 0 ), testval );
    _TEST_EQ_( more_rows.get( row, 1 ), oldval );
    _TEST_EQ_( more_rows.get( row, 2 ), testval );

    row = 3;
    std::vector<test_t> row3( 2, testval );
    oldval = more_rows.get( row, 2 );
    more_rows.set_row( row, row3 );
    _TEST_EQ_( more_rows.get( row, 0 ), testval );
    _TEST_EQ_( more_rows.get( row, 1 ), testval );
    _TEST_EQ_( more_rows.get( row, 2 ), oldval );

    row = 4;
    more_rows.set_row( row, testval );
    _TEST_EQ_( more_rows.get( row, 0 ), testval );
    _TEST_EQ_( more_rows.get( row, 1 ), testval );
    _TEST_EQ_( more_rows.get( row, 2 ), testval );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, SetColumnsMethodsWork ) {
    test_t *col0          = NCPA::arrays::zeros<test_t>( dim1 );
    size_t col0_inds[ 3 ] = { 0, 1, 2 };
    NCPA::arrays::fill( col0, dim1, testval );
    size_t col = 0;
    more_cols.set_column( col, dim1, col0_inds, col0 );
    for ( size_t row = 0; row < dim1; row++ ) {
        _TEST_EQ_( more_cols.get( row, col ), testval );
    }

    col           = 1;
    test_t oldval = more_cols.get( 1, col );
    std::vector<test_t> col1( 2, testval );
    std::vector<size_t> col1_inds( { 0, 2 } );
    more_cols.set_column( col, col1_inds, col1 );
    _TEST_EQ_( more_cols.get( 0, col ), testval );
    _TEST_EQ_( more_cols.get( 1, col ), oldval );
    _TEST_EQ_( more_cols.get( 2, col ), testval );

    col    = 2;
    oldval = more_cols.get( 1, col );
    more_cols.set_column( col, { 0, 2 }, { testval, testval } );
    _TEST_EQ_( more_cols.get( 0, col ), testval );
    _TEST_EQ_( more_cols.get( 1, col ), oldval );
    _TEST_EQ_( more_cols.get( 2, col ), testval );

    col = 3;
    std::vector<test_t> col3( 2, testval );
    oldval = more_cols.get( 2, col );
    more_cols.set_column( col, col3 );
    _TEST_EQ_( more_cols.get( 0, col ), testval );
    _TEST_EQ_( more_cols.get( 1, col ), testval );
    _TEST_EQ_( more_cols.get( 2, col ), oldval );

    col = 4;
    more_cols.set_column( col, testval );
    _TEST_EQ_( more_cols.get( 0, col ), testval );
    _TEST_EQ_( more_cols.get( 1, col ), testval );
    _TEST_EQ_( more_cols.get( 2, col ), testval );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest, TransposeWorksCorrectly ) {
    EXPECT_TRUE( symmetric.transpose()->equals( symmetric ) );
    EXPECT_TRUE( identity.transpose()->equals( identity ) );
    EXPECT_TRUE( more_rows.transpose()->equals( more_cols ) );
}

TEST_F( NCPALinearAlgebraDenseMatrixTest,
        MatrixMatrixMultiplicationIsCorrect ) {
    EXPECT_TRUE( more_rows.multiply( more_cols )->equals( product ) );
    EXPECT_TRUE( more_cols.multiply( more_rows )->equals( product2 ) );
    EXPECT_TRUE( square.multiply( identity )->equals( square ) );
}
