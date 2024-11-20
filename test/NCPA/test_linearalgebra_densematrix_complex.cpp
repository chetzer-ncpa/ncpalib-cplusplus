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

#define _TEST_EQ_       EXPECT_COMPLEX_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_COMPLEX_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraLibraryComplexDenseMatrixTest

typedef complex<double> test_t;
typedef details::dense_matrix<test_t> mat_t;

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            testval = test_t( -4.2, 4.2 );
            square    = mat_t( dim1, dim1 );
            more_rows = mat_t( dim2, dim1 );
            more_cols = mat_t( dim1, dim2 );
            product   = mat_t( dim2, dim2 );
            product2  = mat_t( dim1, dim1 );
            symmetric = mat_t( dim1, dim1 );
            identity  = mat_t( dim1, dim1 );
            zeromat   = mat_t( dim2, dim2 );
            for ( size_t i = 0; i < dim1; i++ ) {
                double die = (double)( i + 1 );
                test_t di( die, die );
                square.set_row( i, { di, di, di } );
                more_rows.set_column( i, { di, di, di, di, di } );
                more_cols.set_row( i, { di, di, di, di, di } );
                symmetric.set_row( i, { 0, 0, 0 } );
                identity.set( i, i, NCPA::math::one<test_t>() );
                for ( size_t j = 0; j < dim1; j++ ) {
                    product2.set( i, j,
                                  test_t( 0.0,  10.0 * (double)( i + 1 ) * (double)( j + 1 ) ) );
                }
            }
            symmetric.set( 1, 1, testval );
            product.set( test_t( 0.0, 28.0) );

        }  // void TearDown() override {}

        // declare stuff here
        details::dense_matrix<test_t> empty, square, more_rows, more_cols,
            product, identity, symmetric, product2, zeromat;
        Matrix<test_t> mat1, mat2;
        const size_t dim1 = 3, dim2 = 5;
        test_t testval;
};

TEST_F( _TEST_TITLE_, DefaultConstructorIsCorrect ) {
    EXPECT_EQ( empty.rows(), 0 );
    EXPECT_EQ( empty.columns(), 0 );
}

TEST_F( _TEST_TITLE_, SizedConstructorIsCorrect ) {
    EXPECT_EQ( square.rows(), dim1 );
    EXPECT_EQ( square.columns(), dim1 );
}

TEST_F( _TEST_TITLE_, EqualsReturnsTrueForEqual ) {
    EXPECT_TRUE( square.equals( square ) );
}

TEST_F( _TEST_TITLE_, EqualsReturnsFalseForUnequal ) {
    EXPECT_FALSE( square.equals( more_rows ) );
}

TEST_F( _TEST_TITLE_, CopyConstructorWorks ) {
    empty = mat_t( square );
    EXPECT_TRUE( empty.equals( square ) );
    square.clear();
    EXPECT_FALSE( empty.equals( square ) );
}

TEST_F( _TEST_TITLE_, AssignmentOperatorWorks ) {
    empty = square;
    EXPECT_TRUE( empty.equals( square ) );
    square.clear();
    EXPECT_FALSE( empty.equals( square ) );
}

TEST_F( _TEST_TITLE_, SizeMethodsReturnExpectedValues ) {
    EXPECT_EQ( square.rows(), dim1 );
    EXPECT_EQ( square.columns(), dim1 );
    EXPECT_EQ( more_rows.rows(), dim2 );
    EXPECT_EQ( more_rows.columns(), dim1 );
    EXPECT_EQ( more_cols.rows(), dim1 );
    EXPECT_EQ( more_cols.columns(), dim2 );
    EXPECT_EQ( product.rows(), dim2 );
    EXPECT_EQ( product.columns(), dim2 );
}

TEST_F( _TEST_TITLE_, ResizeWorks ) {
    square.resize( 2, 7 );
    EXPECT_EQ( square.rows(), 2 );
    EXPECT_EQ( square.columns(), 7 );
}

TEST_F( _TEST_TITLE_, ResizePreservesExistingData ) {
    square.resize( dim1 + 1, dim1 + 1 );
    for ( size_t i = 0; i < dim1 + 1; i++ ) {
        for ( size_t j = 0; j < dim1 + 1; j++ ) {
            if ( i < dim1 && j < dim1 ) {
                _TEST_EQ_( more_cols.get( i, j ), square.get( i, j ) );
            } else {
                _TEST_EQ_( square.get( i, j ), NCPA::math::zero<test_t>() );
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

TEST_F( _TEST_TITLE_, SwapWorks ) {
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

// TEST_F( _TEST_TITLE_, GetMethodAllowsModification ) {
//     square.get( 0, 0 ) = testval;
//     _TEST_EQ_( square.get( 0, 0 ), testval );
// }

TEST_F( _TEST_TITLE_, GetRowReturnsCopy ) {
    for ( size_t i = 0; i < dim1; i++ ) {
        for ( size_t j = 0; j < dim1; j++ ) {
            _TEST_EQ_( test_t( (double)( i + 1 ), (double)( i + 1 )) , square.get_row( i )->get( j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, GetColumnReturnsCopy ) {
    for ( size_t i = 0; i < dim1; i++ ) {
        for ( size_t j = 0; j < dim1; j++ ) {
            _TEST_EQ_( test_t( (double)( j + 1 ), (double)( j + 1 )), square.get_column( i )->get( j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, GetDiagonalReturnsExpectedValues ) {
    std::vector<test_t> diag0( { test_t(1.0,1.0), test_t(2.0,2.0), test_t(3.0,3.0) } );
    std::vector<test_t>    diag1( {test_t(1.0,1.0), test_t(2.0,2.0)} ),

        diag2( { test_t(1.0,1.0) } ), diag_1( {test_t(2.0,2.0), test_t(3.0,3.0)} ), 
        diag_2( {test_t(3.0,3.0)} );
    _TEST_ARRAY_EQ_( 3, square.get_diagonal()->as_std(), diag0 );
    _TEST_ARRAY_EQ_( 2, square.get_diagonal( 1 )->as_std(), diag1 );
    _TEST_ARRAY_EQ_( 1, square.get_diagonal( 2 )->as_std(), diag2 );
    _TEST_ARRAY_EQ_( 2, square.get_diagonal( -1 )->as_std(), diag_1 );
    _TEST_ARRAY_EQ_( 1, square.get_diagonal( -2 )->as_std(), diag_2 );
}

TEST_F( _TEST_TITLE_, IndexingOperatorReturnsVector ) {
    for ( size_t i = 0; i < dim1; i++ ) {
        for ( size_t j = 0; j < dim1; j++ ) {
            _TEST_EQ_( square[i][j], square.get( i, j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, ClearWorksProperly ) {
    square.clear();
    EXPECT_EQ( square.rows(), 0 );
    EXPECT_EQ( square.columns(), 0 );
}

TEST_F( _TEST_TITLE_, CloneMethodsWork ) {
    EXPECT_EQ( square.clone()->rows(), dim1 );
    EXPECT_EQ( square.clone()->columns(), dim1 );
    EXPECT_EQ( square.fresh_clone()->rows(), 0 );
    EXPECT_EQ( square.fresh_clone()->columns(), 0 );
}

TEST_F( _TEST_TITLE_, AsArrayWorksWithNullPointer ) {
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

TEST_F( _TEST_TITLE_, AsArrayWorksWithPreallocatedPointer ) {
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

TEST_F( _TEST_TITLE_, AsArrayThrowsIfDimensionsWrong ) {
    size_t nrows = more_rows.rows(), ncols = more_rows.columns() + 2;
    test_t **testarr = NCPA::arrays::zeros<test_t>( nrows, ncols );
    EXPECT_THROW(
        { more_rows.as_array( nrows, ncols, testarr ); },
        std::invalid_argument );
}

TEST_F( _TEST_TITLE_, SetMethodsWork ) {
    _TEST_EQ_( square.get( 0, 0 ), test_t(1.0,1.0) );
    square.set( 0, 0, testval );
    _TEST_EQ_( square.get( 0, 0 ), testval );
    square.set( testval );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), testval );
        }
    }
}

TEST_F( _TEST_TITLE_, SetRowMethodsWork ) {
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

TEST_F( _TEST_TITLE_, SetColumnsMethodsWork ) {
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

TEST_F( _TEST_TITLE_, TransposeWorksCorrectly ) {
    mat_t tsym = symmetric, tid = identity, tmr = more_rows;
    EXPECT_TRUE( tsym.transpose().equals( symmetric ) );
    EXPECT_TRUE( tid.transpose().equals( identity ) );
    EXPECT_TRUE( tmr.transpose().equals( more_cols ) );
}

TEST_F( _TEST_TITLE_, MatrixMatrixMultiplicationIsCorrect ) {
    EXPECT_TRUE( more_rows.multiply( more_cols )->equals( product ) );
    EXPECT_TRUE( more_cols.multiply( more_rows )->equals( product2 ) );
    EXPECT_TRUE( square.multiply( identity )->equals( square ) );
}

TEST_F( _TEST_TITLE_, ScaleWorksWithScalar ) {
    mat_t original = product;
    product.scale( 0.5 );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( product.get( i, j ), 0.5 * original.get( i, j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, ScaleWorksWithMatrix ) {
    mat_t original = product;
    product.scale( product );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( product.get( i, j ),
                       (original.get( i, j ) * original.get( i, j )) );
        }
    }
}

TEST_F( _TEST_TITLE_, AddWorksWithScalar ) {
    mat_t original = square;
    square.add( testval );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ),
                       (original.get( i, j ) + testval) );
        }
    }
}

TEST_F( _TEST_TITLE_, AddWorksWithMatrix ) {
    mat_t original = square;
    square.add( square );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ),
                       (original.get( i, j ) * 2.0) );
        }
    }
}

TEST_F( _TEST_TITLE_, AddWorksWithScaledMatrix ) {
    mat_t original = square;
    square.add( square, -1.0 );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ),
                       NCPA::math::zero<test_t>() );
        }
    }
}

TEST_F( _TEST_TITLE_, ZeroWorks ) {
    square.set( testval );
    _TEST_EQ_( square.get( 0, 0 ), testval );
    square.zero( 0, 0 );
    _TEST_EQ_( square.get( 0, 0 ), NCPA::math::zero<test_t>() );
}

TEST_F( _TEST_TITLE_, SetDiagonalWorks ) {
    mat_t testmat( 5, 5 ), control( 5, 5 );
    test_t zero = NCPA::math::zero<test_t>();
    
    for (auto i = 0; i < 5; i++) {
        _TEST_EQ_( testmat.get(i,i), zero );
        control.set( i, i, testval );
    }
    // set control off-diagonals
    control.set( 0, 2, testval );
    control.set( 1, 3, testval );
    control.set( 2, 4, testval );
    control.set( 3, 0, testval );
    control.set( 4, 1, testval );

    // first version
    test_t diag[ 5 ] = { testval, testval, testval, testval, testval };
    vector<test_t> vdiag( { testval, testval, testval, testval, testval } );
    testmat.set_diagonal( 5, diag );
    testmat.set_diagonal( 3, diag, 2 );
    testmat.set_diagonal( 2, diag, -3 );
    EXPECT_TRUE( testmat.equals( control ) );

    // second version
    testmat.zero();
    testmat.set_diagonal( vdiag );
    testmat.set_diagonal( vdiag, 2 );
    testmat.set_diagonal( vdiag, -3 );
    EXPECT_TRUE( testmat.equals( control ) );

    // third version
    testmat.zero();
    testmat.set_diagonal( { testval, testval, testval, testval, testval } );
    testmat.set_diagonal( { testval, testval, testval, testval, testval }, 2 );
    testmat.set_diagonal( { testval, testval, testval, testval, testval }, -3 );
    EXPECT_TRUE( testmat.equals( control ) );

    // fourth version
    testmat.zero();
    testmat.set_diagonal( testval );
    testmat.set_diagonal( testval, 2 );
    testmat.set_diagonal( testval, -3 );
    EXPECT_TRUE( testmat.equals( control ) );
}

TEST_F( _TEST_TITLE_, PlusEqualOperatorWorksWithScalar ) {
    mat_t original = square;
    square += testval;
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ),
                       (original.get( i, j ) + testval) );
        }
    }
}

TEST_F( _TEST_TITLE_, PlusEqualOperatorWorksWithMatrix ) {
    mat_t original = square;
    square += square;
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ),
                       (original.get( i, j ) * 2.0) );
        }
    }
}

TEST_F( _TEST_TITLE_, MinusEqualOperatorWorksWithScalar ) {
    mat_t original = square;
    square -= testval;
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ),
                       (original.get( i, j ) - testval) );
        }
    }
}

TEST_F( _TEST_TITLE_, MinusEqualWorksWithMatrix ) {
    mat_t original = square;
    square -= square;
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ),
                       NCPA::math::zero<test_t>() );
        }
    }
}

TEST_F( _TEST_TITLE_, TimesEqualOperatorWorksWithScalar ) {
    mat_t original = product;
    product*= 0.5;
    for ( size_t i = 0; i < product.rows(); i++ ) {
        for ( size_t j = 0; j < product.columns(); j++ ) {
            _TEST_EQ_( product.get( i, j ), 0.5 * original.get( i, j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, TimesEqualOperatorWorksWithMatrix ) {
    mat_t original = product;
    product *= product;
    for ( size_t i = 0; i < product.rows(); i++ ) {
        for ( size_t j = 0; j < product.columns(); j++ ) {
            _TEST_EQ_( product.get( i, j ),
                       (original.get( i, j ) * original.get( i, j )) );
        }
    }
}

TEST_F( _TEST_TITLE_, DivideEqualOperatorWorksWithScalar ) {
    mat_t original = product;
    product /= testval;
    for ( size_t i = 0; i < product.rows(); i++ ) {
        for ( size_t j = 0; j < product.columns(); j++ ) {
            _TEST_EQ_( product.get( i, j ),
                       (original.get( i, j ) / testval) );
        }
    }
}

TEST_F( _TEST_TITLE_, IsDiagonalReturnsCorrectly ) {
    EXPECT_TRUE( identity.is_diagonal() );
    EXPECT_FALSE( square.is_diagonal() );
    EXPECT_FALSE( product.is_diagonal() );
    EXPECT_FALSE( more_rows.is_diagonal() );
    EXPECT_FALSE( more_cols.is_diagonal() );
    EXPECT_TRUE( empty.is_diagonal() );
    EXPECT_TRUE( zeromat.is_diagonal() );
} 

TEST_F( _TEST_TITLE_, IsTridiagonalReturnsCorrectly ) {
    EXPECT_TRUE( identity.is_tridiagonal() );
    EXPECT_FALSE( square.is_tridiagonal() );
    EXPECT_FALSE( product.is_tridiagonal() );
    EXPECT_FALSE( more_rows.is_tridiagonal() );
    EXPECT_FALSE( more_cols.is_tridiagonal() );
    EXPECT_TRUE( empty.is_tridiagonal() );
    EXPECT_TRUE( zeromat.is_tridiagonal() );

    square.zero( 0, dim1-1 );
    square.zero( dim1-1, 0 );
    EXPECT_TRUE( square.is_tridiagonal() );
} 

TEST_F( _TEST_TITLE_, IsLowerTriagonalIsCorrect ) {
    EXPECT_TRUE( identity.is_lower_triangular() );
    EXPECT_FALSE( square.is_lower_triangular() );
    EXPECT_FALSE( product.is_lower_triangular() );
    EXPECT_FALSE( more_rows.is_lower_triangular() );
    EXPECT_FALSE( more_cols.is_lower_triangular() );
    EXPECT_TRUE( empty.is_lower_triangular() );
    EXPECT_TRUE( zeromat.is_lower_triangular() );

    for (size_t i = 0; i < zeromat.rows(); i++) {
        for (size_t j = 0; j < i; j++) {
            zeromat.set(i,j,testval);
        }
    }
    EXPECT_TRUE( zeromat.is_lower_triangular() );
    zeromat.set( 0,1,testval);
    EXPECT_FALSE( zeromat.is_lower_triangular() );
}