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
typedef dense_vector<test_t> vec_t;

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here

            square = mat_t( dim1, dim1 );
            square.set_diagonal( diagval );
            square.set_diagonal( supdiagval, 1 );
            square.set_diagonal( supdiagval, -1 );

            more_rows = mat_t( dim2, dim1 );
            more_rows.set_diagonal( diagval )
                .set_diagonal( supdiagval, 1 )
                .set_diagonal( supdiagval, -1 );

            more_cols = mat_t( dim1, dim2 );
            more_cols.set_diagonal( diagval )
                .set_diagonal( supdiagval, 1 )
                .set_diagonal( supdiagval, -1 );

            product  = mat_t( dim2, dim2 );
            product2 = mat_t( dim1, dim1 );

            symmetric = mat_t( dim1, dim1 );
            symmetric.set_diagonal( diagval )
                .set_diagonal( supdiagval, 1 )
                .set_diagonal( supdiagval, -1 );
            symmetric.set( 1, 1, test_t( 42.0, 0.0 ) );

            identity.identity( dim1, dim1 );

            zeromat = mat_t( dim2, dim2 );

            testvec = vec_t( dim1 );
            testvec.set( diagval );
            rightvec = vec_t( dim1 );
            rightvec.set( 0, test_t(0,-4) );
            rightvec.set( 4, test_t(0,-4) );
            leftvec = rightvec;

        }  // void TearDown() override {}

        // declare stuff here
        mat_t empty, square, more_rows, more_cols, product, identity,
            symmetric, product2, zeromat;
        // Matrix<test_t> mat1, mat2;
        const size_t dim1 = 5, dim2 = 8;
        test_t testval    = test_t( -4.2, 2.1 );
        const test_t zero = NCPA::math::zero<test_t>();
        test_t diagval = test_t( -2.0, 2.0 ), supdiagval = test_t( 1.0, -1.0 );
        vec_t testvec, leftvec, rightvec;
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

TEST_F( _TEST_TITLE_, InternalCoordinateTransformationsWork ) {
    int ind1, ind2;
    for ( int row = 0; row < dim1; row++ ) {
        for ( int col = 0; col < dim1; col++ ) {
            if ( abs( row - col ) > 1 ) {
                EXPECT_FALSE( square.rowcol2internal( row, col, ind1, ind2 ) );
            } else {
                EXPECT_TRUE( square.rowcol2internal( row, col, ind1, ind2 ) );
            }
            EXPECT_EQ( ind1, 1 + col - row );
            EXPECT_EQ( ind2, row );
        }
    }
}

TEST_F( _TEST_TITLE_, ResizeWorks ) {
    square.resize( 2, 7 );
    EXPECT_EQ( square.rows(), 2 );
    EXPECT_EQ( square.columns(), 7 );
}

TEST_F( _TEST_TITLE_, ResizePreservesExistingData ) {
    square.resize( dim1 + 1, dim1 + 1 );
    ASSERT_EQ( square.rows(), dim1 + 1 );
    ASSERT_EQ( square.columns(), dim1 + 1 );
    for ( size_t i = 0; i < dim1 + 1; i++ ) {
        for ( size_t j = 0; j < dim1 + 1; j++ ) {
            if ( i < dim1 && j < dim1 ) {
                _TEST_EQ_( more_cols.get( i, j ), square.get( i, j ) );
            } else {
                _TEST_EQ_( square.get( i, j ), zero );
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

TEST_F( _TEST_TITLE_, GetRowReturnsCopy ) {
    for ( size_t i = 0; i < dim1; i++ ) {
        for ( size_t j = 0; j < dim1; j++ ) {
            _TEST_EQ_( square.get_row( i )->get( j ), square.get( i, j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, GetColumnReturnsCopy ) {
    for ( size_t i = 0; i < dim1; i++ ) {
        for ( size_t j = 0; j < dim1; j++ ) {
            _TEST_EQ_( square.get_column( i )->get( j ), square.get( j, i ) );
        }
    }
}

TEST_F( _TEST_TITLE_, GetDiagonalReturnsExpectedValues ) {
    std::vector<test_t> diag0(
        { diagval, diagval, diagval, diagval, diagval } ),
        diag1( { supdiagval, supdiagval, supdiagval, supdiagval } ),
        diag2( { zero, zero, zero } ),
        diag_1( { supdiagval, supdiagval, supdiagval, supdiagval } ),
        diag_2( { zero, zero, zero } );
    _TEST_ARRAY_EQ_( dim1, square.get_diagonal()->as_std(), diag0 );
    _TEST_ARRAY_EQ_( dim1 - 1, square.get_diagonal( 1 )->as_std(), diag1 );
    _TEST_ARRAY_EQ_( dim1 - 2, square.get_diagonal( 2 )->as_std(), diag2 );
    _TEST_ARRAY_EQ_( dim1 - 1, square.get_diagonal( -1 )->as_std(), diag_1 );
    _TEST_ARRAY_EQ_( dim1 - 2, square.get_diagonal( -2 )->as_std(), diag_2 );
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
    _TEST_EQ_( square.get( 0, 0 ), diagval );
    square.set( 0, 0, testval );
    _TEST_EQ_( square.get( 0, 0 ), testval );
    square.set( testval );
    for ( int i = 0; i < square.rows(); i++ ) {
        for ( int j = 0; j < square.columns(); j++ ) {
            if ( abs( i - j ) > 1 ) {
                _TEST_EQ_( square.get( i, j ), zero );
            } else {
                _TEST_EQ_( square.get( i, j ), testval );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, ArraySetRowMethodWorks ) {
    test_t *row0          = NCPA::arrays::zeros<test_t>( dim1 );
    size_t row0_inds[ 5 ] = { 0, 1, 2, 3, 4 };
    NCPA::arrays::fill( row0, dim1, testval );
    size_t row = 0;
    more_rows.set_row( row, 2, row0_inds, row0 );
    for ( size_t col = 0; col < 2; col++ ) {
        _TEST_EQ_( more_rows.get( row, col ), testval );
    }
}

TEST_F( _TEST_TITLE_, VectorSetRowMethodWorks ) {
    size_t row           = 1;
    test_t oldval = more_rows.get( row, 1 );
    std::vector<test_t> row1( 2, testval );
    std::vector<size_t> row1_inds( { 0, 2 } );
    more_rows.set_row( row, row1_inds, row1 );
    _TEST_EQ_( more_rows.get( row, 0 ), testval );
    _TEST_EQ_( more_rows.get( row, 1 ), oldval );
    _TEST_EQ_( more_rows.get( row, 2 ), testval );
}

TEST_F( _TEST_TITLE_, InitListSetRowMethodWorks ) {
    size_t row    = 2;
    test_t oldval = more_rows.get( row, 2 );
    more_rows.set_row( row, { 1, 3 }, { testval, testval } );
    _TEST_EQ_( more_rows.get( row, 1 ), testval );
    _TEST_EQ_( more_rows.get( row, 2 ), oldval );
    _TEST_EQ_( more_rows.get( row, 3 ), testval );
}

TEST_F( _TEST_TITLE_, BareVectorSetRowMethodWorks ) {
    size_t row = 3;
    std::vector<test_t> row3( 3, testval );
    test_t oldval = more_rows.get( row, 4 );
    more_rows.set_row( row, row3 );
    _TEST_EQ_( more_rows.get( row, 2 ), testval );
    _TEST_EQ_( more_rows.get( row, 3 ), testval );
    _TEST_EQ_( more_rows.get( row, 4 ), testval );
}

TEST_F( _TEST_TITLE_, ConstantSetRowMethodWorks ) {
    size_t row = 4;
    more_rows.set_row( row, testval );
    _TEST_EQ_( more_rows.get( row, 3 ), testval );
    _TEST_EQ_( more_rows.get( row, 4 ), testval );
}

TEST_F( _TEST_TITLE_, ArraySetColumnsMethodWorks ) {
    test_t *col0          = NCPA::arrays::zeros<test_t>( dim1 );
    size_t col0_inds[ 5 ] = { 0, 1, 2, 3, 4 };
    NCPA::arrays::fill( col0, dim1, testval );
    size_t col = 0;
    more_cols.set_column( col, 2, col0_inds, col0 );
    for ( size_t row = 0; row < 2; row++ ) {
        _TEST_EQ_( more_cols.get( row, col ), testval );
    }
}
TEST_F( _TEST_TITLE_, VectorSetColumnsMethodWorks ) {
    size_t col           = 1;
    test_t oldval = more_cols.get( 1, col );
    std::vector<test_t> col1( 2, testval );
    std::vector<size_t> col1_inds( { 0, 2 } );
    more_cols.set_column( col, col1_inds, col1 );
    _TEST_EQ_( more_cols.get( 0, col ), testval );
    _TEST_EQ_( more_cols.get( 1, col ), oldval );
    _TEST_EQ_( more_cols.get( 2, col ), testval );
}
TEST_F( _TEST_TITLE_, InitListSetColumnsMethodWorks ) {
    size_t col    = 2;
    test_t oldval = more_cols.get( 2, col );
    more_cols.set_column( col, { 1, 3 }, { testval, testval } );
    _TEST_EQ_( more_cols.get( 1, col ), testval );
    _TEST_EQ_( more_cols.get( 2, col ), oldval );
    _TEST_EQ_( more_cols.get( 3, col ), testval );
}
TEST_F( _TEST_TITLE_, BareVectorSetColumnsMethodWorks ) {
    size_t col = 3;
    std::vector<test_t> col3( 3, testval );
    more_cols.set_column( col, col3 );
    _TEST_EQ_( more_cols.get( 2, col ), testval );
    _TEST_EQ_( more_cols.get( 3, col ), testval );
    _TEST_EQ_( more_cols.get( 4, col ), testval );
}

TEST_F( _TEST_TITLE_, ConstantSetColumnsMethodWorks ) {
    size_t col = 4;
    more_cols.set_column( col, testval );
    _TEST_EQ_( more_cols.get( 3, col ), testval );
    _TEST_EQ_( more_cols.get( 4, col ), testval );
}


TEST_F( _TEST_TITLE_, TransposeWorksCorrectly ) {
    mat_t tsym = symmetric, tid = identity, tmr = more_rows;
    EXPECT_TRUE( tsym.transpose().equals( symmetric ) );
    EXPECT_TRUE( tid.transpose().equals( identity ) );
    EXPECT_TRUE( tmr.transpose().equals( more_cols ) );
}

TEST_F( _TEST_TITLE_, MatrixVectorRightMultiplicationIsCorrect ) {
    auto product = square.right_multiply( testvec );
    // NCPA_DEBUG << "[ ";
    // for (size_t i = 0; i < dim1; i++) {
    //     if (i != 0) {
    //         NCPA_DEBUG << "  ";
    //     }
    //     NCPA_DEBUG << "[ ";
    //     for (size_t j = 0; j < dim1; j++) {
    //         if (j != 0) {
    //             NCPA_DEBUG << ", ";
    //         }
    //         NCPA_DEBUG << square.get(i,j);
    //     }
    //     NCPA_DEBUG << "][ " << testvec[ i ] << "] = [" << product->get(i) << "]" << endl;
    // }

    for (size_t i = 0; i < dim1; i++) {
        EXPECT_NEAR( product->get(i).real(), rightvec.get(i).real(), 1e-10 );
        EXPECT_NEAR( product->get(i).imag(), rightvec.get(i).imag(), 1e-10 );
    }
}

TEST_F( _TEST_TITLE_, VectorMatrixLeftMultiplicationIsCorrect ) {
    auto product = square.left_multiply( testvec );
    for (size_t i = 0; i < dim1; i++) {
        EXPECT_NEAR( product->get(i).real(), leftvec.get(i).real(), 1e-10 );
        EXPECT_NEAR( product->get(i).imag(), leftvec.get(i).imag(), 1e-10 );
    }
}

TEST_F( _TEST_TITLE_, MatrixMatrixMultiplicationIsCorrect ) {
    EXPECT_TRUE( square.multiply( identity )->equals( square ) );
    mat_t lhs = square;
    lhs.zero()
        .set_diagonal( { 1, 2, 3, 4, 5 } )
        .set_diagonal( { 6, 7, 8, 9 }, 1 );
    mat_t rhs = lhs;
    lhs.set_diagonal( { 10, 11, 12, 13 }, -1 );
    // cout << "LHS: " << endl << lhs << "RHS: " << endl << rhs;
    product2.zero()
        .set_diagonal( { 1, 64, 86, 112, 142 } )
        .set_diagonal( { 18, 35, 56, 81 }, 1 )
        .set_diagonal( { 42, 56, 72 }, 2 )
        .set_diagonal( { 10, 22, 36, 52 }, -1 );
    // cout << "Expect product = " << endl << product2;
    auto product = lhs.multiply( rhs );
    // cout << "Actual product = " << endl << *product;
    EXPECT_TRUE( product2.equals( *product ) );
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
                       ( original.get( i, j ) * original.get( i, j ) ) );
        }
    }
}

TEST_F( _TEST_TITLE_, AddWorksWithScalar ) {
    mat_t original = square;
    square.add( testval );
    for ( int i = 0; i < square.rows(); i++ ) {
        for ( int j = 0; j < square.columns(); j++ ) {
            if ( abs( i - j ) > 1 ) {
                _TEST_EQ_( square.get( i, j ), zero );
            } else {
                _TEST_EQ_( square.get( i, j ),
                           ( original.get( i, j ) + testval ) );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, AddWorksWithMatrix ) {
    mat_t original = square;
    // cout << "Before add:" << endl << square;
    square.add( square );
    // cout << "After add:" << endl << square;
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), ( original.get( i, j ) * 2.0 ) );
        }
    }
}

TEST_F( _TEST_TITLE_, AddWorksWithScaledMatrix ) {
    mat_t original = square;
    square.add( square, -1.0 );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), zero );
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
    for ( auto i = 0; i < 5; i++ ) {
        _TEST_EQ_( testmat.get( i, i ), zero );
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
    // cout << "First test passed" << endl;

    // second version
    testmat.zero();
    // cout << "Zero OK" << endl;
    testmat.set_diagonal( vdiag );
    // cout << "Main diagonal set" << endl;
    testmat.set_diagonal( vdiag, 2 );
    // cout << "Upper second off-diagonal set" << endl;
    testmat.set_diagonal( vdiag, -3 );
    // cout << "Lower third off-diagonal set" << endl;
    EXPECT_TRUE( testmat.equals( control ) );
    // cout << "Second test passed" << endl;

    // third version
    testmat.zero();
    testmat.set_diagonal( { testval, testval, testval, testval, testval } );
    testmat.set_diagonal( { testval, testval, testval, testval, testval }, 2 );
    testmat.set_diagonal( { testval, testval, testval, testval, testval },
                          -3 );
    EXPECT_TRUE( testmat.equals( control ) );

    // fourth version
    testmat.zero();
    testmat.set_diagonal( testval );
    testmat.set_diagonal( testval, 2 );
    testmat.set_diagonal( testval, -3 );
    EXPECT_TRUE( testmat.equals( control ) );
}

TEST_F( _TEST_TITLE_, PlusEqualOperatorWorksWithScalar ) {
    mat_t original  = square;
    square         += testval;
    for ( int i = 0; i < square.rows(); i++ ) {
        for ( int j = 0; j < square.columns(); j++ ) {
            if ( abs( i - j ) > 1 ) {
                _TEST_EQ_( square.get( i, j ), zero );
            } else {
                _TEST_EQ_( square.get( i, j ),
                           ( original.get( i, j ) + testval ) );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, PlusEqualOperatorWorksWithMatrix ) {
    mat_t original  = square;
    square         += square;
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), ( original.get( i, j ) * 2.0 ) );
        }
    }
}

TEST_F( _TEST_TITLE_, MinusEqualOperatorWorksWithScalar ) {
    mat_t original  = square;
    square         -= testval;
    for ( int i = 0; i < square.rows(); i++ ) {
        for ( int j = 0; j < square.columns(); j++ ) {
            if ( abs( i - j ) > 1 ) {
                _TEST_EQ_( square.get( i, j ), zero );
            } else {
                _TEST_EQ_( square.get( i, j ),
                           ( original.get( i, j ) - testval ) );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, MinusEqualWorksWithMatrix ) {
    mat_t original  = square;
    square         -= square;
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), zero );
        }
    }
}

TEST_F( _TEST_TITLE_, TimesEqualOperatorWorksWithScalar ) {
    mat_t original  = product;
    product        *= 0.5;
    for ( size_t i = 0; i < product.rows(); i++ ) {
        for ( size_t j = 0; j < product.columns(); j++ ) {
            _TEST_EQ_( product.get( i, j ), 0.5 * original.get( i, j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, TimesEqualOperatorWorksWithMatrix ) {
    mat_t original  = product;
    product        *= product;
    for ( size_t i = 0; i < product.rows(); i++ ) {
        for ( size_t j = 0; j < product.columns(); j++ ) {
            _TEST_EQ_( product.get( i, j ),
                       ( original.get( i, j ) * original.get( i, j ) ) );
        }
    }
}

TEST_F( _TEST_TITLE_, DivideEqualOperatorWorksWithScalar ) {
    mat_t original  = product;
    product        /= testval;
    for ( size_t i = 0; i < product.rows(); i++ ) {
        for ( size_t j = 0; j < product.columns(); j++ ) {
            _TEST_EQ_( product.get( i, j ),
                       ( original.get( i, j ) / testval ) );
        }
    }
}

TEST_F( _TEST_TITLE_, IsDiagonalReturnsCorrectly ) {
    EXPECT_TRUE( identity.is_diagonal() );
    EXPECT_FALSE( square.is_diagonal() );
    // EXPECT_FALSE( product.is_diagonal() );
    EXPECT_FALSE( more_rows.is_diagonal() );
    EXPECT_FALSE( more_cols.is_diagonal() );
    EXPECT_TRUE( empty.is_diagonal() );
    EXPECT_TRUE( zeromat.is_diagonal() );
}

TEST_F( _TEST_TITLE_, IsTridiagonalReturnsCorrectly ) {
    EXPECT_TRUE( identity.is_tridiagonal() );
    EXPECT_TRUE( square.is_tridiagonal() );
    // EXPECT_FALSE( product.is_tridiagonal() );
    EXPECT_TRUE( more_rows.is_tridiagonal() );
    EXPECT_TRUE( more_cols.is_tridiagonal() );
    EXPECT_TRUE( empty.is_tridiagonal() );
    EXPECT_TRUE( zeromat.is_tridiagonal() );

    ASSERT_EQ( square.get_row( 0 )->count_nonzero_indices(), 2 );
    square.zero( 0, 1 );
    ASSERT_EQ( square.get_row( 0 )->count_nonzero_indices(), 1 );
    ASSERT_EQ( square.get_row( 2 )->count_nonzero_indices(), 3 );
    square.zero( 2, 3 );
    ASSERT_EQ( square.get_row( 2 )->count_nonzero_indices(), 2 );

    EXPECT_TRUE( square.is_tridiagonal() );
}

TEST_F( _TEST_TITLE_, IsLowerTriangularIsCorrect ) {
    EXPECT_TRUE( identity.is_lower_triangular() );
    EXPECT_FALSE( square.is_lower_triangular() );
    // EXPECT_FALSE( product.is_lower_triangular() );
    EXPECT_FALSE( more_rows.is_lower_triangular() );
    EXPECT_FALSE( more_cols.is_lower_triangular() );
    EXPECT_TRUE( empty.is_lower_triangular() );
    EXPECT_TRUE( zeromat.is_lower_triangular() );

    for ( size_t i = 1; i < zeromat.rows(); i++ ) {
        zeromat.set( i, i - 1, testval );
    }
    EXPECT_TRUE( zeromat.is_lower_triangular() );
    zeromat.set( 0, 1, testval );
    EXPECT_FALSE( zeromat.is_lower_triangular() );
}

TEST_F( _TEST_TITLE_, IsUpperTriangularIsCorrect ) {
    EXPECT_TRUE( identity.is_upper_triangular() );
    EXPECT_FALSE( square.is_upper_triangular() );
    // EXPECT_FALSE( product.is_upper_triangular() );
    EXPECT_FALSE( more_rows.is_upper_triangular() );
    EXPECT_FALSE( more_cols.is_upper_triangular() );
    EXPECT_TRUE( empty.is_upper_triangular() );
    EXPECT_TRUE( zeromat.is_upper_triangular() );

    for ( size_t i = 1; i < zeromat.rows(); i++ ) {
        zeromat.set( i - 1, i, testval );
    }
    EXPECT_TRUE( zeromat.is_upper_triangular() );
    zeromat.set( 1, 0, testval );
    EXPECT_FALSE( zeromat.is_upper_triangular() );
}
