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

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraLibrarySymmetricMatrixTest

typedef double test_t;
typedef symmetric_matrix<test_t> mat_t;
typedef dense_vector<test_t> vec_t;

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            int ok = 0;
            square = mat_t( dim1, dim1 );
            square.set_diagonal( -2.0 );
            square.set_diagonal( 1.0, 1 );

            product  = mat_t( dim1, dim1 );
            product2 = mat_t( dim1, dim1 );

            symmetric = mat_t( dim1, dim1 );
            symmetric.set_diagonal( -2.0 ).set_diagonal( 1.0, 1 );
            symmetric.set( 1, 1, 42.0 );

            identity.identity( dim1, dim1 );
            zeromat = mat_t( dim1, dim1 );


            testvec = vec_t( dim1 );
            testvec.set( testval );

            rightvec = vec_t( dim1 );
            rightvec.set( 0, 4.2 ).set( dim1 - 1, 4.2 );

            leftvec = rightvec;


        }  // void TearDown() override {}

        // declare stuff here
        mat_t empty, square, product, identity, symmetric, product2, zeromat;
        // Matrix<test_t> mat1, mat2;
        const size_t dim1 = 5, dim2 = 5;
        test_t testval    = -4.2;
        const test_t zero = NCPA::math::zero<test_t>();
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
    EXPECT_FALSE( square.equals( symmetric ) );
}

TEST_F( _TEST_TITLE_, CopyConstructorWorks ) {
    empty = mat_t( square );
    EXPECT_TRUE( empty.equals( square ) );
    square.clear();
    EXPECT_FALSE( empty.equals( square ) );
}

TEST_F( _TEST_TITLE_, CopyOfDenseMatrixIsBandDiagonal ) {
    dense_matrix<test_t> dmat( 5, 5 );
    dmat.set_diagonal( 1.0 ).set_diagonal( -1.0, 1 ).set_diagonal( -1.0, -1 );
    band_diagonal_matrix<test_t> bdmat( dmat );
    EXPECT_TRUE( bdmat.is_band_diagonal() );
    EXPECT_TRUE( bdmat.get_diagonal()->equals( *dmat.get_diagonal() ) );
    EXPECT_TRUE( bdmat.get_diagonal( 1 )->equals( *dmat.get_diagonal( 1 ) ) );
    EXPECT_TRUE(
        bdmat.get_diagonal( -1 )->equals( *dmat.get_diagonal( -1 ) ) );
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
    EXPECT_EQ( product.rows(), dim2 );
    EXPECT_EQ( product.columns(), dim2 );
}

TEST_F( _TEST_TITLE_, InternalCoordinateTransformationsWork ) {
    int ind1, ind2;
    for ( int row = 0; row < dim1; row++ ) {
        for ( int col = 0; col < dim1; col++ ) {
            if ( abs( row - col ) > 1 ) {
                EXPECT_FALSE( square.rowcol2internal( row, col, ind1, ind2 ) );
                // cout << "[" << row << "," << col << "] -> [" << ind1 << ","
                // << ind2 << "]" << endl;
            } else {
                if ( col < row ) {
                    EXPECT_TRUE(
                        square.rowcol2internal( col, row, ind1, ind2 ) );
                    // cout << "[" << row << "," << col << "] -> [" << ind1 <<
                    // "," << ind2 << "]" << endl;
                } else {
                    EXPECT_TRUE(
                        square.rowcol2internal( row, col, ind1, ind2 ) );
                    // cout << "[" << row << "," << col << "] -> [" << ind1 <<
                    // "," << ind2 << "]" << endl;
                }
            }
            EXPECT_EQ( ind1, abs( col - row ) );
            // EXPECT_EQ( ind2, row );
        }
    }
}

TEST_F( _TEST_TITLE_, ResizeWorks ) {
    square.resize( 7, 7 );
    EXPECT_EQ( square.rows(), 7 );
    EXPECT_EQ( square.columns(), 7 );
}

TEST_F( _TEST_TITLE_, ResizePreservesExistingData ) {
    mat_t orig = square;
    square.resize( dim1 + 1, dim1 + 1 );
    ASSERT_EQ( square.rows(), dim1 + 1 );
    ASSERT_EQ( square.columns(), dim1 + 1 );
    ASSERT_EQ( orig.rows(), dim1 );
    ASSERT_EQ( orig.columns(), dim1 );
    for ( size_t i = 0; i < dim1 + 1; i++ ) {
        for ( size_t j = 0; j < dim1 + 1; j++ ) {
            if ( i < dim1 && j < dim1 ) {
                _TEST_EQ_( orig.get( i, j ), square.get( i, j ) );
            } else {
                _TEST_EQ_( square.get( i, j ), 0.0 );
            }
        }
    }

    square.resize( dim1 - 1, dim1 - 1 );
    for ( size_t i = 0; i < dim1 - 1; i++ ) {
        for ( size_t j = 0; j < dim1 - 1; j++ ) {
            _TEST_EQ_( orig.get( i, j ), square.get( i, j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, SwapWorks ) {
    ASSERT_FALSE( square.is_identity() );
    ASSERT_TRUE( identity.is_identity() );
    swap( square, identity );
    ASSERT_TRUE( square.is_identity() );
    ASSERT_FALSE( identity.is_identity() );
}

TEST_F( _TEST_TITLE_, BandWidthReturnsCorrectly ) {
    EXPECT_EQ( square.bandwidth(), 3 );
    square.set_diagonal( testval, 2 );
    EXPECT_EQ( square.bandwidth(), 5 );
}

TEST_F( _TEST_TITLE_, BandRowIndicesReturnsCorrectly ) {
    std::vector<size_t> inds { 0, 1 };
    EXPECT_EQ( inds.size(), square.band_row_indices( 0 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_row_indices( 0 ) );

    inds = { 0, 1, 2 };
    EXPECT_EQ( inds.size(), square.band_row_indices( 1 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_row_indices( 1 ) );

    inds = { 1, 2, 3 };
    EXPECT_EQ( inds.size(), square.band_row_indices( 2 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_row_indices( 2 ) );

    square.set_diagonal( testval, 2 );
    inds = { 0, 1, 2, 3, 4 };
    EXPECT_EQ( inds.size(), square.band_row_indices( 2 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_row_indices( 2 ) );

    inds = { 1, 2, 3, 4 };
    EXPECT_EQ( inds.size(), square.band_row_indices( 3 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_row_indices( 3 ) );

    inds = { 2, 3, 4 };
    EXPECT_EQ( inds.size(), square.band_row_indices( 4 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_row_indices( 4 ) );
}

TEST_F( _TEST_TITLE_, BandColumnIndicesReturnsCorrectly ) {
    std::vector<size_t> inds { 0, 1 };
    EXPECT_EQ( inds.size(), square.band_column_indices( 0 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_column_indices( 0 ) );

    inds = { 0, 1, 2 };
    EXPECT_EQ( inds.size(), square.band_column_indices( 1 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_column_indices( 1 ) );

    inds = { 1, 2, 3 };
    EXPECT_EQ( inds.size(), square.band_column_indices( 2 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_column_indices( 2 ) );

    square.set_diagonal( testval, 2 );
    inds = { 0, 1, 2, 3, 4 };
    EXPECT_EQ( inds.size(), square.band_column_indices( 2 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_column_indices( 2 ) );

    inds = { 1, 2, 3, 4 };
    EXPECT_EQ( inds.size(), square.band_column_indices( 3 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_column_indices( 3 ) );

    inds = { 2, 3, 4 };
    EXPECT_EQ( inds.size(), square.band_column_indices( 4 ).size() );
    EXPECT_ARRAY_EQ( inds.size(), inds, square.band_column_indices( 4 ) );
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
    std::vector<test_t> diag0( { -2, -2, -2, -2, -2 } ),
        diag1( { 1.0, 1.0, 1.0, 1.0 } ), diag2( { 0, 0, 0 } ),
        diag_1( { 1.0, 1.0, 1.0, 1.0 } ), diag_2( { 0, 0, 0 } );
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
    square.as_array( nrows, ncols, testarr );
    EXPECT_EQ( nrows, square.rows() );
    EXPECT_EQ( ncols, square.columns() );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( testarr[ i ][ j ], square.get( i, j ) );
        }
    }
    NCPA::arrays::free_array( testarr, nrows, ncols );
}

TEST_F( _TEST_TITLE_, AsArrayWorksWithPreallocatedPointer ) {
    size_t nrows = square.rows(), ncols = square.columns();
    test_t **testarr = NCPA::arrays::zeros<test_t>( nrows, ncols );
    square.as_array( nrows, ncols, testarr );
    EXPECT_EQ( nrows, square.rows() );
    EXPECT_EQ( ncols, square.columns() );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( testarr[ i ][ j ], square.get( i, j ) );
        }
    }
    NCPA::arrays::free_array( testarr, nrows, ncols );
}

TEST_F( _TEST_TITLE_, AsArrayThrowsIfDimensionsWrong ) {
    size_t nrows = square.rows() + 1, ncols = square.columns() + 2;
    test_t **testarr = NCPA::arrays::zeros<test_t>( nrows, ncols );
    EXPECT_THROW(
        { square.as_array( nrows, ncols, testarr ); }, std::invalid_argument );
}

TEST_F( _TEST_TITLE_, SetMethodsWork ) {
    _TEST_EQ_( square.get( 0, 0 ), -2.0 );
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
    square.set_row( row, 2, row0_inds, row0 );
    for ( size_t col = 0; col < 2; col++ ) {
        _TEST_EQ_( square.get( row, col ), testval );
    }
}

TEST_F( _TEST_TITLE_, VectorSetRowMethodWorks ) {
    size_t row    = 1;
    test_t oldval = square.get( row, 0 );
    std::vector<test_t> row1( 2, testval );
    std::vector<size_t> row1_inds( { 1, 2 } );
    square.set_row( row, row1_inds, row1 );
    _TEST_EQ_( square.get( row, 0 ), oldval );
    _TEST_EQ_( square.get( row, 1 ), testval );
    _TEST_EQ_( square.get( row, 2 ), testval );
}

TEST_F( _TEST_TITLE_, InitListSetRowMethodWorks ) {
    size_t row    = 2;
    test_t oldval = square.get( row, 1 );
    square.set_row( row, { 2, 3 }, { testval, testval } );
    _TEST_EQ_( square.get( row, 1 ), oldval );
    _TEST_EQ_( square.get( row, 2 ), testval );
    _TEST_EQ_( square.get( row, 3 ), testval );
}

TEST_F( _TEST_TITLE_, BareVectorSetRowMethodThrowsException) {
    size_t row = 3;
    std::vector<test_t> row3( 3, testval );
    test_t oldval = square.get( row, 4 );
    EXPECT_THROW( {square.set_row( row, row3 );}, std::logic_error );
    
}

TEST_F( _TEST_TITLE_, ConstantSetRowMethodWorks ) {
    size_t row = 3;
    square.set_row( row, testval );
    _TEST_EQ_( square.get( row, 3 ), testval );
    _TEST_EQ_( square.get( row, 4 ), testval );
}

TEST_F( _TEST_TITLE_, ArraySetColumnMethodWorks ) {
    test_t *col0          = NCPA::arrays::zeros<test_t>( 2 );
    size_t row0_inds[ 2 ] = { 0, 1 };
    NCPA::arrays::fill( col0, 2, testval );
    size_t col = 1;
    square.set_column( col, 2, row0_inds, col0 );
    for ( size_t row = 0; row < 2; row++ ) {
        _TEST_EQ_( square.get( row, col ), testval );
    }
}

TEST_F( _TEST_TITLE_, VectorSetColumnMethodWorks ) {
    size_t col    = 1;
    test_t oldval = square.get( 2, col );
    std::vector<test_t> col1( 2, testval );
    std::vector<size_t> col1_inds( { 0, 1 } );
    square.set_column( col, col1_inds, col1 );
    _TEST_EQ_( square.get( 0, col ), testval );
    _TEST_EQ_( square.get( 1, col ), testval );
    _TEST_EQ_( square.get( 2, col ), oldval );
}

TEST_F( _TEST_TITLE_, InitListSetColumnMethodWorks ) {
    size_t col    = 2;
    test_t oldval = square.get( 3, col );
    square.set_column( col, { 1, 2 }, { testval, testval } );
    _TEST_EQ_( square.get( 1, col ), testval );
    _TEST_EQ_( square.get( 2, col ), testval );
    _TEST_EQ_( square.get( 3, col ), oldval );
}

TEST_F( _TEST_TITLE_, BareVectorSetColumnMethodThrowsLogicError ) {
    size_t col = 3;
    std::vector<test_t> col3( 3, testval );
    EXPECT_THROW( {square.set_column( col, col3 );}, std::logic_error );
}

TEST_F( _TEST_TITLE_, ConstantSetColumnMethodWorks ) {
    size_t col = 3;
    square.set_column( col, testval );
    _TEST_EQ_( square.get( 2, col ), testval );
    _TEST_EQ_( square.get( 3, col ), testval );
}

TEST_F( _TEST_TITLE_, TransposeWorksCorrectly ) {
    mat_t tsym = symmetric, tid = identity, tmr = square;
    EXPECT_TRUE( tmr.transpose().equals( square ) );
    tmr = square;
    EXPECT_TRUE( tmr.transpose().transpose().equals( square ) );
}

TEST_F( _TEST_TITLE_, MatrixVectorRightMultiplicationIsCorrect ) {
    EXPECT_TRUE( square.right_multiply( testvec )->equals( rightvec ) );
}

TEST_F( _TEST_TITLE_, VectorMatrixLeftMultiplicationIsCorrect ) {
    EXPECT_TRUE( square.left_multiply( testvec )->equals( leftvec ) );
}

TEST_F( _TEST_TITLE_, MatrixMatrixMultiplicationIsCorrect ) {
    // cout << "mat1 = " << endl << *square.multiply( identity ) << endl << "mat2 = " << endl << square << endl;
    EXPECT_TRUE( square.multiply( identity )->equals( square ) );
    mat_t lhs = square;
    lhs.zero()
        .set_diagonal( { 1, 2, 3, 4, 5 } )
        .set_diagonal( { 6, 7, 8, 9 }, 1 );
    mat_t rhs = lhs;

    // cout << "LHS: " << endl << lhs << "RHS: " << endl << rhs;
    product2.zero()
        .set_diagonal( { 37, 89, 122, 161, 106 } )
        .set_diagonal( { 18, 35, 56, 81 }, 1 )
        .set_diagonal( { 42, 56, 72 }, 2 );
    // cout << "Expect product = " << endl << product2;
    auto product = lhs.multiply( rhs );
    // cout << "Actual product = " << endl << *product;
    EXPECT_TRUE( product2.equals( *product ) );
}

TEST_F( _TEST_TITLE_, ScaleWorksWithScalar ) {
    mat_t original = square;
    square.scale( 0.5 );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), 0.5 * original.get( i, j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, ScaleWorksWithMatrix ) {
    mat_t original = square;
    square.scale( square );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ),
                       original.get( i, j ) * original.get( i, j ) );
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
                           original.get( i, j ) + testval );
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
            _TEST_EQ_( square.get( i, j ), original.get( i, j ) * 2.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, AddWorksWithScaledMatrix ) {
    mat_t original = square;
    square.add( square, -1.0 );
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), 0.0 );
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
                           original.get( i, j ) + testval );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, PlusEqualOperatorWorksWithMatrix ) {
    mat_t original  = square;
    square         += square;
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), original.get( i, j ) * 2.0 );
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
                           original.get( i, j ) - testval );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, MinusEqualWorksWithMatrix ) {
    mat_t original  = square;
    square         -= square;
    for ( size_t i = 0; i < square.rows(); i++ ) {
        for ( size_t j = 0; j < square.columns(); j++ ) {
            _TEST_EQ_( square.get( i, j ), 0.0 );
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
                       original.get( i, j ) * original.get( i, j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, DivideEqualOperatorWorksWithScalar ) {
    mat_t original  = product;
    product        /= testval;
    for ( size_t i = 0; i < product.rows(); i++ ) {
        for ( size_t j = 0; j < product.columns(); j++ ) {
            _TEST_EQ_( product.get( i, j ), original.get( i, j ) / testval );
        }
    }
}

TEST_F( _TEST_TITLE_, IsDiagonalReturnsCorrectly ) {
    EXPECT_TRUE( identity.is_diagonal() );
    EXPECT_FALSE( square.is_diagonal() );
    // EXPECT_FALSE( product.is_diagonal() );
    EXPECT_FALSE( square.is_diagonal() );
    EXPECT_FALSE( square.is_diagonal() );
    EXPECT_TRUE( empty.is_diagonal() );
    EXPECT_TRUE( zeromat.is_diagonal() );
}

TEST_F( _TEST_TITLE_, IsTridiagonalReturnsCorrectly ) {
    EXPECT_TRUE( identity.is_tridiagonal() );
    EXPECT_TRUE( square.is_tridiagonal() );
    // EXPECT_FALSE( product.is_tridiagonal() );
    EXPECT_TRUE( square.is_tridiagonal() );
    EXPECT_TRUE( square.is_tridiagonal() );
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
    EXPECT_FALSE( square.is_lower_triangular() );
    EXPECT_FALSE( square.is_lower_triangular() );
    EXPECT_TRUE( empty.is_lower_triangular() );
    EXPECT_TRUE( zeromat.is_lower_triangular() );

    for ( size_t i = 1; i < zeromat.rows(); i++ ) {
        zeromat.set( i, i - 1, testval );
    }
    EXPECT_FALSE( zeromat.is_lower_triangular() );
    zeromat.set( 0, 1, testval );
    EXPECT_FALSE( zeromat.is_lower_triangular() );
}

TEST_F( _TEST_TITLE_, IsUpperTriangularIsCorrect ) {
    EXPECT_TRUE( identity.is_upper_triangular() );
    EXPECT_FALSE( square.is_upper_triangular() );
    // EXPECT_FALSE( product.is_upper_triangular() );
    EXPECT_FALSE( square.is_upper_triangular() );
    EXPECT_FALSE( square.is_upper_triangular() );
    EXPECT_TRUE( empty.is_upper_triangular() );
    EXPECT_TRUE( zeromat.is_upper_triangular() );

    for ( size_t i = 1; i < zeromat.rows(); i++ ) {
        zeromat.set( i - 1, i, testval );
    }
    EXPECT_FALSE( zeromat.is_upper_triangular() );
    zeromat.set( 1, 0, testval );
    EXPECT_FALSE( zeromat.is_upper_triangular() );
}
