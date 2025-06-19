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
typedef Matrix<test_t> Mat_t;
typedef BlockMatrix<test_t> bmat_t;
typedef dense_vector<test_t> vec_t;
typedef Vector<test_t> Vec_t;
typedef std::unique_ptr<Vec_t> Vec_ptr_t;

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here

            blockrows = 3;
            blockcols = 3;
            tilerows  = 5;
            tilecols  = 5;

            blockmat.block_type( matrix_t::BAND_DIAGONAL )
                .resize( blockrows, blockcols, tilerows, tilecols )
                .identity();
            empty.block_type( matrix_t::BAND_DIAGONAL );
            // rect.block_type( matrix_t::DENSE );

            ztile = MatrixFactory<test_t>::build( matrix_t::BAND_DIAGONAL );
            ztile.resize( tilerows, tilecols );
            itile = ztile;
            itile.identity();
            tri_d = itile;
            tri_d.set_diagonal( 2.0, 1 ).set_diagonal( 2.0, -1 );
            Imat = ztile;
            Imat.identity( blockmat.rows(), blockmat.columns() );
            recttile = MatrixFactory<test_t>::build( matrix_t::DENSE );
            recttile.resize( 4, 3 );
            for (size_t i = 0; i < 4; ++i) {
                recttile.set_row( i, { i, i, i } );
            }

        }  // void TearDown() override {}

        // declare stuff here
        bmat_t blockmat, empty, rect;
        Mat_t itile, ztile, tri_d, Imat, recttile;
        size_t blockrows, blockcols, tilerows, tilecols;
};

TEST_F( _TEST_TITLE_, DefaultConstructorIsCorrect ) {
    EXPECT_EQ( empty.rows(), 0 );
    EXPECT_EQ( empty.columns(), 0 );
    EXPECT_THROW( { rect.resize( 2, 2, 2, 2 ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, SizedConstructorIsCorrect ) {
    EXPECT_EQ( blockmat.rows(), blockrows * tilerows );
    EXPECT_EQ( blockmat.columns(), blockcols * tilecols );
    EXPECT_EQ( blockmat.block_rows(), blockrows );
    EXPECT_EQ( blockmat.block_columns(), blockcols );
    EXPECT_EQ( blockmat.rows_per_block(), tilerows );
    EXPECT_EQ( blockmat.columns_per_block(), tilecols );
}

TEST_F( _TEST_TITLE_, EqualsMethodReturnsTrueForEqual ) {
    EXPECT_TRUE( blockmat.equals( blockmat ) );
}

TEST_F( _TEST_TITLE_, EqualsMethodReturnsTrueForStandardMatrix ) {
    EXPECT_TRUE( blockmat.equals( Imat ) );
}

TEST_F( _TEST_TITLE_, EqualsMethodReturnsFalseForUnequal ) {
    EXPECT_FALSE( blockmat.equals( empty ) );
}

TEST_F( _TEST_TITLE_, EqualsMethodReturnsFalseForStandardMatrix ) {
    EXPECT_FALSE( blockmat.equals( tri_d ) );
}

TEST_F( _TEST_TITLE_, EqualsOperatorReturnsTrueForEqual ) {
    EXPECT_TRUE( blockmat == blockmat );
}

TEST_F( _TEST_TITLE_, EqualsOperatorReturnsFalseForUnequal ) {
    EXPECT_FALSE( blockmat == empty );
}

TEST_F( _TEST_TITLE_, NotEqualsOperatorReturnsFalseForEqual ) {
    EXPECT_TRUE( blockmat != empty );
}

TEST_F( _TEST_TITLE_, NotEqualsOperatorReturnsTrueForUnequal ) {
    EXPECT_FALSE( blockmat != blockmat );
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

TEST_F( _TEST_TITLE_, SameMethodWorks ) {
    ASSERT_FALSE( empty.same( blockmat ) );
    empty.resize( blockrows, blockcols, tilerows, tilecols );
    EXPECT_TRUE( empty.same( blockmat ) );
}

TEST_F( _TEST_TITLE_, LikeMethodWorks ) {
    ASSERT_TRUE( empty.is_empty() );
    empty.like( blockmat );
    EXPECT_TRUE( empty.same( blockmat ) );
}

TEST_F( _TEST_TITLE_, GetMethodWorks ) {
    EXPECT_DOUBLE_EQ( blockmat.get( 0, 0 ).real(), 1.0 );
    EXPECT_EQ( blockmat.get( 0, 0 ).imag(), 0.0 );
    EXPECT_EQ( blockmat.get( 1, 0 ).real(), 0.0 );
    EXPECT_EQ( blockmat.get( 1, 0 ).imag(), 0.0 );
}

TEST_F( _TEST_TITLE_, GetBlockMethodWorks ) {
    ASSERT_FALSE( ztile == itile );
    for (size_t r = 0; r < blockrows; ++r) {
        for (size_t c = 0; c < blockcols; ++c) {
            if (r == c) {
                EXPECT_TRUE( blockmat.get_block( r, c ) == itile );
            } else {
                EXPECT_TRUE( blockmat.get_block( r, c ) == ztile );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, SetBlockMethodWorks ) {
    ASSERT_FALSE( ztile == itile );
    ASSERT_TRUE( blockrows == blockcols );
    ASSERT_TRUE( tilerows == tilecols );
    empty.resize( blockrows, blockcols, tilerows, tilecols );
    for (size_t i = 0; i < blockrows; ++i) {
        empty.set_block( i, i, itile );
    }
    EXPECT_TRUE( empty.is_identity() );
}

TEST_F( _TEST_TITLE_, IsSquareMethodWorks ) {
    EXPECT_TRUE( blockmat.is_square() );
    EXPECT_TRUE( empty.is_square() );
    blockmat.resize( 3, 1, 2, 2 );
    EXPECT_FALSE( blockmat.is_square() );
}

TEST_F( _TEST_TITLE_, IsEmptyMethodWorks ) {
    EXPECT_FALSE( blockmat.is_empty() );
    EXPECT_TRUE( empty.is_empty() );
}

TEST_F( _TEST_TITLE_, IsIdentityMethodWorks ) {
    EXPECT_TRUE( blockmat.is_identity() );
    EXPECT_TRUE( empty.is_identity() );
    blockmat.set( 2, 2, 2.0 );
    EXPECT_FALSE( blockmat.is_identity() );
}

TEST_F( _TEST_TITLE_, IsDiagonalMethodWorks ) {
    EXPECT_TRUE( blockmat.is_diagonal() );
    EXPECT_TRUE( empty.is_diagonal() );
    blockmat.set( 2, 2, 2.0 );
    EXPECT_TRUE( blockmat.is_diagonal() );
    blockmat.set( 0, 2, 2.0 );
    EXPECT_FALSE( blockmat.is_diagonal() );
}

TEST_F( _TEST_TITLE_, IsSymmetricMethodWorks ) {
    EXPECT_TRUE( blockmat.is_symmetric() );
    EXPECT_TRUE( empty.is_symmetric() );
    blockmat.set( 2, 2, 2.0 );
    EXPECT_TRUE( blockmat.is_symmetric() );
    blockmat.set( 0, 2, 2.0 );
    EXPECT_FALSE( blockmat.is_symmetric() );
    blockmat.set( 2, 0, 2.0 );
    EXPECT_TRUE( blockmat.is_symmetric() );
}

TEST_F( _TEST_TITLE_, IsTridiagonalMethodWorks ) {
    ASSERT_TRUE( blockrows == blockcols );
    ASSERT_TRUE( tilerows == tilecols );
    ASSERT_TRUE( empty.is_empty() );
    ASSERT_TRUE( tri_d.is_tridiagonal() );
    empty.resize( blockrows, blockcols, tilerows, tilecols );
    ASSERT_TRUE( empty.is_zero() );
    for (size_t i = 0; i < blockrows; ++i) {
        empty.set_block( i, i, tri_d );
    }
    EXPECT_FALSE( empty.is_zero() );
    EXPECT_TRUE( empty.is_tridiagonal() );
}

TEST_F( _TEST_TITLE_, IsBlockDiagonalMethodWorks ) {
    ASSERT_TRUE( blockmat.is_block_diagonal() );
    blockmat.set_block( 1, 0, tri_d );
    EXPECT_FALSE( blockmat.is_block_diagonal() );
}

TEST_F( _TEST_TITLE_, IsBlockTridiagonalMethodWorks ) {
    ASSERT_TRUE( blockrows == blockcols );
    ASSERT_TRUE( blockmat.is_block_tridiagonal() );
    for (size_t i = 0; i < blockrows - 1; ++i) {
        blockmat.set_block( i + 1, i, itile ).set_block( i, i + 1, itile );
    }
    EXPECT_TRUE( blockmat.is_block_tridiagonal() );
    blockmat.set_block( 2, 0, itile );
    blockmat.set_block( 0, 2, itile );
    EXPECT_FALSE( blockmat.is_block_tridiagonal() );
}

TEST_F( _TEST_TITLE_, BooleanOperatorWorks ) {
    EXPECT_TRUE( ztile );
    EXPECT_TRUE( itile );
    EXPECT_FALSE( rect );
}

TEST_F( _TEST_TITLE_, ClearMethodWorks ) {
    ASSERT_FALSE( blockmat == empty );
    blockmat.clear();
    EXPECT_TRUE( blockmat == empty );
    EXPECT_TRUE( blockmat.is_empty() );
    EXPECT_EQ( blockmat.rows(), 0 );
    EXPECT_EQ( blockmat.columns(), 0 );
    EXPECT_EQ( blockmat.block_rows(), 0 );
    EXPECT_EQ( blockmat.block_columns(), 0 );
    EXPECT_EQ( blockmat.rows_per_block(), 0 );
    EXPECT_EQ( blockmat.columns_per_block(), 0 );
}

TEST_F( _TEST_TITLE_, BlockTypeMethodReturnsCorrectType ) {
    EXPECT_TRUE( blockmat.block_type() == matrix_t::BAND_DIAGONAL );
    blockmat.block_type( matrix_t::DENSE );
}

TEST_F( _TEST_TITLE_, ResizeMethodShrinksNumberOfBlocksCorrectly ) {
    ASSERT_EQ( blockmat.block_rows(), blockrows );
    ASSERT_EQ( blockmat.block_columns(), blockcols );
    ASSERT_EQ( blockmat.rows(), blockrows * tilerows );
    ASSERT_EQ( blockmat.columns(), blockcols * tilecols );

    // Change number of blocks
    blockmat.resize( blockrows - 1, blockcols, tilerows, tilecols );
    EXPECT_EQ( blockmat.block_rows(), blockrows - 1 );
    EXPECT_EQ( blockmat.block_columns(), blockcols );
    EXPECT_EQ( blockmat.rows(), ( blockrows - 1 ) * tilerows );
    EXPECT_EQ( blockmat.columns(), blockcols * tilecols );
}

TEST_F( _TEST_TITLE_, ResizeMethodGrowsNumberOfBlocksCorrectly ) {
    ASSERT_EQ( blockmat.block_rows(), blockrows );
    ASSERT_EQ( blockmat.block_columns(), blockcols );
    ASSERT_EQ( blockmat.rows(), blockrows * tilerows );
    ASSERT_EQ( blockmat.columns(), blockcols * tilecols );

    // Change number of blocks
    blockmat.resize( blockrows, blockcols + 1, tilerows, tilecols );
    EXPECT_EQ( blockmat.block_rows(), blockrows );
    EXPECT_EQ( blockmat.block_columns(), blockcols + 1 );
    EXPECT_EQ( blockmat.rows(), blockrows * tilerows );
    EXPECT_EQ( blockmat.columns(), ( blockcols + 1 ) * tilecols );
}

TEST_F( _TEST_TITLE_, ResizeMethodShrinksBlockSizeCorrectly ) {
    ASSERT_EQ( blockmat.block_rows(), blockrows );
    ASSERT_EQ( blockmat.block_columns(), blockcols );
    ASSERT_EQ( blockmat.rows(), blockrows * tilerows );
    ASSERT_EQ( blockmat.columns(), blockcols * tilecols );

    // Change number of blocks
    blockmat.resize( blockrows, blockcols, tilerows - 1, tilecols );
    EXPECT_EQ( blockmat.block_rows(), blockrows );
    EXPECT_EQ( blockmat.block_columns(), blockcols );
    EXPECT_EQ( blockmat.rows(), blockrows * ( tilerows - 1 ) );
    EXPECT_EQ( blockmat.columns(), blockcols * tilecols );
}

TEST_F( _TEST_TITLE_, ResizeMethodGrowsBlockSizeCorrectly ) {
    ASSERT_EQ( blockmat.block_rows(), blockrows );
    ASSERT_EQ( blockmat.block_columns(), blockcols );
    ASSERT_EQ( blockmat.rows(), blockrows * tilerows );
    ASSERT_EQ( blockmat.columns(), blockcols * tilecols );

    // Change number of blocks
    blockmat.resize( blockrows, blockcols, tilerows, tilecols + 1 );
    EXPECT_EQ( blockmat.block_rows(), blockrows );
    EXPECT_EQ( blockmat.block_columns(), blockcols );
    EXPECT_EQ( blockmat.rows(), blockrows * tilerows );
    EXPECT_EQ( blockmat.columns(), blockcols * ( tilecols + 1 ) );
}

TEST_F( _TEST_TITLE_, ResizeMethodKeepsExistingBlocksWhenShrinking ) {
    ASSERT_TRUE( blockmat.is_identity() );
    for (size_t i = 0; i < blockrows; ++i) {
        ASSERT_TRUE( blockmat.get_block( i, i ).is_identity() );
    }
    blockmat.resize( blockrows - 1, blockcols - 1 );
    EXPECT_TRUE( blockmat.is_identity() );
    for (size_t i = 0; i < blockrows - 1; ++i) {
        EXPECT_TRUE( blockmat.get_block( i, i ).is_identity() );
    }
}

TEST_F( _TEST_TITLE_, ResizeMethodKeepsExistingBlocksWhenGrowing ) {
    ASSERT_TRUE( blockmat.is_identity() );
    for (size_t i = 0; i < blockrows; ++i) {
        ASSERT_TRUE( blockmat.get_block( i, i ).is_identity() );
    }
    blockmat.resize( blockrows + 1, blockcols + 1 );
    EXPECT_FALSE( blockmat.is_identity() );
    for (size_t i = 0; i < blockrows; ++i) {
        EXPECT_TRUE( blockmat.get_block( i, i ).is_identity() );
    }
    EXPECT_TRUE( blockmat.get_block( blockrows, blockrows ).is_zero() );
}

TEST_F( _TEST_TITLE_, GetRowMethodWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    Vec_ptr_t firstrow = blockmat.get_row( 0 );
    EXPECT_EQ( firstrow->size(), blockmat.columns() );
    EXPECT_DOUBLE_EQ( firstrow->get( 0 ).real(), 1.0 );
    for (size_t i = 1; i < firstrow->size(); ++i) {
        EXPECT_EQ( firstrow->get( i ), 0.0 );
    }

    Vec_ptr_t row9 = blockmat.get_row( 9 );
    EXPECT_EQ( row9->size(), blockmat.columns() );
    for (size_t i = 1; i < row9->size(); ++i) {
        if (i == 9) {
            EXPECT_DOUBLE_EQ( row9->get( i ).real(), 1.0 );
        } else {
            EXPECT_EQ( row9->get( i ).real(), 0.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, GetColumnMethodWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    Vec_ptr_t col6 = blockmat.get_column( 6 );
    EXPECT_EQ( col6->size(), blockmat.rows() );
    for (size_t i = 0; i < col6->size(); ++i) {
        if (i == 6) {
            EXPECT_DOUBLE_EQ( col6->get( i ).real(), 1.0 );
        } else {
            EXPECT_EQ( col6->get( i ).real(), 0.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, GetDiagonalMethodWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    Vec_ptr_t diag = blockmat.get_diagonal();
    EXPECT_EQ( diag->size(), blockmat.rows() );
    for (size_t i = 0; i < diag->size(); ++i) {
        EXPECT_DOUBLE_EQ( diag->get( i ).real(), 1.0 );
    }
    diag = blockmat.get_diagonal( 1 );
    EXPECT_EQ( diag->size(), blockmat.rows() - 1 );
    for (size_t i = 0; i < diag->size(); ++i) {
        EXPECT_EQ( diag->get( i ).real(), 0.0 );
    }

    blockmat.set_block( 1, 0, tri_d );
    diag = blockmat.get_diagonal( -( (int)tilerows ) );
    EXPECT_EQ( diag->size(), blockmat.rows() - tilerows );
    for (size_t i = 0; i < diag->size(); ++i) {
        if (i < tilerows) {
            EXPECT_DOUBLE_EQ( diag->get( i ).real(), 1.0 );
        } else {
            EXPECT_EQ( diag->get( i ).real(), 0.0 );
        }
    }
    diag = blockmat.get_diagonal( -( (int)tilerows ) - 1 );
    EXPECT_EQ( diag->size(), blockmat.rows() - tilerows - 1 );
    for (size_t i = 0; i < diag->size(); ++i) {
        if (i < tilerows - 1) {
            EXPECT_DOUBLE_EQ( diag->get( i ).real(), 2.0 );
        } else {
            EXPECT_EQ( diag->get( i ).real(), 0.0 );
        }
    }
    diag = blockmat.get_diagonal( -( (int)tilerows ) + 1 );
    EXPECT_EQ( diag->size(), blockmat.rows() - tilerows + 1 );
    for (size_t i = 0; i < diag->size(); ++i) {
        if (i == 0) {
            EXPECT_EQ( diag->get( i ).real(), 0.0 );
        } else if (i < tilerows) {
            EXPECT_DOUBLE_EQ( diag->get( i ).real(), 2.0 );
        } else {
            EXPECT_EQ( diag->get( i ).real(), 0.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, AsMatrixMethodWorks ) {
    Mat_t I = blockmat.as_matrix();
    EXPECT_EQ( I.rows(), blockmat.rows() );
    EXPECT_EQ( I.columns(), blockmat.columns() );
    EXPECT_TRUE( I.is_identity() );
    EXPECT_TRUE( I == blockmat );
}

TEST_F( _TEST_TITLE_, ZeroMethodWorksForIndividualElement ) {
    ASSERT_TRUE( blockmat.is_identity() );
    blockmat.zero( 2 * tilerows, 2 * tilerows );
    EXPECT_FALSE( blockmat.is_identity() );
    EXPECT_TRUE( blockmat.get_block( 0, 0 ).is_identity() );
    EXPECT_TRUE( blockmat.get_block( 1, 1 ).is_identity() );
    EXPECT_FALSE( blockmat.get_block( 2, 2 ).is_identity() );
}

TEST_F( _TEST_TITLE_, ZeroMethodWorksForWholeMatrix ) {
    ASSERT_FALSE( blockmat.is_zero() );
    ASSERT_EQ( blockmat.rows(), blockrows * tilerows );
    ASSERT_EQ( blockmat.columns(), blockcols * tilecols );
    blockmat.zero();
    EXPECT_TRUE( blockmat.is_zero() );
    EXPECT_EQ( blockmat.rows(), blockrows * tilerows );
    EXPECT_EQ( blockmat.columns(), blockcols * tilecols );
}

TEST_F( _TEST_TITLE_, TransposeMethodWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    EXPECT_TRUE( blockmat.transpose().is_identity() );

    // Non-square block matrix
    rect.block_type( matrix_t::DENSE ).resize( 2, 3, 4, 3 );
    EXPECT_EQ( rect.rows(), 8 );
    EXPECT_EQ( rect.columns(), 9 );
    BlockMatrix<test_t> rectT = rect;
    rectT.transpose();
    EXPECT_EQ( rectT.block_rows(), rect.block_columns() );
    EXPECT_EQ( rectT.block_columns(), rect.block_rows() );
    EXPECT_EQ( rectT.rows_per_block(), rect.columns_per_block() );
    EXPECT_EQ( rectT.columns_per_block(), rect.rows_per_block() );
    for (size_t r = 0; r < rect.block_rows(); ++r) {
        for (size_t c = 0; c < rect.block_columns(); ++c) {
            Matrix<test_t> orig  = rect.get_block( r, c );
            Matrix<test_t> trans = orig;
            orig.transpose();
            EXPECT_TRUE( rectT.get_block( c, r ) == orig );
        }
    }
}

TEST_F( _TEST_TITLE_, AddMethodWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    blockmat.add( blockmat );
    for (size_t i = 0; i < blockmat.rows(); ++i) {
        EXPECT_DOUBLE_EQ( blockmat.get( i, i ).real(), 2.0 );
    }
    empty.like( blockmat );
    blockmat.add( empty );
    for (size_t i = 0; i < blockmat.rows(); ++i) {
        EXPECT_DOUBLE_EQ( blockmat.get( i, i ).real(), 2.0 );
    }
}

TEST_F( _TEST_TITLE_, AddEqualsOperatorWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    blockmat += blockmat;
    for (size_t i = 0; i < blockmat.rows(); ++i) {
        EXPECT_DOUBLE_EQ( blockmat.get( i, i ).real(), 2.0 );
    }
    empty.like( blockmat );
    blockmat += empty;
    for (size_t i = 0; i < blockmat.rows(); ++i) {
        EXPECT_DOUBLE_EQ( blockmat.get( i, i ).real(), 2.0 );
    }
}

TEST_F( _TEST_TITLE_, AddOperatorWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    bmat_t newmat = blockmat + blockmat;
    for (size_t i = 0; i < newmat.rows(); ++i) {
        EXPECT_DOUBLE_EQ( newmat.get( i, i ).real(), 2.0 );
    }
    empty.like( blockmat );
    bmat_t newmat2 = newmat + empty;
    for (size_t i = 0; i < newmat2.rows(); ++i) {
        EXPECT_DOUBLE_EQ( newmat2.get( i, i ).real(), 2.0 );
    }
}

TEST_F( _TEST_TITLE_, SubtractMethodWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    empty.like( blockmat );
    ASSERT_TRUE( empty.is_zero() );
    blockmat.subtract( empty );
    EXPECT_TRUE( blockmat.is_identity() );
    blockmat.subtract( blockmat );
    EXPECT_TRUE( blockmat.is_zero() );
}

TEST_F( _TEST_TITLE_, MinusEqualOperatorWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    empty.like( blockmat );
    ASSERT_TRUE( empty.is_zero() );
    blockmat -= empty;
    EXPECT_TRUE( blockmat.is_identity() );
    blockmat -= blockmat;
    EXPECT_TRUE( blockmat.is_zero() );
}

TEST_F( _TEST_TITLE_, MinusOperatorWorks ) {
    ASSERT_TRUE( blockmat.is_identity() );
    empty.like( blockmat );
    ASSERT_TRUE( empty.is_zero() );
    bmat_t diff = blockmat - empty;
    EXPECT_TRUE( diff.is_identity() );
    bmat_t diff2 = blockmat - blockmat;
    EXPECT_TRUE( diff2.is_zero() );
}

TEST_F( _TEST_TITLE_, ScaleMethodWorksWithScalar ) {
    ASSERT_TRUE( blockmat.is_identity() );
    empty = blockmat;
    EXPECT_TRUE( empty.scale( 4.0 )
                 == blockmat.add( blockmat ).add( blockmat ) );
}

TEST_F( _TEST_TITLE_, TimesEqualOperatorWorksWithScalar ) {
    ASSERT_TRUE( blockmat.is_identity() );
    empty = blockmat;
    empty *= 3.0;
    EXPECT_TRUE( empty
                 == blockmat + blockmat + blockmat );
}

TEST_F( _TEST_TITLE_, TimesOperatorWorksWithScalar ) {
    ASSERT_TRUE( blockmat.is_identity() );
    empty = blockmat * 3;
    EXPECT_TRUE( empty
                 == blockmat + blockmat + blockmat );
}

TEST_F( _TEST_TITLE_, ScaleMethodWorksWithMatrix ) {
    ASSERT_TRUE( blockmat.is_identity() );
    blockmat.scale( 5.0 ).scale( blockmat );
    for (size_t i = 0; i < blockmat.rows(); ++i) {
        EXPECT_DOUBLE_EQ( blockmat.get( i, i ).real(), 25.0 );
    }
}

TEST_F( _TEST_TITLE_, MultiplyMethodWorksCorrectly ) {
    ASSERT_TRUE( blockmat.is_identity() );
    ASSERT_EQ( blockmat.block_rows(), 3 );
    BlockMatrix<test_t> A = blockmat;
    for (size_t i = 1; i < 3; i++) {
        A.set_block( i, i - 1, -itile );
    }
    BlockMatrix<test_t> B = -A;
    B.transpose();

    Mat_t Amat = Imat;
    Mat_t Bmat = Imat;
    Amat.set_diagonal( -1.0, -5 );
    Bmat.set_diagonal( -1.0 ).set_diagonal( 1.0, 5 );
    EXPECT_TRUE( A == Amat );
    EXPECT_TRUE( B == Bmat );
    A.multiply( B );
    EXPECT_TRUE( A == ( Amat * Bmat ) );
}

TEST_F( _TEST_TITLE_, TimesEqualOperatorWorksCorrectlyWithMatrix ) {
    ASSERT_TRUE( blockmat.is_identity() );
    ASSERT_EQ( blockmat.block_rows(), 3 );
    BlockMatrix<test_t> A = blockmat;
    for (size_t i = 1; i < 3; i++) {
        A.set_block( i, i - 1, -itile );
    }
    BlockMatrix<test_t> B = -A;
    B.transpose();

    Mat_t Amat = Imat;
    Mat_t Bmat = Imat;
    Amat.set_diagonal( -1.0, -5 );
    Bmat.set_diagonal( -1.0 ).set_diagonal( 1.0, 5 );
    EXPECT_TRUE( A == Amat );
    EXPECT_TRUE( B == Bmat );
    A *= B;
    EXPECT_TRUE( A == ( Amat * Bmat ) );
}

TEST_F( _TEST_TITLE_, TimesOperatorWorksCorrectlyWithBlockMatrix ) {
    ASSERT_TRUE( blockmat.is_identity() );
    ASSERT_EQ( blockmat.block_rows(), 3 );
    BlockMatrix<test_t> A = blockmat;
    for (size_t i = 1; i < 3; i++) {
        A.set_block( i, i - 1, -itile );
    }
    BlockMatrix<test_t> B = -A;
    B.transpose();

    Mat_t Amat = Imat;
    Mat_t Bmat = Imat;
    Amat.set_diagonal( -1.0, -5 );
    Bmat.set_diagonal( -1.0 ).set_diagonal( 1.0, 5 );
    EXPECT_TRUE( A == Amat );
    EXPECT_TRUE( B == Bmat );
    bmat_t C = A * B;
    EXPECT_TRUE( C == ( Amat * Bmat ) );
    C = B * A;
    EXPECT_TRUE( C == ( Bmat * Amat ) );
}

TEST_F( _TEST_TITLE_, TimesOperatorWorksCorrectlyWithStandardMatrix ) {
    ASSERT_TRUE( blockmat.is_identity() );
    ASSERT_EQ( blockmat.block_rows(), 3 );
    BlockMatrix<test_t> A = blockmat;
    for (size_t i = 1; i < 3; i++) {
        A.set_block( i, i - 1, -itile );
    }
    BlockMatrix<test_t> B = -A;
    B.transpose();

    Mat_t Amat = Imat;
    Mat_t Bmat = Imat;
    Amat.set_diagonal( -1.0, -5 );
    Bmat.set_diagonal( -1.0 ).set_diagonal( 1.0, 5 );
    EXPECT_TRUE( A == Amat );
    EXPECT_TRUE( B == Bmat );
    Mat_t C = Amat * B;
    EXPECT_TRUE( C == ( Amat * Bmat ) );
    C = A * Bmat;
    EXPECT_TRUE( C == ( Amat * Bmat ) );
    C = B * Amat;
    EXPECT_TRUE( C == ( Bmat * Amat ) );
    C = Bmat * A;
    EXPECT_TRUE( C == ( Bmat * Amat ) );
}