#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/linearalgebra.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <complex>
#include <memory>
#include <vector>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::linear;

typedef double test_t;
typedef details::dense_matrix<test_t> mat_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraMatrixDenseMatrixTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            dim1   = 3;
            testval = 42.0;
            square = mat_t( dim1, dim1 );
            for ( size_t i = 0; i < dim1; i++ ) {
                double di = (double)( i + 1 );
                square.set_row( i, { di, di, di } );
            }

            wrapper1 = Matrix<test_t>( square.clone() );
            product = wrapper1;
            product *= 6.0;

        }  // void TearDown() override {}

        // declare stuff here
        details::dense_matrix<test_t> square, dproduct;
        std::vector<std::vector<test_t>> svd;

        Matrix<test_t> wrapper1, wrapper2, product;
        size_t matrows, matcols, dim1, dim2;
        const test_t zero = NCPA::math::zero<test_t>(),
                     one  = NCPA::math::one<test_t>();
        test_t testval;
};

TEST_F( _TEST_TITLE_, DefaultConstructorGivesEmptyMatrix ) {
    EXPECT_EQ( wrapper2.rows(), 0 );
    EXPECT_EQ( wrapper2.columns(), 0 );
}

TEST_F( _TEST_TITLE_, EqualsMethodWorks ) {
    EXPECT_TRUE( wrapper1.equals( wrapper1 ) );
    EXPECT_FALSE( wrapper1.equals( wrapper2 ) );
}

TEST_F( _TEST_TITLE_, EqualityOperatorWorks ) {
    EXPECT_TRUE( wrapper1 == wrapper1 );
    EXPECT_FALSE( wrapper1 == wrapper2 );
}

TEST_F( _TEST_TITLE_, InequalityOperatorWorks ) {
    EXPECT_FALSE( wrapper1 != wrapper1 );
    EXPECT_TRUE( wrapper1 != wrapper2 );
}

TEST_F( _TEST_TITLE_, PointerConstructorWorks ) {
    wrapper2 = Matrix<test_t>( std::unique_ptr<mat_t>( new mat_t( square ) ) );
    EXPECT_TRUE( wrapper1 == wrapper2 );
}

TEST_F( _TEST_TITLE_, CopyConstructorWorks ) {
    wrapper2 = Matrix<test_t>( wrapper1 );
    EXPECT_TRUE( wrapper1 == wrapper2 );
}

TEST_F( _TEST_TITLE_, AssignmentOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_TRUE( wrapper1 == wrapper2 );
}

TEST_F( _TEST_TITLE_, IndexingOperatorsWork ) {
    wrapper2 = wrapper1;
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1[i][j], wrapper2.get( i, j ) );
        }
    }
}

TEST_F( _TEST_TITLE_, SwapWorks ) {
    Matrix<test_t> mat3 = wrapper1;
    ASSERT_TRUE( mat3 == wrapper1 );
    ASSERT_FALSE( mat3 == wrapper2 );
    swap( wrapper1, wrapper2 );
    EXPECT_FALSE( mat3 == wrapper1 );
    EXPECT_TRUE( mat3 == wrapper2 );
}

TEST_F( _TEST_TITLE_, ResizeWorks ) {
    wrapper1.resize( 2, 2 );
    EXPECT_EQ( wrapper1.rows(), 2 );
    EXPECT_EQ( wrapper1.columns(), 2 );
    EXPECT_THROW( { wrapper1[2][2]; }, std::out_of_range );
}

TEST_F( _TEST_TITLE_, ResizePreservesExistingValues ) {
    wrapper2 = wrapper1;
    wrapper1.resize( dim1 + 1, dim1 + 1 );
    for ( size_t i = 0; i < dim1; i++ ) {
        for ( size_t j = 0; j < dim1; j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ), wrapper2.get( i, j ) );
        }
    }
    _TEST_EQ_( wrapper1.get( dim1, dim1 ), zero );
}

TEST_F( _TEST_TITLE_, IsRowMatrixIsCorrect ) {
    EXPECT_FALSE( wrapper1.is_row_matrix() );
    wrapper1.resize( 1, 3 );
    EXPECT_TRUE( wrapper1.is_row_matrix() );
    EXPECT_FALSE( wrapper1.is_column_matrix() );
}

TEST_F( _TEST_TITLE_, IsColumnMatrixIsCorrect ) {
    EXPECT_FALSE( wrapper1.is_column_matrix() );
    wrapper1.resize( 3, 1 );
    EXPECT_TRUE( wrapper1.is_column_matrix() );
    EXPECT_FALSE( wrapper1.is_row_matrix() );
}

TEST_F( _TEST_TITLE_, IsEmptyIsCorrect ) {
    EXPECT_TRUE( wrapper2.is_empty() );
    EXPECT_FALSE( wrapper1.is_empty() );
}

TEST_F( _TEST_TITLE_, GetReturnsZeroIfEmpty ) {
    _TEST_EQ_( wrapper2.get( 0, 0 ), zero );
}

TEST_F( _TEST_TITLE_, GetRowReturnsRowMatrix ) {
    auto mat3 = *( wrapper1.get_row( 0 ) );
    EXPECT_EQ( mat3.rows(), 1 );
    EXPECT_EQ( mat3.columns(), 3 );
    EXPECT_TRUE( mat3.is_row_matrix() );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( mat3.get( 0, i ), wrapper1.get( 0, i ) );
        _TEST_EQ_( mat3.get( i ), wrapper1.get( 0, i ) );
    }
}

TEST_F( _TEST_TITLE_, GetColumnReturnsRowMatrix ) {
    auto mat3 = *( wrapper1.get_column( 0 ) );
    EXPECT_EQ( mat3.columns(), 1 );
    EXPECT_EQ( mat3.rows(), 3 );
    EXPECT_TRUE( mat3.is_column_matrix() );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( mat3.get( i, 0 ), wrapper1.get( i, 0 ) );
        _TEST_EQ_( mat3.get( i ), wrapper1.get( i, 0 ) );
    }
}

TEST_F( _TEST_TITLE_, GetRowVectorReturnsRowAsVector ) {
    auto vec = *( wrapper1.get_row_vector( 0 ) );
    EXPECT_EQ( vec.size(), 3 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( vec[ i ], wrapper1.get( 0, i ) );
    }
}

TEST_F( _TEST_TITLE_, GetColumnVectorReturnsColumnAsVector ) {
    auto vec = *( wrapper1.get_column_vector( 0 ) );
    EXPECT_EQ( vec.size(), 3 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( vec[ i ], wrapper1.get( i, 0 ) );
    }
}

TEST_F( _TEST_TITLE_, GetDiagonalReturnsExpectedVector ) {
    auto vec = *( wrapper1.get_diagonal( 0 ) );
    EXPECT_EQ( vec.size(), dim1 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( vec[ i ], wrapper1.get( i, 0 ) );
    }

    vec = *( wrapper1.get_diagonal( 1 ) );
    EXPECT_EQ( vec.size(), dim1 - 1 );
    for ( auto i = 0; i < dim1 - 1; i++ ) {
        _TEST_EQ_( vec[ i ], wrapper1.get( i, i + 1 ) );
    }

    vec = *( wrapper1.get_diagonal( -1 ) );
    EXPECT_EQ( vec.size(), dim1 - 1 );
    for ( auto i = 0; i < dim1 - 1; i++ ) {
        _TEST_EQ_( vec[ i ], wrapper1.get( i + 1, i ) );
    }
}

TEST_F( _TEST_TITLE_, SetRowMethodsWork ) {
    wrapper2 = wrapper1;
    for ( auto i = 0; i < dim1; i++ ) {
        EXPECT_NE( wrapper2.get( 1, i ), one );
    }

    // constant
    wrapper2.set_row( 1, one );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;

    // array
    test_t testrow[ 3 ] = { one, one, one };
    size_t inds[ 3 ]    = { 0, 1, 2 };
    wrapper2.set_row( 1, 3, inds, testrow );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;

    // vector with indices
    std::vector<test_t> testrow_v( testrow, testrow + dim1 );
    std::vector<size_t> inds_v( inds, inds + dim1 );
    wrapper2.set_row( 1, inds_v, testrow_v );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;

    // vector without indices
    wrapper2.set_row( 1, testrow_v );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;

    // initializer list with indices
    wrapper2.set_row( 1, { 0, 1, 2 }, { one, one, one } );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;

    // initializer list without indices
    wrapper2.set_row( 1, { one, one, one } );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;

    // abstract_vector
    auto avec = *( wrapper1.get_row_vector( 0 ) );
    wrapper2.set_row( 1, avec );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;

    // Vector
    Vector<test_t> aVec( *( wrapper1.get_row_vector( 0 ) ) );
    wrapper2.set_row( 1, aVec );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;

    // row from another matrix
    wrapper2.set_row( 1, wrapper1, 0 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;

    // row matrix
    wrapper2.set_row( 1, *wrapper1.get_row( 0 ) );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;
    wrapper2.set_row( 1, wrapper1.get_row( 0 ) );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( 1, i ), one );
    }
    wrapper2 = wrapper1;
}

TEST_F( _TEST_TITLE_, SetColumnMethodsWork ) {
    wrapper2       = wrapper1;
    test_t testval = one * 5.0;
    for ( auto i = 0; i < dim1; i++ ) {
        EXPECT_NE( wrapper2.get( i, 1 ), testval );
    }

    // constant
    wrapper2.set_column( 1, testval );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), testval );
    }
    wrapper2 = wrapper1;

    // array
    test_t testcol[ 3 ] = { testval, testval, testval };
    size_t inds[ 3 ]    = { 0, 1, 2 };
    wrapper2.set_column( 1, 3, inds, testcol );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), testval );
    }
    wrapper2 = wrapper1;

    // vector with indices
    std::vector<test_t> testcol_v( testcol, testcol + dim1 );
    std::vector<size_t> inds_v( inds, inds + dim1 );
    wrapper2.set_column( 1, inds_v, testcol_v );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), testval );
    }
    wrapper2 = wrapper1;

    // vector without indices
    wrapper2.set_column( 1, testcol_v );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), testval );
    }
    wrapper2 = wrapper1;

    // initializer list with indices
    wrapper2.set_column( 1, { 0, 1, 2 }, { testval, testval, testval } );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), testval );
    }
    wrapper2 = wrapper1;

    // initializer list without indices
    wrapper2.set_column( 1, { testval, testval, testval } );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), testval );
    }
    wrapper2 = wrapper1;

    // abstract_vector
    auto avec = *( wrapper1.get_column_vector( 0 ) );
    wrapper2.set_column( 1, avec );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
    }
    wrapper2 = wrapper1;

    // Vector
    Vector<test_t> aVec( *( wrapper1.get_column_vector( 0 ) ) );
    wrapper2.set_column( 1, aVec );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
    }
    wrapper2 = wrapper1;

    // row from another matrix
    wrapper2.set_column( 1, wrapper1, 0 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
    }
    wrapper2 = wrapper1;

    // row matrix
    wrapper2.set_column( 1, *wrapper1.get_column( 0 ) );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
    }
    wrapper2 = wrapper1;
    wrapper2.set_column( 1, wrapper1.get_column( 0 ) );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
    }
    wrapper2 = wrapper1;
}

TEST_F( _TEST_TITLE_,
        NCPALinearAlgebraMatrixDenseMatrixTest_SetDiagonalMethodsWork_Test ) {
    wrapper2 = wrapper1;
    wrapper2.clear().resize( 5, 5 );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper2.get( i, i ), zero );
    }

    // constant
    wrapper2.set_diagonal( one );
    EXPECT_TRUE( wrapper2.is_identity() );
    wrapper2.zero();
    EXPECT_FALSE( wrapper2.is_identity() );

    // array
    test_t testrow[ 5 ] = { one, one, one, one, one };
    wrapper2.set_diagonal( 5, testrow );
    EXPECT_TRUE( wrapper2.is_identity() );
    wrapper2.zero();
    EXPECT_FALSE( wrapper2.is_identity() );

    // vector
    vector<test_t> testrow_v( testrow, testrow + 5 );
    wrapper2.set_diagonal( testrow_v );
    EXPECT_TRUE( wrapper2.is_identity() );
    wrapper2.zero();
    EXPECT_FALSE( wrapper2.is_identity() );

    // init list
    wrapper2.set_diagonal( { one, one, one, one, one } );
    EXPECT_TRUE( wrapper2.is_identity() );
    wrapper2.zero();
    EXPECT_FALSE( wrapper2.is_identity() );

    // abstract_vector
    auto avec = *( wrapper1.get_column_vector( 0 ) );
    wrapper2.set_diagonal( avec, 2 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, i + 2 ), wrapper1.get( i, 0 ) );
    }
    wrapper2.zero();

    // Vector
    Vector<test_t> aVec( *( wrapper1.get_column_vector( 0 ) ) );
    wrapper2.set_diagonal( aVec, -2 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i + 2, i ), wrapper1.get( i, 0 ) );
    }
}

TEST_F( _TEST_TITLE_, TransposeWorks ) {
    wrapper2 = wrapper1;
    wrapper2.transpose();
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ), wrapper2.get( j, i ) );
        }
    }
}

TEST_F( _TEST_TITLE_, AddWorksWithMatrix ) {
    wrapper2 = wrapper1;
    wrapper1.add( wrapper1 );
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ), ( 2.0 * wrapper2.get( i, j ) ) );
        }
    }
}

TEST_F( _TEST_TITLE_, AddWorksWithScalar ) {
    wrapper2 = wrapper1;
    wrapper1.add( testval );
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ), (testval + wrapper2.get( i, j ) ) );
        }
    }
}

TEST_F( _TEST_TITLE_, SubtractWorksWithMatrix ) {
    wrapper2 = wrapper1;
    wrapper1.subtract( wrapper1 );
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ), zero );
        }
    }
}

TEST_F( _TEST_TITLE_, SubtractWorksWithScalar ) {
    wrapper2 = wrapper1;
    wrapper1.subtract( testval );
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ), (wrapper2.get( i, j ) - testval) );
        }
    }
}

TEST_F( _TEST_TITLE_, ScaleWorksWithMatrix ) {
    wrapper2 = wrapper1;
    wrapper1.scale( wrapper1 );
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ), (wrapper2.get( i, j ) * wrapper2.get( i, j )) );
        }
    }
}

TEST_F( _TEST_TITLE_, ScaleWorksWithScalar ) {
    wrapper2 = wrapper1;
    wrapper1.scale( testval );
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ), (wrapper2.get( i, j ) * testval) );
        }
    }
}

TEST_F( _TEST_TITLE_, IdentityWorks ) {
    ASSERT_FALSE( wrapper1.is_identity() );
    wrapper1.identity();
    EXPECT_TRUE( wrapper1.is_identity() );
}

TEST_F( _TEST_TITLE_, IdentityDoesNotWorkForNonSquareMatrix ) {
    wrapper1.resize( 2, 4 );
    EXPECT_THROW( { wrapper1.identity(); }, std::invalid_argument);
}

TEST_F( _TEST_TITLE_, MultiplyWorks ) {
    Matrix<test_t> mat1 = wrapper1;
    mat1.identity();
    EXPECT_TRUE( wrapper1.multiply( mat1 ) == wrapper1 );
    // Matrix<test_t> product2 = wrapper1;
    // product2.multiply( wrapper1 );
    EXPECT_TRUE( wrapper1.multiply(wrapper1) == product );
    // for (size_t i = 0; i < 3; i++) {
    //     for (size_t j = 0; j < 3; j++) {
    //         _TEST_EQ_( product2[i][j], product[i][j] );
    //     }
    // }

    Matrix<test_t> mat2 = mat1, mat3 = mat1;
    mat2.resize( 3, 5 ).set( 1.0 );
    mat3.resize( 5, 3 ).set( 1.0 );
    product.resize( 5, 5 ).set( 3.0 );
    EXPECT_TRUE( mat3.multiply(mat2) == product );
    mat2.resize( 3, 5 ).set( 1.0 );
    mat3.resize( 5, 2 ).set( 1.0 );
    product.resize( 3, 2 ).set( 5.0 );
    EXPECT_TRUE( mat2.multiply(mat3) == product );
}

TEST_F( _TEST_TITLE_, MultiplyOperatorWorks ) {
    Matrix<test_t> mat1 = wrapper1;
    mat1.identity();
    EXPECT_TRUE( wrapper1 * mat1 == wrapper1 );
    EXPECT_TRUE( wrapper1 * wrapper1 == product );

    Matrix<test_t> mat2 = mat1, mat3 = mat1;
    mat2.resize( 3, 5 ).set( 1.0 );
    mat3.resize( 5, 3 ).set( 1.0 );
    product.resize( 5, 5 ).set( 3.0 );
    EXPECT_TRUE( mat3 * mat2 == product );
    mat3.resize( 5, 2 );
    product.resize( 3, 2 ).set( 5.0 );
    EXPECT_TRUE( mat2 * mat3 == product );
}

TEST_F( _TEST_TITLE_, OtherBinaryOperatorsWork ) {
    EXPECT_TRUE( wrapper1 + wrapper1 == wrapper1 * 2.0 );
    EXPECT_TRUE( wrapper1 - wrapper1 == wrapper1 * 0.0 );
    wrapper1.set( 1.0 );
    wrapper2 = wrapper1;
    wrapper1 *= 3.0;
    EXPECT_TRUE( wrapper1 == wrapper2 * 3.0 );
}