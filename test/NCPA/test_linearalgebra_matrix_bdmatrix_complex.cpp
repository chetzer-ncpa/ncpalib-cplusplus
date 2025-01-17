#define NCPA_DEBUG_ON

#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/linearalgebra.hpp"
#include "NCPA/logging.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <complex>
#include <memory>
#include <vector>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::linear;

typedef complex<double> test_t;
typedef band_diagonal_matrix<test_t> mat_t;

#define _TEST_EQ_       EXPECT_COMPLEX_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_COMPLEX_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraMatrixBandDiagonalComplexMatrixTest

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

            wrapper1 = Matrix<test_t>( square.clone() );
            wproduct = wrapper1;

        }  // void TearDown() override {}

        // declare stuff here
        mat_t empty, square, more_rows, more_cols, product, identity,
            symmetric, product2, zeromat;
        // Matrix<test_t> mat1, mat2;
        const size_t dim1 = 5, dim2 = 8;
        test_t testval = test_t( -4.2, 2.1 );
        test_t diagval = test_t( -2.0, 2.0 ), supdiagval = test_t( 1.0, -1.0 );

        Matrix<test_t> wrapper1, wrapper2, wproduct;
        const test_t zero = NCPA::math::zero<test_t>(),
                     one  = NCPA::math::one<test_t>();
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
            _TEST_EQ_( wrapper1.get( i, j ), wrapper2.get( i, j ) );
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
    // cout << "Resizing..." << endl;
    wrapper1.resize( 2, 2 );
    // cout << "Successful!" << endl;
    EXPECT_EQ( wrapper1.rows(), 2 );
    EXPECT_EQ( wrapper1.columns(), 2 );
    // EXPECT_THROW( { wrapper1.get( 2, 2 ); }, std::range_error );
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

TEST_F( _TEST_TITLE_, GetRowReturnsRowVector ) {
    auto mat3 = *( wrapper1.get_row( 0 ) );
    EXPECT_EQ( mat3.size(), dim1 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( mat3.get( i ), wrapper1.get( 0, i ) );
    }
}

TEST_F( _TEST_TITLE_, GetColumnReturnsColumnVector ) {
    auto mat3 = *( wrapper1.get_column( 0 ) );
    EXPECT_EQ( mat3.size(), dim1 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( mat3.get( i ), wrapper1.get( i, 0 ) );
    }
}

TEST_F( _TEST_TITLE_, GetDiagonalReturnsExpectedVector ) {
    auto vec = *( wrapper1.get_diagonal( 0 ) );
    EXPECT_EQ( vec.size(), dim1 );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( vec[ i ], wrapper1.get( i, i ) );
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

// TEST_F( _TEST_TITLE_, SetRowMethodsWork ) {
//     wrapper2 = wrapper1;
//     for ( auto i = 0; i < dim1; i++ ) {
//         EXPECT_NE( wrapper2.get( 1, i ), one );
//     }

//     // constant
//     wrapper2.set_row( 1, one );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;

//     // array
//     test_t testrow[ 3 ] = { one, one, one };
//     size_t inds[ 3 ]    = { 0, 1, 2 };
//     wrapper2.set_row( 1, 3, inds, testrow );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;

//     // vector with indices
//     std::vector<test_t> testrow_v( testrow, testrow + dim1 );
//     std::vector<size_t> inds_v( inds, inds + dim1 );
//     wrapper2.set_row( 1, inds_v, testrow_v );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;

//     // vector without indices
//     wrapper2.set_row( 1, testrow_v );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;

//     // initializer list with indices
//     wrapper2.set_row( 1, { 0, 1, 2 }, { one, one, one } );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;

//     // initializer list without indices
//     wrapper2.set_row( 1, { one, one, one } );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;

//     // abstract_vector
//     auto avec = *( wrapper1.get_row_vector( 0 ) );
//     wrapper2.set_row( 1, avec );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;

//     // Vector
//     Vector<test_t> aVec( *( wrapper1.get_row_vector( 0 ) ) );
//     wrapper2.set_row( 1, aVec );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;

//     // row from another matrix
//     wrapper2.set_row( 1, wrapper1, 0 );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;

//     // row matrix
//     wrapper2.set_row( 1, *wrapper1.get_row( 0 ) );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;
//     wrapper2.set_row( 1, wrapper1.get_row( 0 ) );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( 1, i ), one );
//     }
//     wrapper2 = wrapper1;
// }

// TEST_F( _TEST_TITLE_, SetColumnMethodsWork ) {
//     wrapper2       = wrapper1;
//     test_t testval = one * 5.0;
//     for ( auto i = 0; i < dim1; i++ ) {
//         EXPECT_NE( wrapper2.get( i, 1 ), testval );
//     }

//     // constant
//     wrapper2.set_column( 1, testval );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), testval );
//     }
//     wrapper2 = wrapper1;

//     // array
//     test_t testcol[ 3 ] = { testval, testval, testval };
//     size_t inds[ 3 ]    = { 0, 1, 2 };
//     wrapper2.set_column( 1, 3, inds, testcol );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), testval );
//     }
//     wrapper2 = wrapper1;

//     // vector with indices
//     std::vector<test_t> testcol_v( testcol, testcol + dim1 );
//     std::vector<size_t> inds_v( inds, inds + dim1 );
//     wrapper2.set_column( 1, inds_v, testcol_v );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), testval );
//     }
//     wrapper2 = wrapper1;

//     // vector without indices
//     wrapper2.set_column( 1, testcol_v );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), testval );
//     }
//     wrapper2 = wrapper1;

//     // initializer list with indices
//     wrapper2.set_column( 1, { 0, 1, 2 }, { testval, testval, testval } );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), testval );
//     }
//     wrapper2 = wrapper1;

//     // initializer list without indices
//     wrapper2.set_column( 1, { testval, testval, testval } );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), testval );
//     }
//     wrapper2 = wrapper1;

//     // abstract_vector
//     auto avec = *( wrapper1.get_column_vector( 0 ) );
//     wrapper2.set_column( 1, avec );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
//     }
//     wrapper2 = wrapper1;

//     // Vector
//     Vector<test_t> aVec( *( wrapper1.get_column_vector( 0 ) ) );
//     wrapper2.set_column( 1, aVec );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
//     }
//     wrapper2 = wrapper1;

//     // row from another matrix
//     wrapper2.set_column( 1, wrapper1, 0 );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
//     }
//     wrapper2 = wrapper1;

//     // row matrix
//     wrapper2.set_column( 1, *wrapper1.get_column( 0 ) );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
//     }
//     wrapper2 = wrapper1;
//     wrapper2.set_column( 1, wrapper1.get_column( 0 ) );
//     for ( auto i = 0; i < dim1; i++ ) {
//         _TEST_EQ_( wrapper2.get( i, 1 ), wrapper1.get( i, 0 ) );
//     }
//     wrapper2 = wrapper1;
// }

TEST_F( _TEST_TITLE_, SetDiagonalMethodsWork ) {
    wrapper2 = wrapper1;
    wrapper2.clear().resize( dim1, dim1 );
    for ( size_t i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, i ), zero );
    }
    // NCPA_DEBUG << "Cleared matrix" << endl;

    // constant
    wrapper2.set_diagonal( one );
    EXPECT_TRUE( wrapper2.is_identity() );
    wrapper2.zero();
    EXPECT_FALSE( wrapper2.is_identity() );
    // NCPA_DEBUG << "Constant method passed" << endl;

    // array
    test_t testrow[ 5 ] = { one, one, one, one, one };
    wrapper2.set_diagonal( dim1, testrow );
    EXPECT_TRUE( wrapper2.is_identity() );
    wrapper2.zero();
    EXPECT_FALSE( wrapper2.is_identity() );
    // NCPA_DEBUG << "Array method passed" << endl;

    // vector
    vector<test_t> testrow_v( testrow, testrow + dim1 );
    wrapper2.set_diagonal( testrow_v );
    EXPECT_TRUE( wrapper2.is_identity() );
    wrapper2.zero();
    EXPECT_FALSE( wrapper2.is_identity() );
    // NCPA_DEBUG << "std::vector method passed" << endl;

    // init list
    wrapper2.set_diagonal( { one, one, one, one, one } );
    EXPECT_TRUE( wrapper2.is_identity() );
    wrapper2.zero();
    EXPECT_FALSE( wrapper2.is_identity() );
    // NCPA_DEBUG << "Init list method passed" << endl;

    // abstract_vector
    auto avec = *( wrapper1.get_column( 0 ) );
    wrapper2.set_diagonal( avec );
    for ( auto i = 0; i < dim1; i++ ) {
        _TEST_EQ_( wrapper2.get( i, i ), wrapper1.get( i, 0 ) );
    }
    wrapper2.zero();
    // NCPA_DEBUG << "abstract vector method passed" << endl;

    // Vector
    Vector<test_t> aVec( *( wrapper1.get_column( 0 ) ) );
    aVec.resize( dim1 - 2 );
    wrapper2.set_diagonal( aVec, -2 );
    for ( auto i = 0; i < dim1-2; i++ ) {
        _TEST_EQ_( wrapper2.get( i + 2, i ), wrapper1.get( i, 0 ) );
    }
    // NCPA_DEBUG << "Vector method passed" << endl;
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
    for ( int i = 0; i < wrapper1.rows(); i++ ) {
        for ( int j = 0; j < wrapper1.columns(); j++ ) {
            if ( abs( i - j ) > 1 ) {
                _TEST_EQ_( wrapper1.get( i, j ), zero );
            } else {
                _TEST_EQ_( wrapper1.get( i, j ),
                           ( testval + wrapper2.get( i, j ) ) );
            }
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
    for ( int i = 0; i < wrapper1.rows(); i++ ) {
        for ( int j = 0; j < wrapper1.columns(); j++ ) {
            if ( abs( i - j ) > 1 ) {
                _TEST_EQ_( wrapper1.get( i, j ), zero );
            } else {
                _TEST_EQ_( wrapper1.get( i, j ),
                           ( wrapper2.get( i, j ) - testval ) );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, ScaleWorksWithMatrix ) {
    wrapper2 = wrapper1;
    wrapper1.scale( wrapper1 );
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ),
                       ( wrapper2.get( i, j ) * wrapper2.get( i, j ) ) );
        }
    }
}

TEST_F( _TEST_TITLE_, ScaleWorksWithScalar ) {
    wrapper2 = wrapper1;
    wrapper1.scale( testval );
    for ( size_t i = 0; i < wrapper1.rows(); i++ ) {
        for ( size_t j = 0; j < wrapper1.columns(); j++ ) {
            _TEST_EQ_( wrapper1.get( i, j ),
                       ( wrapper2.get( i, j ) * testval ) );
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
    EXPECT_THROW( { wrapper1.identity(); }, std::invalid_argument );
}

TEST_F( _TEST_TITLE_, MultiplyWorks ) {
    wrapper2 = Matrix<test_t>( identity.clone() );
    EXPECT_TRUE( wrapper1.multiply( wrapper2 ).equals( wrapper1 ) );

    Matrix<test_t> lhs = wrapper1;
    lhs.zero()
        .set_diagonal( { 1, 2, 3, 4, 5 } )
        .set_diagonal( { 6, 7, 8, 9 }, 1 );
    Matrix<test_t> rhs = lhs;
    lhs.set_diagonal( { 10, 11, 12, 13 }, -1 );
    // cout << "LHS: " << endl << lhs << "RHS: " << endl << rhs;
    wproduct.zero()
        .set_diagonal( { 1, 64, 86, 112, 142 } )
        .set_diagonal( { 18, 35, 56, 81 }, 1 )
        .set_diagonal( { 42, 56, 72 }, 2 )
        .set_diagonal( { 10, 22, 36, 52 }, -1 );
    // cout << "Expect product = " << endl << product2;
    auto newproduct = lhs.multiply( rhs );
    // cout << "Actual product = " << endl << *product;
    EXPECT_TRUE( wproduct.equals( newproduct ) );
}

TEST_F( _TEST_TITLE_, MultiplyOperatorWorks ) {
    wrapper2 = Matrix<test_t>( identity.clone() );
    EXPECT_TRUE( wrapper1 * wrapper2 == wrapper1 );

    Matrix<test_t> lhs = wrapper1;
    lhs.zero()
        .set_diagonal( { 1, 2, 3, 4, 5 } )
        .set_diagonal( { 6, 7, 8, 9 }, 1 );
    Matrix<test_t> rhs = lhs;
    lhs.set_diagonal( { 10, 11, 12, 13 }, -1 );
    // cout << "LHS: " << endl << lhs << "RHS: " << endl << rhs;
    wproduct.zero()
        .set_diagonal( { 1, 64, 86, 112, 142 } )
        .set_diagonal( { 18, 35, 56, 81 }, 1 )
        .set_diagonal( { 42, 56, 72 }, 2 )
        .set_diagonal( { 10, 22, 36, 52 }, -1 );
    // cout << "Expect product = " << endl << product2;
    auto newproduct = lhs * rhs;
    // cout << "Actual product = " << endl << *product;
    EXPECT_TRUE( wproduct == newproduct );
}

TEST_F( _TEST_TITLE_, OtherBinaryOperatorsWork ) {
    EXPECT_TRUE( wrapper1 + wrapper1 == wrapper1 * 2.0 );
    EXPECT_TRUE( wrapper1 - wrapper1 == wrapper1 * 0.0 );
    wrapper1.set( 1.0 );
    wrapper2  = wrapper1;
    wrapper1 *= 3.0;
    EXPECT_TRUE( wrapper1 == wrapper2 * 3.0 );
}
