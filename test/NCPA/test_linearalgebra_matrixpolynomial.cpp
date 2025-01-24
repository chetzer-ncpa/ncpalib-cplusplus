#include "NCPA/gtest.hpp"
#include "NCPA/linearalgebra.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <vector>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::linear;

typedef double test_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraLibraryMatrixPolynomialTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            densemat = MatrixFactory<test_t>::build( matrix_t::DENSE );
            bdmat    = MatrixFactory<test_t>::build( matrix_t::BAND_DIAGONAL );
            symmat   = MatrixFactory<test_t>::build( matrix_t::SYMMETRIC );
            findiffmat = MatrixFactory<test_t>::build( matrix_t::SYMMETRIC );
        }  // void TearDown() override {}

        // declare stuff here
        Matrix<test_t> densemat, bdmat, symmat, findiffmat;
        MatrixPolynomial<test_t> poly;
        vector<Matrix<test_t>> control;

        const test_t zero = NCPA::math::zero<test_t>(),
                     one  = NCPA::math::one<test_t>();
        size_t order      = 5;
        test_t testval    = -4.2;
};

void setup_mat( Matrix<test_t>& mat, bool setsub = true ) {
    mat.resize( 8, 8 )
        .set_diagonal( { 1, 2, 3, 4, 5, 6, 7, 8 } )
        .set_diagonal( { 1, -1, 1, -1, 1, -1, 1 }, 1 );
    if ( setsub ) {
        mat.set_diagonal( { 2, 0, -2, 0, 2, 0, -2 }, -1 );
    }
}

void setup_fd_mat( Matrix<test_t>& mat ) {
    mat.resize( 8, 8 )
        .set_diagonal( { 1, 2, 3, 4, 5, 6, 7, 8 } )
        .set_diagonal( NCPA::math::one<test_t>(), 1 );
}

vector<Matrix<test_t>> make_polynomial_vector( Matrix<test_t>& mat,
                                               size_t order ) {
    vector<Matrix<test_t>> polyv( order );
    polyv[ 0 ] = mat;
    for ( size_t i = 1; i < order; i++ ) {
        polyv[ i ] = polyv[ i - 1 ] * polyv[ 0 ];
    }
    return polyv;
}

void print_comparison( const vector<Matrix<test_t>>& a,
                       const vector<Matrix<test_t>>& b ) {
    for ( size_t i = 0; i < a.size(); i++ ) {
        cout << "Power " << i + 1 << ":" << endl
             << a[ i ] << endl
             << "- vs -" << endl
             << b[ i ] << endl
             << endl;
    }
}

TEST_F( _TEST_TITLE_, DenseMatPolyWorks ) {
    setup_mat( densemat );
    control = make_polynomial_vector( densemat, order );
    poly    = MatrixPolynomial<test_t>();
    poly.set_base( densemat )
        .set_order( order )
        .set_algorithm( matrix_polynomial_algorithm_t::MULTIPLY );
    // print_comparison( poly.vector(), control );
    EXPECT_TRUE( poly.compute().vector() == control );
}

TEST_F( _TEST_TITLE_, BandDiagonalMatPolyWorks ) {
    setup_mat( bdmat );
    control = make_polynomial_vector( bdmat, order );
    poly    = MatrixPolynomial<test_t>();
    poly.set_base( bdmat ).set_order( order ).set_algorithm(
        matrix_polynomial_algorithm_t::MULTIPLY );
    // print_comparison( poly.vector(), control );
    EXPECT_TRUE( poly.compute().vector() == control );
}

TEST_F( _TEST_TITLE_, SymmetricMatPolyWorks ) {
    setup_mat( symmat, false );
    control = make_polynomial_vector( symmat, order );
    poly    = MatrixPolynomial<test_t>();
    poly.set_base( symmat ).set_order( order ).set_algorithm(
        matrix_polynomial_algorithm_t::SYMMETRIC_REFLECTED );
    poly.compute();
    EXPECT_TRUE( poly.vector() == control );
}

TEST_F( _TEST_TITLE_, FiniteDifferenceMatPolyWorks ) {
    setup_fd_mat( findiffmat );
    control = make_polynomial_vector( findiffmat, order );
    poly    = MatrixPolynomial<test_t>();
    poly.set_base( findiffmat ).set_order( order ).set_algorithm(
        matrix_polynomial_algorithm_t::FINITE_DIFFERENCE_REFLECTED );
    poly.compute();
    // print_comparison( poly.vector(), control );
    EXPECT_TRUE( poly.vector() == control );
}

TEST_F( _TEST_TITLE_, ScaleWorksWithScalar ) {
    setup_mat( bdmat );
    control = make_polynomial_vector( bdmat, order );
    poly    = MatrixPolynomial<test_t>();
    poly.set_base( bdmat ).set_order( order );
    poly.compute().scale( testval );
    for ( size_t i = 0; i < control.size(); i++ ) {
        EXPECT_TRUE( poly.vector()[ i ] == control[ i ] * testval );
    }
}

TEST_F( _TEST_TITLE_, ScaleWorksWithVector ) {
    setup_mat( bdmat );
    control = make_polynomial_vector( bdmat, order );
    poly    = MatrixPolynomial<test_t>();
    poly.set_base( bdmat ).set_order( order );
    vector<test_t> factors( control.size() );
    factors[ 0 ] = testval;
    for ( size_t i = 1; i < control.size(); i++ ) {
        factors[ i ] = -factors[ i - 1 ];
    }
    poly.compute().scale( factors );
    for ( size_t i = 0; i < control.size(); i++ ) {
        EXPECT_TRUE( poly.vector()[ i ] == control[ i ] * factors[ i ] );
    }
}
