#define NCPA_DEBUG_ON

#include "NCPA/gtest.hpp"
#include "NCPA/linearalgebra.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::linear;

typedef double test_t;
typedef complex<double> ctest_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraBasicBandDiagonalSolverTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            dmat
                = MatrixFactory<test_t>::build( matrix_t::BAND_DIAGONAL );
            invec  = VectorFactory<test_t>::build( vector_t::DENSE );
            outvec = VectorFactory<test_t>::build( vector_t::DENSE );
            cmat   = MatrixFactory<ctest_t>::build(
                matrix_t::BAND_DIAGONAL );
            cinvec  = VectorFactory<ctest_t>::build( vector_t::DENSE );
            coutvec = VectorFactory<ctest_t>::build( vector_t::DENSE );
        }  // void TearDown() override {}

        // declare stuff here
        const test_t zero   = NCPA::math::zero<test_t>(),
                     one    = NCPA::math::one<test_t>();
        const ctest_t czero = NCPA::math::zero<ctest_t>(),
                      cone  = NCPA::math::one<ctest_t>();

        Matrix<test_t> dmat, dmatinv;
        Vector<test_t> invec, outvec;
        Matrix<ctest_t> cmat;
        Vector<ctest_t> cinvec, coutvec;

        basic_band_diagonal_linear_system_solver<test_t> solver;
        basic_band_diagonal_linear_system_solver<ctest_t> csolver;
};

TEST_F( _TEST_TITLE_, SolverIsRigorouslyCorrect2x2 ) {
    size_t n = 2;
    dmat.resize( n, n ).set_diagonal( 3 ).set_diagonal( 1, 1 ).set_diagonal(
        -1, -1 );
    dmatinv = dmat;
    dmatinv.transpose().scale( 0.1 );
    invec.resize( n ).set( 0, 2.0 ).set( 1, 4.0 );

    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    outvec                  = dmatinv * invec;
    for ( size_t i = 0; i < n; i++ ) {
        _TEST_EQ_( solution.get( i ), outvec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, SolverIsRigorouslyCorrect3x3 ) {
    size_t n = 3;
    dmat.resize( n, n ).set_diagonal( -3 ).set_diagonal( 1, 1 ).set_diagonal(
        1, -1 );
    dmatinv
        = MatrixFactory<test_t>::build( matrix_t::DENSE );
    dmatinv.resize( n, n )
        .set_row(
            0, { -0.380952380952381, -0.142857142857143, -0.047619047619048 } )
        .set_row(
            1, { -0.142857142857143, -0.428571428571429, -0.142857142857143 } )
        .set_row( 2, { -0.047619047619048, -0.142857142857143,
                       -0.380952380952381 } );

    invec.resize( n ).set( { 2.0, 4.0, 6.0 } );

    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    outvec                  = dmatinv * invec;
    for ( size_t i = 0; i < n; i++ ) {
        EXPECT_NEAR( solution.get( i ), outvec.get( i ), 1e-12 );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForTrivialCase ) {
    dmat.identity( 4, 4 );
    invec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        invec.set( i, (test_t)( i + 1 ) );
    }
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution.get( i ), invec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForComplexTrivialCase ) {
    cmat.identity( 4, 4 );
    cinvec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinvec.set( i, cone * (test_t)( i + 1 ) );
    }
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution.get( i ), cinvec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForKnownCase ) {
    invec.resize( 4 ).set( { -5.0, 17.0, -35.0, 34.0 } );
    dmat.zero()
        .resize( 4, 4 )
        .set_diagonal( { -1, -2, -3, -4 } )
        .set_diagonal( { 4, 5, 6 }, -1 )
        .set_diagonal( { 2, 3, 4 }, 1 );
    Vector<test_t> expected = invec;
    expected.set( { 1.0, -2.0, 3.0, -4.0 } );
    // std::vector<test_t> expected = { 1.0, -2.0, 3.0, -4.0 };
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    // cout << "dmat = " << endl << dmat << endl;
    // cout << "b = " << invec << endl;
    // cout << "expected = " << expected << endl;
    // cout << "solution = " << solution << endl;
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForComplexKnownCase ) {
    ctest_t I = ctest_t( 0.0, 1.0 );
    cinvec.resize( 4 ).set( { ctest_t( -1.0, -7.0 ), ctest_t( 6.0, 20.0 ),
                              ctest_t( -13.0, -39.0 ),
                              ctest_t( 2.0, 34.0 ) } );
    cmat.zero()
        .resize( 4, 4 )
        .set_diagonal( { -1.0 * I, -2.0 * I, -3.0 * I, -4.0 * I } )
        .set_diagonal( { 4.0, 5.0, 6.0 }, -1 )
        .set_diagonal( { 2.0 + I, 3.0 + I, 4.0 + I }, 1 );
    std::vector<ctest_t> expected
        = { 1.0 + 1.0 * I, -2.0 - 2.0 * I, 3.0 + 3.0 * I, -4.0 - 4.0 * I };
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ).real(), expected[ i ].real(), 1.0e-10 );
        EXPECT_NEAR( solution.get( i ).imag(), expected[ i ].imag(), 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForAsymmetricKnownCase ) {
    invec.resize( 4 ).set( { 11.0, -7.7, -38.5, 37.40 } );
    dmat.zero()
        .resize( 4, 4 )
        .set_diagonal( { -1.1, -2.2, -3.3, -4.4 } )
        .set_diagonal( { 4.4, 5.5, 6.6 }, -1 )
        .set_diagonal( { 2.2, 3.3, 4.4 }, 1 )
        .set_diagonal( { 5.5, 6.6 }, 2 );
    Vector<test_t> expected = invec;
    expected.set( { 1.0, -2.0, 3.0, -4.0 } );
    // std::vector<test_t> expected = { 1.0, -2.0, 3.0, -4.0 };
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    // cout << "dmat = " << endl << dmat << endl;
    // cout << "b = " << invec << endl;
    // cout << "expected = " << expected << endl;
    // cout << "solution = " << solution << endl;
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForRandomCase ) {
    size_t n = 5;
    Vector<test_t> expected
        = VectorFactory<test_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<test_t>( n, -5.0, 5.0 ) );
    invec.resize( n ).zero();
    dmat.resize( n, n ).zero();
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( n, -5.0, 5.0 ) );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( n - 1, -5.0, 5.0 ),
                       1 );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( n - 1, -5.0, 5.0 ),
                       -1 );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( n - 2, -5.0, 5.0 ),
                       -2 );

    invec = dmat * expected;
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    // cout << "solution = " << solution << endl;
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, AlternateSolverIsCorrectForRandomCase ) {
    size_t n           = 5;
    Vector<test_t> rhs = VectorFactory<test_t>::build( vector_t::DENSE );
    rhs.set( NCPA::math::random_numbers<test_t>( n, -5.0, 5.0 ) );
    dmat.resize( n, n ).zero();
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( n, -5.0, 5.0 ) );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( n - 1, -5.0, 5.0 ),
                       1 );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( n - 1, -5.0, 5.0 ),
                       -1 );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( n - 2, -5.0, 5.0 ),
                       -2 );
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( rhs );
    Matrix<test_t> inverse  = dmat.invert();
    Vector<test_t> prod     = inverse * rhs;
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( prod.get( i ), solution[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, AlternateSolverIsCorrectForComplexRandomCase ) {
    size_t n = 5;
    Vector<ctest_t> rhs
        = VectorFactory<ctest_t>::build( vector_t::DENSE );
    rhs.set( NCPA::math::random_numbers<ctest_t>( n, -5.0, 5.0 ) );
    cmat.resize( n, n )
        .zero()
        .set_diagonal( NCPA::math::random_numbers<ctest_t>( n, -5.0, 5.0 ) )
        .set_diagonal( NCPA::math::random_numbers<ctest_t>( n - 1, -5.0, 5.0 ),
                       1 ).set_diagonal( NCPA::math::random_numbers<ctest_t>( n - 1, -5.0, 5.0 ),
                       -1 ).set_diagonal( NCPA::math::random_numbers<ctest_t>( n - 2, -5.0, 5.0 ),
                       -2 );
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( rhs );
    Matrix<ctest_t> inverse  = cmat.invert();
    Vector<ctest_t> prod     = inverse * rhs;
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( prod.get( i ).real(), solution[ i ].real(), 1.0e-10 );
        EXPECT_NEAR( prod.get( i ).imag(), solution[ i ].imag(), 1.0e-10 );
    }
}
