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
            dmat   = MatrixFactory<test_t>::build( family_t::NCPA_BAND_DIAGONAL );
            invec  = VectorFactory<test_t>::build( family_t::NCPA_DENSE );
            outvec  = VectorFactory<test_t>::build( family_t::NCPA_DENSE );
            cmat   = MatrixFactory<ctest_t>::build( family_t::NCPA_BAND_DIAGONAL );
            cinvec = VectorFactory<ctest_t>::build( family_t::NCPA_DENSE );
            coutvec = VectorFactory<ctest_t>::build( family_t::NCPA_DENSE );
        }  // void TearDown() override {}

        // declare stuff here
        const test_t zero   = NCPA::math::zero<test_t>(),
                     one    = NCPA::math::one<test_t>();
        const ctest_t czero = NCPA::math::zero<ctest_t>(),
                      cone  = NCPA::math::one<ctest_t>();

        Matrix<test_t> dmat;
        Vector<test_t> invec, outvec;
        Matrix<ctest_t> cmat;
        Vector<ctest_t> cinvec, coutvec;

        details::basic_band_diagonal_linear_system_solver<test_t> solver;
        details::basic_band_diagonal_linear_system_solver<ctest_t> csolver;
};

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
    invec.resize( 4 ).set( { 10.0, -7.0, -35.0, 34.0 } );
    dmat.zero()
        .resize( 4, 4 )
        .set_diagonal( { -1, -2, -3, -4 } )
        .set_diagonal( { 4, 5, 6 }, -1 )
        .set_diagonal( { 2, 3, 4 }, 1 )
        .set_diagonal( { 5, 6 }, 2 );
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
    dmat.resize( 4, 4 ).zero()
        .set_diagonal( NCPA::math::random_numbers<test_t>( 4, -1.0, 1.0 ) )
        .set_diagonal(  NCPA::math::random_numbers<test_t>( 3, -1.0, 1.0 ), -1 )
        .set_diagonal(  NCPA::math::random_numbers<test_t>( 3, -1.0, 1.0 ), 1 );
    Vector<test_t> x = VectorFactory<test_t>::build( family_t::NCPA_DENSE );
    x.resize(4).zero().set( NCPA::math::random_numbers<test_t>( 4, -1.0, 1.0 ) );
    invec = dmat * x;

    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    cout << "dmat = " << endl << dmat << endl;
    cout << "b = " << invec << endl;
    cout << "expected = " << x << endl;
    cout << "solution = " << solution << endl;
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), x[ i ], 1.0e-10 );
    }
}