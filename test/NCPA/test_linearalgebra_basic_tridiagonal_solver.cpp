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
#define _TEST_TITLE_    NCPALinearAlgebraBasicTridiagonalSolverTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            dmat   = MatrixFactory<test_t>::build( matrix_t::DENSE );
            invec  = VectorFactory<test_t>::build( vector_t::DENSE );
            cmat   = MatrixFactory<ctest_t>::build( matrix_t::DENSE );
            cinvec = VectorFactory<ctest_t>::build( vector_t::DENSE );
        }  // void TearDown() override {}

        // declare stuff here
        const test_t zero   = NCPA::math::zero<test_t>(),
                     one    = NCPA::math::one<test_t>();
        const ctest_t czero = NCPA::math::zero<ctest_t>(),
                      cone  = NCPA::math::one<ctest_t>();

        Matrix<test_t> dmat;
        Vector<test_t> invec;
        Matrix<ctest_t> cmat;
        Vector<ctest_t> cinvec;

        basic_tridiagonal_linear_system_solver<test_t> solver;
        basic_tridiagonal_linear_system_solver<ctest_t> csolver;
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

TEST_F( _TEST_TITLE_, SolverIsCorrectForRandomCase ) {
    Vector<test_t> expected
        = VectorFactory<test_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    invec.resize( 4 ).zero();
    dmat.resize( 4, 4 );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ), 1 );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ),
                       -1 );
    for ( size_t i = 0; i < 4; i++ ) {
        invec.set( i, dmat.get_row( i )->dot( expected ) );
    }
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForRandomComplexCase ) {
    Vector<ctest_t> expected
        = VectorFactory<ctest_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cinvec.resize( 4 ).zero();
    cmat.resize( 4, 4 );
    cmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 3, -5.0, 5.0 ),
                       1 );
    cmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 3, -5.0, 5.0 ),
                       -1 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinvec.set( i, cmat.get_row( i )->dot( expected ) );
    }
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ).real(), expected[ i ].real(),
                     1.0e-10 );
        EXPECT_NEAR( solution.get( i ).imag(), expected[ i ].imag(),
                     1.0e-10 );
    }
}
