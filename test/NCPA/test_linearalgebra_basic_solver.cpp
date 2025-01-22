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
#define _TEST_TITLE_    NCPALinearAlgebraBasicSolverTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            dmat   = MatrixFactory<test_t>::build( matrix_t::DENSE );
            inmat  = MatrixFactory<test_t>::build( matrix_t::DENSE );
            cmat   = MatrixFactory<ctest_t>::build( matrix_t::DENSE );
            cinmat = MatrixFactory<ctest_t>::build( matrix_t::DENSE );
        }  // void TearDown() override {}

        // declare stuff here
        const test_t zero   = NCPA::math::zero<test_t>(),
                     one    = NCPA::math::one<test_t>();
        const ctest_t czero = NCPA::math::zero<ctest_t>(),
                      cone  = NCPA::math::one<ctest_t>();

        Matrix<test_t> dmat, inmat;
        Matrix<ctest_t> cmat, cinmat;

        basic_linear_system_solver<test_t> solver;
        basic_linear_system_solver<ctest_t> csolver;
};

TEST_F( _TEST_TITLE_, SolverIsCorrectForTrivialCase ) {
    // cout << "identity():" << endl;
    dmat.identity( 4, 4 );
    // cout << "resize():" << endl;
    inmat.resize( 4, 1 );
    // cout << "for loop:" << endl;
    for ( size_t i = 0; i < 4; i++ ) {
        inmat.set( i, 0, (test_t)( i + 1 ) );
        // inmat[ i ][ 0 ] = (test_t)( i + 1 );
    }
    // cout << "set_system_matrix()" << endl;
    solver.set_system_matrix( dmat );
    // cout << "solve()" << endl;
    Vector<test_t> solution = solver.solve( inmat );
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution.get( i ), inmat.get( i, 0 ) );
        // _TEST_EQ_( solution[ i ][ 0 ], inmat[ i ][ 0 ] );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForComplexTrivialCase ) {
    cmat.identity( 4, 4 );
    cinmat.resize( 4, 1 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinmat.set( i, 0, cone * (test_t)( i + 1 ) );
    }
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( cinmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution.get( i ), cinmat.get( i, 0 ) );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForRandomCase ) {
    Vector<test_t> expected = VectorFactory<test_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    inmat.resize( 4, 1 ).zero();
    dmat.resize( 4, 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        dmat.set_row( i, NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
        inmat.set( i, 0, dmat.get_row( i )->dot( expected ) );
        // inmat[ i ][ 0 ] = dmat.get_row_vector( i )->dot( expected );
    }
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( inmat );
    // cout << "system = " << dmat << endl
    //      << "inmat = " << inmat << endl
    //      << "expected = " << expected << endl
    //      << "got = " << solution << endl;
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), expected.get( i ), 1.0e-10 );
        // EXPECT_NEAR( solution[ i ][ 0 ], expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, SolverIsCorrectForRandomComplexCase ) {
    Vector<ctest_t> expected
        = VectorFactory<ctest_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cinmat.resize( 4, 1 ).zero();
    cmat.resize( 4, 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        cmat.set_row( i, NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
        cinmat.set( i, 0, cmat.get_row( i )->dot( expected ) );
        // cinmat[ i ][ 0 ] = cmat.get_row_vector( i )->dot( expected );
    }
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( cinmat );
    // cout << "system = " << dmat << endl
    //      << "inmat = " << inmat << endl
    //      << "expected = " << expected << endl
    //      << "got = " << solution << endl;
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ).real(), expected[ i ].real(), 1.0e-10 );
        EXPECT_NEAR( solution.get( i ).imag(), expected[ i ].imag(), 1.0e-10 );

        // EXPECT_NEAR( solution[ i ][ 0 ].real(), expected[ i
        // ].real(), 1.0e-10 ); EXPECT_NEAR( solution[ i ][ 0 ].imag(),
        // expected[ i ].imag(), 1.0e-10 );
    }
}
