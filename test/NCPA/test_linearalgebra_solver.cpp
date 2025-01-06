#define NCPA_DEBUG_ON

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

typedef double test_t;
typedef complex<double> ctest_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraSolverTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            dmat   = MatrixFactory<test_t>::build( matrix_t::DENSE );
            invec  = VectorFactory<test_t>::build( vector_t::DENSE );
            cmat   = MatrixFactory<ctest_t>::build( matrix_t::DENSE );
            cinvec = VectorFactory<ctest_t>::build( vector_t::DENSE );
            dsmat
                = MatrixFactory<test_t>::build( matrix_t::BAND_DIAGONAL );
            csmat = MatrixFactory<ctest_t>::build(
                matrix_t::BAND_DIAGONAL );

        }  // void TearDown() override {}

        // declare stuff here
        const test_t zero   = NCPA::math::zero<test_t>(),
                     one    = NCPA::math::one<test_t>();
        const ctest_t czero = NCPA::math::zero<ctest_t>(),
                      cone  = NCPA::math::one<ctest_t>();

        Matrix<test_t> dmat, dsmat;
        Matrix<ctest_t> cmat, csmat;
        Vector<test_t> invec;
        Vector<ctest_t> cinvec;

        Solver<test_t> solver;
        Solver<ctest_t> csolver;
};

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForTrivialDenseCase ) {
    dmat.identity( 4, 4 );
    invec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        invec.set( i, (test_t)( i + 1 ) );
    }
    solver = SolverFactory<test_t>::build( solver_t::BASIC );
    solver.set_system_matrix( dmat );
    // std::cout << "Solving..." << endl;
    Vector<test_t> solution = solver.solve( invec );
    // cout << "Done" << endl;
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution.get( i ), invec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForComplexTrivialDenseCase ) {
    cmat.identity( 4, 4 );
    cinvec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinvec.set( i, cone * (test_t)( i + 1 ) );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::BASIC );
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution.get( i ), cinvec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForRandomDenseCase ) {
    Vector<test_t> expected
        = VectorFactory<test_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    invec.resize( 4 ).zero();
    dmat.resize( 4, 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        dmat.set_row( i, NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
        invec.set( i, dmat.get_row( i )->dot( expected ) );
    }
    solver = SolverFactory<test_t>::build( solver_t::BASIC );
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForRandomComplexDenseCase ) {
    Vector<ctest_t> expected
        = VectorFactory<ctest_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cinvec.resize( 4 ).zero();
    cmat.resize( 4, 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        cmat.set_row( i, NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
        cinvec.set( i, cmat.get_row( i )->dot( expected ) );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::BASIC );
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ).real(), expected[ i ].real(),
                     1.0e-10 );
        EXPECT_NEAR( solution.get( i ).imag(), expected[ i ].imag(),
                     1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForTrivialBandDiagonalCase ) {
    dsmat.identity( 4, 4 );
    invec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        invec.set( i, (test_t)( i + 1 ) );
    }
    solver = SolverFactory<test_t>::build( solver_t::BASIC );
    solver.set_system_matrix( dsmat );
    Vector<test_t> solution = solver.solve( invec );
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution.get( i ), invec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForComplexTrivialBandDiagonalCase ) {
    csmat.identity( 4, 4 );
    cinvec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinvec.set( i, cone * (test_t)( i + 1 ) );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::BASIC );
    csolver.set_system_matrix( csmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution.get( i ), cinvec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForRandomBandDiagonalCase ) {
    Vector<test_t> expected
        = VectorFactory<test_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    invec.resize( 4 ).zero();
    dsmat.resize( 4, 4 );
    dsmat.set_diagonal( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) )
        .set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ), 1 )
        .set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ),
                       -1 );
    for ( size_t i = 0; i < 4; i++ ) {
        invec.set( i, dsmat.get_row( i )->dot( expected ) );
    }
    solver = SolverFactory<test_t>::build( solver_t::TRIDIAGONAL );
    solver.set_system_matrix( dsmat );
    Vector<test_t> solution = solver.solve( invec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForRandomComplexBandDiagonalCase ) {
    Vector<ctest_t> expected
        = VectorFactory<ctest_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<ctest_t>( 5, -5.0, 5.0 ) );
    cinvec.resize( 5 ).zero();
    csmat.resize( 5, 5 );
    csmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 5, -5.0, 5.0 ) )
        .set_diagonal( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ), 1 )
        .set_diagonal( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ),
                       -1 );
    for ( size_t i = 0; i < 5; i++ ) {
        cinvec.set( i, csmat.get_row( i )->dot( expected ) );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::BASIC );
    csolver.set_system_matrix( csmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_NEAR( solution.get( i ).real(), expected[ i ].real(),
                     1.0e-10 );
        EXPECT_NEAR( solution.get( i ).imag(), expected[ i ].imag(),
                     1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForTrivialDenseCase ) {
    dmat.identity( 4, 4 );
    invec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        invec.set( i, (test_t)( i + 1 ) );
    }
    solver = SolverFactory<test_t>::build( solver_t::TRIDIAGONAL );
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution.get( i ), invec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForComplexTrivialDenseCase ) {
    cmat.identity( 4, 4 );
    cinvec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinvec.set( i, cone * (test_t)( i + 1 ) );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::TRIDIAGONAL );
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution.get( i ), cinvec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForRandomDenseCase ) {
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
    solver = SolverFactory<test_t>::build( solver_t::TRIDIAGONAL );
    solver.set_system_matrix( dmat );
    Vector<test_t> solution = solver.solve( invec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForRandomComplexDenseCase ) {
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
    csolver = SolverFactory<ctest_t>::build( solver_t::TRIDIAGONAL );
    csolver.set_system_matrix( cmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ).real(), expected[ i ].real(),
                     1.0e-10 );
        EXPECT_NEAR( solution.get( i ).imag(), expected[ i ].imag(),
                     1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForTrivialBandDiagonalCase ) {
    dsmat.identity( 4, 4 );
    invec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        invec.set( i, (test_t)( i + 1 ) );
    }
    solver = SolverFactory<test_t>::build( solver_t::TRIDIAGONAL );
    solver.set_system_matrix( dsmat );
    Vector<test_t> solution = solver.solve( invec );
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution.get( i ), invec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_,
        TridiagonalSolverIsCorrectForComplexTrivialBandDiagonalCase ) {
    csmat.identity( 4, 4 );
    cinvec.resize( 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinvec.set( i, cone * (test_t)( i + 1 ) );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::TRIDIAGONAL );
    csolver.set_system_matrix( csmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution.get( i ), cinvec.get( i ) );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForRandomBandDiagonalCase ) {
    Vector<test_t> expected
        = VectorFactory<test_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    invec.resize( 4 ).zero();
    dsmat.resize( 4, 4 );
    dsmat.set_diagonal( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    dsmat.set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ),
                        1 );
    dsmat.set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ),
                        -1 );
    for ( size_t i = 0; i < 4; i++ ) {
        invec.set( i, dsmat.get_row( i )->dot( expected ) );
    }
    solver = SolverFactory<test_t>::build( solver_t::TRIDIAGONAL );
    solver.set_system_matrix( dsmat );
    Vector<test_t> solution = solver.solve( invec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ), expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_,
        TridiagonalSolverIsCorrectForRandomComplexBandDiagonalCase ) {
    Vector<ctest_t> expected
        = VectorFactory<ctest_t>::build( vector_t::DENSE );
    expected.set( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cinvec.resize( 4 ).zero();
    csmat.resize( 4, 4 );
    csmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    csmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 3, -5.0, 5.0 ),
                        1 );
    csmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 3, -5.0, 5.0 ),
                        -1 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinvec.set( i, csmat.get_row( i )->dot( expected ) );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::TRIDIAGONAL );
    csolver.set_system_matrix( csmat );
    Vector<ctest_t> solution = csolver.solve( cinvec );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution.get( i ).real(), expected[ i ].real(),
                     1.0e-10 );
        EXPECT_NEAR( solution.get( i ).imag(), expected[ i ].imag(),
                     1.0e-10 );
    }
}
