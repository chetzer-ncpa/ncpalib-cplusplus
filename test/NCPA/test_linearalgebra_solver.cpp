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
#define _TEST_TITLE_    NCPALinearAlgebraSolverTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            dmat   = MatrixFactory<test_t>::build( family_t::NCPA_DENSE );
            inmat  = MatrixFactory<test_t>::build( family_t::NCPA_DENSE );
            cmat   = MatrixFactory<ctest_t>::build( family_t::NCPA_DENSE );
            cinmat = MatrixFactory<ctest_t>::build( family_t::NCPA_DENSE );
            dsmat  = MatrixFactory<test_t>::build( family_t::NCPA_SPARSE );
            csmat  = MatrixFactory<ctest_t>::build( family_t::NCPA_SPARSE );

        }  // void TearDown() override {}

        // declare stuff here
        const test_t zero   = NCPA::math::zero<test_t>(),
                     one    = NCPA::math::one<test_t>();
        const ctest_t czero = NCPA::math::zero<ctest_t>(),
                      cone  = NCPA::math::one<ctest_t>();

        Matrix<test_t> dmat, inmat, dsmat;
        Matrix<ctest_t> cmat, cinmat, csmat;

        Solver<test_t> solver;
        Solver<ctest_t> csolver;
};

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForTrivialDenseCase ) {
    dmat.identity( 4, 4 );
    inmat.resize( 4, 1 );
    for ( size_t i = 0; i < 4; i++ ) {
        inmat[ i ][ 0 ] = (test_t)( i + 1 );
    }
    solver = SolverFactory<test_t>::build( solver_t::BASIC );
    solver.set_system_matrix( dmat );
    // std::cout << "Solving..." << endl;
    Matrix<test_t> solution = solver.solve( inmat );
    // cout << "Done" << endl;
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution[ i ][ 0 ], inmat[ i ][ 0 ] );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForComplexTrivialDenseCase ) {
    cmat.identity( 4, 4 );
    cinmat.resize( 4, 1 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinmat[ i ][ 0 ] = cone * (test_t)( i + 1 );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::BASIC );
    csolver.set_system_matrix( cmat );
    Matrix<ctest_t> solution = csolver.solve( cinmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution[ i ][ 0 ], cinmat[ i ][ 0 ] );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForRandomDenseCase ) {
    Vector<test_t> expected
        = VectorFactory<test_t>::build( family_t::NCPA_DENSE );
    expected.set( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    inmat.resize( 4, 1 ).zero();
    dmat.resize( 4, 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        dmat.set_row( i, NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
        inmat[ i ][ 0 ] = dmat.get_row_vector( i )->dot( expected );
    }
    solver = SolverFactory<test_t>::build( solver_t::BASIC );
    solver.set_system_matrix( dmat );
    Matrix<test_t> solution = solver.solve( inmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution[ i ][ 0 ], expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForRandomComplexDenseCase ) {
    Vector<ctest_t> expected
        = VectorFactory<ctest_t>::build( family_t::NCPA_DENSE );
    expected.set( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cinmat.resize( 4, 1 ).zero();
    cmat.resize( 4, 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        cmat.set_row( i, NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
        cinmat[ i ][ 0 ] = cmat.get_row_vector( i )->dot( expected );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::BASIC );
    csolver.set_system_matrix( cmat );
    Matrix<ctest_t> solution = csolver.solve( cinmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution[ i ][ 0 ].real(), expected[ i ].real(),
                     1.0e-10 );
        EXPECT_NEAR( solution[ i ][ 0 ].imag(), expected[ i ].imag(),
                     1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForTrivialSparseCase ) {
    dsmat.identity( 4, 4 );
    inmat.resize( 4, 1 );
    for ( size_t i = 0; i < 4; i++ ) {
        inmat[ i ][ 0 ] = (test_t)( i + 1 );
    }
    solver = SolverFactory<test_t>::build( solver_t::BASIC );
    solver.set_system_matrix( dsmat );
    Matrix<test_t> solution = solver.solve( inmat );
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution[ i ][ 0 ], inmat[ i ][ 0 ] );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForComplexTrivialSparseCase ) {
    csmat.identity( 4, 4 );
    cinmat.resize( 4, 1 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinmat[ i ][ 0 ] = cone * (test_t)( i + 1 );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::BASIC );
    csolver.set_system_matrix( csmat );
    Matrix<ctest_t> solution = csolver.solve( cinmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution[ i ][ 0 ], cinmat[ i ][ 0 ] );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForRandomSparseCase ) {
    Vector<test_t> expected
        = VectorFactory<test_t>::build( family_t::NCPA_DENSE );
    expected.set( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    inmat.resize( 4, 1 ).zero();
    dsmat.resize( 4, 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        dsmat.set_row( i, NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
        inmat[ i ][ 0 ] = dsmat.get_row_vector( i )->dot( expected );
    }
    solver = SolverFactory<test_t>::build( solver_t::BASIC );
    solver.set_system_matrix( dsmat );
    Matrix<test_t> solution = solver.solve( inmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution[ i ][ 0 ], expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, BasicSolverIsCorrectForRandomComplexSparseCase ) {
    Vector<ctest_t> expected
        = VectorFactory<ctest_t>::build( family_t::NCPA_DENSE );
    expected.set( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cinmat.resize( 4, 1 ).zero();
    csmat.resize( 4, 4 );
    for ( size_t i = 0; i < 4; i++ ) {
        csmat.set_row( i,
                       NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
        cinmat[ i ][ 0 ] = csmat.get_row_vector( i )->dot( expected );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::BASIC );
    csolver.set_system_matrix( csmat );
    Matrix<ctest_t> solution = csolver.solve( cinmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution[ i ][ 0 ].real(), expected[ i ].real(),
                     1.0e-10 );
        EXPECT_NEAR( solution[ i ][ 0 ].imag(), expected[ i ].imag(),
                     1.0e-10 );
    }
}









TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForTrivialDenseCase ) {
    dmat.identity( 4, 4 );
    inmat.resize( 4, 1 );
    for ( size_t i = 0; i < 4; i++ ) {
        inmat[ i ][ 0 ] = (test_t)( i + 1 );
    }
    solver = SolverFactory<test_t>::build( solver_t::TRIDIAGONAL );
    solver.set_system_matrix( dmat );
    Matrix<test_t> solution = solver.solve( inmat );
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution[ i ][ 0 ], inmat[ i ][ 0 ] );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForComplexTrivialDenseCase ) {
    cmat.identity( 4, 4 );
    cinmat.resize( 4, 1 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinmat[ i ][ 0 ] = cone * (test_t)( i + 1 );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::TRIDIAGONAL );
    csolver.set_system_matrix( cmat );
    Matrix<ctest_t> solution = csolver.solve( cinmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution[ i ][ 0 ], cinmat[ i ][ 0 ] );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForRandomDenseCase ) {
    Vector<test_t> expected = VectorFactory<test_t>::build( family_t::NCPA_DENSE );
    expected.set( 
        NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    inmat.resize( 4, 1 ).zero();
    dmat.resize( 4, 4 );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ), 1 );
    dmat.set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ), -1 );
    for ( size_t i = 0; i < 4; i++ ) {
        inmat[ i ][ 0 ] = dmat.get_row_vector( i )->dot( expected );
    }
    solver = SolverFactory<test_t>::build( solver_t::TRIDIAGONAL );
    solver.set_system_matrix( dmat );
    Matrix<test_t> solution = solver.solve( inmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution[ i ][ 0 ], expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForRandomComplexDenseCase ) {
    Vector<ctest_t> expected = VectorFactory<ctest_t>::build( family_t::NCPA_DENSE );
    expected.set( 
        NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cinmat.resize( 4, 1 ).zero();
    cmat.resize( 4, 4 );
    cmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 3, -5.0, 5.0 ), 1 );
    cmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 3, -5.0, 5.0 ), -1 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinmat[ i ][ 0 ] = cmat.get_row_vector( i )->dot( expected );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::TRIDIAGONAL );
    csolver.set_system_matrix( cmat );
    Matrix<ctest_t> solution = csolver.solve( cinmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution[ i ][ 0 ].real(), expected[ i ].real(), 1.0e-10 );
        EXPECT_NEAR( solution[ i ][ 0 ].imag(), expected[ i ].imag(), 1.0e-10 );
    }
}








TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForTrivialSparseCase ) {
    dsmat.identity( 4, 4 );
    inmat.resize( 4, 1 );
    for ( size_t i = 0; i < 4; i++ ) {
        inmat[ i ][ 0 ] = (test_t)( i + 1 );
    }
    solver = SolverFactory<test_t>::build( solver_t::TRIDIAGONAL );
    solver.set_system_matrix( dsmat );
    Matrix<test_t> solution = solver.solve( inmat );
    for ( size_t i = 0; i < 4; i++ ) {
        _TEST_EQ_( solution[ i ][ 0 ], inmat[ i ][ 0 ] );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForComplexTrivialSparseCase ) {
    csmat.identity( 4, 4 );
    cinmat.resize( 4, 1 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinmat[ i ][ 0 ] = cone * (test_t)( i + 1 );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::TRIDIAGONAL );
    csolver.set_system_matrix( csmat );
    Matrix<ctest_t> solution = csolver.solve( cinmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_COMPLEX_DOUBLE_EQ( solution[ i ][ 0 ], cinmat[ i ][ 0 ] );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForRandomSparseCase ) {
    Vector<test_t> expected = VectorFactory<test_t>::build( family_t::NCPA_DENSE );
    expected.set( 
        NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    inmat.resize( 4, 1 ).zero();
    dsmat.resize( 4, 4 );
    dsmat.set_diagonal( NCPA::math::random_numbers<test_t>( 4, -5.0, 5.0 ) );
    dsmat.set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ), 1 );
    dsmat.set_diagonal( NCPA::math::random_numbers<test_t>( 3, -5.0, 5.0 ), -1 );
    for ( size_t i = 0; i < 4; i++ ) {
        inmat[ i ][ 0 ] = dsmat.get_row_vector( i )->dot( expected );
    }
    solver = SolverFactory<test_t>::build( solver_t::TRIDIAGONAL );
    solver.set_system_matrix( dsmat );
    Matrix<test_t> solution = solver.solve( inmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution[ i ][ 0 ], expected[ i ], 1.0e-10 );
    }
}

TEST_F( _TEST_TITLE_, TridiagonalSolverIsCorrectForRandomComplexSparseCase ) {
    Vector<ctest_t> expected = VectorFactory<ctest_t>::build( family_t::NCPA_DENSE );
    expected.set( 
        NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    cinmat.resize( 4, 1 ).zero();
    csmat.resize( 4, 4 );
    csmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 4, -5.0, 5.0 ) );
    csmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 3, -5.0, 5.0 ), 1 );
    csmat.set_diagonal( NCPA::math::random_numbers<ctest_t>( 3, -5.0, 5.0 ), -1 );
    for ( size_t i = 0; i < 4; i++ ) {
        cinmat[ i ][ 0 ] = csmat.get_row_vector( i )->dot( expected );
    }
    csolver = SolverFactory<ctest_t>::build( solver_t::TRIDIAGONAL );
    csolver.set_system_matrix( csmat );
    Matrix<ctest_t> solution = csolver.solve( cinmat );
    for ( size_t i = 0; i < 4; i++ ) {
        EXPECT_NEAR( solution[ i ][ 0 ].real(), expected[ i ].real(), 1.0e-10 );
        EXPECT_NEAR( solution[ i ][ 0 ].imag(), expected[ i ].imag(), 1.0e-10 );
    }
}