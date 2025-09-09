#define NCPA_DEBUG_ON
#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/linearalgebra.hpp"
#include "NCPA/linearalgebra/tridiagonal_eigensolver.hpp"
#include "NCPA/logging.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::linear;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraLibraryTridiagonalEigensystemTest

typedef double test_t;

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            testmat1 = MatrixFactory<double>::build( matrix_t::BAND_DIAGONAL );
            testmat2 = MatrixFactory<double>::build( matrix_t::BAND_DIAGONAL );
            testmat3 = MatrixFactory<double>::build( matrix_t::BAND_DIAGONAL );

            testmat1.identity( 10, 10 );
            testmat2.resize( 10, 10 ).set_diagonal(
                { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 } );
        }  // void TearDown() override {}

        // declare stuff here
        Matrix<test_t> testmat1, testmat2, testmat3;
        TridiagonalEigensolver<double> solver;
};

TEST_F( _TEST_TITLE_, WorksForIdentity ) {
    solver.set( testmat1 );
    solver.solve();
    auto eigs = solver.eigenvalues();
    for (auto it = eigs.begin(); it != eigs.end(); ++it) {
        cout << *it << endl;
    }
}
