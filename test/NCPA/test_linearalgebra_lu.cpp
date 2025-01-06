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

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraLUDecompositionTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            dmat = MatrixFactory<test_t>::build( family_t::NCPA_DENSE );
            dmat.resize( 4, 4 );
            for ( size_t r = 0; r < 4; r++ ) {
                for ( size_t c = 0; c < 4; c++ ) {
                    dmat.set( r, c,
                              NCPA::math::random_number<test_t>( 0.0, 10.0 ) );
                }
            }
            // lu.set( dmat );

        }  // void TearDown() override {}

        // declare stuff here
        const test_t zero = NCPA::math::zero<test_t>(),
                     one  = NCPA::math::one<test_t>();
        LUDecomposition<test_t> lu;
        Matrix<test_t> dmat;
};

TEST_F( _TEST_TITLE_, LUDecompositionIsCorrect ) {
    lu.decompose( dmat, false );
    Matrix<test_t> left  = lu.permutation() * dmat;
    Matrix<test_t> right = lu.lower() * lu.upper();
    for ( size_t r = 0; r < left.rows(); r++ ) {
        for ( size_t c = 0; c < left.columns(); c++ ) {
            EXPECT_NEAR( left.get( r, c ), right.get( r, c ), 1e-10 );
        }
    }
}

TEST_F( _TEST_TITLE_, LUDecompositionIsCorrectWithPivot ) {
    lu.decompose( dmat, true );
    Matrix<test_t> left  = lu.permutation() * dmat;
    Matrix<test_t> right = lu.lower() * lu.upper();
    for ( size_t r = 0; r < left.rows(); r++ ) {
        for ( size_t c = 0; c < left.columns(); c++ ) {
            EXPECT_NEAR( left.get( r, c ), right.get( r, c ), 1e-10 );
        }
    }
}

TEST_F( _TEST_TITLE_, LUDecompositionIsCorrectInsideMatrix ) {
    // lu = dmat.lu();
    Matrix<test_t> left  = dmat.lu().permutation() * dmat;
    Matrix<test_t> right = dmat.lu().lower() * dmat.lu().upper();
    for ( size_t r = 0; r < left.rows(); r++ ) {
        for ( size_t c = 0; c < left.columns(); c++ ) {
            EXPECT_NEAR( left.get( r, c ), right.get( r, c ), 1e-10 );
        }
    }
}