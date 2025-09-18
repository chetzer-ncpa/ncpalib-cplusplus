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
#define _TEST_TITLE_    TridiagonalEigensolverTest

typedef std::complex<double> test_t;

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            testmat1 = MatrixFactory<test_t>::build( matrix_t::BAND_DIAGONAL );
            testmat2 = MatrixFactory<test_t>::build( matrix_t::BAND_DIAGONAL );

            testmat1.resize( 10, 10 );
            for (size_t i = 0; i < 10; ++i) {
                testmat1.set( i, i, { 1.0 / sqrt( (double)( i + 1 ) ), 0.0 } );
                if (i > 0) {
                    testmat1.set( i, i - 1,
                                  { -1.0 / sqrt( (double)i ), 0.0 } );
                    testmat1.set( i - 1, i,
                                  { -1.0 / sqrt( (double)i ), 0.0 } );
                }
            }
            eigenvalues_expected1
                = { -0.593932543092772, -0.341923610136462, -0.176672870943261,
                    0.046033209866685,  0.305732113229965,  0.573265384622933,
                    0.817886540920661,  1.027106288226187,  1.306694219021785,
                    2.056809167576945 };
            eigenvectors_expected1 = {
                {  0.356887269873318,  0.568854233666613,  0.541947083673384,
                 0.402758309036991,  0.255394253222439,  0.144279124543804,
                 0.074410577474755,  0.035499929828542,  0.015587828046509,
                 0.005708821435262 },
                { -0.190208103035641, -0.255244744302796, -0.109674228881427,
                 0.137983196621047,  0.358983513210436,  0.479179636599229,
                 0.487264235258763,  0.410492149179755,  0.286574196527432,
                 0.145140974610012 },
                { -0.216657013556655, -0.254934430151702, -0.012231318760957,
                 0.296255455424263,  0.415059569467444,  0.247806241942316,
                 -0.099628506762764, -0.413859439662887, -0.514160271900477,
                 -0.347710561619812 },
                {  0.243659360584029,  0.232442938102282, -0.127275653068484,
                 -0.401811025061822, -0.217852457454235,  0.253810216819035,
                 0.463836526328004,  0.133198383848741, -0.380006360388433,
                 -0.468805843466926 },
                {  0.255605143018326,  0.177458442490885, -0.260749613061358,
                 -0.340012653681736,  0.168980639172486,  0.433604792200330,
                 -0.076225508384061, -0.482914074027858,  0.016170022354806,
                 0.513546661751555 },
                {  0.253388229388441,  0.108129528609152, -0.337878278961939,
                 -0.134821651339047,  0.409903750905714,  0.035199566064268,
                 -0.463254996473255,  0.201352153137857,  0.370112268823793,
                 -0.479971596934458 },
                {  0.235708382390554,  0.042925668851139, -0.340066994222292,
                 0.089106073284784,  0.336024298490449, -0.378137325873824,
                 0.011328590977921,  0.395249369576036, -0.531204606873975,
                 0.352965424207465 },
                {  0.230659617851977, -0.006252326083638, -0.323372484426843,
                 0.259564758918453,  0.099761948609050, -0.419561494457656,
                 0.526723789638827, -0.451452943999678,  0.296969458636632,
                 -0.139249979360598 },
                {  0.316015568413053, -0.096920147953167, -0.364730578307326,
                 0.579452287998573, -0.513727026719945,  0.339462580685341,
                 -0.184307110518667,  0.086215716095288, -0.035395443790199,
                 0.011912045307142 },
                { -0.626078656697545,  0.661645664022225, -0.377518759656310,
                 0.157043925234755, -0.053053730012061,  0.015368688164590,
                 -0.003943356836809,  0.000915543811942, -0.000195038429539,
                 0.000037351203332 }
            };

            testmat2.resize( 5, 5 )
                .set_diagonal( test_t( 2.0, 0.0 ), 0 )
                .set_diagonal( test_t( -1.0, 0.0 ), 1 )
                .set_diagonal( test_t( -1.0, 0.0 ), -1 );

            eigenvalues_expected2
                = { 0.267949192431843, 1.0, 2.0, 3.0, 3.732050807569423 };
            eigenvectors_expected2 = {
                { -0.288675134594813, -0.5, -0.577350269189626, -0.5,
                 -0.288675134594813                                      },
                {               -0.5, -0.5,                0.0,  0.5, 0.5 },
                {  0.577350269189626,  0.0, -0.577350269189626,  0.0,
                 0.577350269189626                                       },
                {               -0.5,  0.5,                0.0, -0.5, 0.5 },
                {  0.288675134594813, -0.5,  0.577350269189626, -0.5,
                 0.288675134594813                                       }
            };
        }

        void TearDown() override {}

        // declare stuff here
        Matrix<test_t> testmat1, testmat2;
        TridiagonalEigensolver solver;
        vector<double> eigenvalues_expected1, eigenvalues_expected2;
        vector<vector<complex<double>>> eigenvectors_expected1,
            eigenvectors_expected2;
};

TEST_F( _TEST_TITLE_, EigenvalueCountCorrectForSimpleCase ) {
    solver.set( testmat2 ).set_range( 0.0, 5.0 ).solve();
    EXPECT_EQ( solver.count(), 5 );
}

TEST_F( _TEST_TITLE_, EigenvaluesCorrectForSimpleCase ) {
    solver.set( testmat2 ).set_range( 0.0, 5.0 ).solve();
    ASSERT_EQ( solver.count(), 5 );

    vector<double> eigs = solver.eigenvalues();
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_NEAR( eigs[ i ], eigenvalues_expected2[ i ], 1e-8 );
    }
}

TEST_F( _TEST_TITLE_, EigenvectorsCorrectForSimpleCase ) {
    solver.set( testmat2 ).set_range( 0.0, 5.0 ).solve();
    ASSERT_EQ( solver.count(), 5 );
    std::vector<std::vector<std::complex<double>>> vecs
        = solver.eigenvectors();
    double factor = 1.0;
    for (size_t i = 0; i < 5; ++i) {
        if (vecs[ i ][ 0 ].real() * eigenvectors_expected2[ i ][ 0 ].real()
            < 0.0) {
            factor = -1.0;
        } else {
            factor = 1.0;
        }
        for (size_t j = 0; j < 5; ++j) {
            EXPECT_NEAR( vecs[ i ][ j ].real(),
                         factor * eigenvectors_expected2[ i ][ j ].real(),
                         1e-8 );
        }
    }
}

TEST_F( _TEST_TITLE_, CorrectForNonConstantOffDiagonals ) {
    solver.set( testmat1 ).set_range( -1.0, 3.0 ).solve();
    EXPECT_EQ( solver.count(), 10 );
    vector<double> eigs = solver.eigenvalues();
    for (size_t i = 0; i < eigs.size(); ++i) {
        EXPECT_NEAR( eigs[ i ], eigenvalues_expected1[ i ], 1e-8 );
    }
    std::vector<std::vector<std::complex<double>>> vecs
        = solver.eigenvectors();
    double factor = 1.0;
    for (size_t i = 0; i < eigs.size(); ++i) {
        if (vecs[ i ][ 0 ].real() * eigenvectors_expected1[ i ][ 0 ].real()
            < 0.0) {
            factor = -1.0;
        } else {
            factor = 1.0;
        }
        for (size_t j = 0; j < eigs.size(); ++j) {
            EXPECT_NEAR( vecs[ i ][ j ].real(),
                         factor * eigenvectors_expected1[ i ][ j ].real(),
                         1e-8 );
        }
    }
}