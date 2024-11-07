#include "NCPA/arrays.hpp"
#include "NCPA/math.hpp"
#include "NCPA/interpolation.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cmath>
#include <vector>
#include <complex>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace NCPA::interpolation;
using namespace std;

double eval_cubic( double x, vector<double> coeffs ) {
    double f;
    for ( auto e = 0; e < coeffs.size(); e++ ) {
        f += coeffs[ e ] * std::pow( x, (double)e );
    }
    return f;
}

vector<double> make_cubic( vector<double>& x, vector<double>& coeffs ) {
    vector<double> out( x.size() );
    for ( auto i = 0; i < x.size(); i++ ) {
        out[ i ] = eval_cubic( x[i], coeffs );
    }
    return out;
}

class NCPAInterpolationTest : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            x_lin     = { 0, 1, 2, 3, 4, 5 };
            f_lin = { -1, 1, -1, 1, -1, 1 };
            cf_lin = NCPA::math::real2complex( f_lin, f_lin );

            x_cub = NCPA::math::linspace(-5.0, 5.0, 101);
            cub_coeffs = { 2.5, -3.1, -0.45, 1.239 };
            f_cub      = make_cubic( x_cub, cub_coeffs );
            cub_tol = 0.01;

            lin_spline_1d.prep( 6 );
            lin_spline_1d.fill( x_lin, f_lin );
            lin_spline_1d.set();
            cub_spline_1d.prep( 6 );
            cub_spline_1d.fill( x_cub, f_cub );
            cub_spline_1d.set();
            clin_spline_1d.prep( 6 );
            clin_spline_1d.fill( x_lin, cf_lin );
            clin_spline_1d.set();
        }  // void TearDown() override {}

        // declare stuff here
        vector<double> x_lin, f_lin;
        vector<double> x_cub, f_cub, cub_coeffs;
        vector<complex<double>> cf_lin;
        LANL::linear_spline_1D<double,double> lin_spline_1d;
        LANL::linear_spline_1D<double,complex<double>> clin_spline_1d;
        LANL::natural_cubic_spline_1D<double,double> cub_spline_1d;
        double cub_tol = 0.01;
};

TEST_F( NCPAInterpolationTest, LinearInterpolationIsCorrect ) {
    for ( size_t i = 0; i < 6; i++ ) {
        // std::cout << "Checking for i = " << i << std::endl;
        EXPECT_DOUBLE_EQ( lin_spline_1d.eval_f( x_lin[ i ] ), f_lin[ i ] );
    }
    for ( double d = 0.5; d < 5.0; d += 1.0 ) {
        // std::cout << "Checking for d = " << d << std::endl;
        EXPECT_DOUBLE_EQ( lin_spline_1d.eval_f( d ), 0.0 );
    }
    // std::cout << "Wrapping up..." << std::endl;
}

TEST_F( NCPAInterpolationTest, LinearInterpolationOfComplexValuesIsCorrect ) {
    for ( size_t i = 0; i < 6; i++ ) {
        std::complex<double> cval = clin_spline_1d.eval_f( x_lin[ i ] );
        EXPECT_DOUBLE_EQ( cval.real(), cf_lin[ i ].real() );
        EXPECT_DOUBLE_EQ( cval.imag(), cf_lin[ i ].imag() );
    }
    for ( double d = 0.5; d < 5.0; d += 1.0 ) {
        std::complex<double> cval = clin_spline_1d.eval_f( d );
        EXPECT_DOUBLE_EQ( cval.real(), 0.0 );
        EXPECT_DOUBLE_EQ( cval.imag(), 0.0 );
    }
}

TEST_F( NCPAInterpolationTest, CubicSplineInterpolationAgreesWithinTolerance ) {
    for ( size_t i = 0; i < x_cub.size(); i++ ) {
        EXPECT_DOUBLE_EQ( cub_spline_1d.eval_f( x_cub[ i ] ), f_cub[ i ] );
    }
    for ( double d = NCPA::math::min(x_cub); d <= NCPA::math::max(x_cub); d += 0.1 ) {
        double expected = eval_cubic( d, cub_coeffs );
        EXPECT_NEAR( cub_spline_1d.eval_f( d ),
                          expected,
                          std::abs(expected * cub_tol)
                           );
    }
}

