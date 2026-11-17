#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <vector>

#include <gtest/gtest-spi.h>

#if NCPA_USE_GSL_INTERPOLATION
#  include "gsl/gsl_interp.h"
#  include "gsl/gsl_spline.h"
#  include "gsl/gsl_version.h"
#endif

using namespace testing;
using namespace NCPA::interpolation;
using namespace std;

#define _TEST_TITLE_ NCPAInterpolationLibraryTest

template<typename T>
T eval_cubic( double x, vector<T> coeffs ) {
    T f;
    for (auto e = 0; e < coeffs.size(); e++) {
        f += coeffs[ e ] * std::pow( x, (double)e );
    }
    return f;
}

template<typename T>
vector<T> make_cubic( vector<double>& x, vector<T>& coeffs ) {
    vector<T> out( x.size() );
    for (auto i = 0; i < x.size(); i++) {
        out[ i ] = eval_cubic( x[ i ], coeffs );
    }
    return out;
}

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            x_lin  = { 0, 1, 2, 3, 4, 5 };
            f_lin  = { -1, 1, -1, 1, -1, 1 };
            cf_lin = NCPA::math::real2complex( x_lin, -x_lin );

            x_cub      = NCPA::math::linspace( -5.0, 5.0, 101 );
            cub_coeffs = { 2.5, -3.1, -0.45, 1.239 };
            f_cub      = make_cubic( x_cub, cub_coeffs );
            // ccub_coeffs
            //     = { NCPA::math::random_number<complex<double>>( -2.0, 2.0 ),
            //         NCPA::math::random_number<complex<double>>( -2.0, 2.0 ),
            //         NCPA::math::random_number<complex<double>>( -2.0, 2.0 ),
            //         NCPA::math::random_number<complex<double>>( -2.0, 2.0 )
            //         };
            // cf_cub = make_cubic( x_cub, ccub_coeffs );

            cub_tol = 0.05;

            nearest_neighbor.fill( x_lin, f_lin );
            nearest_neighbor.ready();
            lanl_lin_spline_1d.fill( x_lin, f_lin );
            lanl_lin_spline_1d.ready();
            lanl_cub_spline_1d.fill( x_cub, f_cub );
            lanl_cub_spline_1d.ready();
            // lanl_clin_spline_1d.fill( x_lin, cf_lin );
            // lanl_clin_spline_1d.ready();
            // lanl_ccub_spline_1d.fill( x_cub, cf_cub );
            // lanl_ccub_spline_1d.ready();
#ifdef NCPA_USE_GSL_INTERPOLATION
            gsl_lin_spline_1d = GSL::GSLSpline<double>( gsl_interp_linear );
            gsl_lin_spline_1d.fill( x_lin, f_lin );
            gsl_lin_spline_1d.ready();

            // gsl_clin_spline_1d = GSL::gsl_spline_1d<double,
            // complex<double>>(
            //     gsl_interp_linear );
            // gsl_clin_spline_1d.fill( x_lin, cf_lin );
            // gsl_clin_spline_1d.ready();
#endif

            x2 = { 0, 1, 2, 3, 4 };
            y2 = { 0, 1, 2, 3, 4 };
            z2 = { 0, 1, 2, 3, 4 };
            f2d_y.resize( 5, 5 );
            f2d_x.resize( 5, 5 );
            // f2d_x.resize2d( 5, 5 );
            // f2d_y.resize2d( 5, 5 );
            // f3d.resize3d( 5, 5, 5 );

            for (size_t i = 0; i < 5; i++) {
                double di = (double)i;
                vector<double> xrow { di, di, di, di, di };
                f2d_x[ i ] = xrow;
                f2d_y[ i ] = x2;
                //     for (size_t j = 0; j < 5; j++) {
                //         f3d[ i ][ j ] = z2;
                //     }
            }

            lanl_natural_2d.init( 5, 5 );
            lanl_natural_2d.fill( x2, y2, f2d_x );
            lanl_natural_2d.ready();
            // lanl_natural_2d.fill( x2, y2, f2d_x );
            // lanl_natural_2d.ready();
            lanl_bicubic_2d.init( 5, 5 );
            lanl_bicubic_2d.fill( x2, y2, f2d_x );
            lanl_bicubic_2d.ready();
            // lanl_hybrid_3d.init( 5, 5, 5 );
            // lanl_hybrid_3d.fill( x2, y2, z2, f3d );
            // lanl_hybrid_3d.ready();


        }  // void TearDown() override {}

        // declare stuff here
        vector<double> x_lin, f_lin;
        vector<complex<double>> cf_lin;
        vector<double> x_cub, f_cub, cub_coeffs;
        NearestNeighborSpline<double> nearest_neighbor;
        LANL::LinearSpline1D<double> lanl_lin_spline_1d;
        LANL::NaturalCubicSpline1D<double> lanl_cub_spline_1d;
#ifdef NCPA_USE_GSL_INTERPOLATION
        GSL::GSLSpline<double> gsl_lin_spline_1d;
#endif
        double cub_tol = 0.01, nat_tol = 0.001, hyb_tol = 0.01;

        // NCPA::arrays::vector2d_t<double> f2d_x, f2d_y;
        // NCPA::arrays::vector3d_t<double> f3d;
        NCPA::arrays::SmartArray<2, double> f2d_x, f2d_y;
        vector<double> x2, y2, z2;
        LANL::NaturalSpline2D<double> lanl_natural_2d;
        // LANL::natural_spline_2d<double> lanl_natural_2d;
        LANL::BicubicSpline2D<double> lanl_bicubic_2d;
        // LANL::hybrid_spline_3d<double> lanl_hybrid_3d;

        // stratified_spline_2d<double, double> strat;

        // Interpolator<double> wrapper;
        // Interpolator1D<double, double> wrapper;
        // Interpolator1D<double, complex<double>> cwrapper;
        // Interpolator2D<double, double> wrapper2;
};

TEST_F( _TEST_TITLE_, NearestNeighborIsCorrect ) {
    for (size_t base = 0; base < 6; ++base) {
        for (double mod = -0.4; mod <= 0.4; mod += 0.1) {
            EXPECT_DOUBLE_EQ( nearest_neighbor.eval_f( x_lin[ base ] + mod ),
                              f_lin[ base ] );
        }
    }
}

TEST_F( _TEST_TITLE_, NearestNeighborDerivativesAreCorrect ) {
    EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ 0 ] - 0.5, 0 ), 0.0 );
    EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ 5 ] + 0.5, 0 ), 0.0 );

    for (size_t i = 0; i < 5; ++i) {
        EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ i ] + 0.5, 0 ),
                          ( f_lin[ i ] < 0.0 ? 2.0 : -2.0 ) );
        EXPECT_DOUBLE_EQ( nearest_neighbor.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ),
                          0.0 );
        EXPECT_DOUBLE_EQ(
            nearest_neighbor.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ), 0.0 );
    }
}

TEST_F( _TEST_TITLE_, NearestNeighborDerivativesAssumeWRTCorrectly ) {
    for (size_t base = 0; base < 6; ++base) {
        for (double mod = -0.4; mod <= 0.4; mod += 0.1) {
            EXPECT_DOUBLE_EQ(
                nearest_neighbor.eval_df( x_lin[ base ] + mod ),
                nearest_neighbor.eval_df( x_lin[ base ] + mod, 0 ) );
            EXPECT_DOUBLE_EQ(
                nearest_neighbor.eval_ddf( x_lin[ base ] + mod ),
                nearest_neighbor.eval_ddf( x_lin[ base ] + mod, 0, 0 ) );
            EXPECT_DOUBLE_EQ(
                nearest_neighbor.eval_dddf( x_lin[ base ] + mod ),
                nearest_neighbor.eval_dddf( x_lin[ base ] + mod, 0, 0, 0 ) );
        }
    }
}

TEST_F( _TEST_TITLE_, NearestNeighborThrowsOnBadlyFormedDerivatives ) {
    EXPECT_NO_THROW( { nearest_neighbor.eval_df( 1.5 ); } );
    EXPECT_NO_THROW( { nearest_neighbor.eval_df( 1.5, 0 ); } );
    EXPECT_THROW(
        { nearest_neighbor.eval_df( 1.5, 0, 0 ); }, AmbiguousCoordinates );
    EXPECT_THROW(
        { nearest_neighbor.eval_df( 1.5, 0, 0, 0 ); }, AmbiguousCoordinates );

    EXPECT_NO_THROW( { nearest_neighbor.eval_ddf( 1.5 ); } );
    EXPECT_THROW(
        { nearest_neighbor.eval_ddf( 1.5, 0 ); }, AmbiguousCoordinates );
    EXPECT_NO_THROW( { nearest_neighbor.eval_ddf( 1.5, 0, 0 ); } );
    EXPECT_THROW(
        { nearest_neighbor.eval_ddf( 1.5, 0, 0, 0 ); }, AmbiguousCoordinates );

    EXPECT_NO_THROW( { nearest_neighbor.eval_dddf( 1.5 ); } );
    EXPECT_THROW(
        { nearest_neighbor.eval_dddf( 1.5, 0 ); }, AmbiguousCoordinates );
    EXPECT_THROW(
        { nearest_neighbor.eval_dddf( 1.5, 0, 0 ); }, AmbiguousCoordinates );
    EXPECT_NO_THROW( { nearest_neighbor.eval_dddf( 1.5, 0, 0, 0 ); } );
}

TEST_F( _TEST_TITLE_, InterpolatorWrapperIsCorrectWithNearestNeighbor ) {
    Interpolator<double> wrapper( nearest_neighbor );
    for (size_t base = 0; base < 6; ++base) {
        for (double mod = -0.4; mod <= 0.4; mod += 0.1) {
            EXPECT_DOUBLE_EQ( nearest_neighbor.eval_f( x_lin[ base ] + mod ),
                              wrapper.eval_f( x_lin[ base ] + mod ) );
            EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ base ] + mod ),
                              wrapper.eval_df( x_lin[ base ] + mod ) );
            EXPECT_DOUBLE_EQ( nearest_neighbor.eval_ddf( x_lin[ base ] + mod ),
                              wrapper.eval_ddf( x_lin[ base ] + mod ) );
            EXPECT_DOUBLE_EQ(
                nearest_neighbor.eval_dddf( x_lin[ base ] + mod ),
                wrapper.eval_dddf( x_lin[ base ] + mod ) );
        }
    }
}

TEST_F( _TEST_TITLE_, NearestNeighborWrapperThrowsOnBadlyFormedDerivatives ) {
    Interpolator<double> wrapper( nearest_neighbor );
    EXPECT_NO_THROW( { wrapper.eval_df( 1.5 ); } );
    EXPECT_NO_THROW( { wrapper.eval_df( 1.5, 0 ); } );
    EXPECT_THROW( { wrapper.eval_df( 1.5, 0, 0 ); }, AmbiguousCoordinates );
    EXPECT_THROW( { wrapper.eval_df( 1.5, 0, 0, 0 ); }, AmbiguousCoordinates );

    EXPECT_NO_THROW( { wrapper.eval_ddf( 1.5 ); } );
    EXPECT_THROW( { wrapper.eval_ddf( 1.5, 0 ); }, AmbiguousCoordinates );
    EXPECT_NO_THROW( { wrapper.eval_ddf( 1.5, 0, 0 ); } );
    EXPECT_THROW(
        { wrapper.eval_ddf( 1.5, 0, 0, 0 ); }, AmbiguousCoordinates );

    EXPECT_NO_THROW( { wrapper.eval_dddf( 1.5 ); } );
    EXPECT_THROW( { wrapper.eval_dddf( 1.5, 0 ); }, AmbiguousCoordinates );
    EXPECT_THROW( { wrapper.eval_dddf( 1.5, 0, 0 ); }, AmbiguousCoordinates );
    EXPECT_NO_THROW( { wrapper.eval_dddf( 1.5, 0, 0, 0 ); } );
}

TEST_F( _TEST_TITLE_, ComplexNearestNeighborWrapperAgreesWithComponents ) {
    Interpolator<complex<double>> wrapper( nearest_neighbor,
                                           nearest_neighbor );
    for (size_t base = 0; base < 6; ++base) {
        for (double mod = -0.4; mod <= 0.4; mod += 0.1) {
            EXPECT_DOUBLE_EQ( wrapper.eval_f( x_lin[ base ] + mod ).real(),
                              nearest_neighbor.eval_f( x_lin[ base ] + mod ) );
            EXPECT_DOUBLE_EQ( wrapper.eval_f( x_lin[ base ] + mod ).imag(),
                              nearest_neighbor.eval_f( x_lin[ base ] + mod ) );
        }
    }
}

TEST_F( _TEST_TITLE_,
        ComplexNearestNeighborWrapperAgreesWithComponentDerivatives ) {
    Interpolator<complex<double>> wrapper( nearest_neighbor,
                                           nearest_neighbor );
    EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ 0 ] - 0.5, 0 ),
                      wrapper.eval_df( x_lin[ 0 ] - 0.5, 0 ).real() );
    EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ 0 ] - 0.5, 0 ),
                      wrapper.eval_df( x_lin[ 0 ] - 0.5, 0 ).imag() );
    EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ 5 ] + 0.5, 0 ),
                      wrapper.eval_df( x_lin[ 5 ] + 0.5, 0 ).real() );
    EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ 5 ] + 0.5, 0 ),
                      wrapper.eval_df( x_lin[ 5 ] + 0.5, 0 ).imag() );

    for (size_t i = 0; i < 5; ++i) {
        EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ i ] + 0.5, 0 ),
                          wrapper.eval_df( x_lin[ i ] + 0.5, 0 ).real() );
        EXPECT_DOUBLE_EQ( nearest_neighbor.eval_df( x_lin[ i ] + 0.5, 0 ),
                          wrapper.eval_df( x_lin[ i ] + 0.5, 0 ).imag() );
        EXPECT_DOUBLE_EQ( nearest_neighbor.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ),
                          wrapper.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ).real() );
        EXPECT_DOUBLE_EQ( nearest_neighbor.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ),
                          wrapper.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ).imag() );
        EXPECT_DOUBLE_EQ(
            nearest_neighbor.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ),
            wrapper.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ).real() );
        EXPECT_DOUBLE_EQ(
            nearest_neighbor.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ),
            wrapper.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ).imag() );
    }
}

TEST_F( _TEST_TITLE_,
        ComplexRealOnlyNearestNeighborWrapperAgreesWithComponents ) {
    Interpolator<complex<double>> wrapper( nearest_neighbor );
    for (size_t base = 0; base < 6; ++base) {
        for (double mod = -0.4; mod <= 0.4; mod += 0.1) {
            EXPECT_DOUBLE_EQ( wrapper.eval_f( x_lin[ base ] + mod ).real(),
                              nearest_neighbor.eval_f( x_lin[ base ] + mod ) );
            EXPECT_DOUBLE_EQ( wrapper.eval_f( x_lin[ base ] + mod ).imag(),
                              0.0 );
        }
    }
}

TEST_F( _TEST_TITLE_, LANLLinearLimitsAreCorrect ) {
    auto limits = lanl_lin_spline_1d.limits().at( 0 );
    EXPECT_DOUBLE_EQ( limits.first, x_lin.front() );
    EXPECT_DOUBLE_EQ( limits.second, x_lin.back() );
}

TEST_F( _TEST_TITLE_, LANLLinearInterpolationIsCorrect ) {
    for (size_t i = 0; i < 6; i++) {
        // double d1 = lanl_lin_spline_1d.eval_f( x_lin[ i ] );
        EXPECT_DOUBLE_EQ( lanl_lin_spline_1d.eval_f( x_lin[ i ] ),
                          f_lin[ i ] );
    }
    for (double d = 0.5; d < 5.0; d += 1.0) {
        EXPECT_DOUBLE_EQ( lanl_lin_spline_1d.eval_f( d ), 0.0 );
    }
}

TEST_F( _TEST_TITLE_, LANLLinearDerivativesAreCorrect ) {
    EXPECT_THROW(
        { lanl_lin_spline_1d.eval_df( x_lin[ 0 ] - 0.5, 0 ); },
        RequestedPointOutOfBounds );
    EXPECT_THROW(
        { lanl_lin_spline_1d.eval_df( x_lin[ 5 ] + 0.5, 0 ); },
        RequestedPointOutOfBounds );

    auto limits = lanl_lin_spline_1d.limits().at( 0 );
    for (size_t i = 0; i < 5; ++i) {
        double val = x_lin[ i ] + 0.5;
        if (val >= limits.first && val <= limits.second) {
            EXPECT_DOUBLE_EQ(
                lanl_lin_spline_1d.eval_df( x_lin[ i ] + 0.5, 0 ),
                ( f_lin[ i ] < 0.0 ? 2.0 : -2.0 ) );
            EXPECT_DOUBLE_EQ(
                lanl_lin_spline_1d.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ), 0.0 );
            EXPECT_DOUBLE_EQ(
                lanl_lin_spline_1d.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ),
                0.0 );
        } else {
            EXPECT_THROW(
                { lanl_lin_spline_1d.eval_df( x_lin[ i ] + 0.5, 0 ); },
                RequestedPointOutOfBounds );
            EXPECT_THROW(
                { lanl_lin_spline_1d.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ); },
                RequestedPointOutOfBounds );
            EXPECT_THROW(
                { lanl_lin_spline_1d.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ); },
                RequestedPointOutOfBounds );
        }
    }
}

TEST_F( _TEST_TITLE_, LANLLinearDerivativesAssumeWRTCorrectly ) {
    auto limits = lanl_lin_spline_1d.limits().at( 0 );
    for (size_t base = 0; base < 6; ++base) {
        for (double mod = -0.4; mod <= 0.4; mod += 0.1) {
            double val = x_lin[ base ] + mod;
            if (val >= limits.first && val <= limits.second) {
                EXPECT_DOUBLE_EQ(
                    lanl_lin_spline_1d.eval_df( x_lin[ base ] + mod ),
                    lanl_lin_spline_1d.eval_df( x_lin[ base ] + mod, 0 ) );
                EXPECT_DOUBLE_EQ(
                    lanl_lin_spline_1d.eval_ddf( x_lin[ base ] + mod ),
                    lanl_lin_spline_1d.eval_ddf( x_lin[ base ] + mod, 0, 0 ) );
                EXPECT_DOUBLE_EQ(
                    lanl_lin_spline_1d.eval_dddf( x_lin[ base ] + mod ),
                    lanl_lin_spline_1d.eval_dddf( x_lin[ base ] + mod, 0, 0,
                                                  0 ) );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, LANLLinearThrowsOnBadlyFormedDerivatives ) {
    EXPECT_NO_THROW( { lanl_lin_spline_1d.eval_df( 1.5 ); } );
    EXPECT_NO_THROW( { lanl_lin_spline_1d.eval_df( 1.5, 0 ); } );
    EXPECT_THROW(
        { lanl_lin_spline_1d.eval_df( 1.5, 0, 0 ); }, AmbiguousCoordinates );
    EXPECT_THROW(
        { lanl_lin_spline_1d.eval_df( 1.5, 0, 0, 0 ); },
        AmbiguousCoordinates );

    EXPECT_NO_THROW( { lanl_lin_spline_1d.eval_ddf( 1.5 ); } );
    EXPECT_THROW(
        { lanl_lin_spline_1d.eval_ddf( 1.5, 0 ); }, AmbiguousCoordinates );
    EXPECT_NO_THROW( { lanl_lin_spline_1d.eval_ddf( 1.5, 0, 0 ); } );
    EXPECT_THROW(
        { lanl_lin_spline_1d.eval_ddf( 1.5, 0, 0, 0 ); },
        AmbiguousCoordinates );

    EXPECT_NO_THROW( { lanl_lin_spline_1d.eval_dddf( 1.5 ); } );
    EXPECT_THROW(
        { lanl_lin_spline_1d.eval_dddf( 1.5, 0 ); }, AmbiguousCoordinates );
    EXPECT_THROW(
        { lanl_lin_spline_1d.eval_dddf( 1.5, 0, 0 ); }, AmbiguousCoordinates );
    EXPECT_NO_THROW( { lanl_lin_spline_1d.eval_dddf( 1.5, 0, 0, 0 ); } );
}

TEST_F( _TEST_TITLE_, InterpolatorWrapperIsCorrectWithLANLLinear ) {
    Interpolator<double> wrapper( lanl_lin_spline_1d );
    auto limits = wrapper.limits().at( 0 );
    for (size_t base = 0; base < 6; ++base) {
        for (double mod = -0.4; mod <= 0.4; mod += 0.1) {
            double val = x_lin[ base ] + mod;
            if (val >= limits.first && val <= limits.second) {
                EXPECT_DOUBLE_EQ(
                    lanl_lin_spline_1d.eval_f( x_lin[ base ] + mod ),
                    wrapper.eval_f( x_lin[ base ] + mod ) );
                EXPECT_DOUBLE_EQ(
                    lanl_lin_spline_1d.eval_df( x_lin[ base ] + mod ),
                    wrapper.eval_df( x_lin[ base ] + mod ) );
                EXPECT_DOUBLE_EQ(
                    lanl_lin_spline_1d.eval_ddf( x_lin[ base ] + mod ),
                    wrapper.eval_ddf( x_lin[ base ] + mod ) );
                EXPECT_DOUBLE_EQ(
                    lanl_lin_spline_1d.eval_dddf( x_lin[ base ] + mod ),
                    wrapper.eval_dddf( x_lin[ base ] + mod ) );
            } else {
                EXPECT_THROW(
                    { lanl_lin_spline_1d.eval_f( val ); },
                    RequestedPointOutOfBounds );
                EXPECT_THROW(
                    { lanl_lin_spline_1d.eval_df( val ); },
                    RequestedPointOutOfBounds );
                EXPECT_NO_THROW( { lanl_lin_spline_1d.eval_ddf( val ); } );
                EXPECT_NO_THROW( { lanl_lin_spline_1d.eval_dddf( val ); } );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, LANLLinearWrapperThrowsOnBadlyFormedDerivatives ) {
    Interpolator<double> wrapper( lanl_lin_spline_1d );
    EXPECT_NO_THROW( { wrapper.eval_df( 1.5 ); } );
    EXPECT_NO_THROW( { wrapper.eval_df( 1.5, 0 ); } );
    EXPECT_THROW( { wrapper.eval_df( 1.5, 0, 0 ); }, AmbiguousCoordinates );
    EXPECT_THROW( { wrapper.eval_df( 1.5, 0, 0, 0 ); }, AmbiguousCoordinates );

    EXPECT_NO_THROW( { wrapper.eval_ddf( 1.5 ); } );
    EXPECT_THROW( { wrapper.eval_ddf( 1.5, 0 ); }, AmbiguousCoordinates );
    EXPECT_NO_THROW( { wrapper.eval_ddf( 1.5, 0, 0 ); } );
    EXPECT_THROW(
        { wrapper.eval_ddf( 1.5, 0, 0, 0 ); }, AmbiguousCoordinates );

    EXPECT_NO_THROW( { wrapper.eval_dddf( 1.5 ); } );
    EXPECT_THROW( { wrapper.eval_dddf( 1.5, 0 ); }, AmbiguousCoordinates );
    EXPECT_THROW( { wrapper.eval_dddf( 1.5, 0, 0 ); }, AmbiguousCoordinates );
    EXPECT_NO_THROW( { wrapper.eval_dddf( 1.5, 0, 0, 0 ); } );
}

TEST_F( _TEST_TITLE_, ComplexLANLLinearWrapperAgreesWithComponents ) {
    Interpolator<complex<double>> wrapper( lanl_lin_spline_1d,
                                           lanl_lin_spline_1d );
    auto limits = wrapper.limits().at( 0 );
    for (size_t base = 0; base < 6; ++base) {
        for (double mod = -0.4; mod <= 0.4; mod += 0.1) {
            double val = x_lin[ base ] + mod;
            if (val >= limits.first && val <= limits.second) {
                EXPECT_DOUBLE_EQ(
                    wrapper.eval_f( x_lin[ base ] + mod ).real(),
                    lanl_lin_spline_1d.eval_f( x_lin[ base ] + mod ) );
                EXPECT_DOUBLE_EQ(
                    wrapper.eval_f( x_lin[ base ] + mod ).imag(),
                    lanl_lin_spline_1d.eval_f( x_lin[ base ] + mod ) );
            }
        }
    }
}

TEST_F( _TEST_TITLE_,
        ComplexLANLLinearWrapperAgreesWithComponentDerivatives ) {
    Interpolator<complex<double>> wrapper( lanl_lin_spline_1d,
                                           lanl_lin_spline_1d );
    auto limits = wrapper.limits().at( 0 );
    for (size_t i = 0; i < 5; ++i) {
        double val = x_lin[ i ] + 0.5;
        if (val >= limits.first && val <= limits.second) {
            EXPECT_DOUBLE_EQ(
                lanl_lin_spline_1d.eval_df( x_lin[ i ] + 0.5, 0 ),
                wrapper.eval_df( x_lin[ i ] + 0.5, 0 ).real() );
            EXPECT_DOUBLE_EQ(
                lanl_lin_spline_1d.eval_df( x_lin[ i ] + 0.5, 0 ),
                wrapper.eval_df( x_lin[ i ] + 0.5, 0 ).imag() );
            EXPECT_DOUBLE_EQ(
                lanl_lin_spline_1d.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ),
                wrapper.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ).real() );
            EXPECT_DOUBLE_EQ(
                lanl_lin_spline_1d.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ),
                wrapper.eval_ddf( x_lin[ i ] + 0.5, 0, 0 ).imag() );
            EXPECT_DOUBLE_EQ(
                lanl_lin_spline_1d.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ),
                wrapper.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ).real() );
            EXPECT_DOUBLE_EQ(
                lanl_lin_spline_1d.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ),
                wrapper.eval_dddf( x_lin[ i ] + 0.5, 0, 0, 0 ).imag() );
        }
    }
}

TEST_F( _TEST_TITLE_, ComplexRealOnlyLANLLinearWrapperAgreesWithComponents ) {
    Interpolator<complex<double>> wrapper( lanl_lin_spline_1d );
    auto limits = wrapper.limits().at( 0 );
    for (size_t base = 0; base < 6; ++base) {
        for (double mod = -0.4; mod <= 0.4; mod += 0.1) {
            double val = x_lin[ base ] + mod;
            if (val >= limits.first && val <= limits.second) {
                EXPECT_DOUBLE_EQ(
                    wrapper.eval_f( x_lin[ base ] + mod ).real(),
                    lanl_lin_spline_1d.eval_f( x_lin[ base ] + mod ) );
                EXPECT_DOUBLE_EQ( wrapper.eval_f( x_lin[ base ] + mod ).imag(),
                                  0.0 );
            }
        }
    }
}

TEST_F( _TEST_TITLE_, LANLCubicSplineInterpolationAgreesWithinTolerance ) {
    for (size_t i = 0; i < x_cub.size(); i++) {
        EXPECT_DOUBLE_EQ( lanl_cub_spline_1d.eval_f( x_cub[ i ] ),
                          f_cub[ i ] );
    }
    for (double d  = NCPA::math::min( x_cub ); d <= NCPA::math::max( x_cub );
         d        += 0.1) {
        double expected = eval_cubic( d, cub_coeffs ),
               actual   = lanl_cub_spline_1d.eval_f( d );
        EXPECT_NEAR( actual, expected, std::abs( expected * cub_tol ) );
    }
}

// TEST_F( _TEST_TITLE_,
//         LANLComplexCubicSplineInterpolationAgreesWithinTolerance ) {
//     for (size_t i = 0; i < x_cub.size(); i++) {
//         EXPECT_NEAR( lanl_ccub_spline_1d.eval_f( x_cub[ i ] ).real(),
//                      cf_cub[ i ].real(), 1e-12 );
//         EXPECT_NEAR( lanl_ccub_spline_1d.eval_f( x_cub[ i ] ).imag(),
//                      cf_cub[ i ].imag(), 1e-12 );
//     }
//     // std::ofstream testout( "lanl_cubic_spline_output.dat", ios_base::out
//     ); for (double d  = NCPA::math::min( x_cub ); d <= NCPA::math::max(
//     x_cub );
//          d        += 0.1) {
//         complex<double> expected = eval_cubic( d, ccub_coeffs ),
//                         actual   = lanl_ccub_spline_1d.eval_f( d );
//         // testout << d << " " << expected << " " << actual << endl;
//         EXPECT_NEAR( actual.real(), expected.real(),
//                      std::abs( expected.real() * cub_tol ) );
//         EXPECT_NEAR( actual.imag(), expected.imag(),
//                      std::abs( expected.imag() * cub_tol ) );
//     }
// }

TEST_F( _TEST_TITLE_, LANLNatural2DIsCorrect ) {
    EXPECT_NEAR( lanl_natural_2d.eval_f( 2.5, 2.5 ), 2.5, nat_tol );
    EXPECT_NEAR( lanl_natural_2d.eval_f( 1.5, 3.5 ), 1.5, nat_tol );
    // EXPECT_NEAR( lanl_bicubic_2d.eval_f( 2.5, 2.5 ), 2.5, nat_tol );
    // EXPECT_NEAR( lanl_bicubic_2d.eval_f( 1.5, 3.5 ), 1.5, nat_tol );
    lanl_natural_2d.clear();
    lanl_natural_2d.init( 5, 5 );
    lanl_natural_2d.fill( x2, y2, f2d_y );
    lanl_natural_2d.ready();
    // lanl_bicubic_2d.clear();
    // lanl_bicubic_2d.init( 5, 5 );
    // lanl_bicubic_2d.fill( x2, y2, f2d_y );
    // lanl_bicubic_2d.ready();
    EXPECT_NEAR( lanl_natural_2d.eval_f( 2.5, 2.5 ), 2.5, cub_tol );
    EXPECT_NEAR( lanl_natural_2d.eval_f( 1.5, 3.3 ), 3.3, cub_tol );
    // EXPECT_NEAR( lanl_bicubic_2d.eval_f( 2.5, 2.5 ), 2.5, cub_tol );
    // EXPECT_NEAR( lanl_bicubic_2d.eval_f( 1.5, 3.9 ), 3.9, cub_tol );
}

TEST_F( _TEST_TITLE_, LANLBicubic2DIsCorrect ) {
    EXPECT_NEAR( lanl_bicubic_2d.eval_f( 2.5, 2.5 ), 2.5, nat_tol );
    EXPECT_NEAR( lanl_bicubic_2d.eval_f( 1.5, 3.5 ), 1.5, nat_tol );
    lanl_bicubic_2d.clear();
    lanl_bicubic_2d.init( 5, 5 );
    lanl_bicubic_2d.fill( x2, y2, f2d_y );
    lanl_bicubic_2d.ready();
    EXPECT_NEAR( lanl_bicubic_2d.eval_f( 2.5, 2.5 ), 2.5, cub_tol );
    EXPECT_NEAR( lanl_bicubic_2d.eval_f( 1.5, 3.9 ), 3.9, cub_tol );
}

TEST_F( _TEST_TITLE_, LANLNaturalWorksForParaboloid ) {
    // lanl_bicubic_2d.clear();
    // lanl_bicubic_2d.init( 6, 6 );
    f2d_y.clear().resize( 6, 6 );
    std::vector<double> xvals { 0, 1, 2, 3, 4, 5 }, yvals { 0, 1, 2, 3, 4, 5 };
    for (size_t i = 0; i <= 5; ++i) {
        for (size_t j = 0; j <= 5; ++j) {
            double x = (double)i;
            double y = (double)j;
            double z
                = -0.4
                    * ( ( x - 2.5 ) * ( x - 2.5 ) + ( y - 2.5 ) * ( y - 2.5 ) )
                + 5.0;
            f2d_y.at( i, j ) = z;
        }
    }
    // lanl_bicubic_2d.fill( xvals, yvals, f2d_y );
    // lanl_bicubic_2d.ready();

    std::vector<double> xtest { 0.5, 1.5, 2.5, 3.5, 4.5 },
        ytest { 0, 1, 2, 3, 4, 5 };
    NCPA::arrays::SmartVector<double> ztest( xtest.size(), ytest.size() );
    ztest = std::vector<double> { 0.9000, 2.5000, 3.3000, 3.3000, 2.5000,
                                  0.9000, 2.1000, 3.7000, 4.5000, 4.5000,
                                  3.7000, 2.1000, 2.5000, 4.1000, 4.9000,
                                  4.9000, 4.1000, 2.5000, 2.1000, 3.7000,
                                  4.5000, 4.5000, 3.7000, 2.1000, 0.9000,
                                  2.5000, 3.3000, 3.3000, 2.5000, 0.9000 };

    lanl_natural_2d.clear();
    lanl_natural_2d.init( 6, 6 );
    lanl_natural_2d.fill( xvals, yvals, f2d_y );
    lanl_natural_2d.ready();

    

    for (size_t i = 0; i < xtest.size(); ++i) {
        for (size_t j = 0; j < ytest.size(); ++j) {
            // double d = lanl_natural_2d.eval_f( xtest[ i ], ytest[ j ] );
            // cout << "z(" << xtest[ i ] << "," << ytest[ j ] << ") = " << d
            //      << endl;
            EXPECT_NEAR( lanl_natural_2d.eval_f( xtest[ i ], ytest[ j ] ),
                         ztest.at( i, j ), 0.05 );
        }
    }
}

// TEST_F( _TEST_TITLE_, Stratified2DIsCorrect ) {
//     // strat = InterpolatorFactory<double,double>::build(
//     // interpolator_2d_type_t::LANL_LINEAR_Y );
//     strat = stratified_spline_2d<double, double>(
//         stratified_axis_type_t { 0, interpolator_1d_type_t::LANL_LINEAR } );
//     strat.init( 1, 5 );
//     strat.fill( x2, y2, f2d_y );
//     strat.ready();
//     EXPECT_NEAR( strat.eval_f( 2.5, 2.5 ), 2.5, nat_tol );
//     EXPECT_NEAR( strat.eval_f( 1.5, 3.5 ), 3.5, nat_tol );
//     EXPECT_NEAR( strat.eval_f( -3.0, 3.5 ), 3.5, nat_tol );
//     EXPECT_THROW( { strat.eval_f( 1.5, 6.5 ); }, std::logic_error );
// }

// TEST_F( _TEST_TITLE_, LANL3DIsCorrect ) {
//     EXPECT_NEAR( lanl_hybrid_3d.eval_f( 2.5, 2.5, 2.5 ), 2.5, hyb_tol );
//     EXPECT_NEAR( lanl_hybrid_3d.eval_f( 1.5, 3.5, 1.5 ), 1.5, hyb_tol );
//     EXPECT_NEAR( lanl_hybrid_3d.eval_f( 1.2, 2.3, 3.4 ), 3.4, hyb_tol );
//     EXPECT_NEAR( lanl_hybrid_3d.eval_f( 3.9, 3.9, 3.9 ), 3.9, hyb_tol );
// }

#if NCPA_USE_GSL_INTERPOLATION
TEST_F( _TEST_TITLE_, GSLLinearInterpolationIsCorrect ) {
    for (size_t i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ( gsl_lin_spline_1d.eval_f( x_lin[ i ] ), f_lin[ i ] );
    }
    for (double d = 0.5; d < 5.0; d += 1.0) {
        EXPECT_DOUBLE_EQ( gsl_lin_spline_1d.eval_f( d ), 0.0 );
    }
}

TEST_F( _TEST_TITLE_, GSLLinearInterpolationAgreesWithLANL ) {
    for (size_t i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ( gsl_lin_spline_1d.eval_f( x_lin[ i ] ),
                          lanl_lin_spline_1d.eval_f( x_lin[ i ] ) );
    }
    for (double d = 0.5; d < 5.0; d += 1.0) {
        EXPECT_DOUBLE_EQ( gsl_lin_spline_1d.eval_f( d ),
                          lanl_lin_spline_1d.eval_f( d ) );
    }
}

// TEST_F( _TEST_TITLE_, GSLLinearInterpolationOfComplexValuesIsCorrect ) {
//     for (size_t i = 0; i < 6; i++) {
//         std::complex<double> cval = gsl_clin_spline_1d.eval_f( x_lin[ i ] );
//         EXPECT_DOUBLE_EQ( cval.real(), cf_lin[ i ].real() );
//         EXPECT_DOUBLE_EQ( cval.imag(), cf_lin[ i ].imag() );
//     }
//     for (double d = 0.5; d < 5.0; d += 1.0) {
//         std::complex<double> cval = gsl_clin_spline_1d.eval_f( d );
//         EXPECT_DOUBLE_EQ( cval.real(), 0.0 );
//         EXPECT_DOUBLE_EQ( cval.imag(), 0.0 );
//     }
// }
#endif
