#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/math.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <gtest/gtest-spi.h>

#ifdef HAVE_GSL_INTERPOLATION_LIBRARY
#  include "gsl/gsl_interp.h"
#  include "gsl/gsl_spline.h"
#  include "gsl/gsl_version.h"
#endif

using namespace testing;
using namespace NCPA::interpolation;
using namespace std;

#define _TEST_TITLE_ NCPAInterpolationExtrapolationLibraryTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            interp = InterpolatorFactory<double, double>::build(
                interpolator_1d_type_t::LANL_LINEAR );
            interp.init( 5 );
            interp.fill( { 0.0, 1.0, 2.0, 3.0, 4.0 },
                         { 0.0, 2.0, 4.0, 5.0, 6.0 } );
            interp.ready();

        }  // void TearDown() override {}

        Interpolator1D<double, double> interp;
        _constant_extrapolator_1d<double, double> constant_ext;
        _zero_extrapolator_1d<double, double> zero_ext;
        _forbidden_extrapolator_1d<double, double> forbidden_ext;
        _linear_extrapolator_1d<double, double> linear_ext;
        _periodic_extrapolator_1d<double, double> periodic_ext;

        Extrapolator1D<double,double> wrapper;
};

TEST_F( _TEST_TITLE_, ConstantExtrapolatorIsCorrect ) {
    EXPECT_DOUBLE_EQ( constant_ext.extrapolate( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( constant_ext.extrapolate( 5.0, interp ), 6.0 );
    EXPECT_DOUBLE_EQ( constant_ext.extrapolate_df( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( constant_ext.extrapolate_df( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( constant_ext.extrapolate_ddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( constant_ext.extrapolate_ddf( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( constant_ext.extrapolate_dddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( constant_ext.extrapolate_dddf( 5.0, interp ), 0.0 );
    EXPECT_THROW( { constant_ext.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { constant_ext.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { constant_ext.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { constant_ext.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, ZeroExtrapolatorIsCorrect ) {
    EXPECT_DOUBLE_EQ( zero_ext.extrapolate( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( zero_ext.extrapolate( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( zero_ext.extrapolate_df( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( zero_ext.extrapolate_df( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( zero_ext.extrapolate_ddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( zero_ext.extrapolate_ddf( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( zero_ext.extrapolate_dddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( zero_ext.extrapolate_dddf( 5.0, interp ), 0.0 );
    EXPECT_THROW( { zero_ext.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { zero_ext.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { zero_ext.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { zero_ext.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, LinearExtrapolatorIsCorrect ) {
    EXPECT_DOUBLE_EQ( linear_ext.extrapolate( -1.0, interp ), -2.0 );
    EXPECT_DOUBLE_EQ( linear_ext.extrapolate( 5.0, interp ), 7.0 );
    EXPECT_DOUBLE_EQ( linear_ext.extrapolate_df( -1.0, interp ), 2.0 );
    EXPECT_DOUBLE_EQ( linear_ext.extrapolate_df( 5.0, interp ), 1.0 );
    EXPECT_DOUBLE_EQ( linear_ext.extrapolate_ddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( linear_ext.extrapolate_ddf( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( linear_ext.extrapolate_dddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( linear_ext.extrapolate_dddf( 5.0, interp ), 0.0 );
    EXPECT_THROW( { linear_ext.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { linear_ext.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { linear_ext.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { linear_ext.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, ForbiddenExtrapolatorIsCorrect ) {
    EXPECT_THROW( { forbidden_ext.extrapolate( -1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate_df( -1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate_ddf( -1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate_dddf( -1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate( 10.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate_df( 10.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate_ddf( 10.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate_dddf( 10.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { forbidden_ext.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, PeriodicExtrapolatorIsCorrect ) {
    EXPECT_DOUBLE_EQ( periodic_ext.extrapolate( -1.0, interp ), 5.0 );
    EXPECT_DOUBLE_EQ( periodic_ext.extrapolate( 5.0, interp ), 2.0 );
    EXPECT_DOUBLE_EQ( periodic_ext.extrapolate_df( -1.0, interp ), 1.0 );
    EXPECT_DOUBLE_EQ( periodic_ext.extrapolate_df( 5.0, interp ), 2.0 );
    EXPECT_DOUBLE_EQ( periodic_ext.extrapolate_ddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( periodic_ext.extrapolate_ddf( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( periodic_ext.extrapolate_dddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( periodic_ext.extrapolate_dddf( 5.0, interp ), 0.0 );
    EXPECT_THROW( { periodic_ext.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { periodic_ext.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { periodic_ext.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { periodic_ext.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, ConstantExtrapolatorWrapperIsCorrect ) {
    wrapper = InterpolatorFactory<double,double>::build( extrapolator_1d_type_t::CONSTANT );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate( 5.0, interp ), 6.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_df( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_df( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_ddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_ddf( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_dddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_dddf( 5.0, interp ), 0.0 );
    EXPECT_THROW( { wrapper.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, ZeroExtrapolatorWrapperIsCorrect ) {
    wrapper = InterpolatorFactory<double,double>::build( extrapolator_1d_type_t::ZERO );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_df( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_df( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_ddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_ddf( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_dddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_dddf( 5.0, interp ), 0.0 );
    EXPECT_THROW( { wrapper.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, LinearExtrapolatorWrapperIsCorrect ) {
    wrapper = InterpolatorFactory<double,double>::build( extrapolator_1d_type_t::LINEAR );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate( -1.0, interp ), -2.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate( 5.0, interp ), 7.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_df( -1.0, interp ), 2.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_df( 5.0, interp ), 1.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_ddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_ddf( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_dddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_dddf( 5.0, interp ), 0.0 );
    EXPECT_THROW( { wrapper.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, ForbiddenExtrapolatorWrapperIsCorrect ) {
    wrapper = InterpolatorFactory<double,double>::build( extrapolator_1d_type_t::FORBIDDEN );
    EXPECT_THROW( { wrapper.extrapolate( -1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_df( -1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_ddf( -1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_dddf( -1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate( 10.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_df( 10.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_ddf( 10.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_dddf( 10.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}

TEST_F( _TEST_TITLE_, PeriodicExtrapolatorWrapperIsCorrect ) {
    wrapper = InterpolatorFactory<double,double>::build( extrapolator_1d_type_t::PERIODIC );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate( -1.0, interp ), 5.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate( 5.0, interp ), 2.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_df( -1.0, interp ), 1.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_df( 5.0, interp ), 2.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_ddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_ddf( 5.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_dddf( -1.0, interp ), 0.0 );
    EXPECT_DOUBLE_EQ( wrapper.extrapolate_dddf( 5.0, interp ), 0.0 );
    EXPECT_THROW( { wrapper.extrapolate( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_df( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_ddf( 1.0, interp ); }, std::logic_error );
    EXPECT_THROW( { wrapper.extrapolate_dddf( 1.0, interp ); }, std::logic_error );
}