#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/abstract_spline_1d.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/Extrapolator1D.hpp"
#include "NCPA/interpolation/gsl.hpp"
#include "NCPA/interpolation/Interpolator1D.hpp"
#include "NCPA/interpolation/lanl.hpp"
#include "NCPA/interpolation/nearest_neighbor_spline_1d.hpp"
#include "NCPA/interpolation/stratified_spline_2d.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace interpolation {
        template<typename INDEPTYPE, typename DEPTYPE>
        class InterpolatorFactory {
            public:
                static bool can_build( interpolator_1d_type_t interptype ) {
                    switch ( interptype ) {
                        case interpolator_1d_type_t::NEAREST_NEIGHBOR:
                        case interpolator_1d_type_t::LANL_LINEAR:
                        case interpolator_1d_type_t::LANL_CUBIC:
#ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE
                        case interpolator_1d_type_t::GSL_LINEAR:
                        case interpolator_1d_type_t::GSL_POLYNOMIAL:
                        case interpolator_1d_type_t::GSL_CUBIC:
                        case interpolator_1d_type_t::GSL_CUBIC_PERIODIC:
                        case interpolator_1d_type_t::GSL_AKIMA:
                        case interpolator_1d_type_t::GSL_AKIMA_PERIODIC:
#  if GSL_MAJOR_VERSION >= 2
                        case interpolator_1d_type_t::GSL_STEFFEN:
#  endif
#endif
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static bool can_build( interpolator_2d_type_t interptype ) {
                    switch ( interptype ) {
                        case interpolator_2d_type_t::NEAREST_NEIGHBOR:
                        case interpolator_2d_type_t::LANL_NATURAL:
                        case interpolator_2d_type_t::LANL_BICUBIC:
                        case interpolator_2d_type_t::LANL_LINEAR_X:
                        case interpolator_2d_type_t::LANL_LINEAR_Y:
                        case interpolator_2d_type_t::LANL_CUBIC_X:
                        case interpolator_2d_type_t::LANL_CUBIC_Y:
#ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE
                        case interpolator_2d_type_t::GSL_LINEAR_X:
                        case interpolator_2d_type_t::GSL_LINEAR_Y:
#  if GSL_MAJOR_VERSION >= 2
                        case interpolator_2d_type_t::GSL_STEFFEN_X:
                        case interpolator_2d_type_t::GSL_STEFFEN_Y:
#  endif
#endif
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static bool can_build( interpolator_3d_type_t interptype ) {
                    switch ( interptype ) {
                        case interpolator_3d_type_t::LANL_HYBRID:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                // template<typename INDEPTYPE, typename DEPTYPE>
                static Interpolator1D<INDEPTYPE, DEPTYPE> build(
                    interpolator_1d_type_t interptype ) {
                    Interpolator1D<INDEPTYPE, DEPTYPE> interp;

                    switch ( interptype ) {
                        case interpolator_1d_type_t::NEAREST_NEIGHBOR:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new nearest_neighbor_spline_1d<
                                        INDEPTYPE, DEPTYPE>() ) );
                            break;
                        case interpolator_1d_type_t::LANL_LINEAR:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new LANL::linear_spline_1d<INDEPTYPE,
                                                               DEPTYPE>() ) );
                            break;
                        case interpolator_1d_type_t::LANL_CUBIC:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new LANL::natural_cubic_spline_1d<
                                        INDEPTYPE, DEPTYPE>() ) );
                            break;
#ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE
                        case interpolator_1d_type_t::GSL_LINEAR:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                        gsl_interp_linear ) ) );
                            break;
                        case interpolator_1d_type_t::GSL_POLYNOMIAL:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                        gsl_interp_polynomial ) ) );
                            break;
                        case interpolator_1d_type_t::GSL_CUBIC:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                        gsl_interp_cspline ) ) );
                            break;
                        case interpolator_1d_type_t::GSL_CUBIC_PERIODIC:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                        gsl_interp_cspline_periodic ) ) );
                            break;
                        case interpolator_1d_type_t::GSL_AKIMA:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                        gsl_interp_akima ) ) );
                            break;
                        case interpolator_1d_type_t::GSL_AKIMA_PERIODIC:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                        gsl_interp_akima_periodic ) ) );
                            break;
#  if GSL_MAJOR_VERSION >= 2
                        case interpolator_1d_type_t::GSL_STEFFEN:
                            interp.set_engine(
                                spline_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                        gsl_interp_steffen ) ) );
                            break;
#  endif
#endif
                        default:
                            throw std::range_error(
                                "Requested 1-D interpolator type "
                                "unrecognized; either it is undefined, not "
                                "applicable, or you don't have the libraries "
                                "available." );
                    }
                    return interp;
                }

                // template<typename INDEPTYPE, typename DEPTYPE>
                static Interpolator2D<INDEPTYPE, DEPTYPE> build(
                    interpolator_2d_type_t interptype ) {
                    Interpolator2D<INDEPTYPE, DEPTYPE> interp;

                    switch ( interptype ) {
                        case interpolator_2d_type_t::NEAREST_NEIGHBOR:
                            interp.set_engine(
                                spline_engine_2d_t<INDEPTYPE, DEPTYPE>(
                                    new nearest_neighbor_spline_2d<
                                        INDEPTYPE, DEPTYPE>() ) );
                            break;
                        case interpolator_2d_type_t::LANL_NATURAL:
                            interp.set_engine(
                                spline_engine_2d_t<INDEPTYPE, DEPTYPE>(
                                    new LANL::natural_spline_2d<INDEPTYPE,
                                                                DEPTYPE>() ) );
                            break;
                        case interpolator_2d_type_t::LANL_BICUBIC:
                            interp.set_engine(
                                spline_engine_2d_t<INDEPTYPE, DEPTYPE>(
                                    new LANL::bicubic_spline_2d<INDEPTYPE,
                                                                DEPTYPE>() ) );
                            break;
                        case interpolator_2d_type_t::LANL_LINEAR_X:
                            interp.set_engine( spline_engine_2d_t<INDEPTYPE,
                                                                  DEPTYPE>(
                                new stratified_spline_2d<INDEPTYPE, DEPTYPE>(
                                    std::pair<size_t, interpolator_1d_type_t> {
                                        1, interpolator_1d_type_t::
                                               LANL_LINEAR } ) ) );
                            break;
                        case interpolator_2d_type_t::LANL_LINEAR_Y:
                            interp.set_engine( spline_engine_2d_t<INDEPTYPE,
                                                                  DEPTYPE>(
                                new stratified_spline_2d<INDEPTYPE, DEPTYPE>(
                                    std::pair<size_t, interpolator_1d_type_t> {
                                        0, interpolator_1d_type_t::
                                               LANL_LINEAR } ) ) );
                            break;
                        case interpolator_2d_type_t::LANL_CUBIC_X:
                            interp.set_engine( spline_engine_2d_t<INDEPTYPE,
                                                                  DEPTYPE>(
                                new stratified_spline_2d<INDEPTYPE, DEPTYPE>(
                                    std::pair<size_t, interpolator_1d_type_t> {
                                        1, interpolator_1d_type_t::
                                               LANL_CUBIC } ) ) );
                            break;
                        case interpolator_2d_type_t::LANL_CUBIC_Y:
                            interp.set_engine( spline_engine_2d_t<INDEPTYPE,
                                                                  DEPTYPE>(
                                new stratified_spline_2d<INDEPTYPE, DEPTYPE>(
                                    std::pair<size_t, interpolator_1d_type_t> {
                                        0, interpolator_1d_type_t::
                                               LANL_CUBIC } ) ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 2-D interpolator type "
                                "unrecognized; either it is undefined, not "
                                "applicable, or you don't have the libraries "
                                "available." );
                    }
                    return interp;
                }

                // template<typename INDEPTYPE, typename DEPTYPE>
                static Interpolator3D<INDEPTYPE, DEPTYPE> build(
                    interpolator_3d_type_t interptype ) {
                    Interpolator3D<INDEPTYPE, DEPTYPE> interp;

                    switch ( interptype ) {
                        case interpolator_3d_type_t::LANL_HYBRID:
                            interp.set_engine(
                                spline_engine_3d_t<INDEPTYPE, DEPTYPE>(
                                    new LANL::hybrid_spline_3d<INDEPTYPE,
                                                               DEPTYPE>() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 3-D interpolator type "
                                "unrecognized; either it is undefined, not "
                                "applicable, or you don't have the libraries "
                                "available." );
                    }
                    return interp;
                }

                static Extrapolator1D<INDEPTYPE, DEPTYPE> build(
                    extrapolator_1d_type_t interptype ) {
                    Extrapolator1D<INDEPTYPE, DEPTYPE> extrap;

                    switch ( interptype ) {
                        case extrapolator_1d_type_t::CONSTANT:
                            extrap.set_engine(
                                extrapolator_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new _constant_extrapolator_1d<
                                        INDEPTYPE, DEPTYPE>() ) );
                            break;
                        case extrapolator_1d_type_t::ZERO:
                            extrap.set_engine(
                                extrapolator_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new _zero_extrapolator_1d<INDEPTYPE,
                                                              DEPTYPE>() ) );
                            break;
                        case extrapolator_1d_type_t::LINEAR:
                            extrap.set_engine(
                                extrapolator_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new _linear_extrapolator_1d<INDEPTYPE,
                                                                DEPTYPE>() ) );
                            break;
                        case extrapolator_1d_type_t::FORBIDDEN:
                            extrap.set_engine(
                                extrapolator_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new _forbidden_extrapolator_1d<
                                        INDEPTYPE, DEPTYPE>() ) );
                            break;
                        case extrapolator_1d_type_t::PERIODIC:
                            extrap.set_engine(
                                extrapolator_engine_1d_t<INDEPTYPE, DEPTYPE>(
                                    new _periodic_extrapolator_1d<
                                        INDEPTYPE, DEPTYPE>() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 1-D extrapolator type "
                                "unrecognized; either it is undefined, not "
                                "applicable, or you don't have the libraries "
                                "available." );
                    }
                    return extrap;
                }
        };

    }  // namespace interpolation
}  // namespace NCPA
