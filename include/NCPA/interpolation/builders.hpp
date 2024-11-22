#pragma once


#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/abstract_spline_1d.hpp"
#include "NCPA/interpolation/interpolator_1d.hpp"
#include "NCPA/interpolation/lanl_interpolators.hpp"
#include "NCPA/interpolation/gsl_interpolators.hpp"
#include "NCPA/interpolation/nearest_neighbor_spline_1d.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace interpolation {
        class InterpolatorFactory {
            public:
                static bool can_build( interpolator_type_t interptype ) {
                    switch (interptype) {
                        case interpolator_type_t::NEAREST_NEIGHBOR:
                        case interpolator_type_t::LANL_LINEAR:
                        case interpolator_type_t::LANL_CUBIC:
#ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE
                        case interpolator_type_t::GSL_LINEAR:
                        case interpolator_type_t::GSL_POLYNOMIAL:
                        case interpolator_type_t::GSL_CUBIC:
                        case interpolator_type_t::GSL_CUBIC_PERIODIC:
                        case interpolator_type_t::GSL_AKIMA:
                        case interpolator_type_t::GSL_AKIMA_PERIODIC:
#  if GSL_MAJOR_VERSION >= 2
                        case interpolator_type_t::GSL_STEFFEN:
#endif
#endif
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                template<typename INDEPTYPE, typename DEPTYPE>
                static Interpolator1D<INDEPTYPE, DEPTYPE> build(
                    interpolator_type_t interptype ) {
                    Interpolator1D<INDEPTYPE, DEPTYPE> interp;

                    switch ( interptype ) {
                        case interpolator_type_t::NEAREST_NEIGHBOR:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new nearest_neighbor_spline_1d<INDEPTYPE,
                                                           DEPTYPE>() ) );
                            break;
                        case interpolator_type_t::LANL_LINEAR:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new LANL::linear_spline_1d<INDEPTYPE,
                                                           DEPTYPE>() ) );
                            break;
                        case interpolator_type_t::LANL_CUBIC:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new LANL::natural_cubic_spline_1d<
                                    INDEPTYPE, DEPTYPE>() ) );
                            break;
#ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE
                        case interpolator_type_t::GSL_LINEAR:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                    gsl_interp_linear ) ) );
                            break;
                        case interpolator_type_t::GSL_POLYNOMIAL:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                    gsl_interp_polynomial ) ) );
                            break;
                        case interpolator_type_t::GSL_CUBIC:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                    gsl_interp_cspline ) ) );
                            break;
                        case interpolator_type_t::GSL_CUBIC_PERIODIC:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                    gsl_interp_cspline_periodic ) ) );
                            break;
                        case interpolator_type_t::GSL_AKIMA:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                    gsl_interp_akima ) ) );
                            break;
                        case interpolator_type_t::GSL_AKIMA_PERIODIC:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                    gsl_interp_akima_periodic ) ) );
                            break;
#  if GSL_MAJOR_VERSION >= 2
                        case interpolator_type_t::GSL_STEFFEN:
                            interp.set_engine( std::unique_ptr<
                                               _spline_1d<INDEPTYPE, DEPTYPE>>(
                                new GSL::gsl_spline_1d<INDEPTYPE, DEPTYPE>(
                                    gsl_interp_steffen ) ) );
                            break;
#  endif
#endif
                        default:
                            throw std::range_error(
                                "Requested interpolator type unrecognized; "
                                "either it is undefined or you don't have the "
                                "libraries available." );
                    }
                    return interp;
                }
        };
    }
}