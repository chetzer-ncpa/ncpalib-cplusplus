#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/abstract_spline.hpp"
#include "NCPA/interpolation/exceptions.hpp"
#include "NCPA/interpolation/gsl.hpp"
#include "NCPA/interpolation/lanl.hpp"
#include "NCPA/interpolation/NearestNeighborSpline.hpp"
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
                        case interpolator_type_t::NEAREST_NEIGHBOR_2D:
                        case interpolator_type_t::STRATIFIED_2D:
                        case interpolator_type_t::LANL_LINEAR:
                        case interpolator_type_t::LANL_CUBIC:
                        case interpolator_type_t::LANL_NATURAL_2D:
                        case interpolator_type_t::LANL_BICUBIC_2D:
                        // case interpolator_type_t::LANL_HYBRID_3D:
#ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE
                        case interpolator_type_t::GSL_LINEAR:
                        case interpolator_type_t::GSL_POLYNOMIAL:
                        case interpolator_type_t::GSL_CUBIC:
                        case interpolator_type_t::GSL_CUBIC_PERIODIC:
                        case interpolator_type_t::GSL_AKIMA:
                        case interpolator_type_t::GSL_AKIMA_PERIODIC:
#  if GSL_MAJOR_VERSION >= 2
                        case interpolator_type_t::GSL_STEFFEN:
#  endif
#endif
                            return true;
                        default:
                            return false;
                    }
                }

                // static bool can_build( interpolator_2d_type_t interptype ) {
                //     switch (interptype) {
                //         case interpolator_2d_type_t::NEAREST_NEIGHBOR:
                //         case interpolator_2d_type_t::LANL_NATURAL:
                //         case interpolator_2d_type_t::LANL_BICUBIC:
                //             return true;
                //             break;
                //         default:
                //             return false;
                //     }
                // }

                // static bool can_build( interpolator_3d_type_t interptype ) {
                //     switch (interptype) {
                //         case interpolator_3d_type_t::LANL_HYBRID:
                //             return true;
                //             break;
                //         default:
                //             return false;
                //     }
                // }

                template<typename FloatType = double>
                static spline_ptr_t<FloatType> build(
                    interpolator_type_t interptype ) {
                    switch (interptype) {
                        case interpolator_type_t::NEAREST_NEIGHBOR:
                            return spline_ptr_t<FloatType>(
                                new NearestNeighborSpline<FloatType>() );
                        case interpolator_type_t::LANL_LINEAR:
                            return spline_ptr_t<FloatType>(
                                new LANL::LinearSpline1D<FloatType>() );
                        case interpolator_type_t::LANL_CUBIC:
                            return spline_ptr_t<FloatType>(
                                new LANL::NaturalCubicSpline1D<FloatType>() );
                        case interpolator_type_t::LANL_NATURAL_2D:
                            return spline_ptr_t<FloatType>(
                                new LANL::NaturalSpline2D<FloatType>() );
                        case interpolator_type_t::LANL_BICUBIC_2D:
                            return spline_ptr_t<FloatType>(
                                new LANL::BicubicSpline2D<FloatType>() );

#if NCPA_USE_GSL_INTERPOLATION
                        case interpolator_type_t::GSL_LINEAR:
                        case interpolator_type_t::GSL_AKIMA:
                        case interpolator_type_t::GSL_AKIMA_PERIODIC:
                        case interpolator_type_t::GSL_CUBIC_PERIODIC:
                        case interpolator_type_t::GSL_CUBIC:
                        case interpolator_type_t::GSL_POLYNOMIAL:
                        case interpolator_type_t::GSL_STEFFEN:
                            return spline_ptr_t<FloatType>(
                                new GSL::GSLSpline<FloatType>(
                                    GSL::gsl_interpolator_type(
                                        interptype ) ) );
#endif
                        default:
                            throw UnsupportedInterpolator(
                                "Requested interpolator type "
                                "unrecognized; either it is undefined, not "
                                "applicable, or you don't have the libraries "
                                "available." );
                    }
                }

                // template<typename INDEPTYPE, typename DEPTYPE>
                //                 static Interpolator<FloatType> build(
                //                     interpolator_type_t interptype ) {
                //                     Interpolator<FloatType> interp;

                //                     switch (interptype) {
                //                         case
                //                         interpolator_type_t::NEAREST_NEIGHBOR:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 NearestNeighborSpline_1d<
                //                                     FloatType>() ) );
                //                             break;
                //                         case
                //                         interpolator_type_t::LANL_LINEAR:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 LANL::linear_spline_1d<FloatType>()
                //                                 ) );
                //                             break;
                //                         case
                //                         interpolator_type_t::LANL_CUBIC:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 LANL::natural_cubic_spline_1d<
                //                                     FloatType>() ) );
                //                             break;
                // #ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE
                //                         case
                //                         interpolator_type_t::GSL_LINEAR:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 GSL::gsl_spline_1d<FloatType>(
                //                                     gsl_interp_linear ) ) );
                //                             break;
                //                         case
                //                         interpolator_type_t::GSL_POLYNOMIAL:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 GSL::gsl_spline_1d<FloatType>(
                //                                     gsl_interp_polynomial )
                //                                     ) );
                //                             break;
                //                         case interpolator_type_t::GSL_CUBIC:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 GSL::gsl_spline_1d<FloatType>(
                //                                     gsl_interp_cspline ) )
                //                                     );
                //                             break;
                //                         case
                //                         interpolator_type_t::GSL_CUBIC_PERIODIC:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 GSL::gsl_spline_1d<FloatType>(
                //                                     gsl_interp_cspline_periodic
                //                                     ) ) );
                //                             break;
                //                         case interpolator_type_t::GSL_AKIMA:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 GSL::gsl_spline_1d<FloatType>(
                //                                     gsl_interp_akima ) ) );
                //                             break;
                //                         case
                //                         interpolator_type_t::GSL_AKIMA_PERIODIC:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 GSL::gsl_spline_1d<FloatType>(
                //                                     gsl_interp_akima_periodic
                //                                     ) ) );
                //                             break;
                // #  if GSL_MAJOR_VERSION >= 2
                //                         case
                //                         interpolator_type_t::GSL_STEFFEN:
                //                             interp.set_engine(
                //                             spline_engine_1d_t<FloatType>(
                //                                 new
                //                                 GSL::gsl_spline_1d<FloatType>(
                //                                     gsl_interp_steffen ) )
                //                                     );
                //                             break;
                // #  endif
                // #endif
                //                         default:
                //                             throw std::range_error(
                //                                 "Requested interpolator type
                //                                 " "unrecognized; either it
                //                                 is undefined, not "
                //                                 "applicable, or you don't
                //                                 have the libraries "
                //                                 "available." );
                //                     }
                //                     return interp;
                //                 }

                // template<typename INDEPTYPE, typename DEPTYPE>
                //         static Interpolator2D<FloatType> build(
                //             interpolator_2d_type_t interptype ) {
                //             Interpolator2D<FloatType> interp;

                //             switch (interptype) {
                //                 case
                //                 interpolator_2d_type_t::NEAREST_NEIGHBOR:
                //                     interp.set_engine(
                //                         spline_engine_2d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new NearestNeighborSpline_2d<
                //                                 FloatType>() ) );
                //                     break;
                //                 case interpolator_2d_type_t::LANL_NATURAL:
                //                     interp.set_engine(
                //                         spline_engine_2d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new
                //                             LANL::natural_spline_2d<INDEPTYPE,
                //                                                         DEPTYPE>()
                //                                                         ) );
                //                     break;
                //                 case interpolator_2d_type_t::LANL_BICUBIC:
                //                     interp.set_engine(
                //                         spline_engine_2d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new
                //                             LANL::bicubic_spline_2d<INDEPTYPE,
                //                                                         DEPTYPE>()
                //                                                         ) );
                //                     break;
                //                 default:
                //                     throw std::range_error(
                //                         "Requested 2-D interpolator type "
                //                         "unrecognized; either it is
                //                         undefined, not " "applicable, you
                //                         don't have the libraries "
                //                         "available, or you need to use a
                //                         different " "build() method." );
                //             }
                //             return interp;
                //         }

                //         static Interpolator2D<FloatType> build(
                //             interpolator_2d_type_t interp2type,
                //             interpolator_type_t interp1type,
                //             size_t stratified_axis = 0 ) {
                //             Interpolator2D<FloatType> interp;
                //             stratified_axis_type_t params( stratified_axis,
                //                                            interp1type );

                //             switch (interp2type) {
                //                 case interpolator_2d_type_t::STRATIFIED:
                //                     interp.set_engine(
                //                         spline_engine_2d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new stratified_spline_2d<
                //                                 FloatType,
                //                                 stratified_axis_type_t>(
                //                                 params ) ) );
                //                     break;
                //                 default:
                //                     throw std::range_error(
                //                         "Requested 2-D interpolator type "
                //                         "unrecognized; either it is
                //                         undefined, not " "applicable, you
                //                         don't have the libraries "
                //                         "available, or you need to use a
                //                         different " "build() method."  );
                //             }
                //             return interp;
                //         }

                //         // template<typename INDEPTYPE, typename DEPTYPE>
                //         static Interpolator3D<FloatType> build(
                //             interpolator_3d_type_t interptype ) {
                //             Interpolator3D<FloatType> interp;

                //             switch (interptype) {
                //                 case interpolator_3d_type_t::LANL_HYBRID:
                //                     interp.set_engine(
                //                         spline_engine_3d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new
                //                             LANL::hybrid_spline_3d<INDEPTYPE,
                //                                                        DEPTYPE>()
                //                                                        ) );
                //                     break;
                //                 default:
                //                     std::ostringstream oss;
                //                     oss << "Requested 3-D interpolator type
                //                     "
                //                            "unrecognized; either it is
                //                            undefined, " "not applicable, or
                //                            you don't have the " "libraries
                //                            available.";
                //                     throw std::range_error( oss.str() );
                //             }
                //             return interp;
                //         }

                //         static Extrapolator1D<FloatType> build(
                //             extrapolator_1d_type_t interptype ) {
                //             Extrapolator1D<FloatType> extrap;

                //             switch (interptype) {
                //                 case extrapolator_1d_type_t::CONSTANT:
                //                     extrap.set_engine(
                //                         extrapolator_engine_1d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new _constant_extrapolator_1d<
                //                                 FloatType>() ) );
                //                     break;
                //                 case extrapolator_1d_type_t::ZERO:
                //                     extrap.set_engine(
                //                         extrapolator_engine_1d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new
                //                             _zero_extrapolator_1d<INDEPTYPE,
                //                                                       DEPTYPE>()
                //                                                       ) );
                //                     break;
                //                 case extrapolator_1d_type_t::LINEAR:
                //                     extrap.set_engine(
                //                         extrapolator_engine_1d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new
                //                             _linear_extrapolator_1d<INDEPTYPE,
                //                                                         DEPTYPE>()
                //                                                         ) );
                //                     break;
                //                 case extrapolator_1d_type_t::FORBIDDEN:
                //                     extrap.set_engine(
                //                         extrapolator_engine_1d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new _forbidden_extrapolator_1d<
                //                                 FloatType>() ) );
                //                     break;
                //                 case extrapolator_1d_type_t::PERIODIC:
                //                     extrap.set_engine(
                //                         extrapolator_engine_1d_t<INDEPTYPE,
                //                         DEPTYPE>(
                //                             new _periodic_extrapolator_1d<
                //                                 FloatType>() ) );
                //                     break;
                //                 default:
                //                     throw std::range_error(
                //                         "Requested 1-D extrapolator type "
                //                         "unrecognized; either it is
                //                         undefined, not " "applicable, or you
                //                         don't have the libraries "
                //                         "available." );
                //             }
                //             return extrap;
                //         }
        };

    }  // namespace interpolation
}  // namespace NCPA
