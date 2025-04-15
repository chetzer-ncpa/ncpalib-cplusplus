#pragma once

#include "NCPA/defines.hpp"
#include "NCPA/interpolation/Interpolator1D.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"

#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace interpolation {
        template<typename INDEPTYPE, typename DEPTYPE>
        class _constant_extrapolator_1d;
        // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _constant_extrapolator_1d,
        //                                        _extrapolator_1d );
        // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _linear_extrapolator_1d,
        //                                        _extrapolator_1d );
        // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _zero_extrapolator_1d,
        //                                        _extrapolator_1d );
        // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _periodic_extrapolator_1d,
        //                                        _extrapolator_1d );
        // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _forbidden_extrapolator_1d,
        //                                        _extrapolator_1d );

        template<typename INDEPTYPE, typename DEPTYPE>
        class _linear_extrapolator_1d;
        template<typename INDEPTYPE, typename DEPTYPE>
        class _zero_extrapolator_1d;
        template<typename INDEPTYPE, typename DEPTYPE>
        class _forbidden_extrapolator_1d;
        template<typename INDEPTYPE, typename DEPTYPE>
        class _periodic_extrapolator_1d;
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_abstract_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_abstract_extrapolator_1d<T, U>& b ) noexcept;
// template<typename T, typename U>
// static void swap( NCPA::interpolation::_extrapolator_1d<T, U>& a,
//                   NCPA::interpolation::_extrapolator_1d<T, U>& b ) noexcept;
template<typename T, typename U>
static void swap(
    NCPA::interpolation::_constant_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_constant_extrapolator_1d<T, U>& b ) noexcept;
template<typename T, typename U>
static void swap(
    NCPA::interpolation::_linear_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_linear_extrapolator_1d<T, U>& b ) noexcept;
template<typename T, typename U>
static void swap(
    NCPA::interpolation::_zero_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_zero_extrapolator_1d<T, U>& b ) noexcept;
template<typename T, typename U>
static void swap(
    NCPA::interpolation::_forbidden_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_forbidden_extrapolator_1d<T, U>& b ) noexcept;
template<typename T, typename U>
static void swap(
    NCPA::interpolation::_periodic_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_periodic_extrapolator_1d<T, U>& b ) noexcept;
template<typename T, typename U>
static void swap( NCPA::interpolation::Extrapolator1D<T, U>& a,
                  NCPA::interpolation::Extrapolator1D<T, U>& b ) noexcept;

#define NCPA_EXTRAPOLATION_TERNARY_RETURN( _X_, _LOWLIMIT_, _HIGHLIMIT_,    \
                                           _LOWEXPR_, _HIGHEXPR_ )          \
    if ( _X_ <= _LOWLIMIT_ ) {                                              \
        return _LOWEXPR_;                                                   \
    } else if ( _X_ >= _HIGHLIMIT_ ) {                                      \
        return _HIGHEXPR_;                                                  \
    } else {                                                                \
        throw std::logic_error( "Tried to extrapolate but point is within " \
                                "interpolator limits!" );                   \
    }

namespace NCPA {
    namespace interpolation {
        // template<typename INDEPTYPE, typename DEPTYPE>
        // struct _interpolator_endpoint_t {
        //         INDEPTYPE x  = NCPA::math::zero<INDEPTYPE>();
        //         DEPTYPE f    = NCPA::math::zero<DEPTYPE>(),
        //                 df   = NCPA::math::zero<DEPTYPE>(),
        //                 ddf  = NCPA::math::zero<DEPTYPE>(),
        //                 dddf = NCPA::math::zero<DEPTYPE>();
        // };

        template<typename INDEPTYPE, typename DEPTYPE>
        class _abstract_extrapolator_1d {
            public:
                virtual ~_abstract_extrapolator_1d() {}

                DECLARE_FRIEND_SWAP_METHOD_TEMPLATE2( _abstract_extrapolator_1d, INDEPTYPE, DEPTYPE );

                virtual DEPTYPE extrapolate(
                    INDEPTYPE x, Interpolator1D<INDEPTYPE, DEPTYPE>& interp )
                    = 0;
                virtual DEPTYPE extrapolate_df(
                    INDEPTYPE x, Interpolator1D<INDEPTYPE, DEPTYPE>& interp )
                    = 0;
                virtual DEPTYPE extrapolate_ddf(
                    INDEPTYPE x, Interpolator1D<INDEPTYPE, DEPTYPE>& interp )
                    = 0;
                virtual DEPTYPE extrapolate_dddf(
                    INDEPTYPE x, Interpolator1D<INDEPTYPE, DEPTYPE>& interp )
                    = 0;

            // protected:
            //     virtual void _endpoints(
            //         Interpolator1D<INDEPTYPE, DEPTYPE>& interp,
            //         _interpolator_endpoint_t& minendpoint,
            //         _interpolator_endpoint_t& maxendpoint ) {
            //         std::pair<INDEPTYPE, INDEPTYPE> lim = interp.limits();
            //         minendpoint.x = lim.first;
            //         minendpoint.f = interp.eval_f( lim.first );
            //         minendpoint.df = interp.eval_df( lim.first );
            //         try {
            //             minendpoint.ddf = interp.eval_ddf( lim.first );
            //         } catch (NCPA::NotImplementedError e) {
            //             minendpoint.ddf = NCPA::math::zero<DEPTYPE>();
            //         }
            //         try {
            //             minendpoint.dddf = interp.eval_dddf( lim.first );
            //         } catch (NCPA::NotImplementedError e) {
            //             minendpoint.dddf = NCPA::math::zero<DEPTYPE>();
            //         }
            //         maxendpoint.x = lim.second;
            //         maxendpoint.f = interp.eval_f( lim.second );
            //         maxendpoint.df = interp.eval_df( lim.second );
            //         try {
            //             maxendpoint.ddf = interp.eval_ddf( lim.second );
            //         } catch (NCPA::NotImplementedError e) {
            //             maxendpoint.ddf = NCPA::math::zero<DEPTYPE>();
            //         }
            //         try {
            //             maxendpoint.dddf = interp.eval_dddf( lim.second );
            //         } catch (NCPA::NotImplementedError e) {
            //             maxendpoint.dddf = NCPA::math::zero<DEPTYPE>();
            //         }
            //     }
        };

        // _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
        //     class _extrapolator_1d<INDEPTYPE, DEPTYPE,
        //                            _ENABLE_IF_INDEP_IS_REAL,
        //                            _ENABLE_IF_DEP_IS_REAL>
        //     : public _abstract_extrapolator_1d<INDEPTYPE, DEPTYPE> {
        //     public:
        //         virtual ~_extrapolator_1d() {}
        // };

        // _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
        //     class _extrapolator_1d<INDEPTYPE, DEPTYPE,
        //                            _ENABLE_IF_INDEP_IS_REAL,
        //                            _ENABLE_IF_DEP_IS_COMPLEX>
        //     : public _abstract_extrapolator_1d<INDEPTYPE, DEPTYPE> {
        //     public:
        //         virtual ~_extrapolator_1d() {}

        //         friend void ::swap<INDEPTYPE, DEPTYPE>(
        //             _extrapolator_1d<INDEPTYPE, DEPTYPE>& a,
        //             _extrapolator_1d<INDEPTYPE, DEPTYPE>& b ) noexcept;

        //         virtual DEPTYPE extrapolate(
        //             INDEPTYPE x,
        //             Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
        //             return DEPTYPE( real()->extrapolate( x, interp ),
        //                             imag()->extrapolate( x, interp ) );
        //         }

        //         virtual DEPTYPE extrapolate_df(
        //             INDEPTYPE x,
        //             Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
        //             return DEPTYPE( real()->extrapolate_df( x, interp ),
        //                             imag()->extrapolate_df( x, interp ) );
        //         }

        //         virtual DEPTYPE extrapolate_ddf(
        //             INDEPTYPE x,
        //             Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
        //             return DEPTYPE( real()->extrapolate_ddf( x, interp ),
        //                             imag()->extrapolate_ddf( x, interp ) );
        //         }

        //         virtual DEPTYPE extrapolate_dddf(
        //             INDEPTYPE x,
        //             Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
        //             return DEPTYPE( real()->extrapolate_dddf( x, interp ),
        //                             imag()->extrapolate_dddf( x, interp ) );
        //         }

        //         virtual _SUBEXTRAPOLATOR_PTR_T real()             = 0;
        //         virtual const _SUBEXTRAPOLATOR_PTR_T real() const = 0;
        //         virtual _SUBEXTRAPOLATOR_PTR_T imag()             = 0;
        //         virtual const _SUBEXTRAPOLATOR_PTR_T imag() const = 0;
        // };

        template<typename INDEPTYPE, typename DEPTYPE>
        class _constant_extrapolator_1d
            : public _abstract_extrapolator_1d<INDEPTYPE, DEPTYPE> {
        // template<typename INDEPTYPE, typename DEPTYPE>
        // class _constant_extrapolator_1d<INDEPTYPE, DEPTYPE, void,
        //                                 ENABLE_IF_REAL( INDEPTYPE ),
        //                                 ENABLE_IF_REAL( DEPTYPE )>
        //     : public _extrapolator_1d<INDEPTYPE, DEPTYPE> {
            public:
                _constant_extrapolator_1d() {}

                _constant_extrapolator_1d(
                    const _constant_extrapolator_1d<INDEPTYPE, DEPTYPE>&
                        source ) :
                    _constant_extrapolator_1d<INDEPTYPE, DEPTYPE>() {}

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2(
                    _constant_extrapolator_1d, INDEPTYPE, DEPTYPE )

                virtual DEPTYPE extrapolate(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( interp.eval_f( limits.first ) ),
                        ( interp.eval_f( limits.second ) ) );
                }

                virtual DEPTYPE extrapolate_df(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( NCPA::math::zero<DEPTYPE>() ),
                        ( NCPA::math::zero<DEPTYPE>() ) );
                }

                virtual DEPTYPE extrapolate_ddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( NCPA::math::zero<DEPTYPE>() ),
                        ( NCPA::math::zero<DEPTYPE>() ) );
                }

                virtual DEPTYPE extrapolate_dddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( NCPA::math::zero<DEPTYPE>() ),
                        ( NCPA::math::zero<DEPTYPE>() ) );
                }
        };

        template<typename INDEPTYPE, typename DEPTYPE>
        class _zero_extrapolator_1d
            : public _abstract_extrapolator_1d<INDEPTYPE, DEPTYPE> {
        // template<typename INDEPTYPE, typename DEPTYPE>
        // class _zero_extrapolator_1d<INDEPTYPE, DEPTYPE, void,
        //                             ENABLE_IF_REAL( INDEPTYPE ),
        //                             ENABLE_IF_REAL( DEPTYPE )>
        //     : public _extrapolator_1d<INDEPTYPE, DEPTYPE> {
            public:
                _zero_extrapolator_1d() {}

                _zero_extrapolator_1d(
                    const _zero_extrapolator_1d<INDEPTYPE, DEPTYPE>& source ) :
                    _zero_extrapolator_1d<INDEPTYPE, DEPTYPE>() {}

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2( _zero_extrapolator_1d,
                                                        INDEPTYPE, DEPTYPE )

                virtual DEPTYPE extrapolate(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( NCPA::math::zero<DEPTYPE>() ),
                        ( NCPA::math::zero<DEPTYPE>() ) );
                }

                virtual DEPTYPE extrapolate_df(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( NCPA::math::zero<DEPTYPE>() ),
                        ( NCPA::math::zero<DEPTYPE>() ) );
                }

                virtual DEPTYPE extrapolate_ddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( NCPA::math::zero<DEPTYPE>() ),
                        ( NCPA::math::zero<DEPTYPE>() ) );
                }

                virtual DEPTYPE extrapolate_dddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( NCPA::math::zero<DEPTYPE>() ),
                        ( NCPA::math::zero<DEPTYPE>() ) );
                }
        };

        template<typename INDEPTYPE, typename DEPTYPE>
        class _linear_extrapolator_1d
            : public _abstract_extrapolator_1d<INDEPTYPE, DEPTYPE> {
        // template<typename INDEPTYPE, typename DEPTYPE>
        // class _linear_extrapolator_1d<INDEPTYPE, DEPTYPE, void,
        //                               ENABLE_IF_REAL( INDEPTYPE ),
        //                               ENABLE_IF_REAL( DEPTYPE )>
        //     : public _extrapolator_1d<INDEPTYPE, DEPTYPE> {
            public:
                _linear_extrapolator_1d() {}

                _linear_extrapolator_1d(
                    const _linear_extrapolator_1d<INDEPTYPE, DEPTYPE>&
                        source ) :
                    _linear_extrapolator_1d<INDEPTYPE, DEPTYPE>() {}

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2(
                    _linear_extrapolator_1d, INDEPTYPE, DEPTYPE )

                virtual DEPTYPE extrapolate(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( interp.eval_f( limits.first )
                          + ( x - limits.first )
                                * interp.eval_df( limits.first ) ),
                        ( interp.eval_f( limits.second )
                          + ( x - limits.second )
                                * interp.eval_df( limits.second ) ) );
                }

                virtual DEPTYPE extrapolate_df(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( interp.eval_df( limits.first ) ),
                        ( interp.eval_df( limits.second ) ) );
                }

                virtual DEPTYPE extrapolate_ddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( NCPA::math::zero<DEPTYPE>() ),
                        ( NCPA::math::zero<DEPTYPE>() ) );
                }

                virtual DEPTYPE extrapolate_dddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    NCPA_EXTRAPOLATION_TERNARY_RETURN(
                        x, limits.first, limits.second,
                        ( NCPA::math::zero<DEPTYPE>() ),
                        ( NCPA::math::zero<DEPTYPE>() ) );
                }
        };

        template<typename INDEPTYPE, typename DEPTYPE>
        class _forbidden_extrapolator_1d
            : public _abstract_extrapolator_1d<INDEPTYPE, DEPTYPE> {
        // template<typename INDEPTYPE, typename DEPTYPE>
        // class _forbidden_extrapolator_1d<INDEPTYPE, DEPTYPE, void,
        //                                  ENABLE_IF_REAL( INDEPTYPE ),
        //                                  ENABLE_IF_REAL( DEPTYPE )>
        //     : public _extrapolator_1d<INDEPTYPE, DEPTYPE> {
            public:
                _forbidden_extrapolator_1d() {}

                _forbidden_extrapolator_1d(
                    const _forbidden_extrapolator_1d<INDEPTYPE, DEPTYPE>&
                        source ) :
                    _forbidden_extrapolator_1d<INDEPTYPE, DEPTYPE>() {}

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2(
                    _forbidden_extrapolator_1d, INDEPTYPE, DEPTYPE )

                virtual DEPTYPE extrapolate(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    throw std::logic_error( "Extrapolation forbidden!" );
                }

                virtual DEPTYPE extrapolate_df(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    throw std::logic_error( "Extrapolation forbidden!" );
                }

                virtual DEPTYPE extrapolate_ddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    throw std::logic_error( "Extrapolation forbidden!" );
                }

                virtual DEPTYPE extrapolate_dddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    throw std::logic_error( "Extrapolation forbidden!" );
                }
        };

        template<typename INDEPTYPE, typename DEPTYPE>
        class _periodic_extrapolator_1d
            : public _abstract_extrapolator_1d<INDEPTYPE, DEPTYPE> {
        // template<typename INDEPTYPE, typename DEPTYPE>
        // class _periodic_extrapolator_1d<INDEPTYPE, DEPTYPE, void,
        //                                 ENABLE_IF_REAL( INDEPTYPE ),
        //                                 ENABLE_IF_REAL( DEPTYPE )>
        //     : public _extrapolator_1d<INDEPTYPE, DEPTYPE> {
            public:
                _periodic_extrapolator_1d() {}

                _periodic_extrapolator_1d(
                    const _periodic_extrapolator_1d<INDEPTYPE, DEPTYPE>&
                        source ) :
                    _periodic_extrapolator_1d<INDEPTYPE, DEPTYPE>() {}

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2(
                    _periodic_extrapolator_1d, INDEPTYPE, DEPTYPE )

                virtual DEPTYPE extrapolate(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    if ( x >= limits.first && x <= limits.second ) {
                        throw std::logic_error(
                            "Tried to extrapolate but point is within "
                            "interpolator limits!" );
                    }
                    INDEPTYPE dlimits = limits.second - limits.first;
                    while ( x < limits.first ) {
                        x += dlimits;
                    }
                    while ( x > limits.second ) {
                        x -= dlimits;
                    }
                    return interp.eval_f( x );
                }

                virtual DEPTYPE extrapolate_df(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    if ( x >= limits.first && x <= limits.second ) {
                        throw std::logic_error(
                            "Tried to extrapolate but point is within "
                            "interpolator limits!" );
                    }
                    INDEPTYPE dlimits = limits.second - limits.first;
                    while ( x < limits.first ) {
                        x += dlimits;
                    }
                    while ( x > limits.second ) {
                        x -= dlimits;
                    }
                    return interp.eval_df( x );
                }

                virtual DEPTYPE extrapolate_ddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    if ( x >= limits.first && x <= limits.second ) {
                        throw std::logic_error(
                            "Tried to extrapolate but point is within "
                            "interpolator limits!" );
                    }
                    INDEPTYPE dlimits = limits.second - limits.first;
                    while ( x < limits.first ) {
                        x += dlimits;
                    }
                    while ( x > limits.second ) {
                        x -= dlimits;
                    }
                    return interp.eval_ddf( x );
                }

                virtual DEPTYPE extrapolate_dddf(
                    INDEPTYPE x,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::pair<INDEPTYPE, INDEPTYPE> limits = interp.limits();
                    if ( x >= limits.first && x <= limits.second ) {
                        throw std::logic_error(
                            "Tried to extrapolate but point is within "
                            "interpolator limits!" );
                    }
                    INDEPTYPE dlimits = limits.second - limits.first;
                    while ( x < limits.first ) {
                        x += dlimits;
                    }
                    while ( x > limits.second ) {
                        x -= dlimits;
                    }
                    return interp.eval_dddf( x );
                }
        };

        // DEFINE_COMPLEX_VERSION_OF_EXTRAPOLATOR( _constant_extrapolator_1d,
        //                                         _extrapolator_1d )
        // DEFINE_COMPLEX_VERSION_OF_EXTRAPOLATOR( _zero_extrapolator_1d,
        //                                         _extrapolator_1d )
        // DEFINE_COMPLEX_VERSION_OF_EXTRAPOLATOR( _linear_extrapolator_1d,
        //                                         _extrapolator_1d )
        // DEFINE_COMPLEX_VERSION_OF_EXTRAPOLATOR( _forbidden_extrapolator_1d,
        //                                         _extrapolator_1d )
        // DEFINE_COMPLEX_VERSION_OF_EXTRAPOLATOR( _periodic_extrapolator_1d,
        //                                         _extrapolator_1d )

        template<typename INDEPTYPE, typename DEPTYPE>
        class Extrapolator1D {
            public:
                Extrapolator1D() {}

                Extrapolator1D(
                    extrapolator_engine_1d_t<INDEPTYPE, DEPTYPE> engine ) {
                    set_engine( engine );
                }

                Extrapolator1D(
                    const Extrapolator1D<INDEPTYPE, DEPTYPE>& other ) :
                    Extrapolator1D<INDEPTYPE, DEPTYPE>() {
                    _engine = other._engine;
                }

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2( Extrapolator1D,
                                                        INDEPTYPE, DEPTYPE )

                void set_engine(
                    extrapolator_engine_1d_t<INDEPTYPE, DEPTYPE> engine ) {
                    _engine = std::move( engine );
                }

                virtual DEPTYPE extrapolate(
                    INDEPTYPE x, Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) {
                    check_engine();
                    return _engine->extrapolate( x, interp );
                }

                virtual DEPTYPE extrapolate_df(
                    INDEPTYPE x, Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) {
                    check_engine();
                    return _engine->extrapolate_df( x, interp );
                }

                virtual DEPTYPE extrapolate_ddf(
                    INDEPTYPE x, Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) {
                    check_engine();
                    return _engine->extrapolate_ddf( x, interp );
                }

                virtual DEPTYPE extrapolate_dddf(
                    INDEPTYPE x, Interpolator1D<INDEPTYPE, DEPTYPE>& interp ) {
                    check_engine();
                    return _engine->extrapolate_dddf( x, interp );
                }

                virtual void check_engine() {
                    if ( !_engine ) {
                        throw std::logic_error(
                            "Extrapolator1D: no engine has been set" );
                    }
                }

                explicit operator bool() const {
                    return ( _engine ? true : false );
                }

            protected:
                extrapolator_engine_1d_t<INDEPTYPE, DEPTYPE> _engine;
        };
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_abstract_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_abstract_extrapolator_1d<T, U>& b ) noexcept {}

// template<typename T, typename U>
// static void swap( NCPA::interpolation::_extrapolator_1d<T, U>& a,
//                   NCPA::interpolation::_extrapolator_1d<T, U>& b ) noexcept {
//     ::swap(
//         dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>(
//             a ),
//         dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>(
//             b ) );
// }

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_constant_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_constant_extrapolator_1d<T, U>& b ) noexcept {
    ::swap( dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( a ),
            dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( b ) );
}

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_linear_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_linear_extrapolator_1d<T, U>& b ) noexcept {
    ::swap( dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( a ),
            dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( b ) );
}

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_zero_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_zero_extrapolator_1d<T, U>& b ) noexcept {
    ::swap( dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( a ),
            dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( b ) );
}

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_forbidden_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_forbidden_extrapolator_1d<T, U>& b ) noexcept {
    ::swap( dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( a ),
            dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( b ) );
}

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_periodic_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_periodic_extrapolator_1d<T, U>& b ) noexcept {
    ::swap( dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( a ),
            dynamic_cast<NCPA::interpolation::_abstract_extrapolator_1d<T, U>&>( b ) );
}

template<typename T, typename U>
static void swap( NCPA::interpolation::Extrapolator1D<T, U>& a,
                  NCPA::interpolation::Extrapolator1D<T, U>& b ) noexcept {
    using std::swap;
    swap( a._engine, b._engine );
}
