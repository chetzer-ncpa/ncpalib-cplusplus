#pragma once

#include "NCPA/defines.hpp"
#include "NCPA/interpolation/Interpolator1D.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include <array>
#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace interpolation {
        template<typename INDEPTYPE, typename DEPTYPE>
        class _constant_extrapolator_2d;
        template<typename INDEPTYPE, typename DEPTYPE>
        class _zero_extrapolator_2d;
        template<typename INDEPTYPE, typename DEPTYPE>
        class _linear_extrapolator_2d;
        template<typename INDEPTYPE, typename DEPTYPE>
        class _forbidden_extrapolator_2d;
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_abstract_extrapolator_2d<T, U>& a,
    NCPA::interpolation::_abstract_extrapolator_2d<T, U>& b ) noexcept;
template<typename T, typename U>
static void swap( NCPA::interpolation::Extrapolator2D<T, U>& a,
                  NCPA::interpolation::Extrapolator2D<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        template<typename INDEPTYPE, typename DEPTYPE>
        class _abstract_extrapolator_2d {
            public:
                virtual ~_abstract_extrapolator_2d() {}

                DECLARE_FRIEND_SWAP_METHOD( _abstract_extrapolator_2d );

                virtual DEPTYPE extrapolate(
                    INDEPTYPE x, INDEPTYPE y,
                    Interpolator2D<INDEPTYPE, DEPTYPE>& interp )
                    = 0;
                virtual DEPTYPE extrapolate_df(
                    INDEPTYPE x, INDEPTYPE y, size_t rel1,
                    Interpolator2D<INDEPTYPE, DEPTYPE>& interp )
                    = 0;
                virtual DEPTYPE extrapolate_ddf(
                    INDEPTYPE x, INDEPTYPE y, size_t rel1, size_t rel2,
                    Interpolator2D<INDEPTYPE, DEPTYPE>& interp )
                    = 0;
                virtual DEPTYPE extrapolate_dddf(
                    INDEPTYPE x, INDEPTYPE y, size_t rel1, size_t rel2,
                    size_t rel3, Interpolator2D<INDEPTYPE, DEPTYPE>& interp )
                    = 0;

            protected:
                virtual std::pair<int, int> _check_limits(
                    INDEPTYPE x, INDEPTYPE y,
                    const std::array<INDEPTYPE, 4>& lims ) const {
                        std::pair<int,int> inlimits({0,0});
                        if (x < lims[0]) {
                            inlimits.first = -1;
                        }
                        if (x > lims[1]) {
                            inlimits.first = 1;
                        }
                        if (y < lims[2]) {
                            inlimits.second = -1;
                        }
                        if (y > lims[3]) {
                            inlimits.second = 1;
                        }
                        return inlimits;
                    }
        };

        template<typename INDEPTYPE, typename DEPTYPE>
        class _constant_extrapolator_2d
            : public _abstract_extrapolator_2d<INDEPTYPE, DEPTYPE> {
            public:
                _constant_extrapolator_2d() {}

                _constant_extrapolator_2d(
                    const _constant_extrapolator_2d<INDEPTYPE, DEPTYPE>&
                        source ) :
                    _constant_extrapolator_2d<INDEPTYPE, DEPTYPE>() {}

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2(
                    _constant_extrapolator_2d, INDEPTYPE, DEPTYPE )

                virtual DEPTYPE extrapolate(
                    INDEPTYPE x,
                    INDEPTYPE y,
                    Interpolator2D<INDEPTYPE, DEPTYPE>& interp ) override {
                    std::array<INDEPTYPE, 4> limits = interp.limits();
                    std::pair<int,int> limstatus = this->_check_limits(x, y, limits);
                    if (limstatus.first == 0 && limstatus.second == 0) {
                        throw std::logic_error( "Tried to extrapolate but point is within "  
                            "interpolator limits!" ); 
                    }
                    INDEPTYPE use_coords[2];
                    if (limstatus.first == 0) {
                        use_coords[0] = x;
                        use_coords[1] = (limstatus.second == -1 ? limits[2] : limits[3] );
                    } else  {
                        use_coords[0] = (limstatus.first == -1 ? limits[0] : limits[1]);
                        if (limstatus.second == -1) {
                            use_coords[1] = limits[2];
                        } else if (limstatus.second == 1) {
                            use_coords[1] = limits[3];
                        } else {
                            use_coords[1] = y;
                        }
                    }
                    return interp.eval_f( use_coords[0], use_coords[1] );
                }

                virtual DEPTYPE extrapolate_df(
                    INDEPTYPE x, INDEPTYPE y, size_t rel1,
                    Interpolator2D<INDEPTYPE, DEPTYPE>& interp ) override {
                        std::array<INDEPTYPE, 4> limits = interp.limits();
                        std::pair<int,int> limstatus = this->_check_limits(x, y, limits);
                        if (limstatus.first == 0 && limstatus.second == 0) {
                            throw std::logic_error( "Tried to extrapolate but point is within "  
                                "interpolator limits!" ); 
                        }
                        INDEPTYPE use_coords[2];
                        if (limstatus.first == 0) {
                            use_coords[0] = x;
                            use_coords[1] = (limstatus.second == -1 ? limits[2] : limits[3] );
                        } else  {
                            use_coords[0] = (limstatus.first == -1 ? limits[0] : limits[1]);
                            if (limstatus.second == -1) {
                                use_coords[1] = limits[2];
                            } else if (limstatus.second == 1) {
                                use_coords[1] = limits[3];
                            } else {
                                use_coords[1] = y;
                            }
                        }
                        return interp.eval_df( use_coords[0], use_coords[1], rel1 );
                }

                virtual DEPTYPE extrapolate_ddf(
                    INDEPTYPE x, INDEPTYPE y, size_t rel1, size_t rel2,
                    Interpolator2D<INDEPTYPE, DEPTYPE>& interp ) override {
                        std::array<INDEPTYPE, 4> limits = interp.limits();
                        std::pair<int,int> limstatus = this->_check_limits(x, y, limits);
                        if (limstatus.first == 0 && limstatus.second == 0) {
                            throw std::logic_error( "Tried to extrapolate but point is within "  
                                "interpolator limits!" ); 
                        }
                        INDEPTYPE use_coords[2];
                        if (limstatus.first == 0) {
                            use_coords[0] = x;
                            use_coords[1] = (limstatus.second == -1 ? limits[2] : limits[3] );
                        } else  {
                            use_coords[0] = (limstatus.first == -1 ? limits[0] : limits[1]);
                            if (limstatus.second == -1) {
                                use_coords[1] = limits[2];
                            } else if (limstatus.second == 1) {
                                use_coords[1] = limits[3];
                            } else {
                                use_coords[1] = y;
                            }
                        }
                        return interp.eval_ddf( use_coords[0], use_coords[1], rel1, rel2 );
                }

                virtual DEPTYPE extrapolate_dddf(
                    INDEPTYPE x, INDEPTYPE y, size_t rel1, size_t rel2, size_t rel3,
                    Interpolator2D<INDEPTYPE, DEPTYPE>& interp ) override {
                        std::array<INDEPTYPE, 4> limits = interp.limits();
                        std::pair<int,int> limstatus = this->_check_limits(x, y, limits);
                        if (limstatus.first == 0 && limstatus.second == 0) {
                            throw std::logic_error( "Tried to extrapolate but point is within "  
                                "interpolator limits!" ); 
                        }
                        INDEPTYPE use_coords[2];
                        if (limstatus.first == 0) {
                            use_coords[0] = x;
                            use_coords[1] = (limstatus.second == -1 ? limits[2] : limits[3] );
                        } else  {
                            use_coords[0] = (limstatus.first == -1 ? limits[0] : limits[1]);
                            if (limstatus.second == -1) {
                                use_coords[1] = limits[2];
                            } else if (limstatus.second == 1) {
                                use_coords[1] = limits[3];
                            } else {
                                use_coords[1] = y;
                            }
                        }
                        return interp.eval_dddf( use_coords[0], use_coords[1], rel1, rel2, rel3 );
                }
        };

        template<typename INDEPTYPE, typename DEPTYPE>
        class _zero_extrapolator_2d
            : public _abstract_extrapolator_2d<INDEPTYPE, DEPTYPE> {
            public:
                _zero_extrapolator_2d() {}

                _zero_extrapolator_2d(
                    const _zero_extrapolator_2d<INDEPTYPE, DEPTYPE>& source ) :
                    _zero_extrapolator_2d<INDEPTYPE, DEPTYPE>() {}

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2( _zero_extrapolator_2d,
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
        class _linear_extrapolator_2d
            : public _abstract_extrapolator_2d<INDEPTYPE, DEPTYPE> {
            public:
                _linear_extrapolator_2d() {}

                _linear_extrapolator_2d(
                    const _linear_extrapolator_2d<INDEPTYPE, DEPTYPE>&
                        source ) :
                    _linear_extrapolator_2d<INDEPTYPE, DEPTYPE>() {}

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2(
                    _linear_extrapolator_2d, INDEPTYPE, DEPTYPE )

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
        class _forbidden_extrapolator_2d
            : public _abstract_extrapolator_2d<INDEPTYPE, DEPTYPE> {
            public:
                _forbidden_extrapolator_2d() {}

                _forbidden_extrapolator_2d(
                    const _forbidden_extrapolator_2d<INDEPTYPE, DEPTYPE>&
                        source ) :
                    _forbidden_extrapolator_2d<INDEPTYPE, DEPTYPE>() {}

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2(
                    _forbidden_extrapolator_2d, INDEPTYPE, DEPTYPE )

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
        class Extrapolator2D {
            public:
                Extrapolator2D() {}

                Extrapolator2D(
                    extrapolator_engine_1d_t<INDEPTYPE, DEPTYPE> engine ) {
                    set_engine( engine );
                }

                Extrapolator2D(
                    const Extrapolator2D<INDEPTYPE, DEPTYPE>& other ) :
                    Extrapolator2D<INDEPTYPE, DEPTYPE>() {
                    _engine = other._engine;
                }

                DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2( Extrapolator2D,
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
                            "Extrapolator2D: no engine has been set" );
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
    NCPA::interpolation::_abstract_extrapolator_2d<T, U>& a,
    NCPA::interpolation::_abstract_extrapolator_2d<T, U>& b ) noexcept {}

// template<typename T, typename U>
// static void swap( NCPA::interpolation::_extrapolator_1d<T, U>& a,
//                   NCPA::interpolation::_extrapolator_1d<T, U>& b ) noexcept
//                   {
//     ::swap(
//         dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
//             a ),
//         dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
//             b ) );
// }

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_constant_extrapolator_2d<T, U>& a,
    NCPA::interpolation::_constant_extrapolator_2d<T, U>& b ) noexcept {
    ::swap(
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            a ),
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            b ) );
}

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_linear_extrapolator_2d<T, U>& a,
    NCPA::interpolation::_linear_extrapolator_2d<T, U>& b ) noexcept {
    ::swap(
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            a ),
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            b ) );
}

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_zero_extrapolator_2d<T, U>& a,
    NCPA::interpolation::_zero_extrapolator_2d<T, U>& b ) noexcept {
    ::swap(
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            a ),
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            b ) );
}

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_forbidden_extrapolator_2d<T, U>& a,
    NCPA::interpolation::_forbidden_extrapolator_2d<T, U>& b ) noexcept {
    ::swap(
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            a ),
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            b ) );
}

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_periodic_extrapolator_1d<T, U>& a,
    NCPA::interpolation::_periodic_extrapolator_1d<T, U>& b ) noexcept {
    ::swap(
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            a ),
        dynamic_cast<NCPA::interpolation::_abstract_extrapolator_2d<T, U>&>(
            b ) );
}

template<typename T, typename U>
static void swap( NCPA::interpolation::Extrapolator2D<T, U>& a,
                  NCPA::interpolation::Extrapolator2D<T, U>& b ) noexcept {
    using std::swap;
    swap( a._engine, b._engine );
}
