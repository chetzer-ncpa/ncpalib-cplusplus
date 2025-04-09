#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/abstract_spline_2d.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

template<typename T, typename U>
static void swap(
    NCPA::interpolation::stratified_spline_2d<T, U>& a,
    NCPA::interpolation::stratified_spline_2d<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        // paramtype is a std::pair<size_t,interpolator_1d_type_t>
        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAM  //
            class stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE,
                                       ENABLE_IF_REAL( INDEPTYPE ),
                                       ENABLE_IF_REAL( DEPTYPE )>
            : public NCPA::interpolation::_spline_2d<INDEPTYPE, DEPTYPE> {
            public:
                stratified_spline_2d() {}

                stratified_spline_2d( PARAMTYPE params ) :
                    stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE>() {
                    _stratified_axis = params.first;
                    _interp1d = InterpolatorFactory<INDEPTYPE, DEPTYPE>::build(
                        params.second );
                }

                stratified_spline_2d(
                    const stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE>&
                        other ) :
                    stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE>() {
                    _x        = other._x;
                    _size_x   = other._size_x;
                    _f        = other._f;
                    _interp1d = other._interp1d;
                }

                stratified_spline_2d(
                    stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE>&&
                        source ) noexcept :
                    stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE>() {
                    ::swap( *this, source );
                }

                friend void ::swap<INDEPTYPE, DEPTYPE>(
                    stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE>& a,
                    stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE>&
                        b ) noexcept;

                stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE>& operator=(
                    stratified_spline_2d<INDEPTYPE, DEPTYPE, PARAMTYPE>
                        other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~stratified_spline_2d() {}

                virtual void fill( size_t N1, size_t N2, const INDEPTYPE *x1,
                                   const INDEPTYPE *x2,
                                   const DEPTYPE **f ) override {
                    if ( _stratified_axis == 0 ) {
                        if ( N2 != _x.size() ) {
                            this->init( N1, N2 );
                        }
                        _x.assign( x2, x2 + N2 );
                        _f.assign( f[ 0 ], f[ 0 ] + N2 );
                        _size_x = N2;
                    } else {
                        if ( N1 != _x.size() ) {
                            this->init( N1, N2 );
                        }
                        _x.assign( x1, x1 + N1 );
                        _f.resize( N1, NCPA::math::zero<DEPTYPE>() );
                        for ( size_t i = 0; i < N1; ++i ) {
                            _f[ i ] = f[ i ][ 0 ];
                        }
                        _size_x = N1;
                    }
                }

                virtual void fill(
                    const std::vector<INDEPTYPE>& x1,
                    const std::vector<INDEPTYPE>& x2,
                    const NCPA::arrays::vector2d_t<DEPTYPE>& f ) override {
                    if ( _stratified_axis == 0 ) {
                        if ( x2.size() != f.dim( 1 ) ) {
                            throw std::invalid_argument(
                                "Vector sizes must be equal!" );
                        }
                        _x      = x2;
                        _f      = f[ 0 ];
                        _size_x = x2.size();
                    } else {
                        if ( x1.size() != f.dim( 0 ) ) {
                            throw std::invalid_argument(
                                "Vector sizes must be equal!" );
                        }
                        _x      = x1;
                        _size_x = x1.size();
                        _f.resize( _size_x, NCPA::math::zero<DEPTYPE>() );
                        for ( size_t i = 0; i < _size_x; ++i ) {
                            _f[ i ] = f[ i ][ 0 ];
                        }
                    }
                }

                virtual void init( size_t N1, size_t N2 ) override {
                    clear();
                    if ( _stratified_axis == 0 ) {
                        _size_x = N2;
                    } else {
                        _size_x = N1;
                    }
                }

                virtual void clear() override {
                    _x.clear();
                    _f.clear();
                    _size_x = 0;
                }

                virtual void ready() override {
                    if ( _size_x == 0 || _x.empty() || _f.empty() ) {
                        throw std::logic_error(
                            "Interpolator has not been set up!" );
                    }
                    _interp1d.init( _size_x ).fill( _x, _f ).ready();
                }

                virtual DEPTYPE eval_f( INDEPTYPE x1, INDEPTYPE x2 ) override {
                    INDEPTYPE x = ( _stratified_axis == 0 ? x2 : x1 );
                    return _interp1d.eval_f( x );
                }

                virtual DEPTYPE eval_df( INDEPTYPE x1, INDEPTYPE x2,
                                         size_t wrt ) override {
                    if ( wrt != _stratified_axis ) {
                        return NCPA::math::zero<DEPTYPE>();
                    }
                    INDEPTYPE x = ( _stratified_axis == 0 ? x2 : x1 );
                    return _interp1d.eval_df( x );
                }

                virtual DEPTYPE eval_ddf( INDEPTYPE x1, INDEPTYPE x2,
                                          size_t wrt1, size_t wrt2 ) override {
                    if ( wrt1 != _stratified_axis
                         || wrt2 != _stratified_axis ) {
                        return NCPA::math::zero<DEPTYPE>();
                    }
                    INDEPTYPE x = ( _stratified_axis == 0 ? x2 : x1 );
                    return _interp1d.eval_ddf( x );
                }

                virtual DEPTYPE eval_dddf( INDEPTYPE x1, INDEPTYPE x2,
                                           size_t wrt1, size_t wrt2,
                                           size_t wrt3 ) override {
                    if ( wrt1 != _stratified_axis || wrt2 != _stratified_axis
                         || wrt3 != _stratified_axis ) {
                        return NCPA::math::zero<DEPTYPE>();
                    }
                    INDEPTYPE x = ( _stratified_axis == 0 ? x2 : x1 );
                    return _interp1d.eval_dddf( x );
                }

                virtual std::array<INDEPTYPE, 4> limits() const override {
                    return std::array<INDEPTYPE, 4> { _x.front(), _x.back(),
                                                      _x.front(), _x.back() };
                }

            private:
                size_t _stratified_axis = 0;
                size_t _size_x          = 0;
                std::vector<INDEPTYPE> _x;
                std::vector<DEPTYPE> _f;
                Interpolator1D<INDEPTYPE, DEPTYPE> _interp1d;
        };

        DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAM(
            stratified_spline_2d, _spline_2d, stratified_axis_type_t );
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap(
    NCPA::interpolation::stratified_spline_2d<T, U>& a,
    NCPA::interpolation::stratified_spline_2d<T, U>& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::interpolation::_spline_2d<T, U>&>( a ),
            dynamic_cast<NCPA::interpolation::_spline_2d<T, U>&>( b ) );
    swap( a._x, b._x );
    swap( a._f, b._f );
    swap( a._size_x, b._size_x );
    swap( a._stratified_axis, b._stratified_axis );
    swap( a._interp1d, b._interp1d );
}
