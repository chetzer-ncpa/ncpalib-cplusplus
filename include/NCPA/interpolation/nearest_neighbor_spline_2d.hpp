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
#include <vector>

template<typename T, typename U>
static void swap(
    NCPA::interpolation::nearest_neighbor_spline_2d<T, U>& a,
    NCPA::interpolation::nearest_neighbor_spline_2d<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
            class nearest_neighbor_spline_2d<INDEPTYPE, DEPTYPE, void,
                                             ENABLE_IF_REAL( INDEPTYPE ),
                                             ENABLE_IF_REAL( DEPTYPE )>
            : public NCPA::interpolation::_spline_2d<INDEPTYPE, DEPTYPE> {
            public:
                nearest_neighbor_spline_2d() :
                    _size_x1 { 0 }, _size_x2 { 0 } {}

                nearest_neighbor_spline_2d(
                    const nearest_neighbor_spline_2d<INDEPTYPE, DEPTYPE>&
                        other ) :
                    nearest_neighbor_spline_2d<INDEPTYPE, DEPTYPE>() {
                    _x1      = other._x1;
                    _x2      = other._x2;
                    _size_x1 = other._size_x1;
                    _size_x2 = other._size_x2;
                    _f       = other._f;
                }

                nearest_neighbor_spline_2d(
                    nearest_neighbor_spline_2d<INDEPTYPE, DEPTYPE>&&
                        source ) noexcept :
                    nearest_neighbor_spline_2d<INDEPTYPE, DEPTYPE>() {
                    ::swap( *this, source );
                }

                friend void ::swap<INDEPTYPE, DEPTYPE>(
                    nearest_neighbor_spline_2d<INDEPTYPE, DEPTYPE>& a,
                    nearest_neighbor_spline_2d<INDEPTYPE, DEPTYPE>&
                        b ) noexcept;

                nearest_neighbor_spline_2d<INDEPTYPE, DEPTYPE>& operator=(
                    nearest_neighbor_spline_2d<INDEPTYPE, DEPTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~nearest_neighbor_spline_2d() {}

                virtual void fill( size_t N1, size_t N2, const INDEPTYPE *x1,
                                   const INDEPTYPE *x2,
                                   const DEPTYPE **f ) override {
                    if ( N1 != _x1.size() || N2 != _x2.size() ) {
                        this->init( N1, N2 );
                    }
                    _x1.assign( x1, x1 + N1 );
                    _x2.assign( x2, x2 + N2 );
                    for ( size_t i = 0; i < N1; ++i ) {
                        _f[ i ].assign( f[ i ], f[ i ] + N2 );
                    }
                    _size_x1 = N1;
                    _size_x2 = N2;
                }

                virtual void fill(
                    const std::vector<INDEPTYPE>& x1,
                    const std::vector<INDEPTYPE>& x2,
                    const NCPA::arrays::vector2d_t<DEPTYPE>& f ) override {
                    if ( x1.size() != f.dim( 0 )
                         || x2.size() != f.dim( 1 ) ) {
                        throw std::invalid_argument(
                            "Vector sizes must be equal!" );
                    }
                    _x1      = x1;
                    _x2      = x2;
                    _f       = f;
                    _size_x1 = x1.size();
                    _size_x2 = x2.size();
                }

                virtual void init( size_t N1, size_t N2 ) override {
                    clear();
                    _size_x1 = N1;
                    _size_x2 = N2;
                }

                virtual void clear() override {
                    _x1.clear();
                    _x2.clear();
                    _f.clear();
                    _size_x1 = 0;
                    _size_x2 = 0;
                }

                virtual void ready() override {
                    if ( _size_x1 == 0 || _size_x2 == 0 || _x1.empty()
                         || _x2.empty() || _f.empty() ) {
                        throw std::logic_error(
                            "Interpolator has not been set up!" );
                    }
                }

                virtual DEPTYPE eval_f( INDEPTYPE x1, INDEPTYPE x2 ) override {
                    std::pair<size_t, size_t> coord
                        = NCPA::math::find_closest_point( _x1, _x2, x1, x2 );
                    return _f[ coord.first ][ coord.second ];
                }

                virtual DEPTYPE eval_df( INDEPTYPE x1, INDEPTYPE x2,
                                         size_t wrt ) override {
                    return NCPA::math::zero<DEPTYPE>();
                }

                virtual DEPTYPE eval_ddf( INDEPTYPE x1, INDEPTYPE x2,
                                          size_t wrt1, size_t wrt2 ) override {
                    return NCPA::math::zero<DEPTYPE>();
                }

                virtual DEPTYPE eval_dddf( INDEPTYPE x1, INDEPTYPE x2,
                                           size_t wrt1, size_t wrt2,
                                           size_t wrt3 ) override {
                    return NCPA::math::zero<DEPTYPE>();
                }

                virtual std::array<INDEPTYPE,4> limits() const override {
                    return std::array<INDEPTYPE,4>{
                        _x1.front(), _x1.back(), _x2.front(), _x2.back()
                    };
                }

            private:
                size_t _size_x1 = 0, _size_x2 = 0;

                std::vector<INDEPTYPE> _x1, _x2;
                NCPA::arrays::vector2d_t<DEPTYPE> _f;
        };
        DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR( nearest_neighbor_spline_2d,
                                                _abstract_spline_2d )
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap(
    NCPA::interpolation::nearest_neighbor_spline_2d<T, U>& a,
    NCPA::interpolation::nearest_neighbor_spline_2d<T, U>& b ) noexcept {
    using std::swap;
    ::swap(
        dynamic_cast<NCPA::interpolation::_abstract_spline_2d<T, U>&>( a ),
        dynamic_cast<NCPA::interpolation::_abstract_spline_2d<T, U>&>( b ) );
    swap( a._x1, b._x1 );
    swap( a._x2, b._x2 );
    swap( a._f, b._f );
    swap( a._size_x1, b._size_x1 );
    swap( a._size_x2, b._size_x2 );
}
