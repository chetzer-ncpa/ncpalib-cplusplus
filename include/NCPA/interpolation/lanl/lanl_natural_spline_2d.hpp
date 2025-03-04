/*********************************************************************
Elements of this code are modified from LANL InfraGA under the following
license: Copyright (c) 2014, Triad National Security, LLC

All rights reserved.

Copyright 2014. Triad National Security, LLC. This software was produced under
U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory
(LANL), which is operated by Los Alamos National Security, LLC for the U.S.
Department of Energy. The U.S. Government has rights to use, reproduce, and
distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL
SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative
works, such modified software should be clearly marked, so as not to confuse it
with the version available from LANL.

Additionally, permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including without
limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Original InfraGA software obtained from
https://github.com/LANL-Seismoacoustics/infraGA

Modified for use in libncpa by Claus Hetzer, NCPA University of Mississippi,
claus@olemiss.edu
*/


#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/abstract_spline_2d.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/lanl/lanl_common.hpp"
#include "NCPA/interpolation/lanl/lanl_declarations.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"
#include <array>

template<typename T, typename U>
static void swap(
    NCPA::interpolation::LANL::natural_spline_2d<T, U>& a,
    NCPA::interpolation::LANL::natural_spline_2d<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        namespace LANL {

            /**
             * 2D INTERPOLATORS
             */
            // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
            //     natural_spline_2d, NCPA::interpolation::_spline_2d );

            template<typename INDEPTYPE, typename DEPTYPE>
            class natural_spline_2d<INDEPTYPE, DEPTYPE, void,
                                    ENABLE_IF_REAL( INDEPTYPE ),
                                    ENABLE_IF_REAL( DEPTYPE )>
                : public NCPA::interpolation::_spline_2d<INDEPTYPE, DEPTYPE> {
                public:
                    natural_spline_2d() :
                        NCPA::interpolation::_spline_2d<INDEPTYPE, DEPTYPE>() {
                    }

                    virtual ~natural_spline_2d() { this->clear(); }

                    natural_spline_2d(
                        const natural_spline_2d<INDEPTYPE, DEPTYPE>& other ) :
                        NCPA::interpolation::_spline_2d<INDEPTYPE, DEPTYPE>(
                            other ) {
                        _length_x   = other._length_x;
                        _length_y   = other._length_y;
                        _accel[ 0 ] = other._accel[ 0 ];
                        _accel[ 1 ] = other._accel[ 1 ];
                        NCPA::arrays::copy( other._x_vals, _length_x, _x_vals );
                        NCPA::arrays::copy( other._y_vals, _length_y, _y_vals );
                        NCPA::arrays::copy( other._f_vals, _length_x, _length_y, _f_vals );
                        NCPA::arrays::copy( other._f_slopes, _length_x, _length_y, _f_slopes );
                        NCPA::arrays::copy( other._new_c, _length_y-1, _new_c );
                        NCPA::arrays::copy( other._new_d, _length_y, _new_d );
                        NCPA::arrays::copy( other._node_vals, _length_x, _node_vals );
                        NCPA::arrays::copy( other._node_slopes, _length_x, _node_slopes );
                    }

                    natural_spline_2d( linear_spline_1d<INDEPTYPE, DEPTYPE>&&
                                           source ) noexcept :
                        natural_spline_2d<INDEPTYPE, DEPTYPE>() {
                        ::swap( *this, source );
                    }

                    friend void ::swap<INDEPTYPE, DEPTYPE>(
                        natural_spline_2d<INDEPTYPE, DEPTYPE>& a,
                        natural_spline_2d<INDEPTYPE, DEPTYPE>& b ) noexcept;

                    natural_spline_2d<INDEPTYPE, DEPTYPE>& operator=(
                        natural_spline_2d<INDEPTYPE, DEPTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    void init( size_t nx, size_t ny ) override {
                        this->clear();
                        _length_x   = nx;
                        _accel[ 0 ] = 0;
                        _length_y   = ny;
                        _accel[ 1 ] = 0;
                        _x_vals = NCPA::arrays::zeros<INDEPTYPE>( _length_x );
                        _y_vals = NCPA::arrays::zeros<INDEPTYPE>( _length_y );
                        _f_vals = NCPA::arrays::zeros<DEPTYPE>( _length_x,
                                                                _length_y );
                        _f_slopes = NCPA::arrays::zeros<DEPTYPE>( _length_x,
                                                                  _length_y );
                        _new_c = NCPA::arrays::zeros<DEPTYPE>( _length_y - 1 );
                        _new_d = NCPA::arrays::zeros<DEPTYPE>( _length_y );
                        _node_vals = NCPA::arrays::zeros<DEPTYPE>( _length_x );
                        _node_slopes
                            = NCPA::arrays::zeros<DEPTYPE>( _length_x );
                    }

                    virtual void fill( size_t N1, size_t N2,
                                       const INDEPTYPE *x1,
                                       const INDEPTYPE *x2,
                                       const DEPTYPE **f ) override {
                        if ( !_initialized() || _length_x != N1
                             || _length_y != N2 ) {
                            init( N1, N2 );
                        }
                        std::memcpy( _x_vals, x1, N1 * sizeof( INDEPTYPE ) );
                        std::memcpy( _y_vals, x2, N2 * sizeof( INDEPTYPE ) );
                        for ( size_t i = 0; i < N1; ++i ) {
                            std::memcpy( _f_vals[ i ], f[ i ],
                                         N2 * sizeof( DEPTYPE ) );
                        }
                    }

                    virtual void fill(
                        const std::vector<INDEPTYPE>& x1,
                        const std::vector<INDEPTYPE>& x2,
                        const NCPA::arrays::vector2d_t<DEPTYPE>& f ) override {
                        size_t N1 = x1.size(), N2 = x2.size();
                        if ( N1 != f.dim( 0 ) ) {
                            throw std::invalid_argument(
                                "Vectors must be same size" );
                        }
                        if ( !_initialized() || _length_x != N1
                             || _length_y != N2 ) {
                            init( N1, N2 );
                        }
                        std::memcpy( _x_vals, &x1[ 0 ],
                                     N1 * sizeof( INDEPTYPE ) );
                        std::memcpy( _y_vals, &x2[ 0 ],
                                     N2 * sizeof( INDEPTYPE ) );
                        for ( size_t i = 0; i < N1; ++i ) {
                            std::memcpy( _f_vals[ i ], &f[ i ][ 0 ],
                                         N2 * sizeof( DEPTYPE ) );
                        }
                    }

                    virtual void clear() override {
                        if ( _x_vals != nullptr ) {
                            delete[] _x_vals;
                            _x_vals = nullptr;
                        }
                        if ( _y_vals != nullptr ) {
                            delete[] _y_vals;
                            _y_vals = nullptr;
                        }
                        if ( _f_vals != nullptr ) {
                            NCPA::arrays::free_array( _f_vals, _length_x,
                                                      _length_y );
                            _f_vals = nullptr;
                        }
                        if ( _f_slopes != nullptr ) {
                            NCPA::arrays::free_array( _f_slopes, _length_x,
                                                      _length_y );
                            _f_slopes = nullptr;
                        }
                        if ( _new_c != nullptr ) {
                            delete[] _new_c;
                            _new_c = nullptr;
                        }
                        if ( _new_d != nullptr ) {
                            delete[] _new_d;
                            _new_d = nullptr;
                        }
                        if ( _node_vals != nullptr ) {
                            delete[] _node_vals;
                            _node_vals = nullptr;
                        }
                        if ( _node_slopes != nullptr ) {
                            delete[] _node_slopes;
                            _node_slopes = nullptr;
                        }
                        _length_x   = 0;
                        _length_y   = 0;
                        _accel[ 0 ] = 0;
                        _accel[ 1 ] = 0;
                    }

                    void ready() override {
                        INDEPTYPE dy1, dy2, ai, bi, ci;
                        DEPTYPE di;

                        for ( int nx = 0; nx < _length_x; ++nx ) {
                            bi = 2.0 / ( _y_vals[ 1 ] - _y_vals[ 0 ] );
                            ci = 1.0 / ( _y_vals[ 1 ] - _y_vals[ 0 ] );
                            di = 3.0
                               * ( _f_vals[ nx ][ 1 ] - _f_vals[ nx ][ 0 ] )
                               / std::pow( _y_vals[ 1 ] - _y_vals[ 0 ], 2 );

                            _new_c[ 0 ] = ci / bi;
                            _new_d[ 0 ] = di / bi;

                            for ( int ny = 1; ny < _length_y - 1; ++ny ) {
                                dy1 = _y_vals[ ny ] - _y_vals[ ny - 1 ],
                                dy2 = _y_vals[ ny + 1 ] - _y_vals[ ny ];

                                ai = 1.0 / dy1;
                                bi = 2.0 * ( 1.0 / dy1 + 1.0 / dy2 );
                                ci = 1.0 / dy2;
                                di = 3.0
                                   * ( ( _f_vals[ nx ][ ny ]
                                         - _f_vals[ nx ][ ny - 1 ] )
                                           / std::pow( dy1, 2 )
                                       + ( _f_vals[ nx ][ ny + 1 ]
                                           - _f_vals[ nx ][ ny ] )
                                             / std::pow( dy2, 2 ) );

                                _new_c[ ny ]
                                    = ci / ( bi - _new_c[ ny - 1 ] * ai );
                                _new_d[ ny ] = ( di - _new_d[ ny - 1 ] * ai )
                                             / ( bi - _new_c[ ny - 1 ] * ai );
                            }

                            ai = 1.0
                               / ( _y_vals[ _length_y - 1 ]
                                   - _y_vals[ _length_y - 2 ] );
                            bi = 2.0
                               / ( _y_vals[ _length_y - 1 ]
                                   - _y_vals[ _length_y - 2 ] );
                            di = 3.0
                               * ( _f_vals[ nx ][ _length_y - 1 ]
                                   - _f_vals[ nx ][ _length_y - 2 ] )
                               / std::pow( _y_vals[ _length_y - 1 ]
                                               - _y_vals[ _length_y - 2 ],
                                           2 );

                            _new_d[ _length_y - 1 ]
                                = ( di - _new_d[ _length_y - 2 ] * ai )
                                / ( bi - _new_c[ _length_y - 2 ] * ai );

                            _f_slopes[ nx ][ _length_y - 1 ]
                                = _new_d[ _length_y - 1 ];
                            for ( int ny = _length_y - 2; ny >= 0; ny-- ) {
                                _f_slopes[ nx ][ ny ]
                                    = _new_d[ ny ]
                                    - _new_c[ ny ] * _f_slopes[ nx ][ ny + 1 ];
                            }
                        }
                    }

                    DEPTYPE eval_f( INDEPTYPE x, INDEPTYPE y ) override {
                        int kx, ky;
                        INDEPTYPE dx, X, A, B;
                        DEPTYPE df;

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );

                        for ( size_t nx = 0; nx < _length_x; nx++ ) {
                            _node_vals[ nx ] = _eval_node_f( y, nx, ky );
                        }
                        _set_node_slopes();

                        df = _node_vals[ kx + 1 ] - _node_vals[ kx ];
                        dx = _x_vals[ kx + 1 ] - _x_vals[ kx ];

                        X = ( x - _x_vals[ kx ] ) / dx;
                        A = _node_slopes[ kx ] * dx - df;
                        B = -_node_slopes[ kx + 1 ] * dx + df;

                        return _node_vals[ kx ]
                             + X
                                   * ( df
                                       + ( 1.0 - X )
                                             * ( A * ( 1.0 - X ) + B * X ) );
                    }

                    DEPTYPE eval_df( INDEPTYPE x, INDEPTYPE y,
                                     size_t n ) override {
                        size_t kx, ky;
                        DEPTYPE df, A, B;
                        INDEPTYPE dx, X;

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );

                        for ( size_t nx = 0; nx < _length_x; nx++ ) {
                            if ( n == 0 ) {
                                _node_vals[ nx ] = _eval_node_f( y, nx, ky );
                            } else {
                                _node_vals[ nx ]
                                    = _eval_node_dfdy( y, nx, ky );
                            }
                        }
                        _set_node_slopes();

                        df = _node_vals[ kx + 1 ] - _node_vals[ kx ];
                        dx = _x_vals[ kx + 1 ] - _x_vals[ kx ];

                        X = ( x - _x_vals[ kx ] ) / dx;
                        A = _node_slopes[ kx ] * dx - df;
                        B = -_node_slopes[ kx + 1 ] * dx + df;

                        if ( n == 0 ) {
                            return df / dx
                                 + ( 1.0 - 2.0 * X )
                                       * ( A * ( 1.0 - X ) + B * X ) / dx
                                 + X * ( 1.0 - X ) * ( B - A ) / dx;
                        } else {
                            return _node_vals[ kx ]
                                 + X
                                       * ( df
                                           + ( 1.0 - X )
                                                 * ( A * ( 1.0 - X )
                                                     + B * X ) );
                        }
                    }

                    DEPTYPE eval_ddf( INDEPTYPE x, INDEPTYPE y, size_t n1,
                                      size_t n2 ) override {
                        size_t kx, ky;
                        DEPTYPE df, A, B;
                        INDEPTYPE dx, X;

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );

                        for ( size_t nx = 0; nx < _length_x; nx++ ) {
                            if ( n1 + n2 == 0 ) {
                                _node_vals[ nx ] = _eval_node_f( y, nx, ky );
                            } else if ( n1 + n2 == 1 ) {
                                _node_vals[ nx ]
                                    = _eval_node_dfdy( y, nx, ky );
                            } else {
                                _node_vals[ nx ]
                                    = _eval_node_ddfdydy( y, nx, ky );
                            }
                        }
                        _set_node_slopes();

                        df = _node_vals[ kx + 1 ] - _node_vals[ kx ];
                        dx = _x_vals[ kx + 1 ] - _x_vals[ kx ];

                        X = ( x - _x_vals[ kx ] ) / dx;
                        A = _node_slopes[ kx ] * dx - df;
                        B = -_node_slopes[ kx + 1 ] * dx + df;

                        if ( n1 + n2 == 0 ) {
                            return 2.0 * ( B - 2.0 * A + ( A - B ) * 3.0 * X )
                                 / std::pow( dx, 2 );
                        } else if ( n1 + n2 == 1 ) {
                            return df / dx
                                 + ( 1.0 - 2.0 * X )
                                       * ( A * ( 1.0 - X ) + B * X ) / dx
                                 + X * ( 1.0 - X ) * ( B - A ) / dx;
                        } else {
                            return _node_vals[ kx ]
                                 + X
                                       * ( df
                                           + ( 1.0 - X )
                                                 * ( A * ( 1.0 - X )
                                                     + B * X ) );
                        }
                    }

                    DEPTYPE eval_dddf( INDEPTYPE x, INDEPTYPE y, size_t n1,
                                       size_t n2, size_t n3 ) override {
                        size_t kx, ky;
                        INDEPTYPE X, dx;
                        DEPTYPE df, A, B;

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );

                        size_t dimsum = n1 + n2 + n3;
                        for ( int nx = 0; nx < _length_x; nx++ ) {
                            if ( dimsum == 0 ) {
                                _node_vals[ nx ] = _eval_node_f( y, nx, ky );
                            } else if ( dimsum == 1 ) {
                                _node_vals[ nx ]
                                    = _eval_node_dfdy( y, nx, ky );
                            } else if ( dimsum == 2 ) {
                                _node_vals[ nx ]
                                    = _eval_node_ddfdydy( y, nx, ky );
                            } else {
                                _node_vals[ nx ]
                                    = _eval_node_dddfdydydy( y, nx, ky );
                            }
                        }
                        _set_node_slopes();
                        df = _node_vals[ kx + 1 ] - _node_vals[ kx ];
                        dx = _x_vals[ kx + 1 ] - _x_vals[ kx ];

                        X = ( x - _x_vals[ kx ] ) / dx;
                        A = _node_slopes[ kx ] * dx - df;
                        B = -_node_slopes[ kx + 1 ] * dx + df;

                        if ( dimsum == 0 ) {
                            return 6.0 * ( A - B ) / std::pow( dx, 3 );
                        } else if ( dimsum == 1 ) {
                            return 2.0 * ( B - 2.0 * A + ( A - B ) * 3.0 * X )
                                 / std::pow( dx, 2 );
                        } else if ( dimsum == 2 ) {
                            return df / dx
                                 + ( 1.0 - 2.0 * X )
                                       * ( A * ( 1.0 - X ) + B * X ) / dx
                                 + X * ( 1.0 - X ) * ( B - A ) / dx;
                        } else {
                            return _node_vals[ kx ]
                                 + X
                                       * ( df
                                           + ( 1.0 - X )
                                                 * ( A * ( 1.0 - X )
                                                     + B * X ) );
                        }
                    }

                    virtual std::array<INDEPTYPE,4> limits() const override {
                        return std::array<INDEPTYPE,4>{
                            _x_vals[0], _x_vals[_length_x-1], _y_vals[0], _y_vals[_length_y-1]
                        };
                    }

                protected:
                    bool _initialized() const {
                        return ( _length_x > 0 && _length_y > 0
                                 && _x_vals != nullptr && _y_vals != nullptr
                                 && _f_vals != nullptr
                                 && _f_slopes != nullptr );
                    }

                    DEPTYPE _eval_node_f( INDEPTYPE y, size_t nx, size_t ny ) {
                        DEPTYPE df, A, B;
                        INDEPTYPE dy, X;

                        df = _f_vals[ nx ][ ny + 1 ] - _f_vals[ nx ][ ny ];
                        dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];

                        X = ( y - _y_vals[ ny ] ) / dy;
                        A = _f_slopes[ nx ][ ny ] * dy - df;
                        B = -_f_slopes[ nx ][ ny + 1 ] * dy + df;

                        return _f_vals[ nx ][ ny ]
                             + X
                                   * ( df
                                       + ( 1.0 - X )
                                             * ( A * ( 1.0 - X ) + B * X ) );
                    }

                    DEPTYPE _eval_node_dfdy( INDEPTYPE y, size_t nx,
                                             size_t ny ) {
                        DEPTYPE df, A, B;
                        INDEPTYPE dy, X;

                        df = _f_vals[ nx ][ ny + 1 ] - _f_vals[ nx ][ ny ];
                        dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];

                        X = ( y - _y_vals[ ny ] ) / dy;
                        A = _f_slopes[ nx ][ ny ] * dy - df;
                        B = -_f_slopes[ nx ][ ny + 1 ] * dy + df;

                        return df / dy
                             + ( 1.0 - 2.0 * X ) * ( A * ( 1.0 - X ) + B * X )
                                   / dy
                             + X * ( 1.0 - X ) * ( B - A ) / dy;
                    }

                    DEPTYPE _eval_node_ddfdydy( INDEPTYPE y, size_t nx,
                                                size_t ny ) {
                        DEPTYPE df, A, B;
                        INDEPTYPE dy, X;

                        df = _f_vals[ nx ][ ny + 1 ] - _f_vals[ nx ][ ny ];
                        dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];

                        X = ( y - _y_vals[ ny ] ) / dy;
                        A = _f_slopes[ nx ][ ny ] * dy - df;
                        B = -_f_slopes[ nx ][ ny + 1 ] * dy + df;

                        return 2.0 * ( B - 2.0 * A + ( A - B ) * 3.0 * X )
                             / std::pow( dy, 2 );
                    }

                    DEPTYPE _eval_node_dddfdydydy( INDEPTYPE y, size_t nx,
                                                   size_t ny ) {
                        DEPTYPE df, A, B;
                        INDEPTYPE dy, X;

                        df = _f_vals[ nx ][ ny + 1 ] - _f_vals[ nx ][ ny ];
                        dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];

                        X = ( y - _y_vals[ ny ] ) / dy;
                        A = _f_slopes[ nx ][ ny ] * dy - df;
                        B = -_f_slopes[ nx ][ ny + 1 ] * dy + df;

                        return 6.0 * ( A - B ) / std::pow( dy, 3 );
                    }

                    void _set_node_slopes() {
                        INDEPTYPE dx1, dx2, ai, bi, ci;
                        DEPTYPE di;

                        bi = 2.0 / ( _x_vals[ 1 ] - _x_vals[ 0 ] );
                        ci = 1.0 / ( _x_vals[ 1 ] - _x_vals[ 0 ] );
                        di = 3.0 * ( _node_vals[ 1 ] - _node_vals[ 0 ] )
                           / std::pow( _x_vals[ 1 ] - _x_vals[ 0 ], 2 );

                        _new_c[ 0 ] = ci / bi;
                        _new_d[ 0 ] = di / bi;

                        for ( int nx = 1; nx < _length_x - 1; nx++ ) {
                            dx1 = _x_vals[ nx ] - _x_vals[ nx - 1 ];
                            dx2 = _x_vals[ nx + 1 ] - _x_vals[ nx ];

                            ai = 1.0 / dx1;
                            bi = 2.0 * ( 1.0 / dx1 + 1.0 / dx2 );
                            ci = 1.0 / dx2;
                            di = 3.0
                               * ( ( _node_vals[ nx ] - _node_vals[ nx - 1 ] )
                                       / std::pow( dx1, 2 )
                                   + ( _node_vals[ nx + 1 ]
                                       - _node_vals[ nx ] )
                                         / std::pow( dx2, 2 ) );

                            _new_c[ nx ] = ci / ( bi - _new_c[ nx - 1 ] * ai );
                            _new_d[ nx ] = ( di - _new_d[ nx - 1 ] * ai )
                                         / ( bi - _new_c[ nx - 1 ] * ai );
                        }

                        ai = 1.0
                           / ( _x_vals[ _length_x - 1 ]
                               - _x_vals[ _length_x - 2 ] );
                        bi = 2.0
                           / ( _x_vals[ _length_x - 1 ]
                               - _x_vals[ _length_x - 2 ] );
                        di = 3.0
                           * ( _node_vals[ _length_x - 1 ]
                               - _node_vals[ _length_x - 2 ] )
                           / std::pow( _x_vals[ _length_x - 1 ]
                                           - _x_vals[ _length_x - 2 ],
                                       2 );

                        _new_d[ _length_x - 1 ]
                            = ( di - _new_d[ _length_x - 2 ] * ai )
                            / ( bi - _new_c[ _length_x - 2 ] * ai );

                        _node_slopes[ _length_x - 1 ]
                            = _new_d[ _length_x - 1 ];
                        for ( int nx = _length_x - 2; nx >= 0; nx-- ) {
                            _node_slopes[ nx ]
                                = _new_d[ nx ]
                                - _new_c[ nx ] * _node_slopes[ nx + 1 ];
                        }
                    }

                    size_t _length_x;              // Number of x nodes
                    size_t _length_y;              // Number of y nodes
                    size_t _accel[ 2 ];            // Indices of x,y point

                    INDEPTYPE *_x_vals = nullptr;  // 1D array of x values
                    INDEPTYPE *_y_vals = nullptr;  // 1D array of y values

                    DEPTYPE **_f_vals
                        = nullptr;  // 2D array of f(x,y) values,
                                    // f_vals[i][j] = f(x[i], y[j])
                    DEPTYPE **_f_slopes
                        = nullptr;  // Slopes used to generate natural cubic
                                    // spline solution at each x = constant
                                    // node, describes f(x[i], y)

                    // cache variables
                    DEPTYPE *_new_c       = nullptr;    // size() = _length_y - 1
                    DEPTYPE *_new_d       = nullptr;    // size() = _length_y
                    DEPTYPE *_node_vals   = nullptr;    // size() = _length_x
                    DEPTYPE *_node_slopes = nullptr;    // size() = _length_x

                    
            };

            DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(
                natural_spline_2d, NCPA::interpolation::_spline_2d );

        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap(
    NCPA::interpolation::LANL::natural_spline_2d<T, U>& a,
    NCPA::interpolation::LANL::natural_spline_2d<T, U>& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::interpolation::_spline_2d<T, U>&>( a ),
            dynamic_cast<NCPA::interpolation::_spline_2d<T, U>&>( b ) );
    swap( a._length_x, b._length_x );
    swap( a._length_y, b._length_y );
    swap( a._accel, b._accel );
    swap( a._x_vals, b._x_vals );
    swap( a._y_vals, b._y_vals );
    swap( a._f_vals, b._f_vals );
    swap( a._f_slopes, b._f_slopes );
    swap( a._new_c, b._new_c );
    swap( a._new_d, b._new_d );
    swap( a._node_vals, b._node_vals );
    swap( a._node_slopes, b._node_slopes );
}
