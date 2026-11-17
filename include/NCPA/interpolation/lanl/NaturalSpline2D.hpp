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
#include "NCPA/cloneable.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/lanl/declarations.hpp"
#include "NCPA/interpolation/lanl/lanl_common.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <array>

namespace NCPA {
    namespace interpolation {
        namespace LANL {

            /**
             * 2D INTERPOLATORS
             */
            // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
            //     NaturalSpline2D, NCPA::interpolation::_spline_2d );

            // template<typename FloatType, typename FloatType>
            // class NaturalSpline2D<FloatType, void,
            //                         ENABLE_IF_REAL( FloatType ),
            //                         ENABLE_IF_REAL( FloatType )>
            //     : public NCPA::interpolation::_spline_2d<FloatType>
            //     {
            template<typename FloatType>
            class NaturalSpline2D
                : public Cloneable<NaturalSpline2D<FloatType>,
                                   spline2d<FloatType>> {
                    typedef Cloneable<NaturalSpline2D<FloatType>,
                                      spline2d<FloatType>>
                        parent_t;
                    typedef spline2d<FloatType> clone_t;

                public:
                    using clone_t::init;

                    NaturalSpline2D() = default;

                    virtual ~NaturalSpline2D() { this->clear(); }

                    NaturalSpline2D( const NaturalSpline2D<FloatType>& other )
                        = default;
                    // _length_x   = other._length_x;
                    // _length_y   = other._length_y;
                    //     _accel[ 0 ] = other._accel[ 0 ];
                    //     _accel[ 1 ] = other._accel[ 1 ];
                    //     _x_vals = other._x_vals;
                    //     _y_vals = other._y_vals;
                    //     _f_vals = other._f_vals;
                    //     _f_slopes = other._f_slopes;
                    //     _new_c = other._new_c;
                    //     NCPA::arrays::copy( other._x_vals, _length_x,
                    //                         _x_vals );
                    //     NCPA::arrays::copy( other._y_vals, _length_y,
                    //                         _y_vals );
                    //     NCPA::arrays::copy( other._f_vals, _length_x,
                    //                         _length_y, _f_vals );
                    //     NCPA::arrays::copy( other._f_slopes, _length_x,
                    //                         _length_y, _f_slopes );
                    //     NCPA::arrays::copy( other._new_c, _length_y - 1,
                    //                         _new_c );
                    //     NCPA::arrays::copy( other._new_d, _length_y, _new_d
                    //     ); NCPA::arrays::copy( other._node_vals, _length_x,
                    //                         _node_vals );
                    //     NCPA::arrays::copy( other._node_slopes, _length_x,
                    //                         _node_slopes );
                    // }

                    NaturalSpline2D(
                        NaturalSpline2D<FloatType>&& source ) noexcept
                        = default;

                    friend void swap(
                        NaturalSpline2D<FloatType>& a,
                        NaturalSpline2D<FloatType>& b ) noexcept {
                        using std::swap;
                        swap( dynamic_cast<parent_t&>( a ),
                              dynamic_cast<parent_t&>( b ) );
                        // swap( a._length_x, b._length_x );
                        // swap( a._length_y, b._length_y );
                        swap( a._accel, b._accel );
                        // swap( a._x_vals, b._x_vals );
                        // swap( a._y_vals, b._y_vals );
                        // swap( a._f_vals, b._f_vals );
                        swap( a._f_slopes, b._f_slopes );
                        swap( a._new_c, b._new_c );
                        swap( a._new_d, b._new_d );
                        swap( a._node_vals, b._node_vals );
                        swap( a._node_slopes, b._node_slopes );
                    }

                    NaturalSpline2D<FloatType>& operator=(
                        NaturalSpline2D<FloatType> other ) {
                        swap( *this, other );
                        return *this;
                    }

                    interpolator_type_t interptype() const override {
                        return NCPA::interpolation::interpolator_type_t::
                            LANL_NATURAL_2D;
                    }

                    void init( std::vector<size_t> dims ) override {
                        // void init( size_t nx, size_t ny ) override {
                        this->_check_dimensions( dims );
                        size_t nx = dims.at( 0 );
                        size_t ny = dims.at( 1 );
                        this->clear();
                        clone_t::init( {nx, ny });
                        if (nx > 0 && ny > 0) {
                            // this->axes( 0 ).resize( nx );
                            // this->axes( 1 ).resize( ny );
                            // this->f().resize( nx, ny );
                            _f_slopes.resize( nx, ny );
                            _new_c.resize( ny - 1 );
                            _new_d.resize( ny );
                            _node_vals.resize( nx );
                            _node_slopes.resize( nx );
                        }
                        // _length_x   = nx;
                        // _accel[ 0 ] = 0;
                        // _length_y   = ny;
                        // _accel[ 1 ] = 0;
                        // _x_vals = NCPA::arrays::zeros<FloatType>( _length_x
                        // ); _y_vals = NCPA::arrays::zeros<FloatType>(
                        // _length_y ); _f_vals =
                        // NCPA::arrays::zeros<FloatType>( _length_x,
                        //                                         _length_y );
                        // _f_slopes = NCPA::arrays::zeros<FloatType>(
                        // _length_x,
                        //                                           _length_y
                        //                                           );
                        // _new_c = NCPA::arrays::zeros<FloatType>( _length_y -
                        // 1
                        // ); _new_d = NCPA::arrays::zeros<FloatType>(
                        // _length_y
                        // ); _node_vals = NCPA::arrays::zeros<FloatType>(
                        // _length_x ); _node_slopes
                        //     = NCPA::arrays::zeros<FloatType>( _length_x );
                    }

                    // virtual void fill( size_t N1, size_t N2,
                    //                    const FloatType *x1,
                    //                    const FloatType *x2,
                    //                    const FloatType **f ) override {
                    //     if (!_initialized() || _x_vals.size() != N1
                    //         || _y_vals.size() != N2) {
                    //         init( N1, N2 );
                    //     }
                    //     std::memcpy( _x_vals.data(), x1,
                    //                  N1 * sizeof( FloatType ) );
                    //     std::memcpy( _y_vals.data(), x2,
                    //                  N2 * sizeof( FloatType ) );
                    //     for (size_t i = 0; i < N1; ++i) {
                    //         for (size_t j = 0; j < N2; ++j) {
                    //             _f_vals( i, j ) = f[ i ][ j ];
                    //         }
                    //     }
                    // }

                    // void fill( const std::vector<FloatType>& x1,
                    //            const std::vector<FloatType>& x2,
                    //            const NCPA::arrays::SmartVector<FloatType>& f
                    //            )
                    //     override {
                    //     size_t N1 = x1.size(), N2 = x2.size();
                    //     if (N1 != f.dim( 0 ) || N2 != f.dim( 1 )) {
                    //         throw InputSizesDoNotMatch(
                    //             "Input axis sizes do not match data "
                    //             "dimensions" );
                    //     }
                    //     if (!_initialized() || _x_vals.size() != N1
                    //         || _y_vals.size() != N2) {
                    //         init( N1, N2 );
                    //     }

                    //     _x_vals = x1;
                    //     _y_vals = x2;
                    //     _f_vals = f;
                    //     // std::memcpy( _x_vals, &x1[ 0 ],
                    //     //              N1 * sizeof( FloatType ) );
                    //     // std::memcpy( _y_vals, &x2[ 0 ],
                    //     //              N2 * sizeof( FloatType ) );
                    //     // for (size_t i = 0; i < N1; ++i) {
                    //     //     std::memcpy( _f_vals[ i ], &f[ i ][ 0 ],
                    //     //                  N2 * sizeof( FloatType ) );
                    //     // }
                    // }

                    virtual void clear() override {
                        parent_t::clear();
                        _f_slopes.clear();
                        _new_c.clear();
                        _new_d.clear();
                        _node_slopes.clear();
                        _node_vals.clear();
                        _accel[ 0 ] = 0;
                        _accel[ 1 ] = 0;
                    }

                    void ready() override {
                        FloatType dy1, dy2, ai, bi, ci;
                        FloatType di;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        for (int nx = 0; nx < _length_x; ++nx) {
                            bi = 2.0 / ( _y_vals[ 1 ] - _y_vals[ 0 ] );
                            ci = 1.0 / ( _y_vals[ 1 ] - _y_vals[ 0 ] );
                            di = 3.0
                               * ( this->f().at( nx, 1 )
                                   - this->f().at( nx, 0 ) )
                               //    * ( _f_vals[ nx ][ 1 ] - _f_vals[ nx ][ 0
                               //    ] )
                               / std::pow( _y_vals[ 1 ] - _y_vals[ 0 ], 2 );

                            _new_c[ 0 ] = ci / bi;
                            _new_d[ 0 ] = di / bi;

                            for (int ny = 1; ny < _length_y - 1; ++ny) {
                                dy1 = _y_vals[ ny ] - _y_vals[ ny - 1 ],
                                dy2 = _y_vals[ ny + 1 ] - _y_vals[ ny ];

                                ai = 1.0 / dy1;
                                bi = 2.0 * ( 1.0 / dy1 + 1.0 / dy2 );
                                ci = 1.0 / dy2;
                                di = 3.0
                                   * ( ( this->f().at( nx, ny )
                                         - this->f().at( nx, ny - 1 ) )
                                           / std::pow( dy1, 2 )
                                       + ( this->f().at( nx, ny + 1 )
                                           - this->f().at( nx, ny ) )
                                             / std::pow( dy2, 2 ) );
                                // di = 3.0
                                //        * ( ( _f_vals[ nx ][ ny ]
                                //              - _f_vals[ nx ][ ny - 1 ] )
                                //                / std::pow( dy1, 2 )
                                //            + ( _f_vals[ nx ][ ny + 1 ]
                                //                - _f_vals[ nx ][ ny ] )
                                //                  / std::pow( dy2, 2 ) );

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
                               * ( this->f().at( nx, _length_y - 1 )
                                   - this->f().at( nx, _length_y - 2 ) )
                               / std::pow( _y_vals[ _length_y - 1 ]
                                               - _y_vals[ _length_y - 2 ],
                                           2 );
                            // di = 3.0
                            //    * ( _f_vals[ nx ][ _length_y - 1 ]
                            //        - _f_vals[ nx ][ _length_y - 2 ] )
                            //    / std::pow( _y_vals[ _length_y - 1 ]
                            //                    - _y_vals[ _length_y - 2 ],
                            //                2 );

                            _new_d[ _length_y - 1 ]
                                = ( di - _new_d[ _length_y - 2 ] * ai )
                                / ( bi - _new_c[ _length_y - 2 ] * ai );

                            _f_slopes.at( nx, _length_y - 1 )
                                // _f_slopes[ nx ][ _length_y - 1 ]
                                = _new_d[ _length_y - 1 ];
                            for (int ny = _length_y - 2; ny >= 0; ny--) {
                                // _f_slopes[ nx ][ ny ]
                                _f_slopes.at( nx, ny )
                                    = _new_d[ ny ]
                                    - _new_c[ ny ]
                                          * _f_slopes.at( nx, ny + 1 );
                                // * _f_slopes[ nx ][ ny + 1 ];
                            }
                        }
                    }

                    FloatType eval_f(
                        const std::vector<FloatType>& x ) override {
                        this->_check_dimensions( x );
                        return this->eval_f( x.at( 0 ), x.at( 1 ) );
                    }

                    FloatType eval_df( const std::vector<FloatType>& x,
                                       size_t wrt ) override {
                        this->_check_dimensions( x );
                        this->_check_wrt( wrt );
                        return this->eval_df( x.at( 0 ), x.at( 1 ), wrt );
                    }

                    FloatType eval_ddf( const std::vector<FloatType>& x,
                                        size_t wrt1, size_t wrt2 ) override {
                        this->_check_dimensions( x );
                        this->_check_wrt( wrt1 );
                        this->_check_wrt( wrt2 );
                        return this->eval_ddf( x.at( 0 ), x.at( 1 ), wrt1,
                                               wrt2 );
                    }

                    FloatType eval_dddf( const std::vector<FloatType>& x,
                                         size_t wrt1, size_t wrt2,
                                         size_t wrt3 ) override {
                        this->_check_dimensions( x );
                        this->_check_wrt( wrt1 );
                        this->_check_wrt( wrt2 );
                        this->_check_wrt( wrt3 );
                        return this->eval_dddf( x.at( 0 ), x.at( 1 ), wrt1,
                                                wrt2, wrt3 );
                    }

                    FloatType eval_f( FloatType x, FloatType y ) {
                        int kx, ky;
                        FloatType dx, X, A, B;
                        FloatType df;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );

                        for (size_t nx = 0; nx < _length_x; nx++) {
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

                    FloatType eval_df( FloatType x, FloatType y, size_t n ) {
                        size_t kx, ky;
                        FloatType df, A, B;
                        FloatType dx, X;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );

                        for (size_t nx = 0; nx < _length_x; nx++) {
                            if (n == 0) {
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

                        if (n == 0) {
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

                    FloatType eval_ddf( FloatType x, FloatType y, size_t n1,
                                        size_t n2 ) {
                        size_t kx, ky;
                        FloatType df, A, B;
                        FloatType dx, X;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );

                        for (size_t nx = 0; nx < _length_x; nx++) {
                            if (n1 + n2 == 0) {
                                _node_vals[ nx ] = _eval_node_f( y, nx, ky );
                            } else if (n1 + n2 == 1) {
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

                        if (n1 + n2 == 0) {
                            return 2.0 * ( B - 2.0 * A + ( A - B ) * 3.0 * X )
                                 / std::pow( dx, 2 );
                        } else if (n1 + n2 == 1) {
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

                    FloatType eval_dddf( FloatType x, FloatType y, size_t n1,
                                         size_t n2, size_t n3 ) {
                        size_t kx, ky;
                        FloatType X, dx;
                        FloatType df, A, B;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );

                        size_t dimsum = n1 + n2 + n3;
                        for (int nx = 0; nx < _length_x; nx++) {
                            if (dimsum == 0) {
                                _node_vals[ nx ] = _eval_node_f( y, nx, ky );
                            } else if (dimsum == 1) {
                                _node_vals[ nx ]
                                    = _eval_node_dfdy( y, nx, ky );
                            } else if (dimsum == 2) {
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

                        if (dimsum == 0) {
                            return 6.0 * ( A - B ) / std::pow( dx, 3 );
                        } else if (dimsum == 1) {
                            return 2.0 * ( B - 2.0 * A + ( A - B ) * 3.0 * X )
                                 / std::pow( dx, 2 );
                        } else if (dimsum == 2) {
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

                    // virtual std::array<FloatType, 4> limits() const override {
                    //     return std::array<FloatType, 4> {
                    //         this->axes( 0 ).front(), this->axes( 0 ).back(),
                    //         this->axes( 1 ).front(), this->axes( 1 ).back()
                    //         // _x_vals[ 0 ], _x_vals[ _length_x - 1 ],
                    //         // _y_vals[ 0 ], _y_vals[ _length_y - 1 ]
                    //     };
                    // }

                protected:
                    bool _initialized() const {
                        return !( this->axes( 0 ).empty()
                                  || this->axes( 1 ).empty()
                                  || this->f().empty() || _f_slopes.empty() );
                        // return ( _length_x > 0 && _length_y > 0
                        //          && _x_vals != nullptr && _y_vals != nullptr
                        //          && _f_vals != nullptr
                        //          && _f_slopes != nullptr );
                    }

                    FloatType _eval_node_f( FloatType y, size_t nx,
                                            size_t ny ) {
                        FloatType df, A, B;
                        FloatType dy, X;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        df = this->f().at( nx, ny + 1 )
                           - this->f().at( nx, ny );
                        dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];
                        X  = ( y - _y_vals[ ny ] ) / dy;
                        A  = _f_slopes.at( nx, ny ) * dy - df;
                        B  = -_f_slopes.at( nx, ny + 1 ) * dy + df;
                        // df = _f_vals[ nx ][ ny + 1 ] - _f_vals[ nx ][ ny ];
                        // dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];
                        // X = ( y - _y_vals[ ny ] ) / dy;
                        // A = _f_slopes[ nx ][ ny ] * dy - df;
                        // B = -_f_slopes[ nx ][ ny + 1 ] * dy + df;

                        // return _f_vals[ nx ][ ny ]
                        return this->f().at( nx, ny )
                             + X
                                   * ( df
                                       + ( 1.0 - X )
                                             * ( A * ( 1.0 - X ) + B * X ) );
                    }

                    FloatType _eval_node_dfdy( FloatType y, size_t nx,
                                               size_t ny ) {
                        FloatType df, A, B;
                        FloatType dy, X;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        df = this->f().at( nx, ny + 1 )
                           - this->f().at( nx, ny );
                        dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];
                        X  = ( y - _y_vals[ ny ] ) / dy;
                        A  = _f_slopes.at( nx, ny ) * dy - df;
                        B  = -_f_slopes.at( nx, ny + 1 ) * dy + df;

                        // df = _f_vals[ nx ][ ny + 1 ] - _f_vals[ nx ][ ny ];
                        // dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];
                        // X = ( y - _y_vals[ ny ] ) / dy;
                        // A = _f_slopes[ nx ][ ny ] * dy - df;
                        // B = -_f_slopes[ nx ][ ny + 1 ] * dy + df;

                        return df / dy
                             + ( 1.0 - 2.0 * X ) * ( A * ( 1.0 - X ) + B * X )
                                   / dy
                             + X * ( 1.0 - X ) * ( B - A ) / dy;
                    }

                    FloatType _eval_node_ddfdydy( FloatType y, size_t nx,
                                                  size_t ny ) {
                        FloatType df, A, B;
                        FloatType dy, X;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        df = this->f().at( nx, ny + 1 )
                           - this->f().at( nx, ny );
                        dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];
                        X  = ( y - _y_vals[ ny ] ) / dy;
                        A  = _f_slopes.at( nx, ny ) * dy - df;
                        B  = -_f_slopes.at( nx, ny + 1 ) * dy + df;


                        // df = _f_vals[ nx ][ ny + 1 ] - _f_vals[ nx ][ ny ];
                        // dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];
                        // X = ( y - _y_vals[ ny ] ) / dy;
                        // A = _f_slopes[ nx ][ ny ] * dy - df;
                        // B = -_f_slopes[ nx ][ ny + 1 ] * dy + df;

                        return 2.0 * ( B - 2.0 * A + ( A - B ) * 3.0 * X )
                             / std::pow( dy, 2 );
                    }

                    FloatType _eval_node_dddfdydydy( FloatType y, size_t nx,
                                                     size_t ny ) {
                        FloatType df, A, B;
                        FloatType dy, X;
                        FloatType *_y_vals = this->axes( 1 ).data();

                        df = this->f().at( nx, ny + 1 )
                           - this->f().at( nx, ny );
                        // df = _f_vals[ nx ][ ny + 1 ] - _f_vals[ nx ][ ny ];
                        dy = _y_vals[ ny + 1 ] - _y_vals[ ny ];

                        X = ( y - _y_vals[ ny ] ) / dy;
                        A = _f_slopes.at( nx, ny ) * dy - df;
                        // A = _f_slopes[ nx ][ ny ] * dy - df;
                        B = -_f_slopes.at( nx, ny + 1 ) * dy + df;
                        // B = -_f_slopes[ nx ][ ny + 1 ] * dy + df;

                        return 6.0 * ( A - B ) / std::pow( dy, 3 );
                    }

                    void _set_node_slopes() {
                        FloatType dx1, dx2, ai, bi, ci;
                        FloatType di;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        bi = 2.0 / ( _x_vals[ 1 ] - _x_vals[ 0 ] );
                        ci = 1.0 / ( _x_vals[ 1 ] - _x_vals[ 0 ] );
                        di = 3.0 * ( _node_vals[ 1 ] - _node_vals[ 0 ] )
                           / std::pow( _x_vals[ 1 ] - _x_vals[ 0 ], 2 );

                        _new_c[ 0 ] = ci / bi;
                        _new_d[ 0 ] = di / bi;

                        for (int nx = 1; nx < _length_x - 1; nx++) {
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
                        for (int nx = _length_x - 2; nx >= 0; nx--) {
                            _node_slopes[ nx ]
                                = _new_d[ nx ]
                                - _new_c[ nx ] * _node_slopes[ nx + 1 ];
                        }
                    }

                    // size_t _length_x;              // Number of x nodes
                    // size_t _length_y;              // Number of y nodes
                    int _accel[ 2 ];  // Indices of x,y point

                    // FloatType *_x_vals = nullptr;  // 1D array of x values
                    // FloatType *_y_vals = nullptr;  // 1D array of y values
                    // std::vector<FloatType> _x_vals {}, _y_vals {};

                    // 2D array of f(x,y) values, f_vals(i,j) = f(x[i], y[j])
                    // NCPA::arrays::SmartVector<FloatType> _f_vals;

                    // Slopes used to generate natural cubic spline solution at
                    // each x = constant node, describes f(x[i], y)
                    NCPA::arrays::SmartVector<FloatType> _f_slopes;
                    // FloatType **_f_vals
                    //     = nullptr;  // 2D array of f(x,y) values,
                    //                 // f_vals[i][j] = f(x[i], y[j])
                    // FloatType **_f_slopes
                    //     = nullptr;  // Slopes used to generate natural cubic
                    //                 // spline solution at each x = constant
                    //                 // node, describes f(x[i], y)

                    // cache variables
                    std::vector<FloatType>
                        _new_c {},        // size() = _length_y - 1
                        _new_d {},        // size() = _length_y
                        _node_vals {},    // size() = _length_x
                        _node_slopes {};  // size() = _length_x

                    // FloatType *_new_c       = nullptr;    // size() =
                    // _length_y - 1 FloatType *_new_d       = nullptr;    //
                    // size() = _length_y FloatType *_node_vals   = nullptr; //
                    // size() = _length_x FloatType *_node_slopes = nullptr; //
                    // size() = _length_x
            };

            // DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(
            //     NaturalSpline2D, NCPA::interpolation::_spline_2d );

        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA
