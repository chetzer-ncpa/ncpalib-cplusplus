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
#include "NCPA/exceptions.hpp"
#include "NCPA/interpolation/abstract_spline.hpp"
#include "NCPA/interpolation/declarations.hpp"
#include "NCPA/interpolation/exceptions.hpp"
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
            //     BicubicSpline2D, NCPA::interpolation::_spline_2d );

            template<typename FloatType>
            class BicubicSpline2D
                : public Cloneable<BicubicSpline2D<FloatType>,
                                   spline2d<FloatType>> {
                    using parent_t = Cloneable<BicubicSpline2D<FloatType>,
                                               spline2d<FloatType>>;
                    using clone_t  = spline2d<FloatType>;

                public:
                using clone_t::init;
                    BicubicSpline2D() = default;

                    virtual ~BicubicSpline2D() { this->clear(); }

                    BicubicSpline2D( const BicubicSpline2D<FloatType>& other )
                        = default;
                    //     :
                    //     parent_t(
                    //         other ) {
                    //     // _length_x   = other._length_x;
                    //     // _length_y   = other._length_y;
                    //     _accel[ 0 ] = other._accel[ 0 ];
                    //     _accel[ 1 ] = other._accel[ 1 ];
                    //     NCPA::arrays::copy( other._x_vals, _length_x,
                    //                         _x_vals );
                    //     NCPA::arrays::copy( other._y_vals, _length_y,
                    //                         _y_vals );
                    //     NCPA::arrays::copy( other._f_vals, _length_x,
                    //                         _length_y, _f_vals );
                    //     NCPA::arrays::copy( other._dfdx_vals, _length_x,
                    //                         _length_y, _dfdx_vals );
                    //     NCPA::arrays::copy( other._dfdy_vals, _length_x,
                    //                         _length_y, _dfdy_vals );
                    //     for (size_t i = 0; i < 4; i++) {
                    //         for (size_t j = 0; j < 4; j++) {
                    //             _dfdx_coeffs[ i ][ j ]
                    //                 = other._dfdx_coeffs[ i ][ j ];
                    //             _dfdy_coeffs[ i ][ j ]
                    //                 = other._dfdy_coeffs[ i ][ j ];
                    //             _f_coeffs[ i ][ j ]
                    //                 = other._f_coeffs[ i ][ j ];
                    //         }
                    //     }
                    // }

                    BicubicSpline2D(
                        BicubicSpline2D<FloatType>&& source ) noexcept
                        = default;

                    friend void swap(
                        BicubicSpline2D<FloatType>& a,
                        BicubicSpline2D<FloatType>& b ) noexcept {
                        using std::swap;
                        swap( dynamic_cast<parent_t&>( a ),
                              dynamic_cast<parent_t&>( b ) );
                        // swap( a._length_x, b._length_x );
                        // swap( a._length_y, b._length_y );
                        swap( a._accel, b._accel );
                        // swap( a._x_vals, b._x_vals );
                        // swap( a._y_vals, b._y_vals );
                        // swap( a._f_vals, b._f_vals );
                        swap( a._dfdx_vals, b._dfdx_vals );
                        swap( a._dfdy_vals, b._dfdy_vals );
                        swap( a._f_coeffs, b._f_coeffs );
                        swap( a._dfdx_coeffs, b._dfdx_coeffs );
                        swap( a._dfdy_coeffs, b._dfdy_coeffs );
                    }

                    BicubicSpline2D<FloatType>& operator=(
                        BicubicSpline2D<FloatType> other ) {
                        swap( *this, other );
                        return *this;
                    }

                    virtual void init( std::vector<size_t> dims ) override {
                        this->_check_dimensions( dims );
                        size_t nx = dims.at( 0 );
                        size_t ny = dims.at( 1 );
                        this->clear();
                        clone_t::init( {nx, ny} );
                        // _length_x   = nx;
                        _accel[ 0 ] = 0;
                        // _length_y   = ny;
                        _accel[ 1 ] = 0;
                        // _x_vals = NCPA::arrays::zeros<FloatType>( _length_x
                        // ); _y_vals = NCPA::arrays::zeros<FloatType>(
                        // _length_y ); _f_vals = NCPA::arrays::zeros<FloatType>(
                        // _length_x,
                        //                                         _length_y );
                        _dfdx_vals.clear();
                        _dfdy_vals.clear();
                        if (nx > 0 && ny > 0) {
                            _dfdx_vals.resize( nx, ny );
                            _dfdy_vals.resize( nx, ny );
                        }

                        // _dfdx_vals = NCPA::arrays::zeros<FloatType>(
                        // _length_x,
                        //                                            _length_y
                        //                                            );
                        // _dfdy_vals = NCPA::arrays::zeros<FloatType>(
                        // _length_x,
                        //                                            _length_y
                        //                                            );
                    }

                    interpolator_type_t interptype() const override {
                        return NCPA::interpolation::interpolator_type_t::
                            LANL_BICUBIC_2D;
                    }



















                    // virtual void fill( size_t N1, size_t N2,
                    //                    const FloatType *x1,
                    //                    const FloatType *x2,
                    //                    const FloatType **f ) override {
                    //     if (!_initialized() || _length_x != N1
                    //         || _length_y != N2) {
                    //         init( N1, N2 );
                    //     }
                    //     std::memcpy( _x_vals, x1, N1 * sizeof( FloatType ) );
                    //     std::memcpy( _y_vals, x2, N2 * sizeof( FloatType ) );
                    //     for (size_t i = 0; i < N1; ++i) {
                    //         std::memcpy( _f_vals[ i ], f[ i ],
                    //                      N2 * sizeof( FloatType ) );
                    //     }
                    // }

                    // virtual void fill(
                    //     const std::vector<FloatType>& x1,
                    //     const std::vector<FloatType>& x2,
                    //     const NCPA::arrays::vector2d_t<FloatType>& f ) override {
                    //     size_t N1 = x1.size(), N2 = x2.size();
                    //     if (N1 != f.size()) {
                    //         throw std::invalid_argument(
                    //             "Vectors must be same size" );
                    //     }
                    //     if (!_initialized() || _length_x != N1
                    //         || _length_y != N2) {
                    //         init( N1, N2 );
                    //     }
                    //     std::memcpy( _x_vals, &x1[ 0 ],
                    //                  N1 * sizeof( FloatType ) );
                    //     std::memcpy( _y_vals, &x2[ 0 ],
                    //                  N2 * sizeof( FloatType ) );
                    //     for (size_t i = 0; i < N1; ++i) {
                    //         std::memcpy( _f_vals[ i ], &f[ i ][ 0 ],
                    //                      N2 * sizeof( FloatType ) );
                    //     }
                    // }

                    virtual void clear() override {
                        parent_t::clear();
                        _dfdx_vals.clear();
                        _dfdy_vals.clear();
                        _accel[ 0 ] = 0;
                        _accel[ 1 ] = 0;
                        for (size_t i = 0; i < 4; ++i) {
                            for (size_t j = 0; j < 4; ++j) {
                                _dfdx_coeffs[i][j] = 0;
                                _dfdy_coeffs[i][j] = 0;
                                _f_coeffs[i][j] = 0;
                            }
                        }
                    }

                    void ready() override {
                        _update_spline_coeffs( _accel[ 0 ], _accel[ 1 ] );
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

                    FloatType eval_f( FloatType x, FloatType y )  {
                        size_t nx, ny;
                        FloatType x_scaled, y_scaled;
                        FloatType result;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        if (!_same_cell( x, y )) {
                            nx = find_segment( x, _x_vals, _length_x,
                                               _accel[ 0 ] );
                            ny = find_segment( y, _y_vals, _length_y,
                                               _accel[ 1 ] );
                            _update_spline_coeffs( nx, ny );
                        } else {
                            nx = _accel[ 0 ];
                            ny = _accel[ 1 ];
                        }

                        x_scaled = ( x - _x_vals[ nx ] )
                                 / ( _x_vals[ nx + 1 ] - _x_vals[ nx ] );
                        y_scaled = ( y - _y_vals[ ny ] )
                                 / ( _y_vals[ ny + 1 ] - _y_vals[ ny ] );

                        result = 0.0;
                        for (size_t j = 0; j < 4; j++) {
                            for (size_t k = 0; k < 4; k++) {
                                result += std::pow( x_scaled, (FloatType)j )
                                        * std::pow( y_scaled, (FloatType)k )
                                        * _f_coeffs[ j ][ k ];
                            }
                        }
                        return result;
                    }

                    FloatType eval_df( FloatType x, FloatType y,
                                     size_t n )  {
                        size_t nx, ny;
                        FloatType dx, dy, x_scaled, y_scaled;
                        FloatType result = 0.0;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        if (!_same_cell( x, y )) {
                            nx = find_segment( x, _x_vals, _length_x,
                                               _accel[ 0 ] );
                            ny = find_segment( y, _y_vals, _length_y,
                                               _accel[ 1 ] );
                            _update_spline_coeffs( nx, ny );
                        } else {
                            nx = _accel[ 0 ];
                            ny = _accel[ 1 ];
                        }

                        dx       = _x_vals[ nx + 1 ] - _x_vals[ nx ];
                        x_scaled = ( x - _x_vals[ nx ] ) / dx;
                        dy       = _y_vals[ ny + 1 ] - _y_vals[ ny ];
                        y_scaled = ( y - _y_vals[ ny ] ) / dy;

                        result = 0.0;
                        if (n == 0) {
                            for (size_t j = 1; j < 4; j++) {
                                for (size_t k = 0; k < 4; k++) {
                                    result
                                        += std::pow( x_scaled,
                                                     (FloatType)( j - 1 ) )
                                         * std::pow( y_scaled, (FloatType)k )
                                         * _f_coeffs[ j ][ k ] * (FloatType)j
                                         / dx;
                                }
                            }
                        } else if (n == 1) {
                            for (size_t j = 0; j < 4; j++) {
                                for (size_t k = 1; k < 4; k++) {
                                    result
                                        += std::pow( x_scaled, (FloatType)j )
                                         * std::pow( y_scaled,
                                                     (FloatType)( k - 1 ) )
                                         * _f_coeffs[ j ][ k ] * (FloatType)k
                                         / dy;
                                }
                            }
                        } else {
                            throw std::logic_error(
                                "Invalid dimension requested, must be 0 or "
                                "1" );
                        }

                        return result;
                    }

                    FloatType eval_ddf( FloatType x, FloatType y, size_t n1,
                                      size_t n2 )  {
                        size_t nx, ny;
                        FloatType dx, dy, x_scaled, y_scaled;
                        FloatType result = 0.0;
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();

                        if (!_same_cell( x, y )) {
                            nx = find_segment( x, _x_vals, _length_x,
                                               _accel[ 0 ] );
                            ny = find_segment( y, _y_vals, _length_y,
                                               _accel[ 1 ] );
                            _update_spline_coeffs( nx, ny );
                        } else {
                            nx = _accel[ 0 ];
                            ny = _accel[ 1 ];
                        }

                        dx       = _x_vals[ nx + 1 ] - _x_vals[ nx ];
                        x_scaled = ( x - _x_vals[ nx ] ) / dx;
                        dy       = _y_vals[ ny + 1 ] - _y_vals[ ny ];
                        y_scaled = ( y - _y_vals[ ny ] ) / dy;

                        result = 0.0;
                        if (n1 == 0 && n2 == 0) {
                            for (size_t j = 2; j < 4; j++) {
                                for (size_t k = 0; k < 4; k++) {
                                    result
                                        += std::pow( x_scaled,
                                                     (FloatType)( j - 2 ) )
                                         * std::pow( y_scaled, (FloatType)k )
                                         * _f_coeffs[ j ][ k ]
                                         * (FloatType)( j * ( j - 1 ) )
                                         / ( dx * dx );
                                }
                            }
                        } else if (n1 == 1 && n2 == 1) {
                            for (size_t j = 0; j < 4; j++) {
                                for (size_t k = 2; k < 4; k++) {
                                    result
                                        += std::pow( x_scaled, (FloatType)j )
                                         * std::pow( y_scaled,
                                                     (FloatType)( k - 2 ) )
                                         * _f_coeffs[ j ][ k ]
                                         * (FloatType)( k * ( k - 1 ) )
                                         / ( dy * dy );
                                }
                            }
                        } else if (( n1 == 0 && n2 == 1 )
                                   || ( n1 == 1 && n2 == 0 )) {
                            for (size_t j = 1; j < 4; j++) {
                                for (size_t k = 1; k < 4; k++) {
                                    result += std::pow( x_scaled,
                                                        (FloatType)( j - 1 ) )
                                            * std::pow( y_scaled,
                                                        (FloatType)( k - 1 ) )
                                            * _f_coeffs[ j ][ k ]
                                            * (FloatType)( j * k )
                                            / ( dx * dy );
                                }
                            }
                        } else {
                            throw std::logic_error(
                                "Invalid dimension requested, must be 0 or "
                                "1" );
                        }

                        return result;
                    }

                    FloatType eval_dddf( FloatType x, FloatType y, size_t n1,
                                       size_t n2, size_t n3 )  {
                        throw NCPA::NotImplementedError(
                            "3rd derivatives not defined for bicubic "
                            "interpolators." );
                    }

                    // virtual std::array<FloatType, 4> limits() const override {
                    //     return std::array<FloatType, 4> {
                    //         _x_vals[ 0 ], _x_vals[ _length_x - 1 ],
                    //         _y_vals[ 0 ], _y_vals[ _length_y - 1 ]
                    //     };
                    // }

                    // virtual interpolator_2d_type_t interptype()
                    //     const override {
                    //     return NCPA::interpolation::interpolator_2d_type_t::
                    //         LANL_BICUBIC;
                    // }


                protected:
                    bool _initialized() const {
                        return !(this->axes(0).empty() || this->axes(1).empty() || this->f().empty()
                            || _dfdx_vals.empty() || _dfdy_vals.empty());
                        // return ( _length_x > 0 && _length_y > 0
                        //          && _x_vals != nullptr && _y_vals != nullptr
                        //          && _f_vals != nullptr && _dfdx_vals != nullptr
                        //          && _dfdy_vals != nullptr );
                    }

                    bool _same_cell( FloatType x, FloatType y ) const {
                        const FloatType *_x_vals = this->axes( 0 ).data();
                        const FloatType *_y_vals = this->axes( 1 ).data();
                        if (_x_vals[ _accel[ 0 ] ] <= x
                            && x <= _x_vals[ _accel[ 0 ] + 1 ]) {
                            if (_y_vals[ _accel[ 1 ] ] <= y
                                && y <= _y_vals[ _accel[ 1 ] + 1 ]) {
                                return true;
                            }
                        }
                        return false;
                    }

                    void _update_spline_coeffs( size_t nx1, size_t ny1 ) {
                        FloatType *_x_vals = this->axes( 0 ).data();
                        FloatType *_y_vals = this->axes( 1 ).data();
                        size_t _length_x   = this->axes( 0 ).size();
                        size_t _length_y   = this->axes( 1 ).size();




                        FloatType A[ 16 ], X[ 16 ];
                        size_t nx1_up, nx1_dn, nx2, nx2_up, nx2_dn, ny1_up,
                            ny1_dn, ny2, ny2_up, ny2_dn;

                        nx1_up = nx1 + 1;
                        nx1_dn = nx1 > 0 ? nx1 - 1 : 0;
                        // nx1_dn = std::max( nx1 - 1, 0 );
                        nx2    = nx1 + 1;
                        ny1_up = ny1 + 1;
                        ny1_dn = ny1 > 0 ? ny1 - 1 : 0;
                        // ny1_dn = std::max( ny1 - 1, 0 );
                        ny2    = ny1 + 1;

                        nx2_up = std::min( nx2 + 1, _length_x - 1 );
                        nx2_dn = nx2 - 1;
                        ny2_up = std::min( ny2 + 1, _length_y - 1 );
                        ny2_dn = ny2 - 1;

                        X[ 0 ] = this->f().at( nx1, ny1 );
                        X[ 1 ] = this->f().at( nx2, ny1 );
                        X[ 2 ] = this->f().at( nx1, ny2 );
                        X[ 3 ] = this->f().at( nx2, ny2 );

                        X[ 4 ] = this->f().at( nx1_up, ny1 ) - this->f().at( nx1_dn, ny1 );
                        X[ 5 ] = this->f().at( nx2_up, ny1 ) - this->f().at( nx2_dn, ny1 );
                        X[ 6 ] = this->f().at( nx1_up, ny2 ) - this->f().at( nx1_dn, ny2 );
                        X[ 7 ] = this->f().at( nx2_up, ny2 ) - this->f().at( nx2_dn, ny2 );

                        X[ 8 ] = this->f().at( nx1, ny1_up ) - this->f().at( nx1, ny1_dn );
                        X[ 9 ] = this->f().at( nx2, ny1_up ) - this->f().at( nx2, ny1_dn );
                        X[ 10 ] = this->f().at( nx1, ny2_up ) - this->f().at( nx1, ny2_dn );
                        X[ 11 ] = this->f().at( nx2, ny2_up ) - this->f().at( nx2, ny2_dn );

                        X[ 12 ] = this->f().at( nx1_up, ny1_up ) - this->f().at( nx1_dn, ny1_dn );
                        X[ 13 ] = this->f().at( nx2_up, ny1_up ) - this->f().at( nx2_dn, ny1_dn );
                        X[ 14 ] = this->f().at( nx1_up, ny2_up ) - this->f().at( nx1_dn, ny2_dn );
                        X[ 15 ] = this->f().at( nx2_up, ny2_up ) - this->f().at( nx2_dn, ny2_dn );


                        // X[ 0 ] = _f_vals[ nx1 ][ ny1 ];
                        // X[ 1 ] = _f_vals[ nx2 ][ ny1 ];
                        // X[ 2 ] = _f_vals[ nx1 ][ ny2 ];
                        // X[ 3 ] = _f_vals[ nx2 ][ ny2 ];

                        // X[ 4 ] = _f_vals[ nx1_up ][ ny1 ]
                        //        - _f_vals[ nx1_dn ][ ny1 ];
                        // X[ 5 ] = _f_vals[ nx2_up ][ ny1 ]
                        //        - _f_vals[ nx2_dn ][ ny1 ];
                        // X[ 6 ] = _f_vals[ nx1_up ][ ny2 ]
                        //        - _f_vals[ nx1_dn ][ ny2 ];
                        // X[ 7 ] = _f_vals[ nx2_up ][ ny2 ]
                        //        - _f_vals[ nx2_dn ][ ny2 ];

                        // X[ 8 ]  = _f_vals[ nx1 ][ ny1_up ]
                        //         - _f_vals[ nx1 ][ ny1_dn ];
                        // X[ 9 ]  = _f_vals[ nx2 ][ ny1_up ]
                        //         - _f_vals[ nx2 ][ ny1_dn ];
                        // X[ 10 ] = _f_vals[ nx1 ][ ny2_up ]
                        //         - _f_vals[ nx1 ][ ny2_dn ];
                        // X[ 11 ] = _f_vals[ nx2 ][ ny2_up ]
                        //         - _f_vals[ nx2 ][ ny2_dn ];

                        // X[ 12 ] = _f_vals[ nx1_up ][ ny1_up ]
                        //         - _f_vals[ nx1_dn ][ ny1_dn ];
                        // X[ 13 ] = _f_vals[ nx2_up ][ ny1_up ]
                        //         - _f_vals[ nx2_dn ][ ny1_dn ];
                        // X[ 14 ] = _f_vals[ nx1_up ][ ny2_up ]
                        //         - _f_vals[ nx1_dn ][ ny2_dn ];
                        // X[ 15 ] = _f_vals[ nx2_up ][ ny2_up ]
                        //         - _f_vals[ nx2_dn ][ ny2_dn ];

                        for (int j = 0; j < 16; j++) {
                            A[ j ] = 0.0;
                            for (int k = 0; k < 16; k++) {
                                // A[ j ] +=
                                // BicubicSpline2D<FloatType,FloatType>::_bicubic_mat[
                                // j ][ k ] * X[ k ];
                                A[ j ]
                                    += constants<FloatType>::
                                           bicubic_conversion_matrix[ j ][ k ]
                                     * X[ k ];
                            }
                        }
                        for (int j = 0; j < 4; j++) {
                            for (int k = 0; k < 4; k++) {
                                _f_coeffs[ j ][ k ] = A[ k * 4 + j ];
                            }
                        }
                    }

                    // size_t _length_x;              // Number of x nodes
                    // size_t _length_y;              // Number of y nodes
                    int _accel[ 2 ];  // Indices of x,y point

                    // FloatType *_x_vals = nullptr;  // 1D array of x values
                    // FloatType *_y_vals = nullptr;  // 1D array of y values

                    // FloatType **_f_vals
                    //     = nullptr;  // 2D array of f(x,y) values,
                    //                 // f_vals[i][j] = f(x[i], y[j])
                    typename parent_t::FType _dfdx_vals, _dfdy_vals;
                    // FloatType **_dfdx_vals
                    //     = nullptr;  // 2D array of df/dx(x,y) values,
                    //                 // fdx_vals[i][j] = df/dy(x[i], y[j])
                    // FloatType **_dfdy_vals
                    //     = nullptr;  // 2D array of df/dy(x,y) values,
                    //                 // fdy_vals[i][j] = df/dy(x[i], y[j])

                    // cache variables
                    FloatType _f_coeffs[ 4 ][ 4 ];  // Array of coefficients to
                                                    // compute f
                    FloatType _dfdx_coeffs[ 4 ][ 4 ];  // Array of coefficients
                                                       // to compute df/dx
                    FloatType _dfdy_coeffs[ 4 ][ 4 ];  // Array of coefficients
                                                       // to compute df/dy

                    // static constexpr bicubic_conversion_matrix_t<FloatType>
                    // _bicubic_mat
                    //     = bicubic_conversion_matrix<FloatType>();
            };
        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA
