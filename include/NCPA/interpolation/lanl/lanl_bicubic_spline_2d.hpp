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
    NCPA::interpolation::LANL::bicubic_spline_2d<T, U>& a,
    NCPA::interpolation::LANL::bicubic_spline_2d<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        namespace LANL {

            /**
             * 2D INTERPOLATORS
             */
            // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
            //     bicubic_spline_2d, NCPA::interpolation::_spline_2d );

            template<typename INDEPTYPE, typename DEPTYPE>
            class bicubic_spline_2d<INDEPTYPE, DEPTYPE, void,
                                    ENABLE_IF_REAL( INDEPTYPE ),
                                    ENABLE_IF_REAL( DEPTYPE )>
                : public NCPA::interpolation::_spline_2d<INDEPTYPE, DEPTYPE> {
                public:
                    bicubic_spline_2d() : _length_x { 0 }, _length_y { 0 } {}

                    virtual ~bicubic_spline_2d() { this->clear(); }

                    bicubic_spline_2d(
                        const bicubic_spline_2d<INDEPTYPE, DEPTYPE>& other ) :
                        NCPA::interpolation::_spline_2d<INDEPTYPE, DEPTYPE>(
                            other ) {
                        _length_x   = other._length_x;
                        _length_y   = other._length_y;
                        _accel[ 0 ] = other._accel[ 0 ];
                        _accel[ 1 ] = other._accel[ 1 ];
                        NCPA::arrays::copy( other._x_vals, _length_x,
                                            _x_vals );
                        NCPA::arrays::copy( other._y_vals, _length_y,
                                            _y_vals );
                        NCPA::arrays::copy( other._f_vals, _length_x,
                                            _length_y, _f_vals );
                        NCPA::arrays::copy( other._dfdx_vals, _length_x,
                                            _length_y, _dfdx_vals );
                        NCPA::arrays::copy( other._dfdy_vals, _length_x,
                                            _length_y, _dfdy_vals );
                        for ( size_t i = 0; i < 4; i++ ) {
                            for ( size_t j = 0; j < 4; j++ ) {
                                _dfdx_coeffs[ i ][ j ]
                                    = other._dfdx_coeffs[ i ][ j ];
                                _dfdy_coeffs[ i ][ j ]
                                    = other._dfdy_coeffs[ i ][ j ];
                                _f_coeffs[ i ][ j ]
                                    = other._f_coeffs[ i ][ j ];
                            }
                        }
                    }

                    bicubic_spline_2d( linear_spline_1d<INDEPTYPE, DEPTYPE>&&
                                           source ) noexcept :
                        bicubic_spline_2d<INDEPTYPE, DEPTYPE>() {
                        ::swap( *this, source );
                    }

                    friend void ::swap<INDEPTYPE, DEPTYPE>(
                        bicubic_spline_2d<INDEPTYPE, DEPTYPE>& a,
                        bicubic_spline_2d<INDEPTYPE, DEPTYPE>& b ) noexcept;

                    bicubic_spline_2d<INDEPTYPE, DEPTYPE>& operator=(
                        bicubic_spline_2d<INDEPTYPE, DEPTYPE> other ) {
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
                        _dfdx_vals = NCPA::arrays::zeros<DEPTYPE>( _length_x,
                                                                   _length_y );
                        _dfdy_vals = NCPA::arrays::zeros<DEPTYPE>( _length_x,
                                                                   _length_y );
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
                        if ( N1 != f.size() ) {
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
                        if ( _dfdx_vals != nullptr ) {
                            NCPA::arrays::free_array( _dfdx_vals, _length_x,
                                                      _length_y );
                            _dfdx_vals = nullptr;
                        }
                        if ( _dfdy_vals != nullptr ) {
                            NCPA::arrays::free_array( _dfdy_vals, _length_x,
                                                      _length_y );
                            _dfdy_vals = nullptr;
                        }
                        _length_x   = 0;
                        _length_y   = 0;
                        _accel[ 0 ] = 0;
                        _accel[ 1 ] = 0;
                    }

                    void ready() override {}

                    DEPTYPE eval_f( INDEPTYPE x, INDEPTYPE y ) override {
                        size_t nx, ny;
                        INDEPTYPE x_scaled, y_scaled;
                        DEPTYPE result;

                        if ( !_same_cell( x, y ) ) {
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
                        for ( size_t j = 0; j < 4; j++ ) {
                            for ( size_t k = 0; k < 4; k++ ) {
                                result += std::pow( x_scaled, (INDEPTYPE)j )
                                        * std::pow( y_scaled, (INDEPTYPE)k )
                                        * _f_coeffs[ j ][ k ];
                            }
                        }
                        return result;
                    }

                    DEPTYPE eval_df( INDEPTYPE x, INDEPTYPE y,
                                     size_t n ) override {
                        size_t nx, ny;
                        INDEPTYPE dx, dy, x_scaled, y_scaled;
                        DEPTYPE result = 0.0;

                        if ( !_same_cell( x, y ) ) {
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
                        if ( n == 0 ) {
                            for ( size_t j = 1; j < 4; j++ ) {
                                for ( size_t k = 0; k < 4; k++ ) {
                                    result
                                        += std::pow( x_scaled,
                                                     (INDEPTYPE)( j - 1 ) )
                                         * std::pow( y_scaled, (INDEPTYPE)k )
                                         * _f_coeffs[ j ][ k ] * (INDEPTYPE)j
                                         / dx;
                                }
                            }
                        } else if ( n == 1 ) {
                            for ( size_t j = 0; j < 4; j++ ) {
                                for ( size_t k = 1; k < 4; k++ ) {
                                    result
                                        += std::pow( x_scaled, (INDEPTYPE)j )
                                         * std::pow( y_scaled,
                                                     (INDEPTYPE)( k - 1 ) )
                                         * _f_coeffs[ j ][ k ] * (INDEPTYPE)k
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

                    DEPTYPE eval_ddf( INDEPTYPE x, INDEPTYPE y, size_t n1,
                                      size_t n2 ) override {
                        size_t nx, ny;
                        INDEPTYPE dx, dy, x_scaled, y_scaled;
                        DEPTYPE result = 0.0;

                        if ( !_same_cell( x, y ) ) {
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
                        if ( n1 == 0 && n2 == 0 ) {
                            for ( size_t j = 2; j < 4; j++ ) {
                                for ( size_t k = 0; k < 4; k++ ) {
                                    result
                                        += std::pow( x_scaled,
                                                     (INDEPTYPE)( j - 2 ) )
                                         * std::pow( y_scaled, (INDEPTYPE)k )
                                         * _f_coeffs[ j ][ k ]
                                         * (INDEPTYPE)( j * ( j - 1 ) )
                                         / ( dx * dx );
                                }
                            }
                        } else if ( n1 == 1 && n2 == 1 ) {
                            for ( size_t j = 0; j < 4; j++ ) {
                                for ( size_t k = 2; k < 4; k++ ) {
                                    result
                                        += std::pow( x_scaled, (INDEPTYPE)j )
                                         * std::pow( y_scaled,
                                                     (INDEPTYPE)( k - 2 ) )
                                         * _f_coeffs[ j ][ k ]
                                         * (INDEPTYPE)( k * ( k - 1 ) )
                                         / ( dy * dy );
                                }
                            }
                        } else if ( ( n1 == 0 && n2 == 1 )
                                    || ( n1 == 1 && n2 == 0 ) ) {
                            for ( size_t j = 1; j < 4; j++ ) {
                                for ( size_t k = 1; k < 4; k++ ) {
                                    result += std::pow( x_scaled,
                                                        (INDEPTYPE)( j - 1 ) )
                                            * std::pow( y_scaled,
                                                        (INDEPTYPE)( k - 1 ) )
                                            * _f_coeffs[ j ][ k ]
                                            * (INDEPTYPE)( j * k )
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

                    DEPTYPE eval_dddf( INDEPTYPE x, INDEPTYPE y, size_t n1,
                                       size_t n2, size_t n3 ) override {
                        throw NCPA::NotImplementedError(
                            "3rd derivatives not defined for bicubic "
                            "interpolators." );
                    }

                    virtual std::array<INDEPTYPE,4> limits() const override {
                        return std::array<INDEPTYPE,4>{
                            _x_vals[0], _x_vals[_length_x-1], _y_vals[0], _y_vals[_length_y-1]
                        };
                    }

                    virtual interpolator_2d_type_t interptype() const override {
                        return NCPA::interpolation::interpolator_2d_type_t::LANL_BICUBIC;
                    }


                protected:
                    bool _initialized() const {
                        return ( _length_x > 0 && _length_y > 0
                                 && _x_vals != nullptr && _y_vals != nullptr
                                 && _f_vals != nullptr && _dfdx_vals != nullptr
                                 && _dfdy_vals != nullptr );
                    }

                    bool _same_cell( INDEPTYPE x, INDEPTYPE y ) const {
                        if ( _x_vals[ _accel[ 0 ] ] <= x
                             && x <= _x_vals[ _accel[ 0 ] + 1 ] ) {
                            if ( _y_vals[ _accel[ 1 ] ] <= y
                                 && y <= _y_vals[ _accel[ 1 ] + 1 ] ) {
                                return true;
                            }
                        }
                        return false;
                    }

                    void _update_spline_coeffs( size_t nx1, size_t ny1 ) {
                        DEPTYPE A[ 16 ], X[ 16 ];
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

                        X[ 0 ] = _f_vals[ nx1 ][ ny1 ];
                        X[ 1 ] = _f_vals[ nx2 ][ ny1 ];
                        X[ 2 ] = _f_vals[ nx1 ][ ny2 ];
                        X[ 3 ] = _f_vals[ nx2 ][ ny2 ];

                        X[ 4 ] = _f_vals[ nx1_up ][ ny1 ]
                               - _f_vals[ nx1_dn ][ ny1 ];
                        X[ 5 ] = _f_vals[ nx2_up ][ ny1 ]
                               - _f_vals[ nx2_dn ][ ny1 ];
                        X[ 6 ] = _f_vals[ nx1_up ][ ny2 ]
                               - _f_vals[ nx1_dn ][ ny2 ];
                        X[ 7 ] = _f_vals[ nx2_up ][ ny2 ]
                               - _f_vals[ nx2_dn ][ ny2 ];

                        X[ 8 ] = _f_vals[ nx1 ][ ny1_up ]
                               - _f_vals[ nx1 ][ ny1_dn ];
                        X[ 9 ] = _f_vals[ nx2 ][ ny1_up ]
                               - _f_vals[ nx2 ][ ny1_dn ];
                        X[ 10 ] = _f_vals[ nx1 ][ ny2_up ]
                                - _f_vals[ nx1 ][ ny2_dn ];
                        X[ 11 ] = _f_vals[ nx2 ][ ny2_up ]
                                - _f_vals[ nx2 ][ ny2_dn ];

                        X[ 12 ] = _f_vals[ nx1_up ][ ny1_up ]
                                - _f_vals[ nx1_dn ][ ny1_dn ];
                        X[ 13 ] = _f_vals[ nx2_up ][ ny1_up ]
                                - _f_vals[ nx2_dn ][ ny1_dn ];
                        X[ 14 ] = _f_vals[ nx1_up ][ ny2_up ]
                                - _f_vals[ nx1_dn ][ ny2_dn ];
                        X[ 15 ] = _f_vals[ nx2_up ][ ny2_up ]
                                - _f_vals[ nx2_dn ][ ny2_dn ];

                        for ( int j = 0; j < 16; j++ ) {
                            A[ j ] = 0.0;
                            for ( int k = 0; k < 16; k++ ) {
                                // A[ j ] +=
                                // bicubic_spline_2d<INDEPTYPE,DEPTYPE>::_bicubic_mat[
                                // j ][ k ] * X[ k ];
                                A[ j ]
                                    += constants<DEPTYPE>::
                                           bicubic_conversion_matrix[ j ][ k ]
                                     * X[ k ];
                            }
                        }
                        for ( int j = 0; j < 4; j++ ) {
                            for ( int k = 0; k < 4; k++ ) {
                                _f_coeffs[ j ][ k ] = A[ k * 4 + j ];
                            }
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
                    DEPTYPE **_dfdx_vals
                        = nullptr;  // 2D array of df/dx(x,y) values,
                                    // fdx_vals[i][j] = df/dy(x[i], y[j])
                    DEPTYPE **_dfdy_vals
                        = nullptr;  // 2D array of df/dy(x,y) values,
                                    // fdy_vals[i][j] = df/dy(x[i], y[j])

                    // cache variables
                    DEPTYPE _f_coeffs[ 4 ][ 4 ];  // Array of coefficients to
                                                  // compute f
                    DEPTYPE _dfdx_coeffs[ 4 ][ 4 ];  // Array of coefficients
                                                     // to compute df/dx
                    DEPTYPE _dfdy_coeffs[ 4 ][ 4 ];  // Array of coefficients
                                                     // to compute df/dy

                    // static constexpr bicubic_conversion_matrix_t<DEPTYPE>
                    // _bicubic_mat
                    //     = bicubic_conversion_matrix<DEPTYPE>();
            };

            DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(
                bicubic_spline_2d, NCPA::interpolation::_spline_2d );

        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap(
    NCPA::interpolation::LANL::bicubic_spline_2d<T, U>& a,
    NCPA::interpolation::LANL::bicubic_spline_2d<T, U>& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::interpolation::_spline_2d<T, U>&>( a ),
            dynamic_cast<NCPA::interpolation::_spline_2d<T, U>&>( b ) );
    swap( a._length_x, b._length_x );
    swap( a._length_y, b._length_y );
    swap( a._accel, b._accel );
    swap( a._x_vals, b._x_vals );
    swap( a._y_vals, b._y_vals );
    swap( a._f_vals, b._f_vals );
    swap( a._dfdx_vals, b._dfdx_vals );
    swap( a._dfdy_vals, b._dfdy_vals );
    swap( a._f_coeffs, b._f_coeffs );
    swap( a._dfdx_coeffs, b._dfdx_coeffs );
    swap( a._dfdy_coeffs, b._dfdy_coeffs );
}
