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
#include "NCPA/interpolation/abstract_spline_3d.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/lanl/lanl_common.hpp"
#include "NCPA/interpolation/lanl/lanl_declarations.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

namespace NCPA {
    namespace interpolation {
        namespace LANL {

            /**
             * 3D INTERPOLATORS
             */
            // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
            //     hybrid_spline_3d, NCPA::interpolation::_spline_2d );

            template<typename INDEPTYPE, typename DEPTYPE>
            class hybrid_spline_3d<INDEPTYPE, DEPTYPE, void,
                                   ENABLE_IF_REAL( INDEPTYPE ),
                                   ENABLE_IF_REAL( DEPTYPE )>
                : public NCPA::interpolation::_spline_3d<INDEPTYPE, DEPTYPE> {
                public:
                    virtual ~hybrid_spline_3d() { clear(); }

                    void init( size_t nx, size_t ny, size_t nz ) override {
                        this->clear();
                        _length_x   = nx;
                        _accel[ 0 ] = 0;
                        _length_y   = ny;
                        _accel[ 1 ] = 0;
                        _length_z   = nz;
                        _accel[ 2 ] = 0;
                        _x_vals = NCPA::arrays::zeros<INDEPTYPE>( _length_x );
                        _y_vals = NCPA::arrays::zeros<INDEPTYPE>( _length_y );
                        _z_vals = NCPA::arrays::zeros<INDEPTYPE>( _length_z );
                        _f_vals = NCPA::arrays::zeros<DEPTYPE>(
                            _length_x, _length_y, _length_z );
                        _f_slopes = NCPA::arrays::zeros<DEPTYPE>(
                            _length_x, _length_y, _length_z );
                        _dfdx_slopes = NCPA::arrays::zeros<DEPTYPE>(
                            _length_x, _length_y, _length_z );
                        _dfdy_slopes = NCPA::arrays::zeros<DEPTYPE>(
                            _length_x, _length_y, _length_z );
                    }

                    virtual void fill( size_t N1, size_t N2, size_t N3,
                                       const INDEPTYPE *x1,
                                       const INDEPTYPE *x2,
                                       const INDEPTYPE *x3,
                                       const DEPTYPE ***f ) override {
                        if ( !_initialized() || _length_x != N1
                             || _length_y != N2 || _length_z != N3 ) {
                            init( N1, N2, N3 );
                        }
                        std::memcpy( _x_vals, x1, N1 * sizeof( INDEPTYPE ) );
                        std::memcpy( _y_vals, x2, N2 * sizeof( INDEPTYPE ) );
                        std::memcpy( _z_vals, x3, N3 * sizeof( INDEPTYPE ) );
                        for ( size_t i = 0; i < N1; ++i ) {
                            for ( size_t j = 0; j < N2; ++j ) {
                                std::memcpy( _f_vals[ i ][ j ], f[ i ][ j ],
                                             N3 * sizeof( DEPTYPE ) );
                            }
                        }
                    }

                    virtual void fill(
                        const std::vector<INDEPTYPE>& x1,
                        const std::vector<INDEPTYPE>& x2,
                        const std::vector<INDEPTYPE>& x3,
                        const NCPA::arrays::vector3d_t<DEPTYPE>& f ) override {
                        size_t N1 = x1.size(), N2 = x2.size(), N3 = x3.size();
                        size_t fN1, fN2, fN3;
                        f.size3d( fN1, fN2, fN3 );
                        if ( N1 != fN1 || N2 != fN2 || N3 != fN3 ) {
                            throw std::invalid_argument(
                                "Vector sizes must agree" );
                        }
                        if ( !_initialized() || _length_x != N1
                             || _length_y != N2 || _length_z != N3 ) {
                            init( N1, N2, N3 );
                        }
                        std::memcpy( _x_vals, &x1[ 0 ],
                                     N1 * sizeof( INDEPTYPE ) );
                        std::memcpy( _y_vals, &x2[ 0 ],
                                     N2 * sizeof( INDEPTYPE ) );
                        std::memcpy( _z_vals, &x3[ 0 ],
                                     N3 * sizeof( INDEPTYPE ) );
                        for ( size_t i = 0; i < N1; ++i ) {
                            for ( size_t j = 0; j < N2; ++j ) {
                                std::memcpy( _f_vals[ i ][ j ],
                                             &f[ i ][ j ][ 0 ],
                                             N3 * sizeof( DEPTYPE ) );
                            }
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
                        if ( _z_vals != nullptr ) {
                            delete[] _z_vals;
                            _z_vals = nullptr;
                        }
                        if ( _f_vals != nullptr ) {
                            NCPA::arrays::free_array( _f_vals, _length_x,
                                                      _length_y, _length_z );
                            _f_vals = nullptr;
                        }
                        if ( _f_slopes != nullptr ) {
                            NCPA::arrays::free_array( _f_slopes, _length_x,
                                                      _length_y, _length_z );
                            _f_slopes = nullptr;
                        }
                        if ( _dfdx_slopes != nullptr ) {
                            NCPA::arrays::free_array( _dfdx_slopes, _length_x,
                                                      _length_y, _length_z );
                            _dfdx_slopes = nullptr;
                        }
                        if ( _dfdy_slopes != nullptr ) {
                            NCPA::arrays::free_array( _dfdy_slopes, _length_x,
                                                      _length_y, _length_z );
                            _dfdy_slopes = nullptr;
                        }
                        _length_x   = 0;
                        _length_y   = 0;
                        _length_z   = 0;
                        _accel[ 0 ] = 0;
                        _accel[ 1 ] = 0;
                        _accel[ 2 ] = 0;
                    }

                    void ready() override {
                        INDEPTYPE dz1, dz2, ai, bi, ci;
                        DEPTYPE di;
                        INDEPTYPE *new_c
                            = NCPA::arrays::zeros<INDEPTYPE>( _length_z - 1 );
                        DEPTYPE *new_d
                            = NCPA::arrays::zeros<DEPTYPE>( _length_z );

                        for ( size_t mx = 0; mx < _length_x; mx++ ) {
                            for ( size_t my = 0; my < _length_y; my++ ) {
                                bi = 2.0 / ( _z_vals[ 1 ] - _z_vals[ 0 ] );
                                ci = 1.0 / ( _z_vals[ 1 ] - _z_vals[ 0 ] );
                                di = 3.0
                                   * ( _f_vals[ mx ][ my ][ 1 ]
                                       - _f_vals[ mx ][ my ][ 0 ] )
                                   / std::pow( _z_vals[ 1 ] - _z_vals[ 0 ],
                                               2.0 );

                                new_c[ 0 ] = ci / bi;
                                new_d[ 0 ] = di / bi;

                                for ( size_t i = 1; i < _length_z - 1; i++ ) {
                                    dz1 = _z_vals[ i ] - _z_vals[ i - 1 ];
                                    dz2 = _z_vals[ i + 1 ] - _z_vals[ i ];

                                    ai = 1.0 / dz1;
                                    bi = 2.0 * ( 1.0 / dz1 + 1.0 / dz2 );
                                    ci = 1.0 / dz2;
                                    di = 3.0
                                       * ( ( _f_vals[ mx ][ my ][ i ]
                                             - _f_vals[ mx ][ my ][ i - 1 ] )
                                               / std::pow( dz1, 2.0 )
                                           + ( _f_vals[ mx ][ my ][ i + 1 ]
                                               - _f_vals[ mx ][ my ][ i ] )
                                                 / std::pow( dz2, 2.0 ) );

                                    new_c[ i ]
                                        = ci / ( bi - new_c[ i - 1 ] * ai );
                                    new_d[ i ] = ( di - new_d[ i - 1 ] * ai )
                                               / ( bi - new_c[ i - 1 ] * ai );
                                }

                                ai = 1.0
                                   / ( _z_vals[ _length_z - 1 ]
                                       - _z_vals[ _length_z - 2 ] );
                                bi = 2.0
                                   / ( _z_vals[ _length_z - 1 ]
                                       - _z_vals[ _length_z - 2 ] );
                                di = 3.0
                                   * ( _f_vals[ mx ][ my ][ _length_z - 1 ]
                                       - _f_vals[ mx ][ my ][ _length_z - 2 ] )
                                   / std::pow( _z_vals[ _length_z - 1 ]
                                                   - _z_vals[ _length_z - 2 ],
                                               2.0 );

                                new_d[ _length_z - 1 ]
                                    = ( di - new_d[ _length_z - 2 ] * ai )
                                    / ( bi - new_c[ _length_z - 2 ] * ai );

                                _f_slopes[ mx ][ my ][ _length_z - 1 ]
                                    = new_d[ _length_z - 1 ];
                                for ( int i = _length_z - 2; i >= 0; i-- ) {
                                    _f_slopes[ mx ][ my ][ i ]
                                        = new_d[ i ]
                                        - new_c[ i ]
                                              * _f_slopes[ mx ][ my ][ i + 1 ];
                                }
                            }
                        }

                        size_t mx_up, mx_dn, my_up, my_dn;
                        DEPTYPE dfdx[ _length_x ][ _length_y ][ _length_z ];
                        DEPTYPE dfdy[ _length_x ][ _length_y ][ _length_z ];

                        for ( size_t mx = 0; mx < _length_x; mx++ ) {
                            for ( size_t my = 0; my < _length_y; my++ ) {
                                for ( size_t mz = 0; mz < _length_z; mz++ ) {
                                    mx_up = std::min( mx + 1, _length_x - 1 );
                                    my_up = std::min( my + 1, _length_y - 1 );
                                    mx_dn = ( mx > 1 ? mx - 1 : 0 );
                                    my_dn = ( my > 1 ? my - 1 : 0 );
                                    // mx_dn = std::max( mx - 1, 0 );
                                    // my_dn = std::max( my - 1, 0 );

                                    dfdx[ mx ][ my ][ mz ]
                                        = ( _f_vals[ mx_up ][ my ][ mz ]
                                            - _f_vals[ mx_dn ][ my ][ mz ] )
                                        / ( _x_vals[ mx_up ]
                                            - _x_vals[ mx_dn ] );
                                    dfdy[ mx ][ my ][ mz ]
                                        = ( _f_vals[ mx ][ my_up ][ mz ]
                                            - _f_vals[ mx ][ my_dn ][ mz ] )
                                        / ( _y_vals[ my_up ]
                                            - _y_vals[ my_dn ] );
                                }
                            }
                        }


                        for ( size_t mx = 0; mx < _length_x; mx++ ) {
                            for ( size_t my = 0; my < _length_y; my++ ) {
                                bi = 2.0 / ( _z_vals[ 1 ] - _z_vals[ 0 ] );
                                ci = 1.0 / ( _z_vals[ 1 ] - _z_vals[ 0 ] );
                                di = 3.0
                                   * ( dfdx[ mx ][ my ][ 1 ]
                                       - dfdx[ mx ][ my ][ 0 ] )
                                   / std::pow( _z_vals[ 1 ] - _z_vals[ 0 ],
                                               2 );

                                new_c[ 0 ] = ci / bi;
                                new_d[ 0 ] = di / bi;

                                for ( size_t i = 1; i < _length_z - 1; i++ ) {
                                    dz1 = _z_vals[ i ] - _z_vals[ i - 1 ],
                                    dz2 = _z_vals[ i + 1 ] - _z_vals[ i ];

                                    ai = 1.0 / dz1;
                                    bi = 2.0 * ( 1.0 / dz1 + 1.0 / dz2 );
                                    ci = 1.0 / dz2;
                                    di = 3.0
                                       * ( ( dfdx[ mx ][ my ][ i ]
                                             - dfdx[ mx ][ my ][ i - 1 ] )
                                               / std::pow( dz1, 2 )
                                           + ( dfdx[ mx ][ my ][ i + 1 ]
                                               - dfdx[ mx ][ my ][ i ] )
                                                 / std::pow( dz2, 2 ) );

                                    new_c[ i ]
                                        = ci / ( bi - new_c[ i - 1 ] * ai );
                                    new_d[ i ] = ( di - new_d[ i - 1 ] * ai )
                                               / ( bi - new_c[ i - 1 ] * ai );
                                }

                                ai = 1.0
                                   / ( _z_vals[ _length_z - 1 ]
                                       - _z_vals[ _length_z - 2 ] );
                                bi = 2.0
                                   / ( _z_vals[ _length_z - 1 ]
                                       - _z_vals[ _length_z - 2 ] );
                                di = 3.0
                                   * ( dfdx[ mx ][ my ][ _length_z - 1 ]
                                       - dfdx[ mx ][ my ][ _length_z - 2 ] )
                                   / std::pow( _z_vals[ _length_z - 1 ]
                                                   - _z_vals[ _length_z - 2 ],
                                               2.0 );

                                new_d[ _length_z - 1 ]
                                    = ( di - new_d[ _length_z - 2 ] * ai )
                                    / ( bi - new_c[ _length_z - 2 ] * ai );
                                _dfdx_slopes[ mx ][ my ][ _length_z - 1 ]
                                    = new_d[ _length_z - 1 ];
                                for ( int i = _length_z - 2; i >= 0; i-- ) {
                                    _dfdx_slopes[ mx ][ my ][ i ]
                                        = new_d[ i ]
                                        - new_c[ i ]
                                              * _dfdx_slopes[ mx ][ my ]
                                                            [ i + 1 ];
                                }
                            }
                        }

                        for ( size_t mx = 0; mx < _length_x; mx++ ) {
                            for ( size_t my = 0; my < _length_y; my++ ) {
                                bi = 2.0 / ( _z_vals[ 1 ] - _z_vals[ 0 ] );
                                ci = 1.0 / ( _z_vals[ 1 ] - _z_vals[ 0 ] );
                                di = 3.0
                                   * ( dfdy[ mx ][ my ][ 1 ]
                                       - dfdy[ mx ][ my ][ 0 ] )
                                   / std::pow( _z_vals[ 1 ] - _z_vals[ 0 ],
                                               2 );

                                new_c[ 0 ] = ci / bi;
                                new_d[ 0 ] = di / bi;

                                for ( size_t i = 1; i < _length_z - 1; i++ ) {
                                    dz1 = _z_vals[ i ] - _z_vals[ i - 1 ];
                                    dz2 = _z_vals[ i + 1 ] - _z_vals[ i ];

                                    ai = 1.0 / dz1;
                                    bi = 2.0 * ( 1.0 / dz1 + 1.0 / dz2 );
                                    ci = 1.0 / dz2;
                                    di = 3.0
                                       * ( ( dfdy[ mx ][ my ][ i ]
                                             - dfdy[ mx ][ my ][ i - 1 ] )
                                               / std::pow( dz1, 2 )
                                           + ( dfdy[ mx ][ my ][ i + 1 ]
                                               - dfdy[ mx ][ my ][ i ] )
                                                 / std::pow( dz2, 2 ) );

                                    new_c[ i ]
                                        = ci / ( bi - new_c[ i - 1 ] * ai );
                                    new_d[ i ] = ( di - new_d[ i - 1 ] * ai )
                                               / ( bi - new_c[ i - 1 ] * ai );
                                }

                                ai = 1.0
                                   / ( _z_vals[ _length_z - 1 ]
                                       - _z_vals[ _length_z - 2 ] );
                                bi = 2.0
                                   / ( _z_vals[ _length_z - 1 ]
                                       - _z_vals[ _length_z - 2 ] );
                                di = 3.0
                                   * ( dfdy[ mx ][ my ][ _length_z - 1 ]
                                       - dfdy[ mx ][ my ][ _length_z - 2 ] )
                                   / std::pow( _z_vals[ _length_z - 1 ]
                                                   - _z_vals[ _length_z - 2 ],
                                               2.0 );

                                new_d[ _length_z - 1 ]
                                    = ( di - new_d[ _length_z - 2 ] * ai )
                                    / ( bi - new_c[ _length_z - 2 ] * ai );
                                _dfdy_slopes[ mx ][ my ][ _length_z - 1 ]
                                    = new_d[ _length_z - 1 ];
                                for ( int i = _length_z - 2; i >= 0; i-- ) {
                                    _dfdy_slopes[ mx ][ my ][ i ]
                                        = new_d[ i ]
                                        - new_c[ i ]
                                              * _dfdy_slopes[ mx ][ my ]
                                                            [ i + 1 ];
                                }
                            }
                        }
                        delete[] new_c;
                        delete[] new_d;
                    }

                    DEPTYPE eval_f( INDEPTYPE x, INDEPTYPE y,
                                    INDEPTYPE z ) override {
                        size_t kx, ky, kz;
                        INDEPTYPE dx, dy, x_scaled, y_scaled;
                        DEPTYPE result, X_vec[ 16 ], A_vec[ 16 ];

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );
                        kz = find_segment( z, _z_vals, _length_z,
                                           _accel[ 2 ] );

                        dx       = _x_vals[ kx + 1 ] - _x_vals[ kx ];
                        x_scaled = ( x - _x_vals[ kx ] ) / dx;
                        dy       = _y_vals[ ky + 1 ] - _y_vals[ ky ];
                        y_scaled = ( y - _y_vals[ ky ] ) / dy;

                        X_vec[ 0 ] = _eval_node_f( z, kx, ky, kz );
                        X_vec[ 8 ] = _finite_diff_dfdy( z, kx, ky, kz ) * dy;
                        X_vec[ 1 ] = _eval_node_f( z, kx + 1, ky, kz );
                        X_vec[ 9 ]
                            = _finite_diff_dfdy( z, kx + 1, ky, kz ) * dy;
                        X_vec[ 2 ] = _eval_node_f( z, kx, ky + 1, kz );
                        X_vec[ 10 ]
                            = _finite_diff_dfdy( z, kx, ky + 1, kz ) * dy;
                        X_vec[ 3 ] = _eval_node_f( z, kx + 1, ky + 1, kz );
                        X_vec[ 11 ]
                            = _finite_diff_dfdy( z, kx + 1, ky + 1, kz ) * dy;

                        X_vec[ 4 ] = _finite_diff_dfdx( z, kx, ky, kz ) * dx;
                        X_vec[ 12 ]
                            = _finite_diff_ddfdxdy( z, kx, ky, kz ) * dx * dy;
                        X_vec[ 5 ]
                            = _finite_diff_dfdx( z, kx + 1, ky, kz ) * dx;
                        X_vec[ 13 ] = _finite_diff_ddfdxdy( z, kx + 1, ky, kz )
                                    * dx * dy;
                        X_vec[ 6 ]
                            = _finite_diff_dfdx( z, kx, ky + 1, kz ) * dx;
                        X_vec[ 14 ] = _finite_diff_ddfdxdy( z, kx, ky + 1, kz )
                                    * dx * dy;
                        X_vec[ 7 ]
                            = _finite_diff_dfdx( z, kx + 1, ky + 1, kz ) * dx;
                        X_vec[ 15 ]
                            = _finite_diff_ddfdxdy( z, kx + 1, ky + 1, kz )
                            * dx * dy;

                        for ( size_t j = 0; j < 16; j++ ) {
                            A_vec[ j ] = 0.0;
                            for ( size_t k = 0; k < 16; k++ ) {
                                A_vec[ j ]
                                    += constants<DEPTYPE>::
                                           bicubic_conversion_matrix[ j ][ k ]
                                     * X_vec[ k ];
                            }
                        }

                        result = 0.0;
                        for ( size_t k1 = 0; k1 < 4; k1++ ) {
                            for ( size_t k2 = 0; k2 < 4; k2++ ) {
                                result += A_vec[ k2 * 4 + k1 ]
                                        * std::pow( x_scaled, k1 )
                                        * std::pow( y_scaled, k2 );
                            }
                        }

                        return result;
                    }

                    DEPTYPE eval_df( INDEPTYPE x, INDEPTYPE y, INDEPTYPE z,
                                     size_t n ) override {
                        size_t kx, ky, kz;
                        INDEPTYPE dx, dy, x_scaled, y_scaled;
                        DEPTYPE X_vec[ 16 ], A_vec[ 16 ], result;

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );
                        kz = find_segment( z, _z_vals, _length_z,
                                           _accel[ 2 ] );

                        dx       = _x_vals[ kx + 1 ] - _x_vals[ kx ];
                        x_scaled = ( x - _x_vals[ kx ] ) / dx;
                        dy       = _y_vals[ ky + 1 ] - _y_vals[ ky ];
                        y_scaled = ( y - _y_vals[ ky ] ) / dy;

                        if ( n == 0 || n == 1 ) {
                            // df/dx or df/dy from d/dx or d/dy of bicubic
                            // interpolation of f
                            X_vec[ 0 ] = _eval_node_f( z, kx, ky, kz );
                            X_vec[ 8 ]
                                = _finite_diff_dfdy( z, kx, ky, kz ) * dy;
                            X_vec[ 1 ] = _eval_node_f( z, kx + 1, ky, kz );
                            X_vec[ 9 ]
                                = _finite_diff_dfdy( z, kx + 1, ky, kz ) * dy;
                            X_vec[ 2 ] = _eval_node_f( z, kx, ky + 1, kz );
                            X_vec[ 10 ]
                                = _finite_diff_dfdy( z, kx, ky + 1, kz ) * dy;
                            X_vec[ 3 ] = _eval_node_f( z, kx + 1, ky + 1, kz );
                            X_vec[ 11 ]
                                = _finite_diff_dfdy( z, kx + 1, ky + 1, kz )
                                * dy;

                            X_vec[ 4 ]
                                = _finite_diff_dfdx( z, kx, ky, kz ) * dx;
                            X_vec[ 12 ] = _finite_diff_ddfdxdy( z, kx, ky, kz )
                                        * dx * dy;
                            X_vec[ 5 ]
                                = _finite_diff_dfdx( z, kx + 1, ky, kz ) * dx;
                            X_vec[ 13 ]
                                = _finite_diff_ddfdxdy( z, kx + 1, ky, kz )
                                * dx * dy;
                            X_vec[ 6 ]
                                = _finite_diff_dfdx( z, kx, ky + 1, kz ) * dx;
                            X_vec[ 14 ]
                                = _finite_diff_ddfdxdy( z, kx, ky + 1, kz )
                                * dx * dy;
                            X_vec[ 7 ]
                                = _finite_diff_dfdx( z, kx + 1, ky + 1, kz )
                                * dx;
                            X_vec[ 15 ]
                                = _finite_diff_ddfdxdy( z, kx + 1, ky + 1, kz )
                                * dx * dy;
                        } else {
                            // df/dz from bicubic interpolation of df/dz
                            X_vec[ 0 ] = _eval_node_dfdz( z, kx, ky, kz );
                            X_vec[ 8 ]
                                = _finite_diff_ddfdydz( z, kx, ky, kz ) * dy;
                            X_vec[ 1 ] = _eval_node_dfdz( z, kx + 1, ky, kz );
                            X_vec[ 9 ]
                                = _finite_diff_ddfdydz( z, kx + 1, ky, kz )
                                * dy;
                            X_vec[ 2 ] = _eval_node_dfdz( z, kx, ky + 1, kz );
                            X_vec[ 10 ]
                                = _finite_diff_ddfdydz( z, kx, ky + 1, kz )
                                * dy;
                            X_vec[ 3 ]
                                = _eval_node_dfdz( z, kx + 1, ky + 1, kz );
                            X_vec[ 11 ]
                                = _finite_diff_ddfdydz( z, kx + 1, ky + 1, kz )
                                * dy;

                            X_vec[ 4 ]
                                = _finite_diff_ddfdxdz( z, kx, ky, kz ) * dx;
                            X_vec[ 12 ]
                                = _finite_diff_dddfdxdydz( z, kx, ky, kz ) * dx
                                * dy;
                            X_vec[ 5 ]
                                = _finite_diff_ddfdxdz( z, kx + 1, ky, kz )
                                * dx;
                            X_vec[ 13 ]
                                = _finite_diff_dddfdxdydz( z, kx + 1, ky, kz )
                                * dx * dy;
                            X_vec[ 6 ]
                                = _finite_diff_ddfdxdz( z, kx, ky + 1, kz )
                                * dx;
                            X_vec[ 14 ]
                                = _finite_diff_dddfdxdydz( z, kx, ky + 1, kz )
                                * dx * dy;
                            X_vec[ 7 ]
                                = _finite_diff_ddfdxdz( z, kx + 1, ky + 1, kz )
                                * dx;
                            X_vec[ 15 ] = _finite_diff_dddfdxdydz( z, kx + 1,
                                                                   ky + 1, kz )
                                        * dx * dy;
                        }

                        for ( size_t j = 0; j < 16; j++ ) {
                            A_vec[ j ] = 0;
                            for ( size_t k = 0; k < 16; k++ ) {
                                A_vec[ j ]
                                    += constants<DEPTYPE>::
                                           bicubic_conversion_matrix[ j ][ k ]
                                     * X_vec[ k ];
                            }
                        }

                        result = 0.0;
                        if ( n == 0 ) {
                            for ( size_t k1 = 1; k1 < 4; k1++ ) {
                                for ( size_t k2 = 0; k2 < 4; k2++ ) {
                                    result += A_vec[ k2 * 4 + k1 ]
                                            * std::pow( x_scaled, k1 - 1 )
                                            * std::pow( y_scaled, k2 )
                                            * (DEPTYPE)( k1 ) / dx;
                                }
                            }
                        } else if ( n == 1 ) {
                            for ( size_t k1 = 0; k1 < 4; k1++ ) {
                                for ( size_t k2 = 1; k2 < 4; k2++ ) {
                                    result += A_vec[ k2 * 4 + k1 ]
                                            * std::pow( x_scaled, k1 )
                                            * std::pow( y_scaled, k2 - 1 )
                                            * (DEPTYPE)( k2 ) / dy;
                                }
                            }
                        } else {
                            for ( size_t k1 = 0; k1 < 4; k1++ ) {
                                for ( size_t k2 = 0; k2 < 4; k2++ ) {
                                    result += A_vec[ k2 * 4 + k1 ]
                                            * std::pow( x_scaled, k1 )
                                            * std::pow( y_scaled, k2 );
                                }
                            }
                        }

                        return result;
                    }

                    DEPTYPE eval_ddf( INDEPTYPE x, INDEPTYPE y, INDEPTYPE z,
                                      size_t n1, size_t n2 ) override {
                        size_t kx, ky, kz;
                        INDEPTYPE dx, dy, x_scaled, y_scaled;
                        DEPTYPE X_vec[ 16 ], A_vec[ 16 ], result;

                        kx = find_segment( x, _x_vals, _length_x,
                                           _accel[ 0 ] );
                        ky = find_segment( y, _y_vals, _length_y,
                                           _accel[ 1 ] );
                        kz = find_segment( z, _z_vals, _length_z,
                                           _accel[ 2 ] );

                        dx       = _x_vals[ kx + 1 ] - _x_vals[ kx ];
                        x_scaled = ( x - _x_vals[ kx ] ) / dx;
                        dy       = _y_vals[ ky + 1 ] - _y_vals[ ky ];
                        y_scaled = ( y - _y_vals[ ky ] ) / dy;

                        if ( n1 == 0 && n2 == 0 ) {
                            // ddf/dxdx from d/dx of bicubic interpolation of
                            // df/dx
                            X_vec[ 0 ] = _eval_node_dfdx( z, kx, ky, kz );
                            X_vec[ 8 ]
                                = _finite_diff_ddfdxdy( z, kx, ky, kz ) * dy;
                            X_vec[ 1 ] = _eval_node_dfdx( z, kx + 1, ky, kz );
                            X_vec[ 9 ]
                                = _finite_diff_ddfdxdy( z, kx + 1, ky, kz )
                                * dy;
                            X_vec[ 2 ] = _eval_node_dfdx( z, kx, ky + 1, kz );
                            X_vec[ 10 ]
                                = _finite_diff_ddfdxdy( z, kx, ky + 1, kz )
                                * dy;
                            X_vec[ 3 ]
                                = _eval_node_dfdx( z, kx + 1, ky + 1, kz );
                            X_vec[ 11 ]
                                = _finite_diff_ddfdxdy( z, kx + 1, ky + 1, kz )
                                * dy;

                            X_vec[ 4 ]
                                = _finite_diff_ddfdxdx( z, kx, ky, kz ) * dx;
                            X_vec[ 12 ]
                                = _finite_diff_dddfdxdxdy( z, kx, ky, kz ) * dx
                                * dy;
                            X_vec[ 5 ]
                                = _finite_diff_ddfdxdx( z, kx + 1, ky, kz )
                                * dx;
                            X_vec[ 13 ]
                                = _finite_diff_dddfdxdxdy( z, kx + 1, ky, kz )
                                * dx * dy;
                            X_vec[ 6 ]
                                = _finite_diff_ddfdxdx( z, kx, ky + 1, kz )
                                * dx;
                            X_vec[ 14 ]
                                = _finite_diff_dddfdxdxdy( z, kx, ky + 1, kz )
                                * dx * dy;
                            X_vec[ 7 ]
                                = _finite_diff_ddfdxdx( z, kx + 1, ky + 1, kz )
                                * dx;
                            X_vec[ 15 ] = _finite_diff_dddfdxdxdy( z, kx + 1,
                                                                   ky + 1, kz )
                                        * dx * dy;
                        } else if ( n1 == 1 && n2 == 1 ) {
                            // ddf/dydy from d/dy of bicubic interpolation of
                            // df/dy
                            X_vec[ 0 ] = _eval_node_dfdy( z, kx, ky, kz );
                            X_vec[ 8 ]
                                = _finite_diff_ddfdydy( z, kx, ky, kz ) * dy;
                            X_vec[ 1 ] = _eval_node_dfdy( z, kx + 1, ky, kz );
                            X_vec[ 9 ]
                                = _finite_diff_ddfdydy( z, kx + 1, ky, kz )
                                * dy;
                            X_vec[ 2 ] = _eval_node_dfdy( z, kx, ky + 1, kz );
                            X_vec[ 10 ]
                                = _finite_diff_ddfdydy( z, kx, ky + 1, kz )
                                * dy;
                            X_vec[ 3 ]
                                = _eval_node_dfdy( z, kx + 1, ky + 1, kz );
                            X_vec[ 11 ]
                                = _finite_diff_ddfdydy( z, kx + 1, ky + 1, kz )
                                * dy;

                            X_vec[ 4 ]
                                = _finite_diff_ddfdxdy( z, kx, ky, kz ) * dx;
                            X_vec[ 12 ]
                                = _finite_diff_dddfdxdydy( z, kx, ky, kz ) * dx
                                * dy;
                            X_vec[ 5 ]
                                = _finite_diff_ddfdxdy( z, kx + 1, ky, kz )
                                * dx;
                            X_vec[ 13 ]
                                = _finite_diff_dddfdxdydy( z, kx + 1, ky, kz )
                                * dx * dy;
                            X_vec[ 6 ]
                                = _finite_diff_ddfdxdy( z, kx, ky + 1, kz )
                                * dx;
                            X_vec[ 14 ]
                                = _finite_diff_dddfdxdydy( z, kx, ky + 1, kz )
                                * dx * dy;
                            X_vec[ 7 ]
                                = _finite_diff_ddfdxdy( z, kx + 1, ky + 1, kz )
                                * dx;
                            X_vec[ 15 ] = _finite_diff_dddfdxdydy( z, kx + 1,
                                                                   ky + 1, kz )
                                        * dx * dy;
                        } else if ( n1 == 2 && n2 == 2 ) {
                            // ddf/dzdz from bicubic interpolation of ddf/dzdz
                            X_vec[ 0 ] = _eval_node_ddfdzdz( z, kx, ky, kz );
                            X_vec[ 8 ]
                                = _finite_diff_dddfdydzdz( z, kx, ky, kz )
                                * dy;
                            X_vec[ 1 ]
                                = _eval_node_ddfdzdz( z, kx + 1, ky, kz );
                            X_vec[ 9 ]
                                = _finite_diff_dddfdydzdz( z, kx + 1, ky, kz )
                                * dy;
                            X_vec[ 2 ]
                                = _eval_node_ddfdzdz( z, kx, ky + 1, kz );
                            X_vec[ 10 ]
                                = _finite_diff_dddfdydzdz( z, kx, ky + 1, kz )
                                * dy;
                            X_vec[ 3 ]
                                = _eval_node_ddfdzdz( z, kx + 1, ky + 1, kz );
                            X_vec[ 11 ] = _finite_diff_dddfdydzdz( z, kx + 1,
                                                                   ky + 1, kz )
                                        * dy;

                            X_vec[ 4 ]
                                = _finite_diff_dddfdxdzdz( z, kx, ky, kz )
                                * dx;
                            X_vec[ 12 ]
                                = _finite_diff_ddddfdxdydzdz( z, kx, ky, kz )
                                * dx * dy;
                            X_vec[ 5 ]
                                = _finite_diff_dddfdxdzdz( z, kx + 1, ky, kz )
                                * dx;
                            X_vec[ 13 ] = _finite_diff_ddddfdxdydzdz(
                                              z, kx + 1, ky, kz )
                                        * dx * dy;
                            X_vec[ 6 ]
                                = _finite_diff_dddfdxdzdz( z, kx, ky + 1, kz )
                                * dx;
                            X_vec[ 14 ] = _finite_diff_ddddfdxdydzdz(
                                              z, kx, ky + 1, kz )
                                        * dx * dy;
                            X_vec[ 7 ] = _finite_diff_dddfdxdzdz( z, kx + 1,
                                                                  ky + 1, kz )
                                       * dx;
                            X_vec[ 15 ] = _finite_diff_ddddfdxdydzdz(
                                              z, kx + 1, ky + 1, kz )
                                        * dx * dy;
                        } else if ( ( n1 == 0 && n2 == 1 )
                                    || ( n1 == 1 && n2 == 0 ) ) {
                            // ddf/dxdy from dd/dxdy of bicubic interpolation
                            // of f
                            X_vec[ 0 ] = _eval_node_f( z, kx, ky, kz );
                            X_vec[ 8 ]
                                = _finite_diff_dfdy( z, kx, ky, kz ) * dy;
                            X_vec[ 1 ] = _eval_node_f( z, kx + 1, ky, kz );
                            X_vec[ 9 ]
                                = _finite_diff_dfdy( z, kx + 1, ky, kz ) * dy;
                            X_vec[ 2 ] = _eval_node_f( z, kx, ky + 1, kz );
                            X_vec[ 10 ]
                                = _finite_diff_dfdy( z, kx, ky + 1, kz ) * dy;
                            X_vec[ 3 ] = _eval_node_f( z, kx + 1, ky + 1, kz );
                            X_vec[ 11 ]
                                = _finite_diff_dfdy( z, kx + 1, ky + 1, kz )
                                * dy;

                            X_vec[ 4 ]
                                = _finite_diff_dfdx( z, kx, ky, kz ) * dx;
                            X_vec[ 12 ] = _finite_diff_ddfdxdy( z, kx, ky, kz )
                                        * dx * dy;
                            X_vec[ 5 ]
                                = _finite_diff_dfdx( z, kx + 1, ky, kz ) * dx;
                            X_vec[ 13 ]
                                = _finite_diff_ddfdxdy( z, kx + 1, ky, kz )
                                * dx * dy;
                            X_vec[ 6 ]
                                = _finite_diff_dfdx( z, kx, ky + 1, kz ) * dx;
                            X_vec[ 14 ]
                                = _finite_diff_ddfdxdy( z, kx, ky + 1, kz )
                                * dx * dy;
                            X_vec[ 7 ]
                                = _finite_diff_dfdx( z, kx + 1, ky + 1, kz )
                                * dx;
                            X_vec[ 15 ]
                                = _finite_diff_ddfdxdy( z, kx + 1, ky + 1, kz )
                                * dx * dy;

                        } else if ( ( n1 == 0 && n2 == 2 )
                                    || ( n1 == 2 && n2 == 0 )
                                    || ( n1 == 1 && n2 == 2 )
                                    || ( n1 == 2 && n2 == 1 ) ) {
                            // ddf/dxdz or ddf/dydz from d/dx or d/dy of
                            // bicubic interpolation of df/dz
                            X_vec[ 0 ] = _eval_node_dfdz( z, kx, ky, kz );
                            X_vec[ 8 ]
                                = _finite_diff_ddfdydz( z, kx, ky, kz ) * dy;
                            X_vec[ 1 ] = _eval_node_dfdz( z, kx + 1, ky, kz );
                            X_vec[ 9 ]
                                = _finite_diff_ddfdydz( z, kx + 1, ky, kz )
                                * dy;
                            X_vec[ 2 ] = _eval_node_dfdz( z, kx, ky + 1, kz );
                            X_vec[ 10 ]
                                = _finite_diff_ddfdydz( z, kx, ky + 1, kz )
                                * dy;
                            X_vec[ 3 ]
                                = _eval_node_dfdz( z, kx + 1, ky + 1, kz );
                            X_vec[ 11 ]
                                = _finite_diff_ddfdydz( z, kx + 1, ky + 1, kz )
                                * dy;

                            X_vec[ 4 ]
                                = _finite_diff_ddfdxdz( z, kx, ky, kz ) * dx;
                            X_vec[ 12 ]
                                = _finite_diff_dddfdxdydz( z, kx, ky, kz ) * dx
                                * dy;
                            X_vec[ 5 ]
                                = _finite_diff_ddfdxdz( z, kx + 1, ky, kz )
                                * dx;
                            X_vec[ 13 ]
                                = _finite_diff_dddfdxdydz( z, kx + 1, ky, kz )
                                * dx * dy;
                            X_vec[ 6 ]
                                = _finite_diff_ddfdxdz( z, kx, ky + 1, kz )
                                * dx;
                            X_vec[ 14 ]
                                = _finite_diff_dddfdxdydz( z, kx, ky + 1, kz )
                                * dx * dy;
                            X_vec[ 7 ]
                                = _finite_diff_ddfdxdz( z, kx + 1, ky + 1, kz )
                                * dx;
                            X_vec[ 15 ] = _finite_diff_dddfdxdydz( z, kx + 1,
                                                                   ky + 1, kz )
                                        * dx * dy;
                        }

                        for ( size_t j = 0; j < 16; j++ ) {
                            A_vec[ j ] = 0.0;
                            for ( size_t k = 0; k < 16; k++ ) {
                                A_vec[ j ]
                                    += constants<DEPTYPE>::
                                           bicubic_conversion_matrix[ j ][ k ]
                                     * X_vec[ k ];
                            }
                        }

                        result = 0.0;
                        if ( n1 == 2 && n2 == 2 ) {
                            for ( size_t k1 = 0; k1 < 4; k1++ ) {
                                for ( size_t k2 = 0; k2 < 4; k2++ ) {
                                    result += A_vec[ k2 * 4 + k1 ]
                                            * std::pow( x_scaled, k1 )
                                            * std::pow( y_scaled, k2 );
                                }
                            }
                        } else if ( ( n1 == 0 && n2 == 0 )
                                    || ( n1 == 0 && n2 == 2 )
                                    || ( n1 == 2 && n2 == 0 ) ) {
                            for ( size_t k1 = 1; k1 < 4; k1++ ) {
                                for ( size_t k2 = 0; k2 < 4; k2++ ) {
                                    result += A_vec[ k2 * 4 + k1 ]
                                            * std::pow( x_scaled, k1 - 1 )
                                            * std::pow( y_scaled, k2 )
                                            * (DEPTYPE)( k1 ) / dx;
                                }
                            }
                        } else if ( ( n1 == 1 && n2 == 1 )
                                    || ( n1 == 1 && n2 == 2 )
                                    || ( n1 == 2 && n2 == 1 ) ) {
                            for ( size_t k1 = 0; k1 < 4; k1++ ) {
                                for ( size_t k2 = 1; k2 < 4; k2++ ) {
                                    result += A_vec[ k2 * 4 + k1 ]
                                            * std::pow( x_scaled, k1 )
                                            * std::pow( y_scaled, k2 - 1 )
                                            * (DEPTYPE)( k2 ) / dy;
                                }
                            }
                        } else {
                            for ( size_t k1 = 1; k1 < 4; k1++ ) {
                                for ( size_t k2 = 1; k2 < 4; k2++ ) {
                                    result += A_vec[ k2 * 4 + k1 ]
                                            * std::pow( x_scaled, k1 - 1 )
                                            * std::pow( y_scaled, k2 - 1 )
                                            * (DEPTYPE)( k1 * k2 )
                                            / ( dx * dy );
                                }
                            }
                        }

                        return result;
                    }

                    DEPTYPE eval_dddf( INDEPTYPE x, INDEPTYPE y, INDEPTYPE z,
                                       size_t n1, size_t n2,
                                       size_t n3 ) override {
                        throw NCPA::NotImplementedError(
                            "3rd derivatives not defined for hybrid 3-D "
                            "interpolators." );
                    }

                protected:
                    bool _initialized() const {
                        return ( _length_x > 0 && _length_y > 0
                                 && _length_z > 0 && _x_vals != nullptr
                                 && _y_vals != nullptr && _z_vals != nullptr
                                 && _f_vals != nullptr && _f_slopes != nullptr
                                 && _dfdx_slopes != nullptr
                                 && _dfdy_slopes != nullptr );
                    }

                    DEPTYPE _eval_node_f( double z, size_t kx, size_t ky,
                                          size_t kz ) {
                        INDEPTYPE dz, X;
                        DEPTYPE df, A, B;

                        dz = _z_vals[ kz + 1 ] - _z_vals[ kz ];
                        X  = ( z - _z_vals[ kz ] ) / dz;

                        df = _f_vals[ kx ][ ky ][ kz + 1 ]
                           - _f_vals[ kx ][ ky ][ kz ];
                        A = _f_slopes[ kx ][ ky ][ kz ] * dz - df;
                        B = -_f_slopes[ kx ][ ky ][ kz + 1 ] * dz + df;

                        return _f_vals[ kx ][ ky ][ kz ]
                             + X
                                   * ( df
                                       + ( 1.0 - X )
                                             * ( A * ( 1.0 - X ) + B * X ) );
                    }

                    DEPTYPE _eval_node_dfdx( INDEPTYPE z, size_t kx, size_t ky,
                                             size_t kz ) {
                        size_t kx_up, kx_dn;
                        DEPTYPE dfdx_kz, dfdx_kzp1, df, A, B;
                        INDEPTYPE dz, X;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );
                        // kx_dn = std::max(kx - 1, 0);

                        dfdx_kz = ( _f_vals[ kx_up ][ ky ][ kz ]
                                    - _f_vals[ kx_dn ][ ky ][ kz ] )
                                / ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] );
                        dfdx_kzp1 = ( _f_vals[ kx_up ][ ky ][ kz + 1 ]
                                      - _f_vals[ kx_dn ][ ky ][ kz + 1 ] )
                                  / ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] );

                        dz = _z_vals[ kz + 1 ] - _z_vals[ kz ];
                        X  = ( z - _z_vals[ kz ] ) / dz;

                        df = dfdx_kzp1 - dfdx_kz;
                        A  = _dfdx_slopes[ kx ][ ky ][ kz ] * dz - df;
                        B  = -_dfdx_slopes[ kx ][ ky ][ kz + 1 ] * dz + df;

                        return dfdx_kz
                             + X
                                   * ( df
                                       + ( 1.0 - X )
                                             * ( A * ( 1.0 - X ) + B * X ) );
                    }

                    DEPTYPE _eval_node_dfdy( INDEPTYPE z, size_t kx, size_t ky,
                                             size_t kz ) {
                        size_t ky_up, ky_dn;
                        DEPTYPE dfdy_kz, dfdy_kzp1, df, A, B;
                        INDEPTYPE dz, X;

                        ky_up = std::min( ky + 1, _length_y - 1 );
                        ky_dn = ( ky > 1 ? ky - 1 : 0 );
                        // ky_dn = std::max(ky - 1, 0);

                        dfdy_kz = ( _f_vals[ kx ][ ky_up ][ kz ]
                                    - _f_vals[ kx ][ ky_dn ][ kz ] )
                                / ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] );
                        dfdy_kzp1 = ( _f_vals[ kx ][ ky_up ][ kz + 1 ]
                                      - _f_vals[ kx ][ ky_dn ][ kz + 1 ] )
                                  / ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] );

                        dz = _z_vals[ kz + 1 ] - _z_vals[ kz ];
                        X  = ( z - _z_vals[ kz ] ) / dz;

                        df = dfdy_kzp1 - dfdy_kz;
                        A  = _dfdx_slopes[ kx ][ ky ][ kz ] * dz - df;
                        B  = -_dfdx_slopes[ kx ][ ky ][ kz + 1 ] * dz + df;

                        return dfdy_kz
                             + X
                                   * ( df
                                       + ( 1.0 - X )
                                             * ( A * ( 1.0 - X ) + B * X ) );
                    }

                    DEPTYPE _eval_node_dfdz( INDEPTYPE z, size_t kx, size_t ky,
                                             size_t kz ) {
                        INDEPTYPE dz, X;
                        DEPTYPE df, A, B;

                        dz = _z_vals[ kz + 1 ] - _z_vals[ kz ];
                        X  = ( z - _z_vals[ kz ] ) / dz;

                        df = _f_vals[ kx ][ ky ][ kz + 1 ]
                           - _f_vals[ kx ][ ky ][ kz ];
                        A = _f_slopes[ kx ][ ky ][ kz ] * dz - df;
                        B = -_f_slopes[ kx ][ ky ][ kz + 1 ] * dz + df;

                        return df / dz
                             + ( 1.0 - 2.0 * X ) * ( A * ( 1.0 - X ) + B * X )
                                   / dz
                             + X * ( 1.0 - X ) * ( B - A ) / dz;
                    }

                    DEPTYPE _eval_node_ddfdxdz( INDEPTYPE z, size_t kx,
                                                size_t ky, size_t kz ) {
                        size_t kx_up, kx_dn;
                        DEPTYPE dfdx_kz, dfdx_kzp1, df, A, B;
                        INDEPTYPE dz, X;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );
                        // kx_dn = std::max(kx - 1, 0);

                        dfdx_kz = ( _f_vals[ kx_up ][ ky ][ kz ]
                                    - _f_vals[ kx_dn ][ ky ][ kz ] )
                                / ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] );
                        dfdx_kzp1 = ( _f_vals[ kx_up ][ ky ][ kz + 1 ]
                                      - _f_vals[ kx_dn ][ ky ][ kz + 1 ] )
                                  / ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] );

                        dz = _z_vals[ kz + 1 ] - _z_vals[ kz ];
                        X  = ( z - _z_vals[ kz ] ) / dz;

                        df = dfdx_kzp1 - dfdx_kz;
                        A  = _dfdx_slopes[ kx ][ ky ][ kz ] * dz - df;
                        B  = -_dfdx_slopes[ kx ][ ky ][ kz + 1 ] * dz + df;

                        return df / dz
                             + ( 1.0 - 2.0 * X ) * ( A * ( 1.0 - X ) + B * X )
                                   / dz
                             + X * ( 1.0 - X ) * ( B - A ) / dz;
                    }

                    DEPTYPE _eval_node_ddfdydz( INDEPTYPE z, size_t kx,
                                                size_t ky, size_t kz ) {
                        size_t ky_up, ky_dn;
                        DEPTYPE dfdy_kz, dfdy_kzp1, df, A, B;
                        INDEPTYPE dz, X;

                        ky_up = std::min( ky + 1, _length_y - 1 );
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );
                        // ky_dn = std::max(ky - 1, 0);

                        dfdy_kz = ( _f_vals[ kx ][ ky_up ][ kz ]
                                    - _f_vals[ kx ][ ky_dn ][ kz ] )
                                / ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] );
                        dfdy_kzp1 = ( _f_vals[ kx ][ ky_up ][ kz + 1 ]
                                      - _f_vals[ kx ][ ky_dn ][ kz + 1 ] )
                                  / ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] );

                        dz = _z_vals[ kz + 1 ] - _z_vals[ kz ];
                        X  = ( z - _z_vals[ kz ] ) / dz;

                        df = dfdy_kzp1 - dfdy_kz;
                        A  = _dfdx_slopes[ kx ][ ky ][ kz ] * dz - df;
                        B  = -_dfdx_slopes[ kx ][ ky ][ kz + 1 ] * dz + df;

                        return df / dz
                             + ( 1.0 - 2.0 * X ) * ( A * ( 1.0 - X ) + B * X )
                                   / dz
                             + X * ( 1.0 - X ) * ( B - A ) / dz;
                    }

                    DEPTYPE _eval_node_ddfdzdz( INDEPTYPE z, size_t kx,
                                                size_t ky, size_t kz ) {
                        DEPTYPE df, A, B;
                        INDEPTYPE X, dz;

                        dz = _z_vals[ kz + 1 ] - _z_vals[ kz ];
                        X  = ( z - _z_vals[ kz ] ) / dz;

                        df = _f_vals[ kx ][ ky ][ kz + 1 ]
                           - _f_vals[ kx ][ ky ][ kz ];
                        A = _f_slopes[ kx ][ ky ][ kz ] * dz - df;
                        B = -_f_slopes[ kx ][ ky ][ kz + 1 ] * dz + df;

                        return 2.0 * ( B - 2.0 * A + ( A - B ) * 3.0 * X )
                             / std::pow( dz, 2 );
                    }

                    // Finite difference derivatives for defining the bicubic
                    // spline corners
                    DEPTYPE _finite_diff_dfdx( INDEPTYPE z, size_t kx,
                                               size_t ky, size_t kz ) {
                        size_t kx_up, kx_dn;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        // kx_dn = std::max(kx - 1, 0);
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );

                        return ( _eval_node_f( z, kx_up, ky, kz )
                                 - _eval_node_f( z, kx_dn, ky, kz ) )
                             / ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] );
                    }

                    DEPTYPE _finite_diff_dfdy( INDEPTYPE z, size_t kx,
                                               size_t ky, size_t kz ) {
                        size_t ky_up, ky_dn;

                        ky_up = std::min( ky + 1, _length_y - 1 );
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );
                        // ky_dn = std::max(ky - 1, 0);

                        return ( _eval_node_f( z, kx, ky_up, kz )
                                 - _eval_node_f( z, kx, ky_dn, kz ) )
                             / ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] );
                    }

                    DEPTYPE _finite_diff_ddfdxdx( INDEPTYPE z, size_t kx,
                                                  size_t ky, size_t kz ) {
                        size_t kx_up, kx_dn;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        // kx_dn = std::max(kx - 1, 0);
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );

                        return ( _eval_node_dfdx( z, kx_up, ky, kz )
                                 - _eval_node_dfdx( z, kx_dn, ky, kz ) )
                             / ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] );
                    }

                    DEPTYPE _finite_diff_ddfdydy( INDEPTYPE z, size_t kx,
                                                  size_t ky, size_t kz ) {
                        size_t ky_up, ky_dn;

                        ky_up = std::min( ky + 1, _length_y - 1 );
                        // ky_dn = std::max(ky - 1, 0);
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );

                        return ( _eval_node_dfdy( z, kx, ky_up, kz )
                                 - _eval_node_dfdy( z, kx, ky_dn, kz ) )
                             / ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] );
                    }

                    DEPTYPE _finite_diff_ddfdxdy( INDEPTYPE z, size_t kx,
                                                  size_t ky, size_t kz ) {
                        size_t kx_up, ky_up, kx_dn, ky_dn;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        ky_up = std::min( ky + 1, _length_y - 1 );

                        kx_dn = ( kx > 0 ? kx - 1 : 0 );
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );
                        // kx_dn = std::max(kx - 1, 0);
                        // ky_dn = std::max(ky - 1, 0);

                        return ( _eval_node_f( z, kx_up, ky_up, kz )
                                 - _eval_node_f( z, kx_up, ky_dn, kz )
                                 - _eval_node_f( z, kx_dn, ky_up, kz )
                                 + _eval_node_f( z, kx_dn, ky_dn, kz ) )
                             / ( ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] )
                                 * ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] ) );
                    }

                    DEPTYPE _finite_diff_ddfdxdz( INDEPTYPE z, size_t kx,
                                                  size_t ky, size_t kz ) {
                        size_t kx_up, kx_dn;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        // kx_dn = std::max(kx - 1, 0);
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );

                        return ( _eval_node_dfdz( z, kx_up, ky, kz )
                                 - _eval_node_dfdz( z, kx_dn, ky, kz ) )
                             / ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] );
                    }

                    DEPTYPE _finite_diff_ddfdydz( INDEPTYPE z, size_t kx,
                                                  size_t ky, size_t kz ) {
                        size_t ky_up, ky_dn;

                        ky_up = std::min( ky + 1, _length_y - 1 );
                        // ky_dn = std::max(ky - 1, 0);
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );

                        return ( _eval_node_dfdz( z, kx, ky_up, kz )
                                 - _eval_node_dfdz( z, kx, ky_dn, kz ) )
                             / ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] );
                    }

                    DEPTYPE _finite_diff_dddfdxdxdy( INDEPTYPE z, size_t kx,
                                                     size_t ky, size_t kz ) {
                        size_t kx_up, ky_up, kx_dn, ky_dn;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        ky_up = std::min( ky + 1, _length_y - 1 );

                        // kx_dn = std::max(kx - 1, 0);
                        // ky_dn = std::max(ky - 1, 0);
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );

                        return ( _eval_node_dfdx( z, kx_up, ky_up, kz )
                                 - _eval_node_dfdx( z, kx_up, ky_dn, kz )
                                 - _eval_node_dfdx( z, kx_dn, ky_up, kz )
                                 + _eval_node_dfdx( z, kx_dn, ky_dn, kz ) )
                             / ( ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] )
                                 * ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] ) );
                    }

                    DEPTYPE _finite_diff_dddfdxdydy( INDEPTYPE z, size_t kx,
                                                     size_t ky, size_t kz ) {
                        size_t kx_up, ky_up, kx_dn, ky_dn;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        ky_up = std::min( ky + 1, _length_y - 1 );

                        // kx_dn = std::max(kx - 1, 0);
                        // ky_dn = std::max(ky - 1, 0);
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );

                        return ( _eval_node_dfdy( z, kx_up, ky_up, kz )
                                 - _eval_node_dfdy( z, kx_up, ky_dn, kz )
                                 - _eval_node_dfdy( z, kx_dn, ky_up, kz )
                                 + _eval_node_dfdy( z, kx_dn, ky_dn, kz ) )
                             / ( ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] )
                                 * ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] ) );
                    }

                    DEPTYPE _finite_diff_dddfdxdydz( INDEPTYPE z, size_t kx,
                                                     size_t ky, size_t kz ) {
                        size_t kx_up, ky_up, kx_dn, ky_dn;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        ky_up = std::min( ky + 1, _length_y - 1 );

                        // kx_dn = std::max(kx - 1, 0);
                        // ky_dn = std::max(ky - 1, 0);
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );

                        return ( _eval_node_dfdz( z, kx_up, ky_up, kz )
                                 - _eval_node_dfdz( z, kx_up, ky_dn, kz )
                                 - _eval_node_dfdz( z, kx_dn, ky_up, kz )
                                 + _eval_node_dfdz( z, kx_dn, ky_dn, kz ) )
                             / ( ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] )
                                 * ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] ) );
                    }

                    DEPTYPE _finite_diff_dddfdxdzdz( INDEPTYPE z, size_t kx,
                                                     size_t ky, size_t kz ) {
                        size_t kx_up, kx_dn;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        // kx_dn = std::max(kx - 1, 0);
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );

                        return ( _eval_node_ddfdzdz( z, kx_up, ky, kz )
                                 - _eval_node_ddfdzdz( z, kx_dn, ky, kz ) )
                             / ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] );
                    }

                    DEPTYPE _finite_diff_dddfdydzdz( INDEPTYPE z, size_t kx,
                                                     size_t ky, size_t kz ) {
                        size_t ky_up, ky_dn;

                        ky_up = std::min( ky + 1, _length_y - 1 );
                        // ky_dn = std::max(ky - 1, 0);
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );

                        return ( _eval_node_ddfdzdz( z, kx, ky_up, kz )
                                 - _eval_node_ddfdzdz( z, kx, ky_dn, kz ) )
                             / ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] );
                    }

                    DEPTYPE _finite_diff_ddddfdxdydzdz( INDEPTYPE z, size_t kx,
                                                        size_t ky,
                                                        size_t kz ) {
                        size_t kx_up, ky_up, kx_dn, ky_dn;

                        kx_up = std::min( kx + 1, _length_x - 1 );
                        ky_up = std::min( ky + 1, _length_y - 1 );

                        // kx_dn = std::max(kx - 1, 0);
                        // ky_dn = std::max(ky - 1, 0);
                        kx_dn = ( kx > 0 ? kx - 1 : 0 );
                        ky_dn = ( ky > 0 ? ky - 1 : 0 );

                        return ( _eval_node_ddfdzdz( z, kx_up, ky_up, kz )
                                 - _eval_node_ddfdzdz( z, kx_up, ky_dn, kz )
                                 - _eval_node_ddfdzdz( z, kx_dn, ky_up, kz )
                                 + _eval_node_ddfdzdz( z, kx_dn, ky_dn, kz ) )
                             / ( ( _x_vals[ kx_up ] - _x_vals[ kx_dn ] )
                                 * ( _y_vals[ ky_up ] - _y_vals[ ky_dn ] ) );
                    }

                    size_t _length_x;              // Number of x nodes
                    size_t _length_y;              // Number of y nodes
                    size_t _length_z;              // Number of z nodes
                    size_t _accel[ 3 ];            // Indices of x,y,z point

                    INDEPTYPE *_x_vals = nullptr;  // 1D array of x values
                    INDEPTYPE *_y_vals = nullptr;  // 1D array of y values
                    INDEPTYPE *_z_vals = nullptr;  // 1D array of z values

                    DEPTYPE ***_f_vals
                        = nullptr;  // 3D array of f(x,y,z) values,
                                    // f_vals[i][j][k] = f(x[i], y[j],
                                    // z[k]) f_vals[i][j] = f(x[i], y[j])
                    DEPTYPE ***_f_slopes
                        = nullptr;  // Slopes used to generate natural
                                    // cubic spline solution at each x and
                                    // y = constant node, describes f(x[i],
                                    // y[i], z)
                    DEPTYPE ***_dfdx_slopes
                        = nullptr;  // Slopes used to generate df/dx at
                                    // each node point, describes df/dx @
                                    // x[i], y[j], z
                    DEPTYPE ***_dfdy_slopes
                        = nullptr;  // Slopes used to generate df/dy at
                                    // each node point, describes df/dy @
                                    // x[i], y[j], z


                    // cache variables
            };

            DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(
                hybrid_spline_3d, NCPA::interpolation::_spline_3d );

        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA
