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
#include "NCPA/interpolation/abstract_spline_1d.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/lanl/lanl_common.hpp"
#include "NCPA/interpolation/lanl/lanl_spline_1d.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

namespace NCPA {
    namespace interpolation {
        namespace LANL {
            // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( natural_cubic_spline_1d,
            //                                        _lanl_spline_1d );

            _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
                class natural_cubic_spline_1d<INDEPTYPE, DEPTYPE, void,
                                              ENABLE_IF_REAL( INDEPTYPE ),
                                              ENABLE_IF_REAL( DEPTYPE )>
                : public _lanl_spline_1d<INDEPTYPE, DEPTYPE> {
                public:
                    virtual ~natural_cubic_spline_1d() { this->clear(); }

                    virtual void ready() override {
                        DEPTYPE *_f_vals   = this->_get_f_vals();
                        INDEPTYPE *_x_vals = this->_get_x_vals();
                        DEPTYPE *_slopes   = this->_get_slopes();

                        INDEPTYPE ai, bi, ci;
                        DEPTYPE di;
                        INDEPTYPE new_c[ this->_get_length() - 1 ];
                        DEPTYPE new_d[ this->_get_length() ];

                        bi = 2.0 / ( _x_vals[ 1 ] - _x_vals[ 0 ] );
                        ci = 1.0 / ( _x_vals[ 1 ] - _x_vals[ 0 ] );
                        di = 3.0 * ( _f_vals[ 1 ] - _f_vals[ 0 ] )
                           / std::pow( _x_vals[ 1 ] - _x_vals[ 0 ], 2 );

                        new_c[ 0 ] = ci / bi;
                        new_d[ 0 ] = di / bi;

                        for ( size_t i = 1; i < this->_get_length() - 1;
                              i++ ) {
                            INDEPTYPE dx1 = _x_vals[ i ] - _x_vals[ i - 1 ],
                                      dx2 = _x_vals[ i + 1 ] - _x_vals[ i ];

                            ai = 1.0 / dx1;
                            bi = 2.0 * ( 1.0 / dx1 + 1.0 / dx2 );
                            ci = 1.0 / dx2;
                            di = 3.0
                               * ( ( _f_vals[ i ] - _f_vals[ i - 1 ] )
                                       / std::pow( dx1, 2 )
                                   + ( _f_vals[ i + 1 ] - _f_vals[ i ] )
                                         / std::pow( dx2, 2 ) );

                            new_c[ i ] = ci / ( bi - new_c[ i - 1 ] * ai );
                            new_d[ i ] = ( di - new_d[ i - 1 ] * ai )
                                       / ( bi - new_c[ i - 1 ] * ai );
                        }

                        ai = 1.0
                           / ( _x_vals[ this->_get_length() - 1 ]
                               - _x_vals[ this->_get_length() - 2 ] );
                        bi = 2.0
                           / ( _x_vals[ this->_get_length() - 1 ]
                               - _x_vals[ this->_get_length() - 2 ] );
                        di = 3.0
                           * ( _f_vals[ this->_get_length() - 1 ]
                               - _f_vals[ this->_get_length() - 2 ] )
                           / std::pow(
                                 _x_vals[ this->_get_length() - 1 ]
                                     - _x_vals[ this->_get_length() - 2 ],
                                 2 );

                        new_d[ this->_get_length() - 1 ]
                            = ( di - new_d[ this->_get_length() - 2 ] * ai )
                            / ( bi - new_c[ this->_get_length() - 2 ] * ai );

                        _slopes[ this->_get_length() - 1 ]
                            = new_d[ this->_get_length() - 1 ];
                        for ( int i = this->_get_length() - 2; i > -1; i-- ) {
                            _slopes[ i ]
                                = new_d[ i ] - new_c[ i ] * _slopes[ i + 1 ];
                        }
                    }

                    DEPTYPE eval_f( INDEPTYPE x ) override {
                        DEPTYPE *_f_vals   = this->_get_f_vals();
                        INDEPTYPE *_x_vals = this->_get_x_vals();
                        DEPTYPE *_slopes   = this->_get_slopes();
                        size_t k;
                        INDEPTYPE dx, X;
                        DEPTYPE df, A, B;

                        k = find_segment( x, _x_vals, this->_get_length(),
                                          this->_get_accel() );

                        dx = _x_vals[ k + 1 ] - _x_vals[ k ];
                        X  = ( x - _x_vals[ k ] ) / dx;

                        df = _f_vals[ k + 1 ] - _f_vals[ k ];
                        A  = _slopes[ k ] * dx - df;
                        B  = -_slopes[ k + 1 ] * dx + df;

                        return _f_vals[ k ]
                             + X
                                   * ( df
                                       + ( 1.0 - X )
                                             * ( A * ( 1.0 - X ) + B * X ) );
                    }

                    DEPTYPE eval_df( INDEPTYPE x ) override {
                        DEPTYPE *_f_vals   = this->_get_f_vals();
                        INDEPTYPE *_x_vals = this->_get_x_vals();
                        DEPTYPE *_slopes   = this->_get_slopes();
                        size_t k;
                        INDEPTYPE dx, X;
                        DEPTYPE df, A, B;

                        k = find_segment( x, _x_vals, this->_get_length(),
                                          this->_get_accel() );

                        dx = _x_vals[ k + 1 ] - _x_vals[ k ];
                        X  = ( x - _x_vals[ k ] ) / dx;

                        df = _f_vals[ k + 1 ] - _f_vals[ k ];
                        A  = _slopes[ k ] * dx - df;
                        B  = -_slopes[ k + 1 ] * dx + df;

                        return df / dx
                             + ( 1.0 - 2.0 * X ) * ( A * ( 1.0 - X ) + B * X )
                                   / dx
                             + X * ( 1.0 - X ) * ( B - A ) / dx;
                    }

                    DEPTYPE eval_ddf( INDEPTYPE x ) override {
                        DEPTYPE *_f_vals   = this->_get_f_vals();
                        INDEPTYPE *_x_vals = this->_get_x_vals();
                        DEPTYPE *_slopes   = this->_get_slopes();
                        size_t k;
                        INDEPTYPE dx, X;
                        DEPTYPE df, A, B;

                        k = find_segment( x, _x_vals, this->_get_length(),
                                          this->_get_accel() );

                        dx = _x_vals[ k + 1 ] - _x_vals[ k ];
                        X  = ( x - _x_vals[ k ] ) / dx;

                        df = _f_vals[ k + 1 ] - _f_vals[ k ];
                        A  = _slopes[ k ] * dx - df;
                        B  = -_slopes[ k + 1 ] * dx + df;

                        return 2.0 * ( B - 2.0 * A + ( A - B ) * 3.0 * X )
                             / std::pow( dx, 2 );
                    }

                    DEPTYPE eval_dddf( INDEPTYPE x ) override {
                        DEPTYPE *_f_vals   = this->_get_f_vals();
                        INDEPTYPE *_x_vals = this->_get_x_vals();
                        DEPTYPE *_slopes   = this->_get_slopes();
                        size_t k;
                        INDEPTYPE dx, X;
                        DEPTYPE df, A, B;

                        k = find_segment( x, _x_vals, this->_get_length(),
                                          this->_get_accel() );

                        dx = _x_vals[ k + 1 ] - _x_vals[ k ];
                        X  = ( x - _x_vals[ k ] ) / dx;

                        df = _f_vals[ k + 1 ] - _f_vals[ k ];
                        A  = _slopes[ k ] * dx - df;
                        B  = -_slopes[ k + 1 ] * dx + df;

                        return 6.0 * ( A - B ) / std::pow( dx, 3 );
                    }
            };

            DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR( natural_cubic_spline_1d,
                                                    _lanl_spline_1d );


        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA
