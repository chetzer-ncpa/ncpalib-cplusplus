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
            // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( linear_spline_1d,
            //                                        _lanl_spline_1d );

            _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
                class linear_spline_1d<INDEPTYPE, DEPTYPE, void,
                                       ENABLE_IF_REAL( INDEPTYPE ),
                                       ENABLE_IF_REAL( DEPTYPE )>
                : public _lanl_spline_1d<INDEPTYPE, DEPTYPE> {
                public:
                    virtual ~linear_spline_1d() {}

                    virtual void ready() override {
                        DEPTYPE *slopes   = this->_get_slopes();
                        INDEPTYPE *x_vals = this->_get_x_vals();
                        DEPTYPE *f_vals   = this->_get_f_vals();
                        for ( int i = 0; i < this->_get_length() - 1; i++ ) {
                            slopes[ i ] = ( f_vals[ i + 1 ] - f_vals[ i ] )
                                        / ( x_vals[ i + 1 ] - x_vals[ i ] );
                        }
                        slopes[ this->_get_length() - 1 ]
                            = slopes[ this->_get_length() - 2 ];
                    }

                    virtual DEPTYPE eval_f( INDEPTYPE x ) override {
                        DEPTYPE *f_vals   = this->_get_f_vals();
                        INDEPTYPE *x_vals = this->_get_x_vals();
                        DEPTYPE *slopes   = this->_get_slopes();
                        size_t k
                            = find_segment( x, x_vals, this->_get_length(),
                                            this->_get_accel() );
                        return f_vals[ k ] + ( x - x_vals[ k ] ) * slopes[ k ];
                    }

                    virtual DEPTYPE eval_df( INDEPTYPE x ) override {
                        DEPTYPE *f_vals   = this->_get_f_vals();
                        INDEPTYPE *x_vals = this->_get_x_vals();
                        DEPTYPE *slopes   = this->_get_slopes();
                        size_t k
                            = find_segment( x, x_vals, this->_get_length(),
                                            this->_get_accel() );
                        return slopes[ k ];
                    }

                    virtual DEPTYPE eval_ddf( INDEPTYPE x ) override {
                        return 0.0;
                    }

                    virtual DEPTYPE eval_dddf( INDEPTYPE x ) override {
                        return 0.0;
                    }
            };

            DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR( linear_spline_1d,
                                                    _lanl_spline_1d );
        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA
