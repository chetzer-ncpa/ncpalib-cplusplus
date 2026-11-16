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
#include "NCPA/interpolation/abstract_spline.hpp"
#include "NCPA/interpolation/declarations.hpp"
#include "NCPA/interpolation/lanl/lanl_common.hpp"
#include "NCPA/interpolation/lanl/lanl_spline_1d.hpp"
// #include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

namespace NCPA {
    namespace interpolation {
        namespace LANL {
            template<typename FloatType>
            class LinearSpline1D
                : public Cloneable<LinearSpline1D<FloatType>,
                                   _lanl_spline_1d<FloatType>> {
                    typedef Cloneable<LinearSpline1D<FloatType>,
                                      _lanl_spline_1d<FloatType>>
                        parent_t;

                public:
                    using parent_t::eval_dddf;
                    using parent_t::eval_ddf;
                    using parent_t::eval_df;
                    using parent_t::eval_f;
                    using parent_t::init;

                    LinearSpline1D() : parent_t() {}

                    virtual ~LinearSpline1D() {}

                    LinearSpline1D(
                        const LinearSpline1D<FloatType>& other ) = default;

                    LinearSpline1D(
                        LinearSpline1D<FloatType>&& source ) noexcept
                        = default;

                    friend void swap(
                        LinearSpline1D<FloatType>& a,
                        LinearSpline1D<FloatType>& b ) noexcept {
                        swap( dynamic_cast<parent_t&>( a ),
                              dynamic_cast<parent_t&>( b ) );
                    }

                    LinearSpline1D<FloatType>& operator=(
                        LinearSpline1D<FloatType> other ) {
                        swap( *this, other );
                        return *this;
                    }

                    void ready() override {
                        FloatType *slopes = this->_get_slopes();
                        FloatType *x_vals = this->_get_x_vals();
                        FloatType *f_vals = this->_get_f_vals();
                        for (int i = 0; i < this->_get_length() - 1; i++) {
                            slopes[ i ] = ( f_vals[ i + 1 ] - f_vals[ i ] )
                                        / ( x_vals[ i + 1 ] - x_vals[ i ] );
                        }
                        slopes[ this->_get_length() - 1 ]
                            = slopes[ this->_get_length() - 2 ];
                        _ready = true;
                    }

                    FloatType eval_f(
                        const std::vector<FloatType>& x ) override {
                        this->_check_dimensions( x );
                        _make_ready();
                        FloatType *f_vals = this->_get_f_vals();
                        FloatType *x_vals = this->_get_x_vals();
                        FloatType *slopes = this->_get_slopes();
                        size_t k          = find_segment( x.front(), x_vals,
                                                          this->_get_length(),
                                                          this->_get_accel() );
                        return f_vals[ k ]
                             + ( x.front() - x_vals[ k ] ) * slopes[ k ];
                    }

                    FloatType eval_df( const std::vector<FloatType>& x,
                                       size_t wrt ) override {
                        this->_check_dimensions( x );
                        this->_check_wrt( wrt );
                        _make_ready();
                        FloatType *f_vals = this->_get_f_vals();
                        FloatType *x_vals = this->_get_x_vals();
                        FloatType *slopes = this->_get_slopes();
                        size_t k          = find_segment( x.front(), x_vals,
                                                          this->_get_length(),
                                                          this->_get_accel() );
                        return slopes[ k ];
                    }

                    FloatType eval_ddf( const std::vector<FloatType>& x,
                                        size_t wrt1, size_t wrt2 ) override {
                        return 0.0;
                    }

                    FloatType eval_dddf( const std::vector<FloatType>& x,
                                         size_t wrt1, size_t wrt2,
                                         size_t wrt3 ) override {
                        return 0.0;
                    }

                    virtual interpolator_type_t interptype() const override {
                        return interpolator_type_t::LANL_LINEAR;
                    }

                protected:
                    bool _ready = false;

                    void _make_ready() {
                        if (!_ready) {
                            this->ready();
                        }
                    }
            };

            // DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR( LinearSpline1D,
            //                                         _lanl_spline_1d );
        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA

// template<typename T, typename U>
// static void swap(
//     NCPA::interpolation::LANL::LinearSpline1D<T, U>& a,
//     NCPA::interpolation::LANL::LinearSpline1D<T, U>& b ) noexcept {
//     ::swap(
//         dynamic_cast<NCPA::interpolation::LANL::_lanl_spline_1d<T, U>&>( a
//         ), dynamic_cast<NCPA::interpolation::LANL::_lanl_spline_1d<T, U>&>(
//         b ) );
// }
