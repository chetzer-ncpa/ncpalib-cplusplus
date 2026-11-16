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
#include "NCPA/interpolation/lanl/declarations.hpp"
#include "NCPA/interpolation/lanl/lanl_common.hpp"
#include "NCPA/interpolation/lanl/lanl_spline_1d.hpp"
// #include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

// template<typename T, typename U>
// static void swap(
//     NCPA::interpolation::LANL::NaturalCubicSpline1D<T, U>& a,
//     NCPA::interpolation::LANL::NaturalCubicSpline1D<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        namespace LANL {
            // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( NaturalCubicSpline1D,
            //                                        _lanl_spline_1d );

            template<typename FloatType>
            class NaturalCubicSpline1D
                : public Cloneable<NaturalCubicSpline1D<FloatType>,
                                   _lanl_spline_1d<FloatType>> {
                    typedef Cloneable<NaturalCubicSpline1D<FloatType>,
                                      _lanl_spline_1d<FloatType>>
                        parent_t;

                public:
                using parent_t::eval_f;
                using parent_t::eval_df;
                using parent_t::eval_ddf;
                using parent_t::eval_dddf;
                using parent_t::init;
                    NaturalCubicSpline1D()          = default;
                    virtual ~NaturalCubicSpline1D() = default;
                    NaturalCubicSpline1D(
                        const NaturalCubicSpline1D<FloatType>& other )
                        = default;

                    NaturalCubicSpline1D(
                        NaturalCubicSpline1D<FloatType>&& source ) noexcept
                        :
                        NaturalCubicSpline1D<FloatType>() {
                        swap( *this, source );
                    }

                    friend void swap(
                        NaturalCubicSpline1D<FloatType>& a,
                        NaturalCubicSpline1D<FloatType>& b ) noexcept {
                        swap( dynamic_cast<parent_t&>( a ),
                              dynamic_cast<parent_t&>( b ) );
                    }

                    NaturalCubicSpline1D<FloatType>& operator=(
                        NaturalCubicSpline1D<FloatType> other ) {
                        swap( *this, other );
                        return *this;
                    }

                    virtual void ready() override {
                        FloatType *_f_vals = this->_get_f_vals();
                        FloatType *_x_vals = this->_get_x_vals();
                        FloatType *_slopes = this->_get_slopes();

                        FloatType ai, bi, ci;
                        FloatType di;
                        FloatType new_c[ this->_get_length() - 1 ];
                        FloatType new_d[ this->_get_length() ];

                        bi = 2.0 / ( _x_vals[ 1 ] - _x_vals[ 0 ] );
                        ci = 1.0 / ( _x_vals[ 1 ] - _x_vals[ 0 ] );
                        di = 3.0 * ( _f_vals[ 1 ] - _f_vals[ 0 ] )
                           / std::pow( _x_vals[ 1 ] - _x_vals[ 0 ], 2 );

                        new_c[ 0 ] = ci / bi;
                        new_d[ 0 ] = di / bi;

                        for (size_t i = 1; i < this->_get_length() - 1; i++) {
                            FloatType dx1 = _x_vals[ i ] - _x_vals[ i - 1 ],
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
                        for (int i = this->_get_length() - 2; i > -1; i--) {
                            _slopes[ i ]
                                = new_d[ i ] - new_c[ i ] * _slopes[ i + 1 ];
                        }
                        _ready = true;
                    }

                    FloatType eval_f(
                        const std::vector<FloatType>& x ) override {
                        this->_check_dimensions( x );
                        _make_ready();
                        FloatType *_f_vals = this->_get_f_vals();
                        FloatType *_x_vals = this->_get_x_vals();
                        FloatType *_slopes = this->_get_slopes();
                        size_t k;
                        FloatType dx, X;
                        FloatType df, A, B;

                        k = find_segment( x.front(), _x_vals,
                                          this->_get_length(),
                                          this->_get_accel() );

                        dx = _x_vals[ k + 1 ] - _x_vals[ k ];
                        X  = ( x.front() - _x_vals[ k ] ) / dx;

                        df = _f_vals[ k + 1 ] - _f_vals[ k ];
                        A  = _slopes[ k ] * dx - df;
                        B  = -_slopes[ k + 1 ] * dx + df;

                        return _f_vals[ k ]
                             + X
                                   * ( df
                                       + ( 1.0 - X )
                                             * ( A * ( 1.0 - X ) + B * X ) );
                    }

                    FloatType eval_df( const std::vector<FloatType>& x,
                                       size_t wrt ) override {
                        this->_check_dimensions( x );
                        this->_check_wrt( wrt );
                        _make_ready();
                        FloatType *_f_vals = this->_get_f_vals();
                        FloatType *_x_vals = this->_get_x_vals();
                        FloatType *_slopes = this->_get_slopes();
                        size_t k;
                        FloatType dx, X;
                        FloatType df, A, B;

                        k = find_segment( x.front(), _x_vals,
                                          this->_get_length(),
                                          this->_get_accel() );

                        dx = _x_vals[ k + 1 ] - _x_vals[ k ];
                        X  = ( x.front() - _x_vals[ k ] ) / dx;

                        df = _f_vals[ k + 1 ] - _f_vals[ k ];
                        A  = _slopes[ k ] * dx - df;
                        B  = -_slopes[ k + 1 ] * dx + df;

                        return df / dx
                             + ( 1.0 - 2.0 * X ) * ( A * ( 1.0 - X ) + B * X )
                                   / dx
                             + X * ( 1.0 - X ) * ( B - A ) / dx;
                    }

                    FloatType eval_ddf( const std::vector<double>& x,
                                        size_t wrt1, size_t wrt2 ) override {
                        this->_check_dimensions( x );
                        this->_check_wrt( wrt1 );
                        this->_check_wrt( wrt2 );
                        _make_ready();
                        FloatType *_f_vals = this->_get_f_vals();
                        FloatType *_x_vals = this->_get_x_vals();
                        FloatType *_slopes = this->_get_slopes();
                        size_t k;
                        FloatType dx, X;
                        FloatType df, A, B;

                        k = find_segment( x.front(), _x_vals,
                                          this->_get_length(),
                                          this->_get_accel() );

                        dx = _x_vals[ k + 1 ] - _x_vals[ k ];
                        X  = ( x.front() - _x_vals[ k ] ) / dx;

                        df = _f_vals[ k + 1 ] - _f_vals[ k ];
                        A  = _slopes[ k ] * dx - df;
                        B  = -_slopes[ k + 1 ] * dx + df;

                        return 2.0 * ( B - 2.0 * A + ( A - B ) * 3.0 * X )
                             / std::pow( dx, 2 );
                    }

                    FloatType eval_dddf( const std::vector<double>& x,
                                         size_t wrt1, size_t wrt2,
                                         size_t wrt3 ) override {
                        this->_check_dimensions( x );
                        this->_check_wrt( wrt1 );
                        this->_check_wrt( wrt2 );
                        this->_check_wrt( wrt3 );
                        _make_ready();
                        FloatType *_f_vals = this->_get_f_vals();
                        FloatType *_x_vals = this->_get_x_vals();
                        FloatType *_slopes = this->_get_slopes();
                        size_t k;
                        FloatType dx, X;
                        FloatType df, A, B;

                        k = find_segment( x.front(), _x_vals,
                                          this->_get_length(),
                                          this->_get_accel() );

                        dx = _x_vals[ k + 1 ] - _x_vals[ k ];
                        X  = ( x.front() - _x_vals[ k ] ) / dx;

                        df = _f_vals[ k + 1 ] - _f_vals[ k ];
                        A  = _slopes[ k ] * dx - df;
                        B  = -_slopes[ k + 1 ] * dx + df;

                        return 6.0 * ( A - B ) / std::pow( dx, 3 );
                    }

                    virtual interpolator_type_t interptype()
                        const override {
                        return interpolator_type_t::LANL_CUBIC;
                    }

                protected:
                    bool _ready = false;

                    void _make_ready() {
                        if (!_ready) {
                            this->ready();
                        }
                    }
            };

            // DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR( NaturalCubicSpline1D,
            //                                         _lanl_spline_1d );


        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA

// template<typename T, typename U>
// static void swap(
//     NCPA::interpolation::LANL::NaturalCubicSpline1D<T, U>& a,
//     NCPA::interpolation::LANL::NaturalCubicSpline1D<T, U>& b ) noexcept {
//     ::swap(
//         dynamic_cast<
//             NCPA::interpolation::LANL::NaturalCubicSpline1D<T, U>&>( a ),
//         dynamic_cast<
//             NCPA::interpolation::LANL::NaturalCubicSpline1D<T, U>&>( b )
//             );
// }
