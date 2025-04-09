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
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <array>

namespace NCPA {
    namespace interpolation {
        namespace LANL {

            //----------------------------------------//
            //------------Common Functions------------//
            //--------Used by All Interpolators-------//
            //----------------------------------------//

            // Find k such that x[k] <= x <= x[k+1]
            template<typename T>
            size_t find_segment( T x, T *x_vals, size_t length,
                                 size_t& prev ) {
                bool done = false;

                if ( x > x_vals[ length - 1 ] || x < x_vals[ 0 ] ) {
                    std::ostringstream oss;

                    oss << "Cannot interpolate outside of given bounds.  "
                           "x = "
                        << x << " is outside of bounds (" << x_vals[ 0 ]
                        << ", " << x_vals[ length - 1 ] << ")." << std::endl;
                    throw std::range_error( oss.str() );
                }
                if ( x >= x_vals[ prev ] && x <= x_vals[ prev + 1 ] ) {
                    done = true;
                }
                if ( !done && prev + 2 <= length - 1 ) {
                    if ( x >= x_vals[ prev + 1 ] && x <= x_vals[ prev + 2 ] ) {
                        prev++;
                        done = true;
                    }
                }
                if ( !done && prev - 1 >= 0 ) {
                    if ( x >= x_vals[ prev - 1 ] && x <= x_vals[ prev ] ) {
                        prev--;
                        done = true;
                    }
                }
                if ( !done ) {
                    for ( size_t n = 0; n < length; n++ ) {
                        if ( x >= x_vals[ n ] && x <= x_vals[ n + 1 ] ) {
                            prev = n;
                            break;
                        } else if ( x >= x_vals[ ( length - 1 ) - ( n + 1 ) ]
                                    && x < x_vals[ ( length - 1 ) - n ] ) {
                            prev = ( length - 1 ) - ( n + 1 );
                            break;
                        }
                    }
                }

                return prev;
            }

            // Return x within x_min <= x <= x_max
            template<typename T>
            T in_interval( T x, T x_min, T x_max ) {
                return std::min( std::max( x, x_min ), x_max );
            }

            // ensure type-safe initialization of the bicubic conversion
            // matrix.  This might well be overkill but it should be handled at
            // compile time so no harm no foul
            template<typename T>
            using bicubic_conversion_matrix_t = T[ 16 ][ 16 ];

            template<typename T>
            class constants {
                public:
                    static constexpr bicubic_conversion_matrix_t<T> bicubic_conversion_matrix = {
                        {  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,0.0,
                         0.0,  0.0,  0.0,  0.0,  0.0,  0.0                         },
                        {  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                         0.0,  0.0,  0.0,  0.0,  0.0,  0.0                         },
                        { -3.0,  3.0,  0.0,  0.0, -2.0, -1.0,  0.0,  0.0,  0.0,  0.0,
                         0.0,  0.0,  0.0,  0.0,  0.0,  0.0                         },
                        {  2.0, -2.0,  0.0,  0.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,
                         0.0,  0.0,  0.0,  0.0,  0.0,  0.0                         },
                        {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,
                         0.0,  0.0,  0.0,  0.0,  0.0,  0.0                         },
                        {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                         0.0,  0.0,  1.0,  0.0,  0.0,  0.0                         },
                        {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -3.0,  3.0,
                         0.0,  0.0, -2.0, -1.0,  0.0,  0.0                         },
                        {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  2.0, -2.0,
                         0.0,  0.0,  1.0,  1.0,  0.0,  0.0                         },
                        { -3.0,  0.0,  3.0,  0.0,  0.0,  0.0,  0.0,  0.0, -2.0,  0.0,
                         -1.0,  0.0,  0.0,  0.0,  0.0,  0.0                        },
                        {  0.0,  0.0,  0.0,  0.0, -3.0,  0.0,  3.0,  0.0,  0.0,  0.0,
                         0.0,  0.0, -2.0,  0.0, -1.0,  0.0                         },
                        {  9.0, -9.0, -9.0,  9.0,  6.0,  3.0, -6.0, -3.0,  6.0,
                         -6.0,                    3.0, -3.0,  4.0,  2.0,  2.0,  1.0 },
                        { -6.0,  6.0,  6.0, -6.0, -3.0, -3.0,  3.0,  3.0, -4.0,
                         4.0,                   -2.0,  2.0, -2.0, -2.0, -1.0, -1.0 },
                        {  2.0,  0.0, -2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,
                         1.0,  0.0,  0.0,  0.0,  0.0,  0.0                         },
                        {  0.0,  0.0,  0.0,  0.0,  2.0,  0.0, -2.0,  0.0,  0.0,  0.0,
                         0.0,  0.0,  1.0,  0.0,  1.0,  0.0                         },
                        { -6.0,  6.0,  6.0, -6.0, -4.0, -2.0,  4.0,  2.0, -3.0,
                         3.0,                   -3.0,  3.0, -2.0, -1.0, -2.0, -1.0 },
                        {  4.0, -4.0, -4.0,  4.0,  2.0,  2.0, -2.0, -2.0,  2.0,
                         -2.0,                    2.0, -2.0,  1.0,  1.0,  1.0,  1.0 }
                    };

                    // static const bicubic_conversion_matrix_t<T>& bicubic_conversion_matrix() {
                    //     return constants<T>::_bicubic_conversion_matrix;
                    // }
            };


        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA

template<typename T>
constexpr NCPA::interpolation::LANL::bicubic_conversion_matrix_t<T> 
NCPA::interpolation::LANL::constants<T>::bicubic_conversion_matrix;