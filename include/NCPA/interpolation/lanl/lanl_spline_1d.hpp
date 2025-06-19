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
#include "NCPA/interpolation/lanl/lanl_declarations.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <cstring>

template<typename T, typename U>
static void swap(
    NCPA::interpolation::LANL::_lanl_spline_1d<T, U>& a,
    NCPA::interpolation::LANL::_lanl_spline_1d<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        namespace LANL {

            /**
             * 1D INTERPOLATORS
             */
            // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
            //     _lanl_spline_1d, NCPA::interpolation::_spline_1d );

            template<typename INDEPTYPE, typename DEPTYPE>
            class _lanl_spline_1d<INDEPTYPE, DEPTYPE, void,
                                  ENABLE_IF_REAL( INDEPTYPE ),
                                  ENABLE_IF_REAL( DEPTYPE )>
                : public NCPA::interpolation::_spline_1d<INDEPTYPE, DEPTYPE> {
                public:
                    _lanl_spline_1d() :
                        _length { 0 },
                        _accel { 0 },
                        _x_vals { nullptr },
                        _f_vals { nullptr },
                        _slopes { nullptr } {}

                    virtual ~_lanl_spline_1d() { this->clear(); }

                    _lanl_spline_1d(
                        const _lanl_spline_1d<INDEPTYPE, DEPTYPE>& other ) :
                        _lanl_spline_1d() {
                        // _length { 0 },
                        // _accel { 0 },
                        // _x_vals { nullptr },
                        // _f_vals { nullptr },
                        // _slopes { nullptr } {
                        _length = other._length;
                        _accel  = other._accel;
                        if ( other._x_vals != nullptr ) {
                            std::memcpy( _x_vals, other._x_vals,
                                         _length * sizeof( INDEPTYPE ) );
                        }
                        if ( other._f_vals != nullptr ) {
                            std::memcpy( _f_vals, other._f_vals,
                                         _length * sizeof( DEPTYPE ) );
                        }
                        if ( other._slopes != nullptr ) {
                            std::memcpy( _slopes, other._slopes,
                                         _length * sizeof( DEPTYPE ) );
                        }
                    }

                    _lanl_spline_1d(
                        _lanl_spline_1d<INDEPTYPE, DEPTYPE>&& source ) noexcept
                        :
                        _lanl_spline_1d<INDEPTYPE, DEPTYPE>() {
                        ::swap( *this, source );
                    }

                    friend void ::swap<INDEPTYPE, DEPTYPE>(
                        _lanl_spline_1d<INDEPTYPE, DEPTYPE>& a,
                        _lanl_spline_1d<INDEPTYPE, DEPTYPE>& b ) noexcept;

                    _lanl_spline_1d<INDEPTYPE, DEPTYPE>& operator=(
                        _lanl_spline_1d<INDEPTYPE, DEPTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual void fill( size_t N, const INDEPTYPE *x,
                                       const DEPTYPE *f ) override {
                        if ( !_initialized() || _length != N ) {
                            init( N );
                        }
                        std::memcpy( _x_vals, x, N * sizeof( INDEPTYPE ) );
                        std::memcpy( _f_vals, f, N * sizeof( DEPTYPE ) );
                    }

                    virtual void fill(
                        const std::vector<INDEPTYPE>& x,
                        const std::vector<DEPTYPE>& f ) override {
                        size_t N = x.size();
                        if ( N != f.size() ) {
                            throw std::invalid_argument(
                                "Vectors must be same size" );
                        }
                        fill( N, &x[ 0 ], &f[ 0 ] );
                    }

                    virtual void init( size_t N ) override {
                        clear();
                        _length = N;
                        _accel  = 0;
                        _x_vals = NCPA::arrays::zeros<DEPTYPE>( N );
                        _f_vals = NCPA::arrays::zeros<INDEPTYPE>( N );
                        _slopes = NCPA::arrays::zeros<INDEPTYPE>( N );
                    }

                    virtual void clear() override {
                        if ( _x_vals != nullptr ) {
                            delete[] _x_vals;
                            _x_vals = nullptr;
                        }
                        if ( _f_vals != nullptr ) {
                            delete[] _f_vals;
                            _f_vals = nullptr;
                        }
                        if ( _slopes != nullptr ) {
                            delete[] _slopes;
                            _slopes = nullptr;
                        }
                    }

                    virtual std::pair<INDEPTYPE, INDEPTYPE> limits()
                        const override {
                        return std::pair<INDEPTYPE, INDEPTYPE>(
                            { this->_get_x_vals()[ 0 ],
                              this->_get_x_vals()[ _get_length() - 1 ] } );
                    }

                    virtual const std::vector<INDEPTYPE> source_x()
                        const override {
                        return std::vector<INDEPTYPE>( _x_vals,
                                                       _x_vals + _length );
                    }

                    virtual const std::vector<DEPTYPE> source_f()
                        const override {
                        return std::vector<DEPTYPE>( _f_vals,
                                                     _f_vals + _length );
                    }

                protected:
                    size_t& _get_length() { return _length; }

                    const size_t& _get_length() const { return _length; }

                    size_t& _get_accel() { return _accel; }

                    const size_t& _get_accel() const { return _accel; }

                    INDEPTYPE *_get_x_vals() { return _x_vals; }

                    const INDEPTYPE *_get_x_vals() const { return _x_vals; }

                    DEPTYPE *_get_f_vals() { return _f_vals; }

                    const DEPTYPE *_get_f_vals() const { return _f_vals; }

                    DEPTYPE *_get_slopes() { return _slopes; }

                    const DEPTYPE *_get_slopes() const { return _slopes; }

                    bool _initialized() const {
                        return ( _length > 0 && _x_vals != nullptr
                                 && _f_vals != nullptr && _slopes != nullptr );
                    }

                    size_t _length;  // Number of base points
                    size_t _accel;   // Index of previous table look up;
                                     // used to increase spline speed
                    INDEPTYPE *_x_vals = nullptr;  // 1D array of x values
                    DEPTYPE *_f_vals   = nullptr;  // 1D array of f(x) values,
                                                   // f_vals[i] = f(x[i])
                    DEPTYPE *_slopes
                        = nullptr;                 // Slopes used to generate
                                    // natural cubic spline solution
            };

            DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(
                _lanl_spline_1d, NCPA::interpolation::_spline_1d );

        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap(
    NCPA::interpolation::LANL::_lanl_spline_1d<T, U>& a,
    NCPA::interpolation::LANL::_lanl_spline_1d<T, U>& b ) noexcept {
    using std::swap;
    swap( a._length, b._length );
    swap( a._accel, b._accel );
    swap( a._x_vals, b._x_vals );
    swap( a._f_vals, b._f_vals );
    swap( a._slopes, b._slopes );
}
