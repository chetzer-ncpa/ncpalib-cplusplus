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
// #include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <cstring>

namespace NCPA {
    namespace interpolation {
        namespace LANL {
            template<typename FloatType>
            class _lanl_spline_1d : public spline1d<FloatType> {
                    typedef spline1d<FloatType> parent_t;

                public:
                    using parent_t::eval_dddf;
                    using parent_t::eval_ddf;
                    using parent_t::eval_df;
                    using parent_t::eval_f;
                    using parent_t::init;

                    _lanl_spline_1d() :
                        // _length { 0 },
                        _accel { 0 },
                        // _x_vals {},
                        // _f_vals {},
                        _slopes {} {}

                    virtual ~_lanl_spline_1d() { this->clear(); }

                    _lanl_spline_1d( const _lanl_spline_1d<FloatType>& other )
                        = default;
                    //     :
                    //     _lanl_spline_1d<FloatType>() {
                    //     _length = other._length;
                    //     _accel  = other._accel;
                    //     if (other._x_vals != nullptr) {
                    //         std::memcpy( _x_vals, other._x_vals,
                    //                      _length * sizeof( FloatType ) );
                    //     }
                    //     if (other._f_vals != nullptr) {
                    //         std::memcpy( _f_vals, other._f_vals,
                    //                      _length * sizeof( FloatType ) );
                    //     }
                    //     if (other._slopes != nullptr) {
                    //         std::memcpy( _slopes, other._slopes,
                    //                      _length * sizeof( FloatType ) );
                    //     }
                    // }

                    _lanl_spline_1d(
                        _lanl_spline_1d<FloatType>&& source ) noexcept
                        = default;

                    friend void swap(
                        _lanl_spline_1d<FloatType>& a,
                        _lanl_spline_1d<FloatType>& b ) noexcept {
                        using std::swap;
                        // swap( a._length, b._length );
                        swap( a._accel, b._accel );
                        // swap( a._x_vals, b._x_vals );
                        // swap( a._f_vals, b._f_vals );
                        swap( a._slopes, b._slopes );
                    }

                    // _lanl_spline_1d<FloatType>& operator=(
                    //     _lanl_spline_1d<FloatType> other ) {
                    //     swap( *this, other );
                    //     return *this;
                    // }

                    void clear() override {
                        // this->axes(0).clear();
                        // this->f().clear();
                        _slopes.clear();
                        parent_t::clear();
                        // if (_x_vals != nullptr) {
                        //     delete[] _x_vals;
                        //     _x_vals = nullptr;
                        // }
                        // if (_f_vals != nullptr) {
                        //     delete[] _f_vals;
                        //     _f_vals = nullptr;
                        // }
                        // if (_slopes != nullptr) {
                        //     delete[] _slopes;
                        //     _slopes = nullptr;
                        // }
                    }

                    size_t dimensions() const override { return 1; }

                    // void fill( const NCPA::arrays::SmartVector<FloatType>&
                    // x,
                    //            const NCPA::arrays::SmartVector<FloatType>& f
                    //            )
                    //     override {
                    //     if (x.rank() != this->dimensions()
                    //         || f.rank() != this->dimensions()) {
                    //         throw DimensionsDoNotMatch(
                    //             "fill() expected vectors of rank 1" );
                    //     }
                    //     this->fill( x.flatten(), f.flatten() );
                    // }

                    // void fill( const std::vector<FloatType>& x,
                    //            const std::vector<FloatType>& f ) override {
                    //     if (x.size() != f.size()) {
                    //         throw InputSizesDoNotMatch(
                    //             "fill() expects input vectors of equal size"
                    //             );
                    //     }
                    //     _x_vals = x;
                    //     _f_vals = f;
                    //     _slopes.resize( x.size(), 0 );
                    //     // if (!_initialized() || _length != x.size()) {
                    //     //     init( x.size() );
                    //     // }
                    //     // std::memcpy( _x_vals, x.data(),
                    //     //              x.size() * sizeof( FloatType ) );
                    //     // std::memcpy( _f_vals, f.data(),
                    //     //              f.size() * sizeof( FloatType ) );
                    // }

                    // virtual void fill( size_t N, const FloatType *x,
                    //                    const FloatType *f ) override {
                    //     if (!_initialized() || _length != N) {
                    //         init( N );
                    //     }
                    //     std::memcpy( _x_vals, x, N * sizeof( FloatType ) );
                    //     std::memcpy( _f_vals, f, N * sizeof( FloatType ) );
                    // }

                    // virtual void fill(
                    //     const std::vector<FloatType>& x,
                    //     const std::vector<FloatType>& f ) override {
                    //     size_t N = x.size();
                    //     if (N != f.size()) {
                    //         throw std::invalid_argument(
                    //             "Vectors must be same size" );
                    //     }
                    //     fill( N, &x[ 0 ], &f[ 0 ] );
                    // }

                    virtual void init( std::vector<size_t> N ) override {
                        clear();
                        this->_check_dimensions( N );
                        // _length = N.front();
                        _accel = 0;
                        // this->axes(0).clear();
                        // this->f().clear();
                        _slopes.clear();
                        _slopes.resize(N.front());
                        parent_t::init(N);
                        // this->axes(0).resize( N.front(), 0 );
                        // this->f().resize( N.front(), 0 );
                        // _slopes.resize( N.front(), 0 );
                        // _x_vals = NCPA::arrays::zeros<FloatType>( N.front()
                        // ); _f_vals = NCPA::arrays::zeros<FloatType>(
                        // N.front() ); _slopes =
                        // NCPA::arrays::zeros<FloatType>( N.front() );
                    }

                    // std::vector<std::pair<FloatType, FloatType>> limits()
                    //     const override {
                    //     return std::vector<std::pair<FloatType, FloatType>>
                    //     {
                    //         this->axes(0).size() >= 2
                    //             ? std::make_pair( this->axes(0).front(),
                    //                               this->axes(0).back() )
                    //             : std::make_pair<FloatType, FloatType>( 0.0,
                    //                                                     0.0
                    //                                                     )
                    //         // this->_get_x_vals()[ 0 ],
                    //         // this->_get_x_vals()[ _get_length() - 1 ]

                    //     };
                    //     // (
                    //     //     { this->_get_x_vals()[ 0 ],
                    //     //       this->_get_x_vals()[ _get_length() - 1 ] }
                    //     );
                    // }

                    // virtual const std::vector<FloatType> source_x()
                    //     const override {
                    //     return std::vector<FloatType>( _x_vals,
                    //                                    _x_vals + _length );
                    // }

                    // virtual const std::vector<FloatType> source_f()
                    //     const override {
                    //     return std::vector<FloatType>( _f_vals,
                    //                                  _f_vals + _length );
                    // }

                protected:
                    size_t _get_length() {
                        return this->axes( 0 ).size();
                        // return _length;
                    }

                    // const size_t& _get_length() const {
                    //     return this->axes(0).size();
                    //     // return _length;
                    // }

                    int& _get_accel() { return _accel; }

                    const int& _get_accel() const { return _accel; }

                    FloatType *_get_x_vals() { return this->axes( 0 ).data(); }

                    const FloatType *_get_x_vals() const {
                        return this->axes( 0 ).data();
                    }

                    FloatType *_get_f_vals() { return this->f().data(); }

                    const FloatType *_get_f_vals() const {
                        return this->f().data();
                    }

                    FloatType *_get_slopes() { return _slopes.data(); }

                    const FloatType *_get_slopes() const {
                        return _slopes.data();
                    }

                    bool _initialized() const {
                        return !( this->axes( 0 ).empty() || this->f().empty()
                                  || _slopes.empty() );
                        // return ( _length > 0 && _x_vals != nullptr
                        //          && _f_vals != nullptr && _slopes != nullptr
                        //          );
                    }

                    // size_t _length;  // Number of base points
                    int _accel;  // Index of previous table look up;
                                 // used to increase spline speed

                    // std::vector<FloatType> _x_vals;
                    // std::vector<FloatType> _f_vals;
                    std::vector<FloatType> _slopes;


                    // FloatType *_x_vals = nullptr;  // 1D array of x values
                    // FloatType *_f_vals = nullptr;  // 1D array of f(x)
                    // values,
                    //                                // f_vals[i] = f(x[i])
                    // FloatType *_slopes
                    //     = nullptr;                 // Slopes used to
                    //     generate
                    //                 // natural cubic spline solution
            };

            // DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(
            //     _lanl_spline_1d, NCPA::interpolation::_spline_1d );

        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA
