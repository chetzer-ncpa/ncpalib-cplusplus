#pragma once

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

#include "NCPA/arrays.hpp"
#include "NCPA/types.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

#if __has_include( "gsl/gsl_spline.h" )
#  define HAVE_GSL_INTERPOLATION_LIBRARY
#  include "gsl/gsl_interp.h"
#  include "gsl/gsl_spline.h"
#  include "gsl/gsl_version.h"
#endif

// Convenience #defines
#define _ENABLE_IF_INDEP_IS_REAL \
    typename std::enable_if<std::is_floating_point<INDEPTYPE>::value>::type
#define _ENABLE_IF_DEP_IS_REAL \
    typename std::enable_if<std::is_floating_point<DEPTYPE>::value>::type
#define _ENABLE_IF_DEP_IS_COMPLEX \
    typename std::enable_if<NCPA::types::is_complex<DEPTYPE>::value>::type
#define _SUBSPLINE_PTR_T \
    NCPA::interpolation::_spline_1d<INDEPTYPE, typename DEPTYPE::value_type> *
#define _SUBSPLINE_T( _CLASSNAME_ ) \
    _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type>

#define DECLARE_GENERIC_INTERPOLATOR_TEMPLATE_WITH_PARAMS(            \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                      \
    template<typename INDEPTYPE, typename DEPTYPE,                    \
             typename PARAMTYPE = _PARAMTYPE_, typename = void,       \
             typename = void>                                         \
    class _CLASSNAME_ : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE> { \
            static_assert( false, "Invalid types for interpolator" ); \
    }
#define DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _CLASSNAME_,           \
                                               _SUPERCLASSNAME_ )     \
    template<typename INDEPTYPE, typename DEPTYPE, typename = void,   \
             typename = void, typename = void>                        \
    class _CLASSNAME_ : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE> { \
            static_assert( false, "Invalid types for interpolator" ); \
    }


#define _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAMS \
    template<typename INDEPTYPE, typename DEPTYPE, typename PARAMTYPE>
#define _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION \
    template<typename INDEPTYPE, typename DEPTYPE>

#define DECLARE_REAL_VERSION_OF_INTERPOLATOR( _CLASSNAME_, _SUPERCLASSNAME_ ) \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION                            \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, void, _ENABLE_IF_INDEP_IS_REAL,     \
                      _ENABLE_IF_DEP_IS_REAL>                                 \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE>

#define DECLARE_REAL_VERSION_OF_INTERPOLATOR_WITH_PARAMS(               \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                        \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAMS          \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, _PARAMTYPE_,                  \
                      _ENABLE_IF_INDEP_IS_REAL, _ENABLE_IF_DEP_IS_REAL> \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE, _PARAMTYPE_>

#define DECLARE_COMPLEX_VERSION_OF_INTERPOLATOR( _CLASSNAME_,             \
                                                 _SUPERCLASSNAME_ )       \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION                        \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, void, _ENABLE_IF_INDEP_IS_REAL, \
                      _ENABLE_IF_DEP_IS_COMPLEX>                          \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE>

#define DECLARE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAMS(               \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                           \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAMS             \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, _PARAMTYPE_,                     \
                      _ENABLE_IF_INDEP_IS_REAL, _ENABLE_IF_DEP_IS_COMPLEX> \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE, _PARAMTYPE_>

#define DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAMS( \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                         \
    DECLARE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAMS(                 \
        _CLASSNAME_, _SUPERCLASSNAME_,                                   \
        _PARAMTYPE_ ) { public: ~_CLASSNAME_() {} };

#define DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR( \
    _CLASSNAME_, _SUPERCLASSNAME_ )                          \
    DECLARE_COMPLEX_VERSION_OF_INTERPOLATOR(                 \
        _CLASSNAME_, _SUPERCLASSNAME_ ) { public: ~_CLASSNAME_() {} };


#define DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAMS(                 \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                            \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAMS              \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, PARAMTYPE,                        \
                      _ENABLE_IF_INDEP_IS_REAL, _ENABLE_IF_DEP_IS_COMPLEX>  \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE> {                     \
        public:                                                             \
            _CLASSNAME_() {}                                                \
            _CLASSNAME_( PARAMTYPE param ) {                                \
                _real_spline                                                \
                    = _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type>( \
                        param );                                            \
                _imag_spline                                                \
                    = _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type>( \
                        param );                                            \
            }                                                               \
            virtual ~_CLASSNAME_() {                                        \
                this->clear();                                              \
            }                                                               \
            virtual _SUBSPLINE_PTR_T real() override {                      \
                return static_cast<_SUBSPLINE_PTR_T>( &_real_spline );      \
            }                                                               \
            virtual _SUBSPLINE_PTR_T imag() override {                      \
                return static_cast<_SUBSPLINE_PTR_T>( &_imag_spline );      \
            }                                                               \
                                                                            \
        private:                                                            \
            _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type>            \
                _real_spline, _imag_spline;                                 \
    };
#define DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR( _CLASSNAME_,              \
                                                _SUPERCLASSNAME_ )        \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION                        \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, void, _ENABLE_IF_INDEP_IS_REAL, \
                      _ENABLE_IF_DEP_IS_COMPLEX>                          \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE> {                   \
        public:                                                           \
            virtual ~_CLASSNAME_() {                                      \
                this->clear();                                            \
            }                                                             \
            virtual _SUBSPLINE_PTR_T real() override {                    \
                return static_cast<_SUBSPLINE_PTR_T>( &_real_spline );    \
            }                                                             \
            virtual _SUBSPLINE_PTR_T imag() override {                    \
                return static_cast<_SUBSPLINE_PTR_T>( &_imag_spline );    \
            }                                                             \
                                                                          \
        private:                                                          \
            _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type>          \
                _real_spline, _imag_spline;                               \
    };

// _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type>


/*
Structure for declaring real and complex interpolators:
First, declare the default, invalid template, as:

    DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( CLASSNAME, SUPERCLASS );

Then declare the specialized templates for real and complex interpolants, as:

    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION //
    class CLASSNAME<INDEPTYPE, DEPTYPE,
        _ENABLE_IF_INDEP_IS_REAL,_ENABLE_IF_DEP_IS_REAL>
        : public SUPERCLASS<INDEPTYPE,DEPTYPE> { ... };

    // If a pure virtual class, use this
    DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(CLASSNAME,SUPERCLASSNAME)

    // If concrete, but you're not adding any methods, use this
    DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR(CLASSNAME,SUPERCLASSNAME)

    // If concrete and you're adding methods, use this, and change _CLASSNAME_
    // and _SUPERCLASSNAME_ as appropriate.  The other macros should still
work. _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION class
_CLASSNAME_<INDEPTYPE, DEPTYPE, _ENABLE_IF_INDEP_IS_REAL,
                      _ENABLE_IF_DEP_IS_COMPLEX>
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE> {
        public:
            virtual ~_CLASSNAME_() {
                this->clear();
            }
            virtual _SUBSPLINE_PTR_T real() override {
                return static_cast<_SUBSPLINE_PTR_T>( &_real_spline );
            }
            virtual _SUBSPLINE_PTR_T imag() override {
                return static_cast<_SUBSPLINE_PTR_T>( &_imag_spline );
            }

            // any new methods can be added here
        private:
            _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type>
                _real_spline, _imag_spline;
    };

*/


namespace NCPA {
    namespace interpolation {
        namespace details {
            template<typename INDEPTYPE, typename DEPTYPE>
            class _abstract_spline_1d {
                public:
                    virtual ~_abstract_spline_1d() {}

                    virtual void fill( size_t N, const INDEPTYPE *x,
                                       const DEPTYPE *f )
                        = 0;
                    virtual void fill( const std::vector<INDEPTYPE>& x,
                                       const std::vector<DEPTYPE>& f )
                        = 0;
                    virtual void init( size_t N )            = 0;
                    virtual void clear()                     = 0;
                    virtual void ready()                     = 0;
                    virtual DEPTYPE eval_f( INDEPTYPE x )    = 0;
                    virtual DEPTYPE eval_df( INDEPTYPE x )   = 0;
                    virtual DEPTYPE eval_ddf( INDEPTYPE x )  = 0;
                    virtual DEPTYPE eval_dddf( INDEPTYPE x ) = 0;
            };
        }  // namespace details

        // 1-D interpolation engines
        DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _spline_1d,
                                               details::_abstract_spline_1d );

        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
            class _spline_1d<INDEPTYPE, DEPTYPE, _ENABLE_IF_INDEP_IS_REAL,
                             _ENABLE_IF_DEP_IS_REAL>
            : public details::_abstract_spline_1d<INDEPTYPE, DEPTYPE> {
            public:
                virtual ~_spline_1d() {}
        };

        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
            class _spline_1d<INDEPTYPE, DEPTYPE, _ENABLE_IF_INDEP_IS_REAL,
                             _ENABLE_IF_DEP_IS_COMPLEX>
            : public details::_abstract_spline_1d<INDEPTYPE, DEPTYPE> {
            public:
                virtual ~_spline_1d() {}

                virtual void fill( size_t N, const INDEPTYPE *x,
                                   const DEPTYPE *f ) override {
                    auto *r
                        = NCPA::arrays::zeros<typename DEPTYPE::value_type>(
                            N ),
                        *i = NCPA::arrays::zeros<typename DEPTYPE::value_type>(
                            N );
                    NCPA::math::complex2real( N, f, r, i );
                    real()->fill( N, x, r );
                    imag()->fill( N, x, i );
                    delete[] r;
                    delete[] i;
                }

                virtual void fill( const std::vector<INDEPTYPE>& x,
                                   const std::vector<DEPTYPE>& f ) override {
                    size_t N = x.size();
                    if ( N != f.size() ) {
                        throw std::invalid_argument(
                            "Vectors must be same size" );
                    }
                    std::vector<typename DEPTYPE::value_type> r, i;
                    NCPA::math::complex2real( f, r, i );
                    real()->fill( x, r );
                    imag()->fill( x, i );
                }

                virtual void init( size_t N ) override {
                    real()->init( N );
                    imag()->init( N );
                }

                virtual void clear() override {
                    real()->clear();
                    imag()->clear();
                }

                virtual void ready() override {
                    real()->ready();
                    imag()->ready();
                }

                virtual DEPTYPE eval_f( INDEPTYPE x ) override {
                    return DEPTYPE( real()->eval_f( x ), imag()->eval_f( x ) );
                }

                virtual DEPTYPE eval_df( INDEPTYPE x ) override {
                    return DEPTYPE( real()->eval_df( x ),
                                    imag()->eval_df( x ) );
                }

                virtual DEPTYPE eval_ddf( INDEPTYPE x ) override {
                    return DEPTYPE( real()->eval_ddf( x ),
                                    imag()->eval_ddf( x ) );
                }

                virtual DEPTYPE eval_dddf( INDEPTYPE x ) override {
                    return DEPTYPE( real()->eval_dddf( x ),
                                    imag()->eval_dddf( x ) );
                }

                virtual _SUBSPLINE_PTR_T real() = 0;
                virtual _SUBSPLINE_PTR_T imag() = 0;
        };

        namespace LANL {
            namespace details {
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
                            << ", " << x_vals[ length - 1 ] << ")."
                            << std::endl;
                        throw std::range_error( oss.str() );
                    }
                    if ( x >= x_vals[ prev ] && x <= x_vals[ prev + 1 ] ) {
                        done = true;
                    }
                    if ( !done && prev + 2 <= length - 1 ) {
                        if ( x >= x_vals[ prev + 1 ]
                             && x <= x_vals[ prev + 2 ] ) {
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
                            } else if ( x >= x_vals[ ( length - 1 )
                                                     - ( n + 1 ) ]
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

                DECLARE_GENERIC_INTERPOLATOR_TEMPLATE(
                    _lanl_spline_1d, NCPA::interpolation::_spline_1d );

                _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
                    class _lanl_spline_1d<INDEPTYPE, DEPTYPE, void,
                                          _ENABLE_IF_INDEP_IS_REAL,
                                          _ENABLE_IF_DEP_IS_REAL>
                    : public NCPA::interpolation::_spline_1d<INDEPTYPE,
                                                             DEPTYPE> {
                    public:
                        virtual ~_lanl_spline_1d() { clear(); }

                        virtual void fill( size_t N, const INDEPTYPE *x,
                                           const DEPTYPE *f ) override {
                            if ( _length != N ) {
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


                    protected:
                        size_t& _get_length() { return _length; }

                        size_t& _get_accel() { return _accel; }

                        INDEPTYPE *_get_x_vals() { return _x_vals; }

                        DEPTYPE *_get_f_vals() { return _f_vals; }

                        DEPTYPE *_get_slopes() { return _slopes; }

                    private:
                        size_t _length;  // Length of input files (x and f(x))
                        size_t _accel;   // Index of previous table look up;
                                         // used to increase spline speed
                        INDEPTYPE *_x_vals = nullptr;  // 1D array of x values
                        DEPTYPE *_f_vals
                            = nullptr;  // 1D array of f(x) values,
                                        // f_vals[i] = f(x[i])
                        DEPTYPE *_slopes
                            = nullptr;  // Slopes used to generate
                                        // natural cubic spline solution
                };

                DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(
                    _lanl_spline_1d, NCPA::interpolation::_spline_1d );

            }  // namespace details

            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( linear_spline_1d,
                                                   details::_lanl_spline_1d );

            _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
                class linear_spline_1d<INDEPTYPE, DEPTYPE, void,
                                       _ENABLE_IF_INDEP_IS_REAL,
                                       _ENABLE_IF_DEP_IS_REAL>
                : public details::_lanl_spline_1d<INDEPTYPE, DEPTYPE> {
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
                        size_t k          = details::find_segment( x, x_vals,
                                                                   this->_get_length(),
                                                                   this->_get_accel() );
                        return f_vals[ k ] + ( x - x_vals[ k ] ) * slopes[ k ];
                    }

                    virtual DEPTYPE eval_df( INDEPTYPE x ) override {
                        DEPTYPE *f_vals   = this->_get_f_vals();
                        INDEPTYPE *x_vals = this->_get_x_vals();
                        DEPTYPE *slopes   = this->_get_slopes();
                        size_t k          = details::find_segment( x, x_vals,
                                                                   this->_get_length(),
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
                                                    details::_lanl_spline_1d );


            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( natural_cubic_spline_1d,
                                                   details::_lanl_spline_1d );

            _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
                class natural_cubic_spline_1d<INDEPTYPE, DEPTYPE, void,
                                              _ENABLE_IF_INDEP_IS_REAL,
                                              _ENABLE_IF_DEP_IS_REAL>
                : public details::_lanl_spline_1d<INDEPTYPE, DEPTYPE> {
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

                        k = details::find_segment( x, _x_vals,
                                                   this->_get_length(),
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

                        k = details::find_segment( x, _x_vals,
                                                   this->_get_length(),
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

                        k = details::find_segment( x, _x_vals,
                                                   this->_get_length(),
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

                        k = details::find_segment( x, _x_vals,
                                                   this->_get_length(),
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
                                                    details::_lanl_spline_1d );


        }  // namespace LANL

#ifdef HAVE_GSL_INTERPOLATION_LIBRARY
        namespace GSL {

            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE_WITH_PARAMS(
                gsl_spline_1d, NCPA::interpolation::_spline_1d,
                const gsl_interp_type * );

            _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAMS  //
                class gsl_spline_1d<INDEPTYPE, DEPTYPE, PARAMTYPE,
                                    _ENABLE_IF_INDEP_IS_REAL,
                                    _ENABLE_IF_DEP_IS_REAL>
                : public NCPA::interpolation::_spline_1d<INDEPTYPE, DEPTYPE> {
                public:
                    gsl_spline_1d() :
                        _interptype { nullptr },
                        _spline { nullptr },
                        _accel { nullptr },
                        _ready { false },
                        _xmin { 0.0 },
                        _xmax { 0.0 } {}

                    gsl_spline_1d( PARAMTYPE interptype ) : gsl_spline_1d() {
                        set_interpolation_type( interptype );
                    }

                    virtual ~gsl_spline_1d() { clear(); }

                    virtual void set_interpolation_type(
                        PARAMTYPE interptype ) {
                        _interptype = interptype;
                    }

                    virtual void fill( size_t N, const INDEPTYPE *x,
                                       const DEPTYPE *f ) override {
                        if ( N != _size ) {
                            init( N );
                        }
                        _set_spline( N, x, f );
                        _ready = true;
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
                        _check_interptype();
                        _allocate_spline( N );
                    }

                    virtual void clear() override {
                        _free_spline();
                        _free_accel();
                    }

                    virtual void ready() override {
                        _check_spline();
                        _check_accel();
                    }

                    DEPTYPE eval_f( INDEPTYPE x ) override {
                        _check_ready();
                        return (DEPTYPE)gsl_spline_eval( _spline, (double)x,
                                                         _accel );
                    }

                    DEPTYPE eval_df( INDEPTYPE x ) override {
                        _check_ready();
                        return (DEPTYPE)gsl_spline_eval_deriv(
                            _spline, (double)x, _accel );
                    }

                    DEPTYPE eval_ddf( INDEPTYPE x ) override {
                        _check_ready();
                        return (DEPTYPE)gsl_spline_eval_deriv2(
                            _spline, (double)x, _accel );
                    }

                    DEPTYPE eval_dddf( INDEPTYPE x ) override {
                        throw std::out_of_range( "GSL splines have no third "
                                                 "derivative capability." );
                    }

                protected:
                    void _allocate_spline( size_t n ) {
                        if ( n < _interptype->min_size ) {
                            std::ostringstream oss;
                            oss << "GSL interpolation type "
                                << _interptype->name << " requires at least "
                                << _interptype->min_size << " points.";
                            throw std::logic_error( oss.str() );
                        }
                        _spline = gsl_spline_alloc( _interptype, n );
                        _accel  = gsl_interp_accel_alloc();
                        _size   = n;
                    }

                    void _set_spline( size_t n, const INDEPTYPE *x,
                                      const DEPTYPE *y ) {
                        double *xd = NCPA::arrays::zeros<double>( n );
                        double *yd = NCPA::arrays::zeros<double>( n );
                        _to_double<INDEPTYPE>( n, x, xd );
                        _to_double<DEPTYPE>( n, y, yd );
                        gsl_spline_init( _spline, xd, yd, n );
                        _xmin = (INDEPTYPE)x[ 0 ];
                        _xmax = (INDEPTYPE)x[ n - 1 ];
                        delete[] xd;
                        delete[] yd;
                    }

                    template<typename T>
                    void _to_double( size_t n, const T *input,
                                     double *& output ) {
                        for ( auto i = 0; i < n; i++ ) {
                            output[ i ] = (double)( input[ i ] );
                        }
                    }

                    template<>
                    void _to_double( size_t n, const double *input,
                                     double *& output ) {
                        std::memcpy( output, input, n * sizeof( double ) );
                    }

                    void _free_interptype() { _interptype = nullptr; }

                    void _free_spline() {
                        if ( _spline != nullptr ) {
                            gsl_spline_free( _spline );
                            _spline = nullptr;
                        }
                        _free_accel();
                        _ready = false;
                    }

                    void _free_accel() {
                        if ( _accel != nullptr ) {
                            gsl_interp_accel_free( _accel );
                            _accel = nullptr;
                        }
                        _ready = false;
                    }

                    void _check_ready() const {
                        if ( !_ready ) {
                            throw std::logic_error( "Spline not ready" );
                        }
                    }

                    void _check_interptype() const {
                        if ( _interptype == nullptr ) {
                            throw std::logic_error(
                                "GSL Spline type has not been set!" );
                        }
                    }

                    void _check_spline() const {
                        if ( _spline == nullptr ) {
                            throw std::logic_error(
                                "Spline has not been initialized!" );
                        }
                    }

                    void _check_accel() const {
                        if ( _accel == nullptr ) {
                            throw std::logic_error(
                                "Spline accelerator has not been "
                                "initialized" );
                        }
                    }


                private:
                    PARAMTYPE _interptype    = nullptr;
                    gsl_spline *_spline      = nullptr;
                    gsl_interp_accel *_accel = nullptr;
                    bool _ready              = false;
                    INDEPTYPE _xmin = 0.0, _xmax = 0.0;
                    size_t _size = 0;
            };

            DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAMS(
                gsl_spline_1d, NCPA::interpolation::_spline_1d,
                const gsl_interp_type * );

        }  // namespace GSL
#endif

    }  // namespace interpolation
}  // namespace NCPA
