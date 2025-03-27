#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

/*
Structure for declaring real and complex interpolators:
First, declare the default, invalid template, as:

    DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( CLASSNAME, SUPERCLASS );

Then declare the specialized templates for real and complex interpolants.  This
will look slightly different depending on whether the new interpolator class's
constructor requires an additional parameter.

NO ADDITIONAL CONSTRUCTOR PARAMETER
-----------------------------------
If there are no additional parameters (i.e. all of the specialization is baked
into the class already), declare and define the real and complex versions.  We
can omit the default constructor and let the compiler do it for us.  For a
concrete class, override all of the pure virtual functions in the real class,
and use the DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR macro:

    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION //
    class CLASSNAME<INDEPTYPE, DEPTYPE, void,
        ENABLE_IF_REAL( INDEPTYPE ),ENABLE_IF_REAL( DEPTYPE )>
        : public SUPERCLASS<INDEPTYPE,DEPTYPE> {
        public:
            virtual ~CLASSNAME() { ... }
            virtual void fill( size_t N, const INDEPTYPE *x,
                               const DEPTYPE *f ) override { ... }
            virtual void fill( const std::vector<INDEPTYPE>& x,
                               const std::vector<DEPTYPE>& f )
                               override { ... }
            virtual void init( size_t N ) override { ... }
            virtual void clear() override { ... }
            virtual void ready() override { ... }
            virtual DEPTYPE eval_f( INDEPTYPE x ) override { ... }
            virtual DEPTYPE eval_df( INDEPTYPE x ) override { ... }
            virtual DEPTYPE eval_ddf( INDEPTYPE x ) override { ... }
            virtual DEPTYPE eval_dddf( INDEPTYPE x ) override { ... }
            // etc
    };
    DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR(CLASSNAME,SUPERCLASSNAME)

If the new class is still pure virtual, override the destructor (you can leave
it empty if it doesn't need to do anything), and use the
DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR macro:

    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION //
    class CLASSNAME<INDEPTYPE, DEPTYPE, void,
        ENABLE_IF_REAL( INDEPTYPE ),ENABLE_IF_REAL( DEPTYPE )>
        : public SUPERCLASS<INDEPTYPE,DEPTYPE> {
        public:
            virtual ~CLASSNAME() {}
    };
    DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(CLASSNAME,SUPERCLASSNAME)

Note for the keen-eyed: the // after the declaration macro is just there to
tell the clang formatter not to put a line break there.  The line break doesn't
affect the functionality, but it makes it less clear to read.

ADDITIONAL CONSTRUCTOR PARAMETER
--------------------------------
If there as an additional parameter, first determine the type of the declare
and define the real and complex versions. For a concrete class, declare default
and parameterized constructores, override all of the pure virtual functions in
the real class, and use the _WITH_PARAM versions of the macros:

    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAM //
    class CLASSNAME<INDEPTYPE, DEPTYPE, PARAMTYPE,
        ENABLE_IF_REAL( INDEPTYPE ),ENABLE_IF_REAL( DEPTYPE )>
        : public SUPERCLASS<INDEPTYPE,DEPTYPE> {
        public:
            CLASSNAME() : SUPERCLASS<INDEPTYPE,DEPTYPE>() {}
            CLASSNAME( PARAMTYPE param )
                : CLASSNAME<INDEPTYPE,DEPTYPE,PARAMTYPE>() { ... }
            virtual ~CLASSNAME() { ... }
            virtual void fill( size_t N, const INDEPTYPE *x,
                               const DEPTYPE *f ) override { ... }
            virtual void fill( const std::vector<INDEPTYPE>& x,
                               const std::vector<DEPTYPE>& f )
                               override { ... }
            virtual void init( size_t N ) override { ... }
            virtual void clear() override { ... }
            virtual void ready() override { ... }
            virtual DEPTYPE eval_f( INDEPTYPE x ) override { ... }
            virtual DEPTYPE eval_df( INDEPTYPE x ) override { ... }
            virtual DEPTYPE eval_ddf( INDEPTYPE x ) override { ... }
            virtual DEPTYPE eval_dddf( INDEPTYPE x ) override { ... }
            // etc
    };
    DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAM(CLASSNAME,SUPERCLASSNAME,
        PARAMTYPE);

If the new class is still pure virtual, override any methods as appropriate,
but at minimum override the destructor (you can leave it empty if it doesn't
need to do anything), and use the
DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR macro:

    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION //
    class CLASSNAME<INDEPTYPE, DEPTYPE, void,
        ENABLE_IF_REAL( INDEPTYPE ),ENABLE_IF_REAL( DEPTYPE )>
        : public SUPERCLASS<INDEPTYPE,DEPTYPE> {
        public:
            virtual ~CLASSNAME() {}
    };
    DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR(CLASSNAME,SUPERCLASSNAME);

*/

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_spline_2d<T, U>& a,
    NCPA::interpolation::_spline_2d<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        template<typename INDEPTYPE, typename DEPTYPE>
        class _abstract_spline_2d {
            public:
                virtual ~_abstract_spline_2d() {}

                virtual void fill( size_t N1, size_t N2, const INDEPTYPE *x1,
                                   const INDEPTYPE *x2, const DEPTYPE **f1 )
                    = 0;
                virtual void fill( const std::vector<INDEPTYPE>& x1,
                                   const std::vector<INDEPTYPE>& x2,
                                   const NCPA::arrays::vector2d_t<DEPTYPE>& f )
                    = 0;
                virtual void init( size_t N1, size_t N2 )            = 0;
                virtual void clear()                                 = 0;
                virtual void ready()                                 = 0;
                virtual DEPTYPE eval_f( INDEPTYPE x1, INDEPTYPE x2 ) = 0;
                virtual DEPTYPE eval_df( INDEPTYPE x1, INDEPTYPE x2,
                                         size_t wrt )
                    = 0;
                virtual DEPTYPE eval_ddf( INDEPTYPE x1, INDEPTYPE x2,
                                          size_t wrt1, size_t wrt2 )
                    = 0;
                virtual DEPTYPE eval_dddf( INDEPTYPE x1, INDEPTYPE x2,
                                           size_t wrt1, size_t wrt2,
                                           size_t wrt3 )
                    = 0;
                virtual std::array<INDEPTYPE,4> limits() const = 0;
                virtual interpolator_2d_type_t interptype() const  = 0;
        };

        // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _spline_2d,
        //                                        _abstract_spline_2d );

        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
            class _spline_2d<INDEPTYPE, DEPTYPE, _ENABLE_IF_INDEP_IS_REAL,
                             _ENABLE_IF_DEP_IS_REAL>
            : public _abstract_spline_2d<INDEPTYPE, DEPTYPE> {
            public:
                virtual ~_spline_2d() {}
        };

        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
            class _spline_2d<INDEPTYPE, DEPTYPE, _ENABLE_IF_INDEP_IS_REAL,
                             _ENABLE_IF_DEP_IS_COMPLEX>
            : public _abstract_spline_2d<INDEPTYPE, DEPTYPE> {
            public:
                virtual ~_spline_2d() {}

                friend void ::swap<INDEPTYPE, DEPTYPE>(
                    _spline_2d<INDEPTYPE, DEPTYPE>& a,
                    _spline_2d<INDEPTYPE, DEPTYPE>& b ) noexcept;

                virtual void fill( size_t N1, size_t N2, const INDEPTYPE *x1,
                                   const INDEPTYPE *x2,
                                   const DEPTYPE **f ) override {
                    auto **rl
                        = NCPA::arrays::zeros<typename DEPTYPE::value_type>(
                            N1, N2 ),
                        **im
                        = NCPA::arrays::zeros<typename DEPTYPE::value_type>(
                            N1, N2 );
                    for ( size_t i = 0; i < N1; ++i ) {
                        NCPA::math::complex2real( N2, f[ i ], rl[ i ],
                                                  im[ i ] );
                    }
                    this->real()->fill( N1, N2, x1, x2, rl );
                    this->imag()->fill( N1, N2, x1, x2, im );
                    NCPA::arrays::free_array( rl, N1, N2 );
                    NCPA::arrays::free_array( im, N1, N2 );
                }

                virtual void fill( const std::vector<INDEPTYPE>& x1,
                                   const std::vector<INDEPTYPE>& x2,
                                   const NCPA::arrays::vector2d_t<DEPTYPE>& f ) override {
                    size_t N1 = x1.size();
                    if ( N1 != f.size() ) {
                        throw std::invalid_argument(
                            "Vectors must be same size" );
                    }
                    NCPA::arrays::vector2d_t<typename DEPTYPE::value_type> rl, im;
                    for ( size_t i = 0; i < N1; ++i ) {
                        NCPA::math::complex2real( f[ i ], rl[ i ], im[ i ] );
                    }
                    real()->fill( x1, x2, rl );
                    imag()->fill( x1, x2, im );
                }

                virtual void init( size_t N1, size_t N2 ) override {
                    real()->init( N1, N2 );
                    imag()->init( N1, N2 );
                }

                virtual void clear() override {
                    real()->clear();
                    imag()->clear();
                }

                virtual void ready() override {
                    real()->ready();
                    imag()->ready();
                }

                virtual DEPTYPE eval_f( INDEPTYPE x1, INDEPTYPE x2 ) override {
                    return DEPTYPE( real()->eval_f( x1, x2 ),
                                    imag()->eval_f( x1, x2 ) );
                }

                virtual DEPTYPE eval_df( INDEPTYPE x1, INDEPTYPE x2, size_t wrt ) override {
                    return DEPTYPE( real()->eval_df( x1, x2, wrt ),
                                    imag()->eval_df( x1, x2, wrt ) );
                }

                virtual DEPTYPE eval_ddf( INDEPTYPE x1, INDEPTYPE x2,
                                          size_t wrt1, size_t wrt2) override {
                    return DEPTYPE( real()->eval_ddf( x1, x2, wrt1, wrt2 ),
                                    imag()->eval_ddf( x1, x2, wrt1, wrt2 ) );
                }

                virtual DEPTYPE eval_dddf( INDEPTYPE x1, INDEPTYPE x2,
                                          size_t wrt1, size_t wrt2, size_t wrt3 ) override {
                    return DEPTYPE( real()->eval_dddf( x1, x2, wrt1, wrt2, wrt3 ),
                                    imag()->eval_dddf( x1, x2, wrt1, wrt2, wrt3 ) );
                }

                virtual std::array<INDEPTYPE,4> limits() const override {
                    return this->real()->limits();
                }

                virtual _SUBSPLINE_2D_PTR_T real() = 0;
                virtual _SUBSPLINE_2D_PTR_T imag() = 0;
        };
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap(
    NCPA::interpolation::_spline_2d<T, U>& a,
    NCPA::interpolation::_spline_2d<T, U>& b ) noexcept {}