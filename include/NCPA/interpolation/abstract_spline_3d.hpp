#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

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
static void swap( NCPA::interpolation::_spline_3d<T, U>& a,
                  NCPA::interpolation::_spline_3d<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        template<typename INDEPTYPE, typename DEPTYPE>
        class _abstract_spline_3d {
            public:
                virtual ~_abstract_spline_3d() {}

                virtual void fill( size_t N1, size_t N2, size_t N3,
                                   const INDEPTYPE *x1, const INDEPTYPE *x2,
                                   const INDEPTYPE *x3, const DEPTYPE ***f1 )
                    = 0;
                virtual void fill( const std::vector<INDEPTYPE>& x1,
                                   const std::vector<INDEPTYPE>& x2,
                                   const std::vector<INDEPTYPE>& x3,
                                   const NCPA::arrays::vector3d_t<DEPTYPE>& f )
                    = 0;
                virtual void init( size_t N1, size_t N2, size_t N3 ) = 0;
                virtual void clear()                                 = 0;
                virtual void ready()                                 = 0;
                virtual DEPTYPE eval_f( INDEPTYPE x1, INDEPTYPE x2,
                                        INDEPTYPE x3 )
                    = 0;
                virtual DEPTYPE eval_df( INDEPTYPE x1, INDEPTYPE x2,
                                         INDEPTYPE x3, size_t dim )
                    = 0;
                virtual DEPTYPE eval_ddf( INDEPTYPE x1, INDEPTYPE x2,
                                          INDEPTYPE x3, size_t dim1,
                                          size_t dim2 )
                    = 0;
                virtual DEPTYPE eval_dddf( INDEPTYPE x1, INDEPTYPE x2,
                                           INDEPTYPE x3, size_t dim1,
                                           size_t dim2, size_t dim3 )
                    = 0;
                virtual interpolator_3d_type_t interptype() const = 0;
        };

        // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _spline_3d,
        //                                        _abstract_spline_3d );

        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
            class _spline_3d<INDEPTYPE, DEPTYPE, _ENABLE_IF_INDEP_IS_REAL,
                             _ENABLE_IF_DEP_IS_REAL>
            : public _abstract_spline_3d<INDEPTYPE, DEPTYPE> {
            public:
                virtual ~_spline_3d() {}
        };

        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
            class _spline_3d<INDEPTYPE, DEPTYPE, _ENABLE_IF_INDEP_IS_REAL,
                             _ENABLE_IF_DEP_IS_COMPLEX>
            : public _abstract_spline_3d<INDEPTYPE, DEPTYPE> {
            public:
                virtual ~_spline_3d() {}

                friend void ::swap<INDEPTYPE, DEPTYPE>(
                    _spline_3d<INDEPTYPE, DEPTYPE>& a,
                    _spline_3d<INDEPTYPE, DEPTYPE>& b ) noexcept;

                virtual void fill( size_t N1, size_t N2, size_t N3,
                                   const INDEPTYPE *x1, const INDEPTYPE *x2,
                                   const INDEPTYPE *x3,
                                   const DEPTYPE ***f ) override {
                    auto ***rl
                        = NCPA::arrays::zeros<typename DEPTYPE::value_type>(
                            N1, N2, N3 ),
                        ***im
                        = NCPA::arrays::zeros<typename DEPTYPE::value_type>(
                            N1, N2, N3 );
                    for ( size_t i = 0; i < N1; ++i ) {
                        for ( size_t j = 0; j < N2; ++j ) {
                            NCPA::math::complex2real(
                                N3, f[ i ][ j ], rl[ i ][ j ], im[ i ][ j ] );
                        }
                    }
                    this->real()->fill( N1, N2, N3, x1, x2, x3, rl );
                    this->imag()->fill( N1, N2, N3, x1, x2, x3, rl );
                    NCPA::arrays::free_array( rl, N1, N2, N3 );
                    NCPA::arrays::free_array( im, N1, N2, N3 );
                }

                virtual void fill(
                    const std::vector<INDEPTYPE>& x1,
                    const std::vector<INDEPTYPE>& x2,
                    const std::vector<INDEPTYPE>& x3,
                    const NCPA::arrays::vector3d_t<DEPTYPE>& f ) override {
                    size_t N1 = x1.size(), N2 = x2.size(), N3 = x3.size();
                    if ( N1 != f.size() ) {
                        throw std::invalid_argument(
                            "Vectors must be same size" );
                    }
                    NCPA::arrays::vector3d_t<typename DEPTYPE::value_type> rl,
                        im;
                    for ( size_t i = 0; i < N1; ++i ) {
                        for ( size_t j = 0; j < N2; ++j ) {
                            NCPA::math::complex2real(
                                f[ i ][ j ], rl[ i ][ j ], im[ i ][ j ] );
                        }
                    }
                    real()->fill( x1, x2, x3, rl );
                    imag()->fill( x1, x2, x3, im );
                }

                virtual void init( size_t N1, size_t N2, size_t N3 ) override {
                    real()->init( N1, N2, N3 );
                    imag()->init( N1, N2, N3 );
                }

                virtual void clear() override {
                    real()->clear();
                    imag()->clear();
                }

                virtual void ready() override {
                    real()->ready();
                    imag()->ready();
                }

                virtual DEPTYPE eval_f( INDEPTYPE x1, INDEPTYPE x2,
                                        INDEPTYPE x3 ) override {
                    return DEPTYPE( real()->eval_f( x1, x2, x3 ),
                                    imag()->eval_f( x1, x2, x3 ) );
                }

                virtual DEPTYPE eval_df( INDEPTYPE x1, INDEPTYPE x2,
                                         INDEPTYPE x3, size_t dim ) override {
                    return DEPTYPE( real()->eval_df( x1, x2, x3, dim ),
                                    imag()->eval_df( x1, x2, x3, dim ) );
                }

                virtual DEPTYPE eval_ddf( INDEPTYPE x1, INDEPTYPE x2,
                                          INDEPTYPE x3, size_t dim1,
                                          size_t dim2 ) override {
                    return DEPTYPE(
                        real()->eval_ddf( x1, x2, x3, dim1, dim2 ),
                        imag()->eval_ddf( x1, x2, x3, dim1, dim2 ) );
                }

                virtual DEPTYPE eval_dddf( INDEPTYPE x1, INDEPTYPE x2,
                                           INDEPTYPE x3, size_t dim1,
                                           size_t dim2,
                                           size_t dim3 ) override {
                    return DEPTYPE(
                        real()->eval_dddf( x1, x2, x3, dim1, dim2, dim3 ),
                        imag()->eval_dddf( x1, x2, x3, dim1, dim2, dim3 ) );
                }

                virtual _SUBSPLINE_2D_PTR_T real() = 0;
                virtual _SUBSPLINE_2D_PTR_T imag() = 0;
        };
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap( NCPA::interpolation::_spline_3d<T, U>& a,
                  NCPA::interpolation::_spline_3d<T, U>& b ) noexcept {}
