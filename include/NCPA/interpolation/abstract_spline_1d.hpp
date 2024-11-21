#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
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
    }  // namespace interpolation
}  // namespace NCPA
