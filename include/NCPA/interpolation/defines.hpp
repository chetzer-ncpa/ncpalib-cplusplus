#pragma once

#include "NCPA/defines.hpp"

#if __has_include( "gsl/gsl_spline.h" )
#  define NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE true
#endif

// Convenience #defines
#define _ENABLE_IF_INDEP_IS_REAL ENABLE_IF_REAL( INDEPTYPE )
    // typename std::enable_if<std::is_floating_point<INDEPTYPE>::value>::type
#define _ENABLE_IF_DEP_IS_REAL ENABLE_IF_REAL( DEPTYPE )
    // typename std::enable_if<std::is_floating_point<DEPTYPE>::value>::type
#define _ENABLE_IF_DEP_IS_COMPLEX ENABLE_IF_COMPLEX( DEPTYPE )
    // typename std::enable_if<NCPA::types::is_complex<DEPTYPE>::value>::type
#define _ENABLE_IF_DEP_IS_NUMERIC ENABLE_IF_NUMERIC( DEPTYPE )
    // typename std::enable_if<NCPA::types::is_numeric<DEPTYPE>::value>::type
#define _SUBSPLINE_PTR_T \
    NCPA::interpolation::_spline_1d<INDEPTYPE, typename DEPTYPE::value_type> *
#define _SUBSPLINE_2D_PTR_T \
    NCPA::interpolation::_spline_2d<INDEPTYPE, typename DEPTYPE::value_type> *
#define _SUBSPLINE_T( _CLASSNAME_ ) \
    _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type>

#define DECLARE_GENERIC_INTERPOLATOR_TEMPLATE_WITH_PARAM(                 \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                          \
    template<typename INDEPTYPE, typename DEPTYPE,                        \
             typename PARAMTYPE = _PARAMTYPE_, typename = void,           \
             typename = void>                                             \
    class _CLASSNAME_ : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE> {     \
            /*static_assert( false, "Invalid types for interpolator" );*/ \
    }
#define DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _CLASSNAME_,               \
                                               _SUPERCLASSNAME_ )         \
    template<typename INDEPTYPE, typename DEPTYPE, typename = void,       \
             typename = void, typename = void>                            \
    class _CLASSNAME_ : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE> {     \
            /*static_assert( false, "Invalid types for interpolator" );*/ \
    }


#define _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAM \
    template<typename INDEPTYPE, typename DEPTYPE, typename PARAMTYPE>
#define _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION \
    template<typename INDEPTYPE, typename DEPTYPE>

#define DECLARE_REAL_VERSION_OF_INTERPOLATOR( _CLASSNAME_, _SUPERCLASSNAME_ ) \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION                            \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, void, _ENABLE_IF_INDEP_IS_REAL,     \
                      _ENABLE_IF_DEP_IS_REAL>                                 \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE>

#define DECLARE_REAL_VERSION_OF_INTERPOLATOR_WITH_PARAM(                \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                        \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAM           \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, _PARAMTYPE_,                  \
                      _ENABLE_IF_INDEP_IS_REAL, _ENABLE_IF_DEP_IS_REAL> \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE, _PARAMTYPE_>

#define DECLARE_COMPLEX_VERSION_OF_INTERPOLATOR( _CLASSNAME_,             \
                                                 _SUPERCLASSNAME_ )       \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION                        \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, void, _ENABLE_IF_INDEP_IS_REAL, \
                      _ENABLE_IF_DEP_IS_COMPLEX>                          \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE>

#define DECLARE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAM(                \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                           \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAM              \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, _PARAMTYPE_,                     \
                      _ENABLE_IF_INDEP_IS_REAL, _ENABLE_IF_DEP_IS_COMPLEX> \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE, _PARAMTYPE_>

#define DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAM( \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                        \
    DECLARE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAM(                 \
        _CLASSNAME_, _SUPERCLASSNAME_,                                  \
        _PARAMTYPE_ ) { public: ~_CLASSNAME_() {} };

#define DEFINE_PURE_VIRTUAL_COMPLEX_VERSION_OF_INTERPOLATOR( \
    _CLASSNAME_, _SUPERCLASSNAME_ )                          \
    DECLARE_COMPLEX_VERSION_OF_INTERPOLATOR(                 \
        _CLASSNAME_, _SUPERCLASSNAME_ ) { public: ~_CLASSNAME_() {} };


#define DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAM(                  \
    _CLASSNAME_, _SUPERCLASSNAME_, _PARAMTYPE_ )                            \
    _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAM               \
    class _CLASSNAME_<INDEPTYPE, DEPTYPE, PARAMTYPE,                        \
                      _ENABLE_IF_INDEP_IS_REAL, _ENABLE_IF_DEP_IS_COMPLEX>  \
        : public _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE> {                     \
        public:                                                             \
            _CLASSNAME_() : _SUPERCLASSNAME_<INDEPTYPE, DEPTYPE>() {}       \
            _CLASSNAME_( PARAMTYPE param ) :                                \
                _CLASSNAME_<INDEPTYPE, DEPTYPE, PARAMTYPE>() {              \
                _real_spline                                                \
                    = _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type,  \
                                  PARAMTYPE>( param );                      \
                _imag_spline                                                \
                    = _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type,  \
                                  PARAMTYPE>( param );                      \
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
            _CLASSNAME_<INDEPTYPE, typename DEPTYPE::value_type, PARAMTYPE> \
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
