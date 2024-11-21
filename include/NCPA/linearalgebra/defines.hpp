#pragma once

#include "NCPA/types.hpp"

#define _ENABLE_IF_ELEMENTTYPE_IS_REAL \
    typename std::enable_if<std::is_floating_point<ELEMENTTYPE>::value>::type

#define _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC \
    typename std::enable_if<NCPA::types::is_numeric<ELEMENTTYPE>::value>::type

#define NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( \
    _CLASSNAME_ )                                                  \
    template<typename ELEMENTTYPE, typename = void>                \
    class _CLASSNAME_ {}

#define NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( _CLASSNAME_,       \
                                                     _SUPERCLASSNAME_ ) \
    template<typename ELEMENTTYPE, typename = void>                     \
    class _CLASSNAME_ : public _SUPERCLASSNAME_<ELEMENTTYPE> {}

#define NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE \
    template<typename ELEMENTTYPE>

#define NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( _CLASSNAME_, \
                                                     _TYPENAME_ ) \
    template<typename _TYPENAME_>                                 \
    static void swap( _CLASSNAME_<_TYPENAME_>& a,                 \
                      _CLASSNAME_<_TYPENAME_>& b ) noexcept;      \
    template<typename _TYPENAME_>                                 \
    bool operator==( const _CLASSNAME_<_TYPENAME_>& a,            \
                     const _CLASSNAME_<_TYPENAME_>& b );          \
    template<typename _TYPENAME_>                                 \
    bool operator!=( const _CLASSNAME_<_TYPENAME_>& a,            \
                     const _CLASSNAME_<_TYPENAME_>& b );

#define NCPA_LINEARALGEBRA_DECLARE_FRIEND_BINARY_OPERATORS( _CLASSNAME_,    \
                                                            _TYPENAME_ )    \
    template<typename _TYPENAME_>                                           \
    _CLASSNAME_<_TYPENAME_> operator+( const _CLASSNAME_<_TYPENAME_>& c1,   \
                                       const _CLASSNAME_<_TYPENAME_>& c2 ); \
    template<typename _TYPENAME_>                                           \
    _CLASSNAME_<_TYPENAME_> operator+( const _CLASSNAME_<_TYPENAME_>& c1,   \
                                       _TYPENAME_ c2 );                     \
    template<typename _TYPENAME_>                                           \
    _CLASSNAME_<_TYPENAME_> operator+( _TYPENAME_ c2,                       \
                                       const _CLASSNAME_<_TYPENAME_>& c1 ); \
    template<typename _TYPENAME_>                                           \
    _CLASSNAME_<_TYPENAME_> operator-( const _CLASSNAME_<_TYPENAME_>& c1,   \
                                       const _CLASSNAME_<_TYPENAME_>& c2 ); \
    template<typename _TYPENAME_>                                           \
    _CLASSNAME_<_TYPENAME_> operator-( const _CLASSNAME_<_TYPENAME_>& c1,   \
                                       _TYPENAME_ c2 );                     \
    template<typename ELEMENTTYPE>                                          \
    _CLASSNAME_<_TYPENAME_> operator*( const _CLASSNAME_<_TYPENAME_>& c1,   \
                                       const _CLASSNAME_<_TYPENAME_>& c2 ); \
    template<typename _TYPENAME_>                                           \
    _CLASSNAME_<_TYPENAME_> operator*( const _CLASSNAME_<_TYPENAME_>& c1,   \
                                       _TYPENAME_ c2 );                     \
    template<typename _TYPENAME_>                                           \
    _CLASSNAME_<_TYPENAME_> operator*( _TYPENAME_ c2,                       \
                                       const _CLASSNAME_<_TYPENAME_>& c1 );

