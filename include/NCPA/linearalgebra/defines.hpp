#pragma once

#include "NCPA/types.hpp"

#define _ENABLE_IF_ELEMENTTYPE_IS_REAL \
    typename std::enable_if<std::is_floating_point<ELEMENTTYPE>::value>::type

#define _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC \
    typename std::enable_if<NCPA::types::is_numeric<ELEMENTTYPE>::value>::type

#define NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( _CLASSNAME_,       \
                                                     _SUPERCLASSNAME_ ) \
    template<typename ELEMENTTYPE, typename = void>                     \
    class _CLASSNAME_ : public _SUPERCLASSNAME_<ELEMENTTYPE> {          \
            static_assert( false, "Invalid type for _CLASSNAME_" );     \
    }

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

#define DECLARE_BINARY_OPERATORS( _CLASSNAME_ )                               \
    template<typename ELEMENTTYPE>                                            \
    _CLASSNAME_<ELEMENTTYPE> operator+( const _CLASSNAME_<ELEMENTTYPE>& c1,   \
                                        const _CLASSNAME_<ELEMENTTYPE>& c2 ); \
    template<typename ELEMENTTYPE>                                            \
    _CLASSNAME_<ELEMENTTYPE> operator-( const _CLASSNAME_<ELEMENTTYPE>& c1,   \
                                        const _CLASSNAME_<ELEMENTTYPE>& c2 ); \
    template<typename ELEMENTTYPE>                                            \
    _CLASSNAME_<ELEMENTTYPE> operator*( const _CLASSNAME_<ELEMENTTYPE>& c1,   \
                                        const _CLASSNAME_<ELEMENTTYPE>& c2 );

