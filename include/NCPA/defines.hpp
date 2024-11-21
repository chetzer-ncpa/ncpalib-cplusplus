#pragma once

#include "NCPA/types.hpp"

// #define ENABLE_IF( CONDITION ) \
//     typename std::enable_if<CONDITION::value, int>::type ENABLER = 0
// #define ENABLE_AND( CONDITION1, CONDITION2 ) \
//     typename std::enable_if<CONDITION1::value, int>::type ENABLER1 = 0, \
//     typename std::enable_if<CONDITION2::value, int>::type ENABLER2 = 0
// #define NO_DEFAULT_ENABLE_IF( CONDITION ) \
//     typename std::enable_if<CONDITION::value, int>::type ENABLER
// #define ENABLE_IF_T( CONDITION, T ) \
//     typename std::enable_if<CONDITION<T>::value, int>::type ENABLER = 0
#define ENABLE_IF_TU( CONDITION, T, U ) \
    typename std::enable_if<CONDITION<T, U>::value, int>::type ENABLER = 0
// #define ENABLE_IF_REAL( T ) ENABLE_IF( std::is_floating_point<T> )

#define ENABLE_IF_REAL( _TYPE_ ) \
    typename std::enable_if<std::is_floating_point<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_REAL( _TYPE_ ) \
    typename std::enable_if<std::is_floating_point<_TYPE_>::value, int>::type ENABLER = 0

#define ENABLE_IF_INTEGRAL( _TYPE_ ) \
    typename std::enable_if<std::is_integral<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_INTEGRAL( _TYPE_ ) \
    typename std::enable_if<std::is_integral<_TYPE_>::value, int>::type ENABLER = 0

#define ENABLE_IF_NUMERIC( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_numeric<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_NUMERIC( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_numeric<_TYPE_>::value, int>::type ENABLER = 0

#define ENABLE_IF_ARITHMETIC( _TYPE_ ) \
    typename std::enable_if<std::is_arithmetic<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_ARITHMETIC( _TYPE_ ) \
    typename std::enable_if<std::is_arithmetic<_TYPE_>::value, int>::type ENABLER = 0

#define ENABLE_IF_COMPLEX( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_complex<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_COMPLEX( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_complex<_TYPE_>::value, int>::type ENABLER = 0

#define ENABLE_IF_DELETEABLE( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_deleteable<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_DELETEABLE( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_deleteable<_TYPE_>::value, int>::type ENABLER = 0

#define ENABLE_FUNCTION_IF_NOT_DELETEABLE( _TYPE_ ) \
    typename std::enable_if<!NCPA::types::is_deleteable<_TYPE_>::value, int>::type ENABLER = 0

#define ENABLE_IF_ITERABLE( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_iterable<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_ITERABLE( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_iterable<_TYPE_>::value, int>::type ENABLER = 0


#define DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( _CLASSNAME_ ) \
    template<typename ELEMENTTYPE, typename = void>           \
    class _CLASSNAME_ {}

#define DECLARE_GENERIC_TEMPLATE( _CLASSNAME_, _SUPERCLASSNAME_ ) \
    template<typename ELEMENTTYPE, typename = void>               \
    class _CLASSNAME_ : public _SUPERCLASSNAME_<ELEMENTTYPE> {}

#define DECLARE_SPECIALIZED_TEMPLATE( _TYPE_ ) template<typename _TYPE_>

#define DECLARE_FRIEND_FUNCTIONS( _CLASSNAME_, _TYPENAME_ )  \
    template<typename _TYPENAME_>                            \
    static void swap( _CLASSNAME_<_TYPENAME_>& a,            \
                      _CLASSNAME_<_TYPENAME_>& b ) noexcept; \
    template<typename _TYPENAME_>                            \
    bool operator==( const _CLASSNAME_<_TYPENAME_>& a,       \
                     const _CLASSNAME_<_TYPENAME_>& b );     \
    template<typename _TYPENAME_>                            \
    bool operator!=( const _CLASSNAME_<_TYPENAME_>& a,       \
                     const _CLASSNAME_<_TYPENAME_>& b );

#define DECLARE_FRIEND_BINARY_OPERATORS( _CLASSNAME_, _TYPENAME_ )          \
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

// #define DEFINE_AS_BOOL( _CONDITION_ )       \
//     explicit operator bool() const { \
//         return _CONDITION_;          \
//     }
