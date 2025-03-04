#pragma once

#include "NCPA/types.hpp"

#include <stdexcept>
#include <string>

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

#define ENABLE_FUNCTION_IF_REAL( _TYPE_ )                                     \
    typename std::enable_if<std::is_floating_point<_TYPE_>::value, int>::type \
        ENABLER                                                               \
        = 0

#define ENABLE_METHOD_IF_REAL( _TYPE_ )                            \
    typename std::enable_if<std::is_floating_point<_TYPE_>::value, \
                            std::nullptr_t>::type                  \
        = nullptr

#define ENABLE_IF_INTEGRAL( _TYPE_ ) \
    typename std::enable_if<std::is_integral<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_INTEGRAL( _TYPE_ )                           \
    typename std::enable_if<std::is_integral<_TYPE_>::value, int>::type \
        ENABLER                                                         \
        = 0

#define ENABLE_METHOD_IF_INTEGRAL( _TYPE_ )                  \
    typename std::enable_if<std::is_integral<_TYPE_>::value, \
                            std::nullptr_t>::type            \
        = nullptr

#define ENABLE_IF_STRING( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_std_string<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_STRING( _TYPE_ )                            \
    typename std::enable_if<NCPA::types::is_std_string<_TYPE_>::value, \
                            int>::type ENABLER                         \
        = 0

#define ENABLE_METHOD_IF_STRING( _TYPE_ )                              \
    typename std::enable_if<NCPA::types::is_std_string<_TYPE_>::value, \
                            std::nullptr_t>::type                      \
        = nullptr

#define ENABLE_IF_NUMERIC( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_numeric<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_NUMERIC( _TYPE_ )                        \
    typename std::enable_if<NCPA::types::is_numeric<_TYPE_>::value, \
                            int>::type ENABLER                      \
        = 0

#define ENABLE_FUNCTION_IF_NOT_NUMERIC( _TYPE_ )                     \
    typename std::enable_if<!NCPA::types::is_numeric<_TYPE_>::value, \
                            int>::type ENABLER                       \
        = 0

#define ENABLE_IF_ARITHMETIC( _TYPE_ ) \
    typename std::enable_if<std::is_arithmetic<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_ARITHMETIC( _TYPE_ )                           \
    typename std::enable_if<std::is_arithmetic<_TYPE_>::value, int>::type \
        ENABLER                                                           \
        = 0

#define ENABLE_IF_COMPLEX( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_complex<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_COMPLEX( _TYPE_ )                        \
    typename std::enable_if<NCPA::types::is_complex<_TYPE_>::value, \
                            int>::type ENABLER                      \
        = 0

#define ENABLE_IF_DELETEABLE( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_deleteable<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_DELETEABLE( _TYPE_ )                        \
    typename std::enable_if<NCPA::types::is_deleteable<_TYPE_>::value, \
                            int>::type ENABLER                         \
        = 0

#define ENABLE_FUNCTION_IF_NOT_DELETEABLE( _TYPE_ )                     \
    typename std::enable_if<!NCPA::types::is_deleteable<_TYPE_>::value, \
                            int>::type ENABLER                          \
        = 0

#define ENABLE_IF_ITERABLE( _TYPE_ ) \
    typename std::enable_if<NCPA::types::is_iterable<_TYPE_>::value>::type

#define ENABLE_FUNCTION_IF_ITERABLE( _TYPE_ )                        \
    typename std::enable_if<NCPA::types::is_iterable<_TYPE_>::value, \
                            int>::type ENABLER                       \
        = 0


#define DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( _CLASSNAME_ ) \
    template<typename ELEMENTTYPE, typename = void>           \
    class _CLASSNAME_ {}

#define DECLARE_GENERIC_TEMPLATE( _CLASSNAME_, _SUPERCLASSNAME_ ) \
    template<typename ELEMENTTYPE, typename = void>               \
    class _CLASSNAME_ : public _SUPERCLASSNAME_<ELEMENTTYPE> {}

#define DECLARE_SPECIALIZED_TEMPLATE( _TYPE_ ) template<typename _TYPE_>

#define DECLARE_SWAP_FUNCTION( _CLASSNAME_ )      \
    template<typename _TYPENAME_>                 \
    static void swap( _CLASSNAME_<_TYPENAME_>& a, \
                      _CLASSNAME_<_TYPENAME_>& b ) noexcept;

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

#define DECLARE_METHOD_ENABLED_IF_STRING( _METHODNAME_, _RETURNTYPE_, \
                                          _TYPENAME_ )                \
    template<typename T_ = _TYPENAME_>                                \
    _RETURNTYPE_ _METHODNAME_( ENABLE_METHOD_IF_STRING( T_ ) )

#define DECLARE_METHOD_ENABLED_IF_INTEGRAL( _METHODNAME_, _RETURNTYPE_, \
                                            _TYPENAME_ )                \
    template<typename T_ = _TYPENAME_>                                  \
    _RETURNTYPE_ _METHODNAME_( ENABLE_METHOD_IF_INTEGRAL( T_ ) )

#define DECLARE_METHOD_ENABLED_IF_REAL( _METHODNAME_, _RETURNTYPE_, \
                                        _TYPENAME_ )                \
    template<typename T_ = _TYPENAME_>                              \
    _RETURNTYPE_ _METHODNAME_( ENABLE_METHOD_IF_REAL( T_ ) )

#define DECLARE_CLONE_FUNCTION( _THISCLASS_, _SUPERCLASS_ )               \
    virtual std::unique_ptr<_SUPERCLASS_> clone() const override {        \
        return std::unique_ptr<_SUPERCLASS_>( new _THISCLASS_( *this ) ); \
    }

#define DECLARE_CLONE_TEMPLATE( _THISCLASS_, _SUPERCLASS_, _TYPE_ )        \
    virtual std::unique_ptr<_SUPERCLASS_<_TYPE_>> clone() const override { \
        return std::unique_ptr<_SUPERCLASS_<_TYPE_>>(                      \
            new _THISCLASS_( *this ) );                                    \
    }

#define DECLARE_CLONE_TEMPLATE2( _THISCLASS_, _SUPERCLASS_, _TYPE1_, \
                                 _TYPE2_ )                           \
    virtual std::unique_ptr<_SUPERCLASS_<_TYPE1_, _TYPE2_>> clone()  \
        const override {                                             \
        return std::unique_ptr<_SUPERCLASS_<_TYPE1_, _TYPE2_>>(      \
            new _THISCLASS_( *this ) );                              \
    }

#define DECLARE_MOVE_CONSTRUCTOR( _THISCLASS_ )                    \
    _THISCLASS_( _THISCLASS_&& source ) noexcept : _THISCLASS_() { \
        ::swap( *this, source );                                   \
    }

#define DECLARE_MOVE_CONSTRUCTOR_TEMPLATE( _THISCLASS_, _TYPE1_ ) \
    _THISCLASS_( _THISCLASS_<_TYPE1_>&& source ) noexcept :       \
        _THISCLASS_<_TYPE1_>() {                                  \
        ::swap( *this, source );                                  \
    }

#define DECLARE_MOVE_CONSTRUCTOR_TEMPLATE2( _THISCLASS_, _TYPE1_, _TYPE2_ ) \
    _THISCLASS_( _THISCLASS_<_TYPE1_, _TYPE2_>&& source ) noexcept :        \
        _THISCLASS_<_TYPE1_, _TYPE2_>() {                                   \
        ::swap( *this, source );                                            \
    }

#define DECLARE_EMPTY_DESTRUCTOR( _THISCLASS_ ) \
    virtual ~_THISCLASS_() {}

#define DECLARE_ASSIGNMENT_OPERATOR( _THISCLASS_ ) \
    _THISCLASS_& operator=( _THISCLASS_ other ) {  \
        ::swap( *this, other );                    \
        return *this;                              \
    }

#define DECLARE_ASSIGNMENT_OPERATOR_TEMPLATE( _THISCLASS_, _TYPE1_ ) \
    _THISCLASS_<_TYPE1_>& operator=( _THISCLASS_<_TYPE1_> other ) {  \
        ::swap( *this, other );                                      \
        return *this;                                                \
    }

#define DECLARE_ASSIGNMENT_OPERATOR_TEMPLATE2( _THISCLASS_, _TYPE1_, \
                                               _TYPE2_ )             \
    _THISCLASS_<_TYPE1_, _TYPE2_>& operator=(                        \
        _THISCLASS_<_TYPE1_, _TYPE2_> other ) {                      \
        ::swap( *this, other );                                      \
        return *this;                                                \
    }

#define DECLARE_FRIEND_SWAP_METHOD( _THISCLASS_ ) \
    friend void ::swap( _THISCLASS_& a, _THISCLASS_& b ) noexcept;

#define DECLARE_FRIEND_SWAP_METHOD_TEMPLATE( _THISCLASS_, _TYPE1_ ) \
    friend void ::swap( _THISCLASS_<_TYPE1_>& a,                    \
                        _THISCLASS_<_TYPE1_>& b ) noexcept;

#define DECLARE_FRIEND_SWAP_METHOD_TEMPLATE2( _THISCLASS_, _TYPE1_, _TYPE2_ ) \
    friend void ::swap( _THISCLASS_<_TYPE1_, _TYPE2_>& a,                     \
                        _THISCLASS_<_TYPE1_, _TYPE2_>& b ) noexcept;

#define DECLARE_WRAPPER_BOILERPLATE_METHODS( _THISCLASS_ ) \
    DECLARE_MOVE_CONSTRUCTOR( _THISCLASS_ )                \
    DECLARE_EMPTY_DESTRUCTOR( _THISCLASS_ )                \
    DECLARE_ASSIGNMENT_OPERATOR( _THISCLASS_ )             \
    DECLARE_FRIEND_SWAP_METHOD( _THISCLASS_ )

#define DECLARE_BOILERPLATE_METHODS( _THISCLASS_, _SUPERCLASS_ ) \
    DECLARE_WRAPPER_BOILERPLATE_METHODS( _THISCLASS_ )           \
    DECLARE_CLONE_FUNCTION( _THISCLASS_, _SUPERCLASS_ )

#define DECLARE_BOILERPLATE_TEMPLATES( _THISCLASS_, _SUPERCLASS_, _TYPE1_ ) \
    DECLARE_WRAPPER_BOILERPLATE_METHODS( _THISCLASS_ )                      \
    DECLARE_CLONE_TEMPLATE( _THISCLASS_, _SUPERCLASS_, _TYPE1_ )

#define DECLARE_BOILERPLATE_TEMPLATES2( _THISCLASS_, _SUPERCLASS_, _TYPE1_, \
                                        _TYPE2_ )                           \
    DECLARE_WRAPPER_BOILERPLATE_METHODS( _THISCLASS_ )                      \
    DECLARE_CLONE_TEMPLATE2( _THISCLASS_, _SUPERCLASS_, _TYPE1_, _TYPE2_ )

#define DECLARE_WRAPPER_BOILERPLATE_TEMPLATES( _THISCLASS_, _TYPE1_ ) \
    DECLARE_MOVE_CONSTRUCTOR_TEMPLATE( _THISCLASS_, _TYPE1_ )         \
    DECLARE_EMPTY_DESTRUCTOR( _THISCLASS_ )                           \
    DECLARE_ASSIGNMENT_OPERATOR_TEMPLATE( _THISCLASS_, _TYPE1_ )      \
    DECLARE_FRIEND_SWAP_METHOD_TEMPLATE( _THISCLASS_, _TYPE1_ )

#define DECLARE_WRAPPER_BOILERPLATE_TEMPLATES2( _THISCLASS_, _TYPE1_,      \
                                                _TYPE2_ )                  \
    DECLARE_MOVE_CONSTRUCTOR_TEMPLATE2( _THISCLASS_, _TYPE1_, _TYPE2_ )    \
    DECLARE_EMPTY_DESTRUCTOR( _THISCLASS_ )                                \
    DECLARE_ASSIGNMENT_OPERATOR_TEMPLATE2( _THISCLASS_, _TYPE1_, _TYPE2_ ) \
    DECLARE_FRIEND_SWAP_METHOD_TEMPLATE2( _THISCLASS_, _TYPE1_, _TYPE2_ )
