#pragma once

#include <stdexcept>
#include <string>

#define NCPA_DEFINE_EXCEPTION( _EXCEPTION_NAME_, _EXCEPTION_BASE_ )       \
    class _EXCEPTION_NAME_ : public _EXCEPTION_BASE_ {                    \
        public:                                                           \
            _EXCEPTION_NAME_( const std::string& msg = "" ) :             \
                _EXCEPTION_BASE_( std::string( #_EXCEPTION_NAME_ ) + ": " \
                                  + msg ) {}                              \
    };

namespace NCPA {
    NCPA_DEFINE_EXCEPTION( NotImplementedError, std::logic_error )
    // class NotImplementedError : public std::logic_error {
    //     public:
    //         NotImplementedError( const std::string& message
    //                              = "Not implemented" ) :
    //             std::logic_error( message ) {}
    // };
}  // namespace NCPA

#ifndef THROW_METHOD_NOT_IMPLEMENTED
#  define THROW_METHOD_NOT_IMPLEMENTED                         \
      throw NCPA::NotImplementedError( std::string( __func__ ) \
                                       + ": Not implemented" );
#endif
