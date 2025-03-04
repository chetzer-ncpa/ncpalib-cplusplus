#pragma once

#include <stdexcept>
#include <string>

#ifndef THROW_METHOD_NOT_IMPLEMENTED
#define THROW_METHOD_NOT_IMPLEMENTED throw NCPA::NotImplementedError( std::string(__func__) + ": Not implemented" );
#endif

namespace NCPA {
    class NotImplementedError : public std::logic_error {
        public:
            NotImplementedError( const std::string& message
                                 = "Not implemented" ) :
                std::logic_error( message ) {}
    };
}  // namespace NCPA