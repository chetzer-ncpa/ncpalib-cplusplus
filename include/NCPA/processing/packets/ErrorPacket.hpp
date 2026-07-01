#pragma once

#include "NCPA/processing/declarations.hpp"
// #include "NCPA/processing/packets.hpp"
#include "NCPA/processing/packets/ResponsePacket.hpp"

#include <memory>

namespace NCPA {
    namespace processing {

        class ErrorPacket : public ResponsePacket,
                            public std::enable_shared_from_this<ErrorPacket> {
            public:
                ErrorPacket( bool stop = false ) :
                    ResponsePacket(( stop ? response_id_t::ERROR_STOP
                                          : response_id_t::ERROR )) {}

                ErrorPacket( const std::string& tag,
                             const std::string& message ) :
                    ResponsePacket( response_id_t::WARNING, tag, message ) {}

                ErrorPacket( const ErrorPacket& other ) :
                    ResponsePacket( other ), _ptr { other._ptr } {}

                ErrorPacket( ErrorPacket&& input ) noexcept : ErrorPacket() {
                    swap( *this, input );
                }

                virtual ~ErrorPacket() {}

                friend void swap( ErrorPacket& a, ErrorPacket& b ) noexcept {
                    using std::swap;
                    swap( dynamic_cast<ResponsePacket&>( a ),
                          dynamic_cast<ResponsePacket&>( b ) );
                }

                ErrorPacket& operator=( ErrorPacket other ) {
                    swap( *this, other );
                    return *this;
                }

            private:
                const AbstractProcessingStep *_ptr;
        };

        // void swap( ErrorPacket& a, ErrorPacket& b ) noexcept {
        //     using std::swap;
        //     swap( dynamic_cast<ResponsePacket&>( a ),
        //             dynamic_cast<ResponsePacket&>( b ) );
        // }
    }  // namespace processing
}  // namespace NCPA
