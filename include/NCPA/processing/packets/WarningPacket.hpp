#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets.hpp"
#include "NCPA/processing/packets/ResponsePacket.hpp"

#include <memory>

namespace NCPA {
    namespace processing {

        class WarningPacket
            : public ResponsePacket,
              public std::enable_shared_from_this<WarningPacket> {
            public:
                WarningPacket() : ResponsePacket( response_id_t::WARNING ) {}

                WarningPacket( const std::string& tag,
                               const std::string& message ) :
                    ResponsePacket( response_id_t::WARNING, tag, message ) {}

                WarningPacket( const WarningPacket& other ) :
                    ResponsePacket( other ), _ptr { other._ptr } {}

                WarningPacket( WarningPacket&& input ) noexcept :
                    WarningPacket() {
                    swap( *this, input );
                }

                virtual ~WarningPacket() {}

                friend void swap( WarningPacket& a,
                                  WarningPacket& b ) noexcept {
                    using std::swap;
                    swap( dynamic_cast<ResponsePacket&>( a ),
                          dynamic_cast<ResponsePacket&>( b ) );
                }

                WarningPacket& operator=( WarningPacket other ) {
                    swap( *this, other );
                    return *this;
                }

            private:
                const AbstractProcessingStep *_ptr;
        };

        // void swap( WarningPacket& a, WarningPacket& b ) noexcept {
        //     using std::swap;
        //     swap( dynamic_cast<ResponsePacket&>( a ),
        //             dynamic_cast<ResponsePacket&>( b ) );
        // }
    }  // namespace processing
}  // namespace NCPA
