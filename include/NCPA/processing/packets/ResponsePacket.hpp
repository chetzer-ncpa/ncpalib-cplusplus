#pragma once

#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets/Packet.hpp"

#include <memory>

void swap( NCPA::processing::ResponsePacket& a,
           NCPA::processing::ResponsePacket& b ) noexcept;

namespace NCPA::processing {
    class ResponsePacket : public Packet {
        public:
            ResponsePacket( response_id_t id ) : Packet(), _ID { id } {}

            ResponsePacket( response_id_t id, const std::string& tag ) :
                Packet( tag ), _ID { id }, _message { "" } {}

            ResponsePacket( response_id_t id, const std::string& tag,
                            const std::string& message ) :
                Packet( tag ), _ID { id }, _message { message } {}

            ResponsePacket( const ResponsePacket& other ) :
                ResponsePacket( other._ID, other._message ) {}

            virtual ~ResponsePacket() {}

            friend void ::swap( ResponsePacket& a,
                                ResponsePacket& b ) noexcept;

            response_id_t ID() const { return _ID; }

            std::string message() const { return _message; }

        private:
            response_id_t _ID;
            std::string _message;
    };
}  // namespace NCPA::processing

inline void swap( NCPA::processing::ResponsePacket& a,
                  NCPA::processing::ResponsePacket& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::processing::Packet&>( a ),
            dynamic_cast<NCPA::processing::Packet&>( b ) );
    swap( a._ID, b._ID );
    swap( a._message, b._message );
}
