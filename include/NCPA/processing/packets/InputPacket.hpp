#pragma once

#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets/Packet.hpp"

#include <memory>

void swap( NCPA::processing::InputPacket& a,
           NCPA::processing::InputPacket& b ) noexcept;

namespace NCPA::processing {
    class InputPacket : public Packet {
        public:
            InputPacket() : Packet(), _ID { input_id_t::INVALID } {}

            InputPacket( input_id_t id ) : Packet(), _ID { id } {}

            InputPacket( input_id_t id, const std::string& tag ) : Packet( tag ), _ID { id } {}            

            InputPacket( const InputPacket& other ) : Packet(other), _ID{ other._ID } {} 

            InputPacket( InputPacket&& other ) noexcept {
                ::swap( *this, other );
            }

            virtual ~InputPacket() {}

            friend void ::swap( InputPacket& a, InputPacket& b ) noexcept;

            virtual input_id_t ID() const { return _ID; }

        private:
            input_id_t _ID;
    };
}  // namespace NCPA::processing

inline void swap( NCPA::processing::InputPacket& a,
           NCPA::processing::InputPacket& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::processing::Packet&>( a ),
          dynamic_cast<NCPA::processing::Packet&>( b ) );
    swap( a._ID, b._ID );
}
