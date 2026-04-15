#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"
#include "NCPA/processing/packets/Packet.hpp"
#include "NCPA/processing/Parameter.hpp"

void swap( NCPA::processing::DataRequestPacket& a,
           NCPA::processing::DataRequestPacket& b ) noexcept;

namespace NCPA::processing {
    class DataRequestPacket : public InputPacket {
        public:
            DataRequestPacket() :
                InputPacket( input_id_t::DATA_REQUEST ) {}

            DataRequestPacket( const std::string& tag ) :
                InputPacket( input_id_t::DATA_REQUEST, tag ) {}

            DataRequestPacket(
                const DataRequestPacket& other ) :
                InputPacket( other ) {}

            DataRequestPacket(
                DataRequestPacket&& other ) noexcept {
                ::swap( *this, other );
            }

            DataRequestPacket& operator=(
                DataRequestPacket other ) {
                ::swap( *this, other );
                return *this;
            }

            virtual ~DataRequestPacket() {}

            friend void ::swap( DataRequestPacket& a,
                                DataRequestPacket& b ) noexcept;

            static std::unique_ptr<InputPacket> build() {
                return std::unique_ptr<InputPacket>(
                    new DataRequestPacket() );
            }

            static std::unique_ptr<InputPacket> build(
                const std::string& tag ) {
                return std::unique_ptr<InputPacket>(
                    new DataRequestPacket( tag ) );
            }

            static std::unique_ptr<InputPacket> build(
                const DataRequestPacket& other ) {
                return std::unique_ptr<InputPacket>(
                    new DataRequestPacket( other ) );
            }
    };
}  // namespace NCPA::processing

void swap( NCPA::processing::DataRequestPacket& a,
           NCPA::processing::DataRequestPacket& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::processing::InputPacket&>( a ),
            dynamic_cast<NCPA::processing::InputPacket&>( b ) );
}
