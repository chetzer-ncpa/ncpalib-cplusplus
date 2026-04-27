#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"
#include "NCPA/processing/packets/Packet.hpp"
#include "NCPA/processing/Parameter.hpp"

// void swap( NCPA::processing::StateRequestPacket& a,
//            NCPA::processing::StateRequestPacket& b ) noexcept;

namespace NCPA {
    namespace processing {
        class StateRequestPacket : public InputPacket {
            public:
                StateRequestPacket() :
                    InputPacket( input_id_t::STATE_REQUEST ) {}

                StateRequestPacket( const std::string& tag ) :
                    InputPacket( input_id_t::STATE_REQUEST, tag ) {}

                StateRequestPacket( const StateRequestPacket& other ) :
                    InputPacket( other ) {}

                StateRequestPacket( StateRequestPacket&& other ) noexcept {
                    swap( *this, other );
                }

                StateRequestPacket& operator=( StateRequestPacket other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual ~StateRequestPacket() {}

                friend void swap( StateRequestPacket& a,
                                  StateRequestPacket& b ) noexcept;

                static std::unique_ptr<InputPacket> build() {
                    return std::unique_ptr<InputPacket>(
                        new StateRequestPacket() );
                }

                static std::unique_ptr<InputPacket> build(
                    const std::string& tag ) {
                    return std::unique_ptr<InputPacket>(
                        new StateRequestPacket( tag ) );
                }

                static std::unique_ptr<InputPacket> build(
                    const StateRequestPacket& other ) {
                    return std::unique_ptr<InputPacket>(
                        new StateRequestPacket( other ) );
                }
        };

        void swap( NCPA::processing::StateRequestPacket& a,
                   NCPA::processing::StateRequestPacket& b ) noexcept {
            using std::swap;
            swap( dynamic_cast<NCPA::processing::InputPacket&>( a ),
                  dynamic_cast<NCPA::processing::InputPacket&>( b ) );
        }

    }  // namespace processing
}  // namespace NCPA
