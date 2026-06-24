#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"
#include "NCPA/processing/packets/Packet.hpp"
#include "NCPA/processing/Parameter.hpp"


namespace NCPA {
    namespace processing {
        class ResetPacket : public InputPacket {
            public:
                ResetPacket() :
                    InputPacket( input_id_t::RESET ) {}

                ResetPacket( const std::string& tag ) :
                    InputPacket( input_id_t::RESET, tag ) {}

                ResetPacket( const ResetPacket& other ) :
                    InputPacket( other ) {}

                ResetPacket( ResetPacket&& other ) noexcept {
                    swap( *this, other );
                }

                ResetPacket& operator=( ResetPacket other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual ~ResetPacket() {}

                friend void swap( ResetPacket& a,
                                  ResetPacket& b ) noexcept;

                static std::unique_ptr<InputPacket> build() {
                    return std::unique_ptr<InputPacket>(
                        new ResetPacket() );
                }

                static std::unique_ptr<InputPacket> build(
                    const std::string& tag ) {
                    return std::unique_ptr<InputPacket>(
                        new ResetPacket( tag ) );
                }

                static std::unique_ptr<InputPacket> build(
                    const ResetPacket& other ) {
                    return std::unique_ptr<InputPacket>(
                        new ResetPacket( other ) );
                }
        };

        void swap( NCPA::processing::ResetPacket& a,
                   NCPA::processing::ResetPacket& b ) noexcept {
            using std::swap;
            ::swap( dynamic_cast<NCPA::processing::InputPacket&>( a ),
                    dynamic_cast<NCPA::processing::InputPacket&>( b ) );
        }

    }  // namespace processing
}  // namespace NCPA
