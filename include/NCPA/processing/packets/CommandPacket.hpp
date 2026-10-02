#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"
#include "NCPA/processing/packets/Packet.hpp"

#include <iostream>

namespace NCPA::processing {
    class CommandPacket : public InputPacket {
        public:
            CommandPacket() : InputPacket( input_id_t::COMMAND ) {}

            CommandPacket( const std::string& tag ) :
                InputPacket( input_id_t::COMMAND, tag ) {}

            CommandPacket( const std::string& tag, int code = 0 ) :
                InputPacket( input_id_t::COMMAND, tag ), _code { code } {}

            CommandPacket( const CommandPacket& other ) :
                InputPacket( other ) {}

            CommandPacket( CommandPacket&& other ) noexcept {
                swap( *this, other );
            }

            CommandPacket& operator=( CommandPacket other ) {
                swap( *this, other );
                return *this;
            }

            virtual ~CommandPacket() {}

            friend void swap( CommandPacket& a, CommandPacket& b ) noexcept {
                using std::swap;
                swap( dynamic_cast<InputPacket&>( a ),
                      dynamic_cast<InputPacket&>( b ) );
            }

            static std::unique_ptr<InputPacket> build() {
                return std::unique_ptr<InputPacket>( new CommandPacket() );
            }

            static std::unique_ptr<InputPacket> build( const std::string& tag,
                                                       int code = 0 ) {
                return std::unique_ptr<InputPacket>(
                    new CommandPacket( tag, code ) );
            }

            static std::unique_ptr<InputPacket> build(
                const CommandPacket& other ) {
                return std::unique_ptr<InputPacket>(
                    new CommandPacket( other ) );
            }

            int code() const { return _code; }

            CommandPacket& set_code( int code ) {
                _code = code;
                return *this;
            }

        protected:
            int _code;
    };
}  // namespace NCPA::processing
