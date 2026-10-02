#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"
#include "NCPA/processing/packets/Packet.hpp"

#include <iostream>

namespace NCPA::processing {
    class RunTestsPacket : public InputPacket {
        public:
            RunTestsPacket() : InputPacket( input_id_t::RUN_TESTS ) {}

            RunTestsPacket( const std::string& tag ) :
                InputPacket( input_id_t::RUN_TESTS, tag ) {}

            RunTestsPacket( const RunTestsPacket& other ) :
                InputPacket( other ) {}

            RunTestsPacket( RunTestsPacket&& other ) noexcept {
                swap( *this, other );
            }

            RunTestsPacket& operator=( RunTestsPacket other ) {
                swap( *this, other );
                return *this;
            }

            virtual ~RunTestsPacket() {}

            friend void swap( RunTestsPacket& a, RunTestsPacket& b ) noexcept {
                using std::swap;
                swap( dynamic_cast<InputPacket&>( a ),
                      dynamic_cast<InputPacket&>( b ) );
            }

            static std::unique_ptr<InputPacket> build() {
                return std::unique_ptr<InputPacket>( new RunTestsPacket() );
            }

            static std::unique_ptr<InputPacket> build( const std::string& tag ) {
                return std::unique_ptr<InputPacket>(
                    new RunTestsPacket( tag ) );
            }

            static std::unique_ptr<InputPacket> build(
                const RunTestsPacket& other ) {
                return std::unique_ptr<InputPacket>(
                    new RunTestsPacket( other ) );
            }
    };
}  // namespace NCPA::processing
