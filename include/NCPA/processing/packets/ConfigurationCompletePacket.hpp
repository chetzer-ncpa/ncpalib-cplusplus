#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"
#include "NCPA/processing/packets/Packet.hpp"
#include "NCPA/processing/Parameter.hpp"

void swap( NCPA::processing::ConfigurationCompletePacket& a,
           NCPA::processing::ConfigurationCompletePacket& b ) noexcept;

namespace NCPA::processing {
    class ConfigurationCompletePacket : public InputPacket {
        public:
            ConfigurationCompletePacket() :
                InputPacket( input_id_t::CONFIGURATION_COMPLETE ) {}

            ConfigurationCompletePacket( const std::string& tag ) :
                InputPacket( input_id_t::CONFIGURATION_COMPLETE, tag ) {}

            ConfigurationCompletePacket(
                const ConfigurationCompletePacket& other ) :
                InputPacket( other ) {}

            ConfigurationCompletePacket(
                ConfigurationCompletePacket&& other ) noexcept {
                ::swap( *this, other );
            }

            ConfigurationCompletePacket& operator=(
                ConfigurationCompletePacket other ) {
                ::swap( *this, other );
                return *this;
            }

            virtual ~ConfigurationCompletePacket() {}

            friend void ::swap( ConfigurationCompletePacket& a,
                                ConfigurationCompletePacket& b ) noexcept;

            static std::unique_ptr<InputPacket> build() {
                return std::unique_ptr<InputPacket>(
                    new ConfigurationCompletePacket() );
            }

            static std::unique_ptr<InputPacket> build(
                const std::string& tag ) {
                return std::unique_ptr<InputPacket>(
                    new ConfigurationCompletePacket( tag ) );
            }

            static std::unique_ptr<InputPacket> build(
                const ConfigurationCompletePacket& other ) {
                return std::unique_ptr<InputPacket>(
                    new ConfigurationCompletePacket( other ) );
            }
    };
}  // namespace NCPA::processing

void swap( NCPA::processing::ConfigurationCompletePacket& a,
           NCPA::processing::ConfigurationCompletePacket& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::processing::InputPacket&>( a ),
            dynamic_cast<NCPA::processing::InputPacket&>( b ) );
}
