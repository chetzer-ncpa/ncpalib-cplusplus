#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"
#include "NCPA/processing/packets/Packet.hpp"
#include "NCPA/processing/Parameter.hpp"
#include "NCPA/processing/ParameterTree.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

void swap( NCPA::processing::ConfigurationQueryPacket& a,
           NCPA::processing::ConfigurationQueryPacket& b ) noexcept;

namespace NCPA::processing {
    class ConfigurationQueryPacket : public InputPacket {
        public:
            ConfigurationQueryPacket() :
                InputPacket( input_id_t::CONFIGURATION_QUERY ) {}

            ConfigurationQueryPacket( const ConfigurationQueryPacket& other ) :
                InputPacket( other ), _parameters { other._parameters } {}

            virtual ~ConfigurationQueryPacket() {}

            friend void ::swap( ConfigurationQueryPacket& a,
                                ConfigurationQueryPacket& b ) noexcept;

            virtual ParameterTree& parameters() { return _parameters; }

            virtual const ParameterTree& parameters() const {
                return _parameters;
            }

            virtual void add_parameter( const std::string& tag,
                                        const Parameter& param ) {
                _parameters.add( tag, param );
            }

            static std::unique_ptr<InputPacket> build() {
                return std::unique_ptr<InputPacket>(
                    new ConfigurationQueryPacket() );
            }

            static std::unique_ptr<InputPacket> build(
                const ConfigurationQueryPacket& other ) {
                return std::unique_ptr<InputPacket>(
                    new ConfigurationQueryPacket( other ) );
            }

        private:
            ParameterTree _parameters;
    };
}  // namespace NCPA::processing

void swap( NCPA::processing::ConfigurationQueryPacket& a,
           NCPA::processing::ConfigurationQueryPacket& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::processing::InputPacket&>( a ),
            dynamic_cast<NCPA::processing::InputPacket&>( b ) );
    swap( a._parameters, b._parameters );
}
