#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"
#include "NCPA/processing/packets/Packet.hpp"
#include "NCPA/processing/Parameter.hpp"

// void swap( NCPA::processing::ConfigurationPacket& a,
//            NCPA::processing::ConfigurationPacket& b ) noexcept;

namespace NCPA {
    namespace processing {
        class ConfigurationPacket : public InputPacket {
            public:
                ConfigurationPacket() :
                    InputPacket( input_id_t::CONFIGURATION ) {}

                ConfigurationPacket( const Parameter& param ) :
                    InputPacket( input_id_t::CONFIGURATION ),
                    _parameter { param.clone() } {}

                ConfigurationPacket( const std::string& tag,
                                     const Parameter& param ) :
                    InputPacket( input_id_t::CONFIGURATION, tag ),
                    _parameter { param.clone() } {}

                ConfigurationPacket( const parameter_ptr_t& param ) :
                    ConfigurationPacket( *param ) {}

                ConfigurationPacket( const std::string& tag,
                                     const parameter_ptr_t& param ) :
                    ConfigurationPacket( tag, *param ) {}

                ConfigurationPacket( const ConfigurationPacket& other ) :
                    InputPacket( other ),
                    _parameter { other._parameter->clone() } {}

                ConfigurationPacket( ConfigurationPacket&& other ) noexcept :
                    ConfigurationPacket() {
                    swap( *this, other );
                }

                virtual ~ConfigurationPacket() {}

                ConfigurationPacket& operator=( ConfigurationPacket other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( ConfigurationPacket& a,
                                  ConfigurationPacket& b ) noexcept;

                virtual const Parameter& parameter() const {
                    return *_parameter;
                }

                // one static build method per constructor
                static std::unique_ptr<InputPacket> build() {
                    return std::unique_ptr<InputPacket>(
                        new ConfigurationPacket() );
                }

                static std::unique_ptr<InputPacket> build(
                    const Parameter& param ) {
                    return std::unique_ptr<InputPacket>(
                        new ConfigurationPacket( param ) );
                }

                static std::unique_ptr<InputPacket> build(
                    const std::string& tag, const Parameter& param ) {
                    return std::unique_ptr<InputPacket>(
                        new ConfigurationPacket( tag, param ) );
                }

                static std::unique_ptr<InputPacket> build(
                    const parameter_ptr_t& param ) {
                    return std::unique_ptr<InputPacket>(
                        new ConfigurationPacket( param ) );
                }

                static std::unique_ptr<InputPacket> build(
                    const std::string& tag, const parameter_ptr_t& param ) {
                    return std::unique_ptr<InputPacket>(
                        new ConfigurationPacket( tag, param ) );
                }

                static std::unique_ptr<InputPacket> build(
                    const ConfigurationPacket& other ) {
                    return std::unique_ptr<InputPacket>(
                        new ConfigurationPacket( other ) );
                }

            private:
                std::unique_ptr<Parameter> _parameter;
        };

        void swap( ConfigurationPacket& a, ConfigurationPacket& b ) noexcept {
            using std::swap;
            swap( dynamic_cast<InputPacket&>( a ),
                  dynamic_cast<InputPacket&>( b ) );
            swap( a._parameter, b._parameter );
        }

    }  // namespace processing
}  // namespace NCPA
