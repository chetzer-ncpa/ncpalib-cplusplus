#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets/ResponsePacket.hpp"
#include "NCPA/processing/Parameter.hpp"
#include "NCPA/processing/ParameterTree.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

// void swap( NCPA::processing::DummyConfigurationPacket& a,
//            NCPA::processing::DummyConfigurationPacket& b ) noexcept;

namespace NCPA {
    namespace processing {
        class DummyConfigurationPacket : public ResponsePacket {
            public:
                DummyConfigurationPacket() :
                    ResponsePacket( response_id_t::DUMMY_CONFIGURATION ) {}

                DummyConfigurationPacket(
                    const ConfigurationQueryPacket& packet ) :
                    ResponsePacket( response_id_t::DUMMY_CONFIGURATION ),
                    _parameters { packet.parameters() } {}

                DummyConfigurationPacket( const InputPacket& packet ) :
                    ResponsePacket( response_id_t::DUMMY_CONFIGURATION ) {
                    if (auto packet_ptr
                        = dynamic_cast<const ConfigurationQueryPacket *>(
                            &packet )) {
                        _parameters = packet_ptr->parameters();
                    } else {
                        throw std::runtime_error(
                            "Can't initialize DummyConfigurationPacket with "
                            "InputPacket that is not a "
                            "ConfigurationQueryPacket!" );
                    }
                }

                DummyConfigurationPacket(
                    DummyConfigurationPacket&& other ) noexcept : DummyConfigurationPacket() {
                    swap( *this, other );
                }

                DummyConfigurationPacket& operator=(
                    DummyConfigurationPacket other ) {
                    swap( *this, other );
                    return *this;
                }

                DummyConfigurationPacket( const ParameterTree& params ) :
                    ResponsePacket( response_id_t::DUMMY_CONFIGURATION ),
                    _parameters { params } {}

                virtual ~DummyConfigurationPacket() {}

                friend void swap( DummyConfigurationPacket& a,
                                  DummyConfigurationPacket& b ) noexcept;

                const ParameterTree& parameters() const { return _parameters; }

                ParameterTree& parameters() { return _parameters; }

                DummyConfigurationPacket& add_parameter(
                    const std::string& tag,
                    const std::unique_ptr<Parameter>& param ) {
                    _parameters.add( tag, param );
                    return *this;
                }

            private:
                ParameterTree _parameters;
        };

        void swap( NCPA::processing::DummyConfigurationPacket& a,
                   NCPA::processing::DummyConfigurationPacket& b ) noexcept {
            using std::swap;
            swap( dynamic_cast<NCPA::processing::ResponsePacket&>( a ),
                  dynamic_cast<NCPA::processing::ResponsePacket&>( b ) );
            swap( a._parameters, b._parameters );
        }
    }  // namespace processing
}  // namespace NCPA
