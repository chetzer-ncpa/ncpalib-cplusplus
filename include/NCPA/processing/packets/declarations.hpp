#pragma once

#include <memory>

namespace NCPA {
    namespace processing {
        class Packet;

        class InputPacket;
        class ConfigurationPacket;
        class ConfigurationCompletePacket;
        class ConfigurationQueryPacket;
        class DataRequestPacket;
        class StateRequestPacket;
        template<typename T>
        class DataPacket;
        template<typename T>
        class GenericPacket;

        class ResponsePacket;
        class DummyConfigurationPacket;
        class StatePacket;
        template<typename T>
        class ProductPacket;


        typedef std::unique_ptr<InputPacket> input_ptr_t;
        typedef std::unique_ptr<ResponsePacket> response_ptr_t;
    }  // namespace processing
}  // namespace NCPA
