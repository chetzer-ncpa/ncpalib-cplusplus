#pragma once

#include <memory>

namespace NCPA {
    namespace processing {
        class Packet;
        class InputPacket;
        class ConfigurationPacket;
        class ConfigurationCompletePacket;
        class ConfigurationQueryPacket;
        class ResponsePacket;
        class DummyConfigurationPacket;
        class DataRequestPacket;

        template<typename T>
        class DataPacket;
        template<typename T>
        class ProductPacket;
        template<typename T>
        class GenericPacket;

        typedef std::unique_ptr<InputPacket> input_ptr_t;
        typedef std::unique_ptr<ResponsePacket> response_ptr_t;
    }
}