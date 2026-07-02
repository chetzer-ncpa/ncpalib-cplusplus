#pragma once

#include <memory>
#include <vector>

namespace NCPA {
    namespace processing {

        enum class input_id_t {
            INVALID,
            OTHER,
            DATA,
            CONFIGURATION,
            CONFIGURATION_COMPLETE,
            CONFIGURATION_QUERY,
            COMMAND,
            DATA_REQUEST,
            STATE_REQUEST,
            RESET
        };

         enum class response_id_t {
            OTHER,
            NO_RESPONSE,
            ACKNOWLEDGE,
            SUCCESS_NO_PRODUCT,
            SUCCESS_PRODUCT,
            WARNING,
            ERROR,
            ERROR_STOP,
            RECONFIGURATION_REQUESTED,
            DUMMY_CONFIGURATION,
            CONFIGURATION_SUCCESS,
            CONFIGURATION_FAILURE,
            STATE
        };

        class Packet;

        class InputPacket;
        class ConfigurationPacket;
        class ConfigurationCompletePacket;
        class ConfigurationQueryPacket;
        class DataRequestPacket;
        class ResetPacket;
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
        class WarningPacket;
        class ErrorPacket;


        typedef std::unique_ptr<InputPacket> input_ptr_t;
        typedef std::unique_ptr<ResponsePacket> response_ptr_t;
    }  // namespace processing
}  // namespace NCPA
