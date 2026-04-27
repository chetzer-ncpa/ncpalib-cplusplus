#pragma once

#include <memory>
#include "NCPA/processing/AbstractProcessingStep.hpp"

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
            STATE_REQUEST
        };

        enum class response_id_t {
            OTHER,
            NO_PRODUCT,
            PRODUCT,
            EARLY_PRODUCT,
            WARNING,
            ERROR,
            ERROR_STOP,
            RECONFIGURATION_REQUESTED,
            DUMMY_CONFIGURATION,
            CONFIGURATION_SUCCESS,
            CONFIGURATION_FAILURE,
            UNRECOGNIZED_REQUEST,
            STATE
        };

        class Packet;
        class InputPacket;
        class ConfigurationPacket;
        class ConfigurationCompletePacket;
        class ConfigurationQueryPacket;
        class ResponsePacket;
        class DummyConfigurationPacket;
        class DataRequestPacket;
        class StateRequestPacket;

        template<typename T>
        class DataPacket;
        template<typename T>
        class ProductPacket;
        template<typename T>
        class GenericPacket;
        // template<class STEPTYPE>
        class StatePacket;

        typedef std::unique_ptr<InputPacket> input_ptr_t;
        typedef std::unique_ptr<ResponsePacket> response_ptr_t;


    }  // namespace processing
}  // namespace NCPA
