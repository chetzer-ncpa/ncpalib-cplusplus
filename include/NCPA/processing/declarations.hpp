#pragma once

#include <chrono>
#include <memory>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#if __has_include( "nlohmann/json.hpp" )
#  define HAVE_NLOHMANN_JSON_HPP 1
#else
#  define HAVE_NLOHMANN_JSON_HPP 0
#endif

namespace NCPA {
    namespace processing {

        typedef std::chrono::time_point<std::chrono::system_clock>
            time_point_t;
        typedef std::chrono::seconds duration_t;

        typedef struct {
                time_point_t time;
                duration_t duration = duration_t::zero();
        } time_interval_t;

        enum class input_id_t {
            INVALID,
            OTHER,
            DATA,
            CONFIGURATION,
            CONFIGURATION_COMPLETE,
            CONFIGURATION_QUERY,
            COMMAND,
            DATA_REQUEST
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
            UNRECOGNIZED_REQUEST
        };


        enum class packet_processing_result_t {
            NOT_APPLICABLE,
            PRODUCT,
            NO_PRODUCT,
            PASS,
            SUCCESS,
            FAILURE,
            NOT_FOUND,
            INVALID_PACKET,
            ERROR
        };

        enum class parameter_type_t { INTEGER, FLOAT, STRING, BOOLEAN, ENUM };

        enum class parameter_form_t { SCALAR, ARRAY };

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

        class Parameter;
        template<typename T>
        class ScalarParameter;
        template<typename T>
        class VectorParameter;

        class IntegerParameter;
        class DoubleParameter;
        class StringParameter;
        class BooleanParameter;
        class IntegerVectorParameter;
        class DoubleVectorParameter;
        class StringVectorParameter;
        class BooleanVectorParameter;

        class ParameterTree;

        typedef std::unique_ptr<Parameter> parameter_ptr_t;
        typedef std::unordered_map<std::string, std::vector<parameter_ptr_t>>
            parameter_tree_t;

        class AbstractDataWrapper;
        template<typename T>
        class DataWrapper;

        class AbstractProcessingStep;
        template<typename intype, typename outtype>
        class ProcessingStep;
        template<typename intype, typename outtype>
        class ProcessingChain;

        std::string to_string( response_id_t t ) {
            switch (t) {
                case response_id_t::OTHER:
                    return "OTHER";
                    break;
                case response_id_t::NO_PRODUCT:
                    return "NO_PRODUCT";
                    break;
                case response_id_t::PRODUCT:
                    return "PRODUCT";
                    break;
                case response_id_t::EARLY_PRODUCT:
                    return "EARLY_PRODUCT";
                    break;
                case response_id_t::WARNING:
                    return "WARNING";
                    break;
                case response_id_t::ERROR:
                    return "ERROR";
                    break;
                case response_id_t::ERROR_STOP:
                    return "ERROR_STOP";
                    break;
                case response_id_t::RECONFIGURATION_REQUESTED:
                    return "RECONFIGURATION_REQUESTED";
                    break;
                case response_id_t::DUMMY_CONFIGURATION:
                    return "DUMMY_CONFIGURATION";
                    break;
                case response_id_t::CONFIGURATION_SUCCESS:
                    return "CONFIGURATION_SUCCESS";
                    break;
                case response_id_t::CONFIGURATION_FAILURE:
                    return "CONFIGURATION_FAILURE";
                    break;
                case response_id_t::UNRECOGNIZED_REQUEST:
                    return "UNRECOGNIZED_REQUEST";
                    break;
                default:
                    throw std::out_of_range(
                        "Unrecognized or unsupported response_id_t value!" );
            }
        }
    }  // namespace processing
}  // namespace NCPA
