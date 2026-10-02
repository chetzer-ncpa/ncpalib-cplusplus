#pragma once

#if __has_include( "nlohmann/json.hpp" )
#  define HAVE_NLOHMANN_JSON_HPP 1
#else
#  define HAVE_NLOHMANN_JSON_HPP 0
#endif

#include "NCPA/processing/packets/declarations.hpp"
#include "NCPA/processing/parameters/declarations.hpp"

#include <chrono>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#define CHECK_PACKET_POINTER_NOT_NULL( _PTR_ )             \
    if (_PTR_ == nullptr) {                                \
        return packet_processing_result_t::PACKET_INVALID; \
    }

namespace NCPA {
    namespace processing {

        typedef std::chrono::time_point<std::chrono::system_clock>
            time_point_t;
        typedef std::chrono::seconds duration_t;

        class TimeInterval {
            public:
                TimeInterval() {}

                TimeInterval( time_point_t t, duration_t d ) :
                    time { t }, duration { d } {}

                virtual ~TimeInterval() {}

                time_point_t time;
                duration_t duration = duration_t::zero();
        };

        typedef TimeInterval time_interval_t;

        enum class packet_processing_result_t {
            ERROR,
            PACKET_INVALID,
            PACKET_NOT_APPLICABLE,
            SUCCESS,
            SUCCESS_PRODUCT,  // success, product generated, continue if
                              // possible
            SUCCESS_NO_PRODUCT,      // success, no product generated, return
            FAILURE,
            FAILURE_PRODUCT,  // failure, product generated anyway, continue if
                              // possible
            FAILURE_NO_PRODUCT  // failure, no product generated, return
        };

        class AbstractDataWrapper;
        template<typename T>
        class DataWrapper;

        class AbstractProcessingStep;
        template<typename intype, typename outtype>
        class ProcessingStep;
        template<typename inouttype, typename state_t>
        class StatefulStep;
        template<typename inouttype>
        class BufferStep;
        template<typename inouttype>
        class PassThroughStep;
        template<typename intype, typename outtype>
        class ProcessingChain;

        typedef std::unique_ptr<AbstractProcessingStep> processing_step_ptr_t;

        std::string to_string( response_id_t t ) {
            switch (t) {
                case response_id_t::OTHER:
                    return "OTHER";
                    break;
                case response_id_t::SUCCESS_NO_PRODUCT:
                    return "SUCCESS_NO_PRODUCT";
                    break;
                case response_id_t::SUCCESS_PRODUCT:
                    return "SUCCESS_PRODUCT";
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
                default:
                    throw std::out_of_range(
                        "Unrecognized or unsupported response_id_t value!" );
            }
        }
    }  // namespace processing
}  // namespace NCPA
