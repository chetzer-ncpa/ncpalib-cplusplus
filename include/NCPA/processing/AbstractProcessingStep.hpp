#pragma once

/*
Normal Workflow:
ConfigurationPacket
ConfigurationPacket
...
ConfigurationCompletePacket
DataPacket
DataPacket
...
*/

#include "NCPA/processing/AbstractDataWrapper.hpp"
#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets.hpp"
#if HAVE_NLOHMANN_JSON_HPP
#  include "nlohmann/json.hpp"
#endif

#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#pragma push_macro( "CHECK_PACKET_POINTER_NOT_NULL" )
#undef CHECK_PACKET_POINTER_NOT_NULL
#define CHECK_PACKET_POINTER_NOT_NULL( _PTR_ )             \
    if (_PTR_ == nullptr) {                                \
        return packet_processing_result_t::PACKET_INVALID; \
    }

void swap( NCPA::processing::AbstractProcessingStep& a,
           NCPA::processing::AbstractProcessingStep& b ) noexcept;

namespace NCPA {
    namespace processing {
        class AbstractProcessingStep {
            public:
                // must be implemented
                virtual bool apply_configuration( std::string& message ) = 0;
                virtual bool product_available() const                   = 0;

            protected:
                // must be implemented

                /*
                Example:
                IntegerParameter dft_points( "dft_points", 0,
                    "Number of points for DFT" );
                add_parameter( dft_points );
                */
                virtual void _define_parameters()                   = 0;
                virtual packet_processing_result_t _process_input() = 0;

            public:
                AbstractProcessingStep() : AbstractProcessingStep( "" ) {}

                AbstractProcessingStep( const std::string& tag ) :
                    _tag { tag } {}

                virtual ~AbstractProcessingStep() {}

                AbstractProcessingStep( const AbstractProcessingStep& other ) :
                    AbstractProcessingStep() {
                    _next                  = other._next;
                    _tag                   = other._tag;
                    _configuration_changed = other._configuration_changed;
                    for (auto it = other._parameters.begin();
                         it != other._parameters.end(); ++it) {
                        _parameters.push_back( ( *it )->clone() );
                    }
                }

                AbstractProcessingStep(
                    AbstractProcessingStep&& other ) noexcept :
                    AbstractProcessingStep() {
                    ::swap( *this, other );
                }

                friend void ::swap( AbstractProcessingStep& a,
                                    AbstractProcessingStep& b ) noexcept;

                virtual AbstractProcessingStep& add_flag(
                    const std::string& flag ) {
                    _flags.push_back( flag );
                    return *this;
                }

                virtual AbstractProcessingStep& add_parameter(
                    Parameter& param ) {
                    _parameters.emplace_back( param.clone() );
                    // _parameters.emplace(param.key(), param.clone());
                    // _parameters[ param.key() ] = std::move(param.clone());
                    return *this;
                }

                virtual bool apply_configuration() {
                    std::string msg;
                    return apply_configuration( msg );
                }

                virtual bool apply_parameter(
                    const NCPA::processing::Parameter& param ) {
                    if (this->parameters().empty()) {
                        this->_define_parameters();
                    }

                    for (auto it = this->parameters().begin();
                         it != this->parameters().end(); ++it) {
                        if (( *it )->key() == param.key()) {
                            *it                    = param.clone();
                            _configuration_changed = true;
                            return true;
                        }
                    }
                    return false;
                }

                virtual bool has_next() const {
                    return ( this->next() != nullptr );
                }

                virtual AbstractProcessingStep *last() {
                    return _next == nullptr ? this : _next->last();
                }

                virtual AbstractProcessingStep *next() { return _next; }

                virtual const AbstractProcessingStep *next() const {
                    return _next;
                }

                // virtual std::unordered_map<std::string, parameter_ptr_t>&
                virtual std::vector<parameter_ptr_t>& parameters() {
                    return _parameters;
                }

                // virtual const std::unordered_map<std::string,
                // parameter_ptr_t>&
                virtual const std::vector<parameter_ptr_t>& parameters()
                    const {
                    return _parameters;
                }

                virtual Parameter& parameter( const std::string& key ) {
                    for (auto it = _parameters.begin();
                         it != _parameters.end(); ++it) {
                        if (( *it )->key() == key) {
                            return **it;
                        }
                    }
                    throw std::out_of_range( "Key " + key + " not found!" );
                }

                virtual response_ptr_t pass_to_next(
                    InputPacket& packet, response_ptr_t response_if_no_next ) {
                    return this->has_next() ? this->next()->process( packet )
                                            : std::move( response_if_no_next );
                }

                virtual response_ptr_t pass_to_next( InputPacket& packet ) {
                    return this->pass_to_next( packet,
                                               _build_product_packet() );
                }

                virtual response_ptr_t process( InputPacket& input ) {
                    if (this->parameters().empty()) {
                        this->_define_parameters();
                    }
                    response_ptr_t resp;
                    switch (input.ID()) {
                        case input_id_t::DATA:
                            resp = this->process_data_packet( input );
                            break;
                        case input_id_t::CONFIGURATION:
                            resp = this->process_configuration_packet( input );
                            break;
                        case input_id_t::CONFIGURATION_QUERY:
                            resp = this->process_configuration_query_packet(
                                input );
                            break;
                        case input_id_t::CONFIGURATION_COMPLETE:
                            resp = this->process_configuration_complete_packet(
                                input );
                            break;
                        case input_id_t::DATA_REQUEST:
                            resp = this->process_data_request_packet( input );
                            break;
                        case input_id_t::RESET:
                            resp = this->process_reset_packet( input );
                            break;
                        default:
                            resp = this->process_other_packet( input );
                    }
                    this->_add_flags_to_packet( resp );
                    return std::move( resp );
                }

                virtual response_ptr_t process(
                    std::unique_ptr<InputPacket>& packet_ptr ) {
                    return this->process( *packet_ptr );
                }

                virtual response_ptr_t process( InputPacket *packet_ptr ) {
                    return this->process( *packet_ptr );
                }

                virtual response_ptr_t process_configuration_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    if (this->parameters().size() == 0) {
                        this->_define_parameters();
                    }
                    switch (this->_process_configuration_packet(
                        _cast_packet<ConfigurationPacket>( packet ), msg )) {
                        case packet_processing_result_t::PACKET_NOT_APPLICABLE:
                            return this->pass_to_next(
                                packet,
                                this->response(
                                    response_id_t::CONFIGURATION_FAILURE,
                                    "No matching tag found for "
                                        + packet.tag() ) );
                            break;
                        case packet_processing_result_t::PACKET_INVALID:
                            return packet_invalid(
                                "process_configuration_packet" );
                            break;
                        case packet_processing_result_t::SUCCESS_NO_PRODUCT:
                        case packet_processing_result_t::SUCCESS_PRODUCT:
                            return this->response(
                                response_id_t::CONFIGURATION_SUCCESS );
                            break;
                        case packet_processing_result_t::FAILURE_NO_PRODUCT:
                        case packet_processing_result_t::FAILURE_PRODUCT:
                            return this->response(
                                response_id_t::CONFIGURATION_FAILURE, msg );
                            break;
                        default:
                            return return_code_unsupported(
                                "process_configuration_packet" );
                    }
                }

                virtual response_ptr_t response( response_id_t resptype,
                                                 const std::string& msg
                                                 = "" ) const {
                    return response_ptr_t(
                        new ResponsePacket( resptype, this->tag(), msg ) );
                }

                virtual response_ptr_t packet_invalid(
                    const std::string& method ) const {
                    return this->response(
                        response_id_t::ERROR,
                        "Invalid packet passed to " + method
                            + ". This should never happen and "
                              "indicates a coding error." );
                }

                // virtual void throw_packet_invalid(
                //     const std::string& method ) const {
                //     throw std::out_of_range( "Invalid packet passed to "
                //                             + method
                //                             + ". This should never happen
                //                             and "
                //                               "indicates a coding error." );
                // }

                virtual response_ptr_t return_code_unsupported(
                    const std::string& method ) const {
                    return this->response(
                        response_id_t::ERROR,
                        "Unrecognized or unsupported return code returned to "
                            + method
                            + ". This should never happen and "
                              "indicates a coding error." );
                }

                // virtual void throw_unsupported_return_code(
                //     const std::string& method ) const {
                //     throw std::out_of_range(
                //         "Unrecognized or unsupported return code returned to
                //         "
                //         + method
                //         + ". This should never happen and "
                //           "indicates a coding error." );
                // }

                virtual response_ptr_t process_configuration_complete_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    switch (this->_process_configuration_complete_packet(
                        _cast_packet<ConfigurationCompletePacket>( packet ),
                        msg )) {
                        case packet_processing_result_t::SUCCESS_NO_PRODUCT:
                        case packet_processing_result_t::SUCCESS_PRODUCT:
                            this->_configuration_changed = false;
                            return this->pass_to_next(
                                packet,
                                this->response(
                                    response_id_t::CONFIGURATION_SUCCESS ) );
                            break;
                        case packet_processing_result_t::FAILURE_NO_PRODUCT:
                        case packet_processing_result_t::FAILURE_PRODUCT:
                            return this->response(
                                response_id_t::CONFIGURATION_FAILURE, msg );
                            break;
                        case packet_processing_result_t::PACKET_INVALID:
                            return packet_invalid(
                                "process_configuration_complete_packet" );
                            break;
                        default:
                            return this->response(
                                response_id_t::ERROR,
                                "_process_configuration_packet() returns "
                                "unsupported result code." );
                    }
                }

                virtual response_ptr_t process_reset_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    switch (this->_process_reset_packet(
                        _cast_packet<ResetPacket>( packet ), msg )) {
                        case packet_processing_result_t::ERROR:
                            return this->response(
                                response_id_t::ERROR,
                                "Error processing reset packet" );
                            break;
                        case packet_processing_result_t::PACKET_INVALID:
                            return packet_invalid( "process_reset_packet" );
                            break;
                        default:
                            return this->pass_to_next(
                                packet,
                                this->response( response_id_t::ACKNOWLEDGE ) );
                    }
                }

                virtual response_ptr_t process_configuration_query_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    switch (this->_process_configuration_query_packet(
                        _cast_packet<ConfigurationQueryPacket>( packet ),
                        msg )) {
                        case packet_processing_result_t::SUCCESS_PRODUCT:
                        case packet_processing_result_t::SUCCESS_NO_PRODUCT:
                        case packet_processing_result_t::PACKET_NOT_APPLICABLE:
                            return this->pass_to_next(
                                packet,
                                response_ptr_t(
                                    new DummyConfigurationPacket( packet ) ) );
                            break;
                        case packet_processing_result_t::FAILURE_NO_PRODUCT:
                        case packet_processing_result_t::FAILURE_PRODUCT:
                            return this->response(
                                response_id_t::ERROR,
                                "Error processing configuration query "
                                "packet" );
                            break;
                        case packet_processing_result_t::PACKET_INVALID:
                            return packet_invalid(
                                "process_configuration_query_packet" );
                            break;
                        default:
                            return this->response(
                                response_id_t::ERROR,
                                "_process_configuration_query_packet() "
                                "returns unsupported result code." );
                    }
                }

                virtual response_ptr_t process_data_request_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    packet_processing_result_t result
                        = this->_process_data_request_packet(
                            _cast_packet<DataRequestPacket>( packet ), msg );
                    switch (result) {
                        case packet_processing_result_t::SUCCESS_PRODUCT:
                        case packet_processing_result_t::FAILURE_PRODUCT:

                            return this->_build_product_packet();
                            break;
                        case packet_processing_result_t::FAILURE_NO_PRODUCT:
                        case packet_processing_result_t::SUCCESS_NO_PRODUCT:
                            return this->response( response_id_t::ERROR,
                                                   "No product available" );
                        case packet_processing_result_t::PACKET_NOT_APPLICABLE:
                            return this->pass_to_next(
                                packet,
                                this->response(
                                    response_id_t::ERROR,
                                    "Data request packet reached last step "
                                    "without being handled!" ) );
                            break;
                        case packet_processing_result_t::PACKET_INVALID:
                            return packet_invalid(
                                "process_data_request_packet" );
                            break;
                        default:
                            return this->response(
                                response_id_t::ERROR,
                                "_process_data_request_packet() "
                                "returns unsupported result code." );
                    }
                }

                virtual response_ptr_t process_state_request_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    switch (this->_process_state_request_packet(
                        _cast_packet<StateRequestPacket>( packet ), msg )) {
                        case packet_processing_result_t::SUCCESS_PRODUCT:
                            return this->_build_state_packet();
                            break;
                        case packet_processing_result_t::PACKET_NOT_APPLICABLE:
                            return this->pass_to_next(
                                packet,
                                this->response(
                                    response_id_t::ERROR,
                                    "State request packet reached last step "
                                    "without being handled!" ) );
                            break;
                        case packet_processing_result_t::SUCCESS_NO_PRODUCT:
                            return this->response(
                                response_id_t::ERROR,
                                "State request returned successful but no "
                                "state to return!" );
                        case packet_processing_result_t::PACKET_INVALID:
                            return packet_invalid(
                                "process_state_request_packet" );
                            break;
                        default:
                            return this->response(
                                response_id_t::ERROR,
                                "_process_state_request_packet() "
                                "returns "
                                "unsupported result code." );
                    }
                }

                // DataPacket is type-specific so we don't define the
                // processing here, but it should call _process_data_packet()
                // and return
                virtual response_ptr_t process_data_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    switch (this->_process_data_packet( packet, msg )) {
                        case packet_processing_result_t::PACKET_NOT_APPLICABLE:
                            return this->pass_to_next(
                                packet, this->response(
                                            response_id_t::ERROR,
                                            "Data packet reached last step "
                                            "without being handled!" ) );
                            break;
                        case packet_processing_result_t::ERROR:
                        case packet_processing_result_t::FAILURE_NO_PRODUCT:
                        case packet_processing_result_t::FAILURE_PRODUCT:
                            return this->response( response_id_t::ERROR, msg );
                            break;
                        case packet_processing_result_t::SUCCESS_NO_PRODUCT:
                            return this->response(
                                response_id_t::SUCCESS_NO_PRODUCT );
                            break;
                        case packet_processing_result_t::SUCCESS_PRODUCT:
                            return this->pass_to_next(
                                *this->_build_next_input_packet(),
                                this->_build_product_packet() );
                            break;
                        default:
                            return return_code_unsupported(
                                "process_data_packet" );
                    }
                }

                virtual response_ptr_t process_other_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    switch (this->_process_other_packet( packet, msg )) {
                        case packet_processing_result_t::PACKET_NOT_APPLICABLE:
                            return this->pass_to_next(
                                packet,
                                this->response( response_id_t::ERROR,
                                                "Packet reached last step "
                                                "without being handled!" ) );
                            break;
                        case packet_processing_result_t::SUCCESS_NO_PRODUCT:
                            return this->response(
                                response_id_t::SUCCESS_NO_PRODUCT );
                            break;
                        case packet_processing_result_t::SUCCESS_PRODUCT:
                            return this->_build_product_packet();
                            break;
                        case packet_processing_result_t::FAILURE_NO_PRODUCT:
                        case packet_processing_result_t::FAILURE_PRODUCT:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR, msg ) );
                            break;
                        default:
                            return return_code_unsupported(
                                "process_other_packet" );
                    }
                }

                virtual AbstractProcessingStep& set_next(
                    AbstractProcessingStep *nextstep ) {
                    _next = nextstep;
                    return *this;
                }

                virtual AbstractProcessingStep& set_tag(
                    const std::string& tag ) {
                    _tag = tag;
                    return *this;
                }

                virtual std::string tag() const { return _tag; }

                // implemented in ProcessingStep
                virtual AbstractDataWrapper& product()  = 0;
                virtual AbstractProcessingStep& reset() = 0;

            protected:
                virtual void _add_flags_to_packet(
                    const response_ptr_t& resp ) {
                    for (auto flag : _flags) {
                        resp->prepend_flag( flag );
                    }
                }

                template<typename T>
                const T *_cast_packet( const InputPacket& packet ) const {
                    if (const T *packet_ptr
                        = dynamic_cast<const T *>( &packet )) {
                        return packet_ptr;
                    } else {
                        return nullptr;
                    }
                }

                virtual packet_processing_result_t
                    _process_configuration_complete_packet(
                        const ConfigurationCompletePacket *packet_ptr,
                        std::string& message ) {
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    if (this->apply_configuration()) {
                        return packet_processing_result_t::SUCCESS_NO_PRODUCT;
                    } else {
                        std::ostringstream oss;
                        oss << "Configuration failure in " << this->tag();
                        message = oss.str();
                        return packet_processing_result_t::FAILURE_NO_PRODUCT;
                    }
                }

                virtual packet_processing_result_t
                    _process_data_request_packet(
                        const DataRequestPacket *packet_ptr,
                        std::string& message ) {
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    if (packet_ptr->tag() == this->tag()) {
                        return (
                            this->product_available()
                                ? packet_processing_result_t::SUCCESS_PRODUCT
                                : packet_processing_result_t::
                                      FAILURE_NO_PRODUCT );
                    } else {
                        return packet_processing_result_t::
                            PACKET_NOT_APPLICABLE;
                    }
                }

                virtual packet_processing_result_t
                    _process_state_request_packet(
                        const StateRequestPacket *packet_ptr,
                        std::string& message ) {
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    if (packet_ptr->tag() == this->tag()) {
                        return packet_processing_result_t::SUCCESS_PRODUCT;
                    } else {
                        return packet_processing_result_t::
                            PACKET_NOT_APPLICABLE;
                    }
                }

                virtual packet_processing_result_t
                    _process_configuration_packet(
                        const ConfigurationPacket *packet_ptr,
                        std::string& message ) {
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    if (packet_ptr->tag() != this->tag()) {
                        return packet_processing_result_t::
                            PACKET_NOT_APPLICABLE;
                    }
                    for (auto it = this->parameters().begin();
                         it != this->parameters().end(); ++it) {
                        if (( *it )->key() == packet_ptr->parameter().key()) {
                            *it = packet_ptr->parameter().clone();
                            _configuration_changed = true;
                            return packet_processing_result_t::
                                SUCCESS_NO_PRODUCT;
                        }
                    }

                    // no matching parameter found, return error
                    std::ostringstream oss;
                    oss << "No parameter " << packet_ptr->parameter().key()
                        << " found in " << packet_ptr->tag() << "!";
                    message = oss.str();
                    return packet_processing_result_t::FAILURE_NO_PRODUCT;
                }

                virtual packet_processing_result_t _process_other_packet(
                    InputPacket& packet, std::string& message ) {
                    return packet_processing_result_t::PACKET_NOT_APPLICABLE;
                }

                virtual packet_processing_result_t _process_reset_packet(
                    const ResetPacket *packet_ptr, std::string& message ) {
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr );
                    if (packet_ptr->tag() == ""
                        || packet_ptr->tag() == this->tag()) {
                        this->reset();
                        return packet_processing_result_t::SUCCESS;
                    } else {
                        return packet_processing_result_t::
                            PACKET_NOT_APPLICABLE;
                    }
                }

                virtual packet_processing_result_t
                    _process_configuration_query_packet(
                        const ConfigurationQueryPacket *packet_ptr,
                        std::string& message ) {
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    ParameterTree tree = packet_ptr->parameters();
                    std::vector<parameter_ptr_t> paramset;
                    if (!this->parameters().empty()) {
                        for (auto it = this->parameters().begin();
                             it != this->parameters().end(); ++it) {
                            tree.add( this->tag(), *it );
                        }
                    }
                    return packet_processing_result_t::SUCCESS_PRODUCT;
                }

                virtual response_ptr_t _build_error_packet(
                    const std::string& msg ) const {
                    return response_ptr_t(
                        new ErrorPacket( this->tag(), msg ) );
                }

                virtual response_ptr_t _build_state_packet() const {
                    return response_ptr_t( new StatePacket( this ) );
                }

                // implemented in ProcessingStep
                virtual response_ptr_t _build_product_packet() const = 0;
                virtual input_ptr_t _build_next_input_packet() const = 0;
                virtual packet_processing_result_t _process_data_packet(
                    InputPacket& packet, std::string& message ) = 0;

                std::vector<parameter_ptr_t> _parameters;
                std::string _tag;
                AbstractProcessingStep *_next = nullptr;
                bool _configuration_changed   = true;
                bool _treat_as_last           = false;
                std::vector<std::string> _flags;
        };
    }  // namespace processing
}  // namespace NCPA

void swap( NCPA::processing::AbstractProcessingStep& a,
           NCPA::processing::AbstractProcessingStep& b ) noexcept {
    using std::swap;
    swap( a._next, b._next );
    swap( a._parameters, b._parameters );
    swap( a._tag, b._tag );
    swap( a._configuration_changed, b._configuration_changed );
    swap( a._treat_as_last, b._treat_as_last );
    swap( a._flags, b._flags );
}

#pragma pop_macro( "CHECK_PACKET_POINTER_NOT_NULL" )
