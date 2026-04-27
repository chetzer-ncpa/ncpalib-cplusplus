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

// #include "NCPA/processing/_step_state.hpp"
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
        return packet_processing_result_t::INVALID_PACKET; \
    }

void swap( NCPA::processing::AbstractProcessingStep& a,
           NCPA::processing::AbstractProcessingStep& b ) noexcept;

namespace NCPA {
    namespace processing {

        // class AbstractProcessingStepState : public _step_state {
        //     public:
        //         AbstractProcessingStepState() : _step_state() {}
        //         virtual ~AbstractProcessingStepState() {}
        //         AbstractProcessingStepState( const AbstractProcessingStepState& other ) : _step_state( other ),
        //             _parameters{other._parameters},
        //             _tag{ other._tag },
        //             _next{ other._next },
        //             _configuration_changed{ other._configuration_changed } {}
                
            
            
        //         std::vector<parameter_ptr_t> _parameters;
        //         std::string _tag;
        //         AbstractProcessingStep *_next = nullptr;
        //         bool _configuration_changed   = true;
        // };

        class AbstractProcessingStep {
            public:
                // must be implemented
                virtual bool apply_configuration( std::string& message ) = 0;

                // virtual const _step_state& state() const { return _state; }


            protected:
                // virtual _step_state& state() { return _state; }

                // must be implemented
                /*
                Example:
                IntegerParameter dft_points( "dft_points", 0,
                    "Number of points for DFT" );
                add_parameter( dft_points );
                */
                virtual void _define_parameters()                   = 0;
                virtual packet_processing_result_t _process_input() = 0;

            // private:
            //     AbstractProcessingStepState _state;

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

                virtual response_ptr_t pass_on( InputPacket& packet,
                                                const std::string& error_msg
                                                = "" ) {
                    return this->has_next()
                             ? this->next()->process( packet )
                             : response_ptr_t( new ResponsePacket(
                                   response_id_t::ERROR, error_msg ) );
                }

                virtual response_ptr_t process( InputPacket& input ) {
                    std::string error_msg;
                    if (this->parameters().empty()) {
                        this->_define_parameters();
                    }
                    switch (input.ID()) {
                        case input_id_t::DATA:
                            return this->process_data_packet( input );
                            break;
                        case input_id_t::CONFIGURATION:
                            return this->process_configuration_packet( input );
                            break;
                        case input_id_t::CONFIGURATION_QUERY:
                            return this->process_configuration_query_packet(
                                input );
                            break;
                        case input_id_t::CONFIGURATION_COMPLETE:
                            return this->process_configuration_complete_packet(
                                input );
                            break;
                        default:
                            return this->process_other_packet( input );
                    }
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
                    switch (
                        // this->_process_configuration_packet( packet, msg ))
                        // {
                        this->_process_configuration_packet(
                            _cast_packet<ConfigurationPacket>( packet ),
                            msg )) {
                        case packet_processing_result_t::NOT_APPLICABLE:
                            if (this->has_next()) {
                                return this->next()->process( packet );
                            } else if (auto packet_ptr = dynamic_cast<
                                           const ConfigurationPacket *>(
                                           &packet )) {
                                std::ostringstream oss;
                                oss << "Configuration packet for "
                                    << packet_ptr->tag()
                                    << "::" << packet_ptr->parameter().key()
                                    << " not found!";
                                return response_ptr_t( new ResponsePacket(
                                    response_id_t::CONFIGURATION_FAILURE,
                                    oss.str() ) );
                            }
                            // fall through, just in case
                        case packet_processing_result_t::INVALID_PACKET:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR,
                                "Invalid packet passed to "
                                "process_configuration_packet().  This should "
                                "never happen and indicates a coding "
                                "error." ) );
                            break;
                        case packet_processing_result_t::SUCCESS:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::CONFIGURATION_SUCCESS ) );
                            break;
                        case packet_processing_result_t::NOT_FOUND:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::CONFIGURATION_FAILURE, msg ) );

                            break;
                        default:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR,
                                "_process_configuration_packet() returns "
                                "unsupported result code." ) );
                    }
                }

                virtual response_ptr_t process_configuration_complete_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    // switch (this->_process_configuration_complete_packet(
                    //     packet, msg )) {
                    switch (this->_process_configuration_complete_packet(
                        _cast_packet<ConfigurationCompletePacket>( packet ),
                        msg )) {
                        case packet_processing_result_t::SUCCESS:
                            this->_configuration_changed = false;
                            return this->has_next()
                                     ? this->next()->process( packet )
                                     : response_ptr_t( new ResponsePacket(
                                           response_id_t::
                                               CONFIGURATION_SUCCESS ) );
                            break;
                        case packet_processing_result_t::FAILURE:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::CONFIGURATION_FAILURE,
                                this->tag(), msg ) );
                            break;
                        case packet_processing_result_t::INVALID_PACKET:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR,
                                "Wrong kind of packet passed to "
                                "process_configuration_complete_packet()!" ) );
                            break;
                        default:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR,
                                "_process_configuration_packet() returns "
                                "unsupported result code." ) );
                    }
                }

                virtual response_ptr_t process_configuration_query_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    // switch (this->_process_configuration_query_packet(
                    // packet, msg )) {
                    switch (this->_process_configuration_query_packet(
                        _cast_packet<ConfigurationQueryPacket>( packet ),
                        msg )) {
                        case packet_processing_result_t::SUCCESS:
                        case packet_processing_result_t::NOT_APPLICABLE:
                            return this->has_next()
                                     ? this->next()->process( packet )
                                     : response_ptr_t(
                                           new DummyConfigurationPacket(
                                               packet ) );
                            break;
                        case packet_processing_result_t::INVALID_PACKET:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR,
                                "Wrong kind of packet passed to "
                                "process_configuration_query_packet()!" ) );
                            break;
                        default:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR,
                                "_process_configuration_query_packet() "
                                "returns "
                                "unsupported result code." ) );
                    }
                }

                virtual response_ptr_t process_data_request_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    switch (this->_process_data_request_packet(
                        _cast_packet<DataRequestPacket>( packet ), msg )) {
                        case packet_processing_result_t::SUCCESS:
                            return this->_build_product_packet();
                            break;
                        case packet_processing_result_t::NOT_APPLICABLE:
                            return this->pass_on(
                                packet,
                                "Data request packet reached last step "
                                "without being handled!" );
                            break;
                        case packet_processing_result_t::INVALID_PACKET:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR,
                                "Wrong kind of packet passed to "
                                "process_data_request_packet()!" ) );
                            break;
                        default:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR,
                                "_process_data_request_packet() "
                                "returns "
                                "unsupported result code." ) );
                    }
                }

                // DataPacket is type-specific so we don't define the
                // processing here, but it should call _process_data_packet()
                // and return
                virtual response_ptr_t process_data_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    switch (this->_process_data_packet( packet, msg )) {
                        case packet_processing_result_t::NO_PRODUCT:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::NO_PRODUCT ) );
                            break;
                        case packet_processing_result_t::PRODUCT:
                            return this->_build_product_packet();
                            break;
                        case packet_processing_result_t::NOT_APPLICABLE:
                            return this->pass_on(
                                packet, "Data packet reached last step "
                                        "without being handled!" );
                            break;
                            // return this->has_next()
                            //          ? this->next()->process( packet )
                            //          : response_ptr_t( new ResponsePacket(
                            //                response_id_t::ERROR,
                            //                "Generic packet reached last step
                            //                " "without being handled!" ) );
                        case packet_processing_result_t::PASS:
                            return this->pass_on(
                                *this->_build_next_input_packet(),
                                "Received continue signal but "
                                "there is no following step!" );
                            // return this->has_next()
                            //          ? this->next()->process(
                            //                *this->_build_next_input_packet()
                            //                )
                            //          : response_ptr_t( new ResponsePacket(
                            //                response_id_t::ERROR,
                            //                "Received continue signal but "
                            //                "there is no following step!" )
                            //                );
                            break;
                        default:
                            throw std::out_of_range(
                                "Unrecognized result returned from "
                                "_process_data_packet()!" );
                    }
                }

                virtual response_ptr_t process_other_packet(
                    InputPacket& packet ) {
                    std::string msg;
                    switch (this->_process_other_packet( packet, msg )) {
                        case packet_processing_result_t::SUCCESS:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::NO_PRODUCT ) );
                            break;
                        case packet_processing_result_t::FAILURE:
                            return response_ptr_t( new ResponsePacket(
                                response_id_t::ERROR, msg ) );
                            break;
                        case packet_processing_result_t::NOT_APPLICABLE:
                            return this->pass_on(
                                packet, "Other packet reached last step "
                                        "without being handled!" );
                            break;
                        default:
                            return this->has_next()
                                     ? this->next()->process( packet )
                                     : response_ptr_t( new ResponsePacket(
                                           response_id_t::
                                               UNRECOGNIZED_REQUEST ) );
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
                virtual AbstractDataWrapper& product() = 0;

            protected:
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
                    // if (packet_ptr == nullptr) {
                    //     return packet_processing_result_t::INVALID_PACKET;
                    // }
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    if (this->apply_configuration()) {
                        return packet_processing_result_t::SUCCESS;
                    } else {
                        std::ostringstream oss;
                        oss << "Configuration failure in " << this->tag();
                        message = oss.str();
                        return packet_processing_result_t::FAILURE;
                    }
                }

                virtual packet_processing_result_t
                    _process_data_request_packet(
                        const DataRequestPacket *packet_ptr,
                        std::string& message ) {
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    // if (packet_ptr == nullptr) {
                    //     return packet_processing_result_t::INVALID_PACKET;
                    // }
                    if (packet_ptr->tag() == this->tag()) {
                        return packet_processing_result_t::SUCCESS;
                    } else {
                        return packet_processing_result_t::NOT_APPLICABLE;
                    }
                }

                // virtual packet_processing_result_t
                //     _process_configuration_complete_packet(
                //         const InputPacket& packet, std::string& message ) {
                //     if (auto packet_ptr
                //         = _cast_packet<ConfigurationCompletePacket>(
                //             packet )) {
                //         if (this->apply_configuration()) {
                //             return packet_processing_result_t::SUCCESS;
                //         } else {
                //             std::ostringstream oss;
                //             oss << "Configuration failure in " <<
                //             this->tag(); message = oss.str(); return
                //             packet_processing_result_t::FAILURE;
                //         }
                //     } else {
                //         return packet_processing_result_t::INVALID_PACKET;
                //     }
                // }

                virtual packet_processing_result_t
                    _process_configuration_packet(
                        const ConfigurationPacket *packet_ptr,
                        std::string& message ) {
                    // if (packet_ptr == nullptr) {
                    //     return packet_processing_result_t::INVALID_PACKET;
                    // }
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    if (packet_ptr->tag() != this->tag()) {
                        return packet_processing_result_t::NOT_APPLICABLE;
                    }
                    for (auto it = this->parameters().begin();
                         it != this->parameters().end(); ++it) {
                        if (( *it )->key() == packet_ptr->parameter().key()) {
                            *it = packet_ptr->parameter().clone();
                            _configuration_changed = true;
                            return packet_processing_result_t::SUCCESS;
                        }
                    }

                    // no matching parameter found, return error
                    std::ostringstream oss;
                    oss << "No parameter " << packet_ptr->parameter().key()
                        << " found in " << packet_ptr->tag() << "!";
                    message = oss.str();
                    return packet_processing_result_t::NOT_FOUND;
                }

                // virtual packet_processing_result_t
                //     _process_configuration_packet( const InputPacket&
                //     packet,
                //                                    std::string& message ) {
                //     if (auto packet_ptr
                //         = dynamic_cast<const ConfigurationPacket *>(
                //             &packet )) {
                //         if (packet_ptr->tag() != this->tag()) {
                //             return
                //             packet_processing_result_t::NOT_APPLICABLE;
                //         }
                //         for (auto it = this->parameters().begin();
                //              it != this->parameters().end(); ++it) {
                //             if (( *it )->key()
                //                 == packet_ptr->parameter().key()) {
                //                 *it = packet_ptr->parameter().clone();
                //                 _configuration_changed = true;
                //                 return packet_processing_result_t::SUCCESS;
                //             }
                //         }

                //         // no matching parameter found, return error
                //         std::ostringstream oss;
                //         oss << "No parameter " <<
                //         packet_ptr->parameter().key()
                //             << " found in " << packet_ptr->tag() << "!";
                //         message = oss.str();
                //         return packet_processing_result_t::NOT_FOUND;
                //     } else {
                //         return packet_processing_result_t::INVALID_PACKET;
                //     }
                // }

                virtual packet_processing_result_t _process_other_packet(
                    InputPacket& packet, std::string& message ) {
                    return packet_processing_result_t::NOT_APPLICABLE;
                }

                virtual packet_processing_result_t
                    _process_configuration_query_packet(
                        const ConfigurationQueryPacket *packet_ptr,
                        std::string& message ) {
                    // if (packet_ptr == nullptr) {
                    //     return packet_processing_result_t::INVALID_PACKET;
                    // }
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    ParameterTree tree = packet_ptr->parameters();
                    std::vector<parameter_ptr_t> paramset;
                    if (!this->parameters().empty()) {
                        for (auto it = this->parameters().begin();
                             it != this->parameters().end(); ++it) {
                            tree.add( this->tag(), *it );
                        }
                    }
                    return packet_processing_result_t::SUCCESS;
                }

                // virtual packet_processing_result_t
                //     _process_configuration_query_packet(
                //         InputPacket& packet, std::string& message ) {
                //     if (auto packet_ptr
                //         = dynamic_cast<ConfigurationQueryPacket *>(
                //             &packet )) {
                //         ParameterTree tree = packet_ptr->parameters();
                //         std::vector<parameter_ptr_t> paramset;
                //         if (!this->parameters().empty()) {
                //             // tree[ this->tag() ] = std::vector<
                //             //     NCPA::processing::parameter_ptr_t>();
                //             for (auto it = this->parameters().begin();
                //                  it != this->parameters().end(); ++it) {
                //                 tree.add( this->tag(), *it );
                //                 // tree[ this->tag() ].push_back(
                //                 //     ( *it )->clone() );
                //             }
                //         }
                //         return packet_processing_result_t::SUCCESS;
                //     } else {
                //         return packet_processing_result_t::INVALID_PACKET;
                //     }
                // }

                // implemented in ProcessingStep
                virtual response_ptr_t _build_product_packet() const = 0;
                virtual input_ptr_t _build_next_input_packet() const = 0;
                virtual packet_processing_result_t _process_data_packet(
                    InputPacket& packet, std::string& message ) = 0;

                std::vector<parameter_ptr_t> _parameters;
                std::string _tag;
                AbstractProcessingStep *_next = nullptr;
                bool _configuration_changed   = true;
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
}

#pragma pop_macro( "CHECK_PACKET_POINTER_NOT_NULL" )
