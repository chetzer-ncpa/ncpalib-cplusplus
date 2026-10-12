#pragma once

#include "NCPA/extern.hpp"
#include "NCPA/processing/AbstractProcessingStep.hpp"
#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets/RunTestsPacket.hpp"

namespace NCPA::processing {
    template<typename intype, typename outtype>
    class ProcessingChain {
        public:
            virtual ProcessingChain<intype, outtype>& define_parameters() = 0;
#if !HAVE_NLOHMANN_JSON_HPP
            virtual ProcessingChain<intype, outtype>& from_json(
                const std::string& json ) = 0;
#endif
            virtual outtype& product() = 0;

        public:
            using input_t   = intype;
            using output_t  = outtype;
            using product_t = outtype;

            ProcessingChain() {}

            ~ProcessingChain() {}

            ProcessingChain( const ProcessingChain& other ) :
                ProcessingChain() {
                _firstlink = other._firstlink;
                _input     = other._input;
            }

            ProcessingChain( ProcessingChain&& other ) noexcept :
                ProcessingChain() {
                swap( *this, other );
            }

            friend void swap( ProcessingChain<intype, outtype>& a,
                              ProcessingChain<intype, outtype>& b ) noexcept {
                using std::swap;
                swap( a._firstlink, b._firstlink );
                swap( a._input, b._input );
            }

            ProcessingChain& add_link( AbstractProcessingStep *link,
                                       const std::string& tag = "" ) {
                if (!tag.empty()) {
                    link->set_tag( tag );
                }
                if (_firstlink == nullptr) {
                    _firstlink = link;
                } else {
                    _firstlink->last()->set_next( link );
                }
                _links[ tag ] = link;
                return *this;
            }

            ProcessingChain& add_link( AbstractProcessingStep& link,
                                       const std::string& tag = "" ) {
                return this->add_link( &link, tag );
            }

            virtual std::string as_json( bool pretty   = false,
                                         size_t indent = 1, char tab = '\t' ) {
                if (pretty) {
                    return this->as_json_str_pretty( indent, tab );
                } else {
                    return this->as_json_str_plain();
                }
            }


#if HAVE_NLOHMANN_JSON_HPP
            virtual std::string as_json_str_plain() const {
                return this->build_json().dump( -1, 0, true );
            }
#else
            virtual std::string as_json_str_plain() const {
                std::ostringstream json;
                json << "{ ";
                for (auto it = this->parameters().begin();
                     it != this->parameters().end(); ++it) {
                    if (it != this->parameters().begin()) {
                        json << ", ";
                    }
                    json << ( *it )->as_json( false );
                }
                json << " }";
                return json.str();
            }
#endif

#if HAVE_NLOHMANN_JSON_HPP
            virtual std::string as_json_str_pretty( int n_indent = 1,
                                                    char tab = '\t' ) const {
                return this->build_json().dump( n_indent, tab, true );
            }
#else
            virtual std::string as_json_str_pretty( size_t n_indent = 0,
                                                    char tab = '\t' ) const {
                std::ostringstream json;
                std::string tabstr;
                if (tab == '\t') {
                    tabstr = "\t";
                } else {
                    tabstr = "    ";
                }
                std::ostringstream indentstr;
                for (size_t i = 0; i < n_indent; ++i) {
                    indentstr << tabstr;
                }
                std::string baseindent = indentstr.str();

                // n_indent = ( tab == '\t' ? 1 : 4 );
                json << baseindent << "{ " << std::endl;
                for (auto it = this->parameters().begin();
                     it != this->parameters().end(); ++it) {
                    if (it != this->parameters().begin()) {
                        json << ",\n";
                    }
                    json << ( *it )->as_json( true, n_indent + 1, tab );
                }
                json << "\n" << baseindent << "}";
                return json.str();
            }
#endif

#if HAVE_NLOHMANN_JSON_HPP
            virtual nlohmann::json build_json() const {
                nlohmann::json out;
                for (auto it = this->parameters().begin();
                     it != this->parameters().end(); ++it) {
                    out[ ( *it )->key() ] = ( *it )->as_json();
                }
                for (auto& it : _links) {
                    out[ it.first ] = it.second->as_json();
                }
                return out;
            }
#endif

            response_id_t finalize_configuration() {
                return this->process( *ConfigurationCompletePacket::build() )
                    ->ID();
            }


#if HAVE_NLOHMANN_JSON_HPP
            virtual ProcessingChain<intype, outtype>& from_json(
                nlohmann::json& json ) {
                // return this->from_json( json.dump( -1, 0, true ) );
                for (auto& linkpair : _links) {
                    linkpair.second->from_json( json.at( linkpair.first ) );
                }
                for (auto& param : _parameters) {
                    param->from_json( json.at( param->key() ) );
                }
                return *this;
            }
#endif

#if HAVE_NLOHMANN_JSON_HPP
            virtual ProcessingChain<intype, outtype>& from_json(
                const std::string& json ) {
                return from_json( nlohmann::json::parse( json ) );
            }
#else
            virtual ProcessingChain<intype, outtype>& from_json(
                const std::string& json ) = 0;
#endif

            template<typename T>
            response_ptr_t generic( T& input ) {
                GenericPacket<T> packet( input );
                return this->process( packet );
            }

            virtual Parameter *parameter( const std::string& key,
                                          bool throw_on_error = true ) {
                for (auto it = this->parameters().begin();
                     it != this->parameters().end(); ++it) {
                    if (( *it )->key() == key) {
                        return it->get();
                    }
                }
                if (throw_on_error) {
                    throw( std::out_of_range( "Parameter " + key
                                              + " not found" ) );
                } else {
                    return nullptr;
                }
            }

            virtual const Parameter *parameter( const std::string& key,
                                                bool throw_on_error
                                                = true ) const {
                for (auto it = this->parameters().begin();
                     it != this->parameters().end(); ++it) {
                    if (( *it )->key() == key) {
                        return it->get();
                    }
                }
                if (throw_on_error) {
                    throw( std::out_of_range( "Parameter " + key
                                              + " not found" ) );
                } else {
                    return nullptr;
                }
            }

            virtual std::vector<parameter_ptr_t>& parameters() {
                return _parameters;
            }

            virtual const std::vector<parameter_ptr_t>& parameters() const {
                return _parameters;
            }

            response_ptr_t process( intype& input ) {
                DataPacket<intype> packet( input );
                return this->process( packet );
            }

            virtual response_ptr_t process( InputPacket& packet ) {
                return _firstlink->process( packet );
            }

            virtual response_id_t reset() {
                return this->process( *ResetPacket::build() )->ID();
            }

            virtual bool run_tests( std::string& message,
                                    bool message_if_passed = true ) {
                bool passed = true;
                std::ostringstream oss;
                for (auto& linkpair : _links) {
                    response_ptr_t response = this->process(
                        *RunTestsPacket::build( linkpair.first ) );
                    if (response->ID() != response_id_t::TESTS_PASSED) {
                        oss << linkpair.first << ": " << response->message()
                            << "\n";
                        passed = false;
                    } else if (message_if_passed) {
                        oss << linkpair.first << ": All tests passed.\n";
                    }
                }
                message = oss.str();
                return passed;
            }

            response_id_t send_configuration( const std::string& tag,
                                              const parameter_ptr_t& param,
                                              bool throw_on_error = false ) {
                response_id_t id = this->process( *ConfigurationPacket::build(
                                                      tag, param ) )
                                       ->ID();
                if (throw_on_error
                    && id != response_id_t::CONFIGURATION_SUCCESS) {
                    throw std::runtime_error( "Configuration error for " + tag
                                              + ":" + param->key() );
                } else {
                    return id;
                }
            }

            response_id_t send_configuration( const std::string& tag,
                                              const Parameter *param,
                                              bool throw_on_error = false ) {
                return this->send_configuration( tag, param->clone(),
                                                 throw_on_error );
            }

            // virtual ProcessingChain<intype, outtype>&
            //     define_passthrough_parameter(
            //         const std::string& module_key,
            //         const std::string& param_key_in,
            //         const std::string& param_key_out ) {
            //     _passthrough.emplace_back( module_key, param_key_in,
            //                                param_key_out );
            //     return *this;
            // }

            // virtual ProcessingChain<intype, outtype>&
            //     define_passthrough_parameter( const std::string& module_key,
            //                                   const std::string& param_key ) {
            //     _passthrough.emplace_back( module_key, param_key, param_key );
            //     return *this;
            // }

            // virtual ProcessingChain<intype, outtype>& pass_parameters_through(
            //     bool throw_on_error = false ) {
            //     for (passthrough_parameter_t& param : _passthrough) {
            //         _pass_parameter_through( param, throw_on_error );
            //     }
            //     return *this;
            // }

        protected:
            virtual void _pass_parameter_through(
                const passthrough_parameter_t& param,
                bool throw_on_error = false ) {
                this->send_configuration(
                    param.tag,
                    this->parameter( param.in_key )->clone_as( param.out_key ),
                    throw_on_error );
            }


        protected:
            AbstractProcessingStep *_firstlink = nullptr;
            DataWrapper<intype> _input;
            std::vector<parameter_ptr_t> _parameters;
            std::vector<passthrough_parameter_t> _passthrough;
            std::unordered_map<std::string, AbstractProcessingStep *> _links;
    };
}  // namespace NCPA::processing
