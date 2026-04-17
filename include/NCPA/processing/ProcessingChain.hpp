#pragma once

#include "NCPA/processing/AbstractProcessingStep.hpp"
#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"

template<typename intype, typename outtype>
void swap( NCPA::processing::ProcessingChain<intype, outtype>& a,
           NCPA::processing::ProcessingChain<intype, outtype>& b ) noexcept;

namespace NCPA::processing {
    template<typename intype, typename outtype>
    class ProcessingChain {
        public:
            ProcessingChain() {}

            ~ProcessingChain() {}

            ProcessingChain( const ProcessingChain& other ) :
                ProcessingChain() {
                _firstlink = other._firstlink;
                _input     = other._input;
            }

            ProcessingChain( ProcessingChain&& other ) noexcept :
                ProcessingChain() {
                ::swap( *this, other );
            }

            friend void ::swap<>(
                ProcessingChain<intype, outtype>& a,
                ProcessingChain<intype, outtype>& b ) noexcept;

            ProcessingChain& operator=( ProcessingChain other ) {
                ::swap( *this, other );
                return *this;
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
                return *this;
            }

            ProcessingChain& add_link( AbstractProcessingStep& link,
                                       const std::string& tag = "" ) {
                return this->add_link( &link, tag );
            }

            response_ptr_t process( intype& input ) {
                // _input.set( input );
                DataPacket<intype> packet( input );
                return this->process( packet );
            }

            virtual response_ptr_t process( InputPacket& packet ) {
                return _firstlink->process( packet );
            }

            template<typename T>
            response_ptr_t generic( T& input ) {
                GenericPacket<T> packet( input );
                return this->process( packet );
            }

            response_id_t finalize_configuration() {
                return this->process( *ConfigurationCompletePacket::build() )
                    ->ID();
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

            virtual std::string as_json_str( bool pretty  = false,
                                             int n_indent = 1,
                                             char tab     = '\t' ) {
                if (pretty) {
                    return this->as_json_str_pretty( n_indent, tab );
                } else {
                    return this->as_json_str_plain();
                }
            }

            virtual outtype& product()                                    = 0;
            virtual ProcessingChain<intype, outtype>& define_parameters() = 0;
#if HAVE_NLOHMANN_JSON_HPP
            virtual nlohmann::json as_json() const {
                nlohmann::json out;
                for (auto it = this->parameters().begin();
                     it != this->parameters().end(); ++it) {
                    out[ ( *it )->key() ] = ( *it )->as_json();
                }
                return out;
            }

            virtual std::string as_json_str_pretty( int n_indent = 1,
                                                    char tab = '\t' ) const {
                return this->as_json().dump( n_indent, tab, true );
            }

            virtual std::string as_json_str_plain() const {
                return this->as_json().dump( -1, 0, true );
            }

            // template<typename T>
            // void set_from_json( nlohmann::json& json, const std::string&
            // key,
            //                     bool is_vector = false ) {
            //     set_from_json<T>( json, key, key, is_vector );
            // }

            // template<typename T>
            // void set_from_json( nlohmann::json& json,
            //                     const std::string& jsonkey,
            //                     const std::string& paramkey,
            //                     bool is_vector = false ) {
            //     if (json.at( jsonkey ).is_null()
            //         || json.at( jsonkey ).at( "value" ).is_null()) {
            //         throw std::out_of_range(
            //             "Parameter " + jsonkey
            //             + " not found or badly structured in JSON!" );
            //     }
            //     if (is_vector) {
            //         this->parameter( paramkey )
            //             ->as_vector<T>()
            //             .set( json.at( jsonkey ).at( "value" ).get<T>() );
            //     } else {
            //         this->parameter( paramkey )
            //             ->as_scalar<T>()
            //             .set( json.at( jsonkey ).at( "value" ).get<T>() );
            //     }
            // }

            virtual ProcessingChain<intype, outtype>& from_json(
                nlohmann::json& json ) = 0;
#else
            virtual std::string as_json( bool pretty            = false,
                                         size_t indent          = 1,
                                         char tab = '\t' ) {
                if (pretty) {
                    return this->as_json_str_pretty( indent, tab );
                } else {
                    return this->as_json_str_plain();
                }
            }

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

            virtual ProcessingChain<intype, outtype>& from_json(
                const std::string& json ) = 0;
#endif


        private:
            AbstractProcessingStep *_firstlink = nullptr;
            DataWrapper<intype> _input;
            std::vector<parameter_ptr_t> _parameters;
    };
}  // namespace NCPA::processing

template<typename intype, typename outtype>
void swap( NCPA::processing::ProcessingChain<intype, outtype>& a,
           NCPA::processing::ProcessingChain<intype, outtype>& b ) noexcept {
    using std::swap;
    swap( a._firstlink, b._firstlink );
    swap( a._input, b._input );
}
