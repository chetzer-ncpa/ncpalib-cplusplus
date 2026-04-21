#pragma once

#include "NCPA/configuration/Argument.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {

        template<typename INTYPE>
        class TypedArgument : public Argument {
            public:
                TypedArgument() : Argument() {}

                TypedArgument( const std::string& tag ) : Argument( tag ) {}

                TypedArgument( const std::string& tag,
                               const INTYPE& defaultval ) :
                    Argument( tag ),
                    _default { defaultval },
                    _value { defaultval } {}

                TypedArgument(
                    const std::string& tag,
                    const mapping_ptr_t<INTYPE, std::string>& mapping ) :
                    TypedArgument<INTYPE>( tag, *mapping ) {}

                TypedArgument(
                    const std::string& tag, const INTYPE& defaultval,
                    const mapping_ptr_t<INTYPE, std::string>& mapping ) :
                    TypedArgument<INTYPE>( tag, *mapping ),
                    _default { defaultval },
                    _value { defaultval } {}

                TypedArgument( const std::string& tag,
                               const Mapping<INTYPE, std::string>& mapping ) :
                    Argument( tag ) {
                    _mappings.push_back( mapping.clone() );
                }

                TypedArgument( const std::string& tag,
                               const INTYPE& defaultval,
                               const Mapping<INTYPE, std::string>& mapping ) :
                    Argument( tag ),
                    _default { defaultval },
                    _value { defaultval } {
                    _mappings.push_back( mapping.clone() );
                }

                virtual ~TypedArgument() {}

                virtual void apply() override {
                    for (auto ptr : _mappings) {
                        ptr->apply( _value );
                    }
                }

                virtual parse_status_t parse(
                    std::vector<std::string>& args ) override {
                    for (auto it = args.begin(); it != args.end(); ++it) {
                        std::string item
                            = it->substr( it->find_first_not_of( "-" ) );
                        if (item.size() > 0 && item == this->tag()) {
                            if (this->expects_value()) {
                                if (( it + 1 ) == args.end()) {
                                    throw std::out_of_range(
                                        "Argument " + item
                                        + " expects a value but none was "
                                          "found (reached end of argument "
                                          "list)!" );
                                } else if (( it + 1 )->at( 0 ) == '-') {
                                    throw std::out_of_range(
                                        "Argument " + item
                                        + " expects a value but none was "
                                          "found (next argument "
                                        + *( it + 1 )
                                        + " appears to be a flag)!" );
                                }
                                this->set( parse_value( *( it + 1 ) ) );
                                args.erase( it, it + 2 );
                            } else {
                                this->set();
                                args.erase( it );
                            }
                            return parse_status_t::SUCCESS;
                        }
                    }
                    if (this->required()) {
                        return parse_status_t::FAILURE;
                    } else {
                        return parse_status_t::NOT_FOUND;
                    }
                }

                virtual static INTYPE parse_value( const std::string& input ) {
                    INTYPE val;
                    parse_string( input, val );
                    return val;
                }

                virtual TypedArgument& set() {
                    throw std::out_of_range( "Cannot set argument "
                                             + this->tag()
                                             + " to an undefined value!" );
                }

                virtual TypedArgument& set( const INTYPE& val ) {
                    _value   = val;
                    _was_set = true;
                    return *this;
                }

                virtual bool was_set() const override { return _was_set; }

            protected:
                std::vector<mapping_ptr_t<INTYPE, std::string>> _mappings;
                INTYPE _default;
                INTYPE _value;

                bool _was_set = false;
        };
    }  // namespace config
}  // namespace NCPA
