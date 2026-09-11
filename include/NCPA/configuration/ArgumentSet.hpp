#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/Argument.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/logging.hpp"
#include "NCPA/strings.hpp"

#include <iostream>
#include <string>
#include <unordered_map>

namespace NCPA {
    namespace config {
        class ArgumentSet : public Cloneable<ArgumentSet>,
                            public Validated<ArgumentSet> {
            public:
                ArgumentSet() {}

                ArgumentSet( const ArgumentSet& other ) {
                    for (auto& argpair : other._args) {
                        _args[ argpair.first ] = argpair.second->clone();
                    }
                }

                ArgumentSet( ArgumentSet&& other ) noexcept {
                    swap( *this, other );
                }

                ArgumentSet& operator=( ArgumentSet other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual ~ArgumentSet() {}

                void swap( ArgumentSet& a, ArgumentSet& b ) noexcept {
                    using std::swap;
                    swap( a._args, b._args );
                }

                NCPA_CLONE_METHOD( ArgumentSet, ArgumentSet )

                ArgumentSet& add( const Argument& a ) {
                    _args[ a.tag() ] = a.clone();
                    return *this;
                }

                ArgumentSet& parse( const std::vector<std::string>& args ) {
                    _expecting_t expecting    = _expecting_t::TAG;
                    Argument *current_arg_ptr = nullptr;
                    int value_count           = 0;

                    for (auto token : args) {
                        switch (_identify_token( token )) {
                            case _parsed_token_t::NONE:
                                break;

                            case _parsed_token_t::TAG:
                                // found a tag.  Are we expecting one?
                                if (expecting == _expecting_t::TAG
                                    || expecting == _expecting_t::EITHER) {
                                    std::string tagname = _strip_tag( token );
                                    auto it = _args.find( tagname );
                                    if (it != _args.end()) {
                                        current_arg_ptr = it->second.get();
                                        current_arg_ptr->mark_set();
                                        if (current_arg_ptr->expects_value()) {
                                            expecting   = _expecting_t::VALUE;
                                            value_count = 0;
                                        }
                                    } else {
                                        throw std::out_of_range(
                                            "No matching arguments for tag "
                                            + token );
                                    }
                                } else {
                                    throw ExpectedValueNotFound((
                                        current_arg_ptr == nullptr
                                            ? "Expected a value but got "
                                              "option "
                                                  + token
                                            : "Expected value for option "
                                                  + current_arg_ptr->tag()
                                                  + " but instead got option "
                                                  + token ));
                                }
                                break;

                            case _parsed_token_t::VALUE:
                                // found a value.  Are we expecting one?
                                if (expecting == _expecting_t::VALUE
                                    || expecting == _expecting_t::EITHER) {
                                    if (current_arg_ptr != nullptr) {
                                        current_arg_ptr->add_value_string(
                                            token );
                                        expecting
                                            = ( current_arg_ptr
                                                        ->can_accept_value()
                                                    ? _expecting_t::EITHER
                                                    : _expecting_t::TAG );
                                    }
                                }
                                break;

                            default:
                                throw std::logic_error(
                                    "Unrecognized or unsupported expecting "
                                    "value" );
                        }
                    }
                    return *this;
                }

                ArgumentSet& parse( int argc, char **argv,
                                    bool skipfirst = true ) {
                    int ind = ( skipfirst ? 1 : 0 );
                    if (argc <= ind) {
                        return *this;
                    }
                    std::vector<std::string> args( argc - ind );
                    for (int i = ind; i < argc; ++i) {
                        args[ i - ind ] = argv[ i ];
                    }
                    return this->parse( args );
                }

                void print_help() const {}

                virtual validation_status_t validate() const override {
                    validation_status_t status;
                    for (auto& argpair : _args) {
                        status = merge( status, argpair.second->validate() );
                    }
                    return status;
                }

                std::string value_string( const std::string& tagname,
                                          size_t n = 0 ) const {
                    auto it = _args.find( tagname );
                    if (it == _args.end()) {
                        throw ArgumentNotFound(
                            "Can't get value string from argument " + tagname
                            + ", argument not defined" );
                    } else {
                        return it->second->value_string( n );
                    }
                }

                std::vector<std::string> value_strings(
                    const std::string& tagname ) const {
                    auto it = _args.find( tagname );
                    if (it == _args.end()) {
                        throw ArgumentNotFound(
                            "Can't get value strings from argument " + tagname
                            + ", argument not defined" );
                    } else {
                        return it->second->value_strings();
                    }
                }

                bool was_set( const std::string& tagname ) const {
                    auto it = _args.find( tagname );
                    if (it == _args.end()) {
                        throw ArgumentNotFound(
                            "Can't check was_found from argument " + tagname
                            + ", argument not defined" );
                    } else {
                        return it->second->was_set();
                    }
                }


            protected:
                std::unordered_map<std::string, argument_ptr_t> _args;
                std::ostream& _outputstr = std::cout;

                enum class _parsed_token_t { NONE, TAG, VALUE };

                enum class _expecting_t { TAG, VALUE, EITHER };

                static _parsed_token_t _identify_token( std::string& token ) {
                    token = NCPA::strings::deblank( token );
                    if (token.size() == 0) {
                        return _parsed_token_t::NONE;
                    }

                    if (token[ 0 ] == '-') {
                        return _parsed_token_t::TAG;
                    }

                    return _parsed_token_t::VALUE;
                }

                static std::string _strip_tag( const std::string& token ) {
                    return NCPA::strings::deblank(
                        token, NCPA_DEFAULT_WHITESPACE "-" );
                }
        };
    }  // namespace config
}  // namespace NCPA
