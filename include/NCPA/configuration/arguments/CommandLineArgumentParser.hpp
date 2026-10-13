#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/arguments/Argument.hpp"
#include "NCPA/configuration/arguments/ArgumentParser.hpp"
#include "NCPA/configuration/arguments/ArgumentParseResult.hpp"
#include "NCPA/configuration/arguments/ArgumentSet.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/exceptions.hpp"
#include "NCPA/strings.hpp"

#include <deque>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        class CommandLineArgumentParser
            : public Cloneable<CommandLineArgumentParser, ArgumentParser> {
                enum class parser_expecting_t { TAG, VALUE, EITHER };
                enum class token_t { NONE, TAG, VALUE };

                struct Token {
                        token_t token_type;
                        std::string text;
                };

                typedef std::pair<std::string, std::vector<std::string>>
                    arg_value_set_t;

            public:
                using parent_t
                    = Cloneable<CommandLineArgumentParser, ArgumentParser>;
                CommandLineArgumentParser() = default;

                CommandLineArgumentParser(
                    const std::vector<std::string>& tokens ) :
                    CommandLineArgumentParser() {
                    this->load( tokens );
                }

                CommandLineArgumentParser(
                    const std::deque<std::string>& tokens ) :
                    CommandLineArgumentParser() {
                    this->load( tokens );
                }

                CommandLineArgumentParser( int argc, char **argv ) :
                    CommandLineArgumentParser() {
                    this->load( argc, argv );
                }

                virtual ~CommandLineArgumentParser() = default;
                CommandLineArgumentParser(
                    const CommandLineArgumentParser& other ) = default;
                CommandLineArgumentParser(
                    CommandLineArgumentParser&& other ) noexcept = default;

                friend void swap( CommandLineArgumentParser& a,
                                  CommandLineArgumentParser& b ) noexcept {
                    using std::swap;
                    swap( static_cast<parent_t&>( a ),
                          static_cast<parent_t&>( b ) );
                }

                virtual CommandLineArgumentParser& load(
                    const std::vector<std::string>& tokens ) {
                    _tokens.assign( tokens.begin(), tokens.end() );
                    return *this;
                }

                virtual CommandLineArgumentParser& load(
                    const std::deque<std::string>& tokens ) {
                    _tokens = tokens;
                    return *this;
                }

                virtual CommandLineArgumentParser& load( int argc,
                                                         char **argv ) {
                    _tokens.clear();
                    for (int i = 1; i < argc; ++i) {
                        _tokens.push_back( argv[ i ] );
                    }
                    return *this;
                }

                virtual ArgumentParseResult parse(
                    ArgumentSet& argset ) const override {
                    if (_tokens.empty()) {
                        return ArgumentParseResult(
                            parse_result_t::NO_ARGUMENTS,
                            "No arguments found" );
                    }

                    parser_expecting_t expecting = parser_expecting_t::EITHER;
                    bool firstargs               = true;
                    bool tag_active              = false;
                    std::string active_tag;
                    std::vector<std::string> active_values;
                    Argument *active_arg = nullptr;

                    std::deque<std::string> tokens = _tokens;
                    arg_value_set_t valset = _parse_next_token( tokens );
                    if (valset.first.empty()) {
                        if (valset.second.empty()) {
                            return ArgumentParseResult(
                                parse_result_t::NO_ARGUMENTS );
                        } else {
                            return ArgumentParseResult(
                                parse_result_t::BARE_VALUES, valset.second );
                        }
                    } else if (valset.first == "h" || valset.first == "help"
                               || valset.first == "usage") {
                        return ArgumentParseResult(
                            parse_result_t::HELP_REQUESTED );
                    } else {
                        argset.argument( valset.first )
                            .from_string( valset.second );
                    }

                    while (!tokens.empty()) {
                        arg_value_set_t valset = _parse_next_token( tokens );
                        argset.argument( valset.first )
                            .from_string( valset.second );
                    }
                    return ArgumentParseResult( parse_result_t::SUCCESS );
                }


            protected:
                std::deque<std::string> _tokens;

                void _finalize_argument(
                    Argument *arg,
                    const std::vector<std::string>& values ) const {
                    if (arg != nullptr) {
                        arg->from_string( values );
                    }
                }

                arg_value_set_t _parse_next_token(
                    std::deque<std::string>& tokens ) const {
                    std::string tagstr = "";
                    std::vector<std::string> valuestrs;
                    if (tokens.empty()) {
                        return std::make_pair( tagstr, valuestrs );
                    }

                    Token t = _parse_token( tokens.front() );
                    // first time through, tag and value get saved
                    switch (t.token_type) {
                        case token_t::TAG:
                            tagstr = t.text;
                            break;
                        case token_t::VALUE:
                            valuestrs.push_back( t.text );
                            break;
                        default:
                            break;
                    }
                    tokens.pop_front();

                    // now values get saved but tags end the parse process
                    while (!tokens.empty()) {
                        t = _parse_token( tokens.front() );
                        switch (t.token_type) {
                            case token_t::TAG:
                                return std::make_pair( tagstr, valuestrs );
                            case token_t::VALUE:
                                valuestrs.push_back( t.text );
                                break;
                            default:
                                break;
                        }
                        tokens.pop_front();
                    }
                    return std::make_pair( tagstr, valuestrs );
                }

                Token _parse_token( const std::string& token ) const {
                    Token t;
                    t.text = NCPA::strings::deblank( token );
                    if (t.text.empty()) {
                        t.token_type = token_t::NONE;
                        return t;
                    } else if (t.text.size() >= 2
                               && t.text.compare( 0, 2, "--" ) == 0) {
                        t.token_type = token_t::TAG;
                        t.text       = NCPA::strings::deblank( t.text, "-" );
                    } else {
                        t.token_type = token_t::VALUE;
                    }
                    return t;
                }
        };
    }  // namespace config
}  // namespace NCPA
