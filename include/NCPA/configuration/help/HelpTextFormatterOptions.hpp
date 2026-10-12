#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/help/HelpTextMarkers.hpp"
#include "NCPA/strings.hpp"

#include <regex>
#include <sstream>
#include <string>

namespace NCPA {
    namespace config {
        using NCPA::strings::ends_with;
        using NCPA::strings::starts_with;

        struct HelpTextFormatterOptions {
                bool break_after_tag                  = false;
                std::string default_value_suffix      = "default=";
                bool double_explicit_newlines         = true;
                size_t first_line_indent              = 0;
                std::string flag_argument_indicator   = "";
                int hanging_indent                    = -1;
                size_t indent_spaces                  = 4;
                bool indent_subsections               = true;
                size_t max_tag_length_without_newline = 30;
                size_t max_width                      = 75;
                bool newline_after_title              = false;
                bool newline_before_title             = true;
                std::string newline_regex
                    = NCPA::strings::NEWLINE_REGEX_STRING;
                std::string optional_argument_indicator = "[val]";
                std::string optional_multiple_argument_indicator
                    = "[val1 ...]";
                bool print_argument_hyphens             = true;
                std::string required_argument_indicator = "val";
                std::string required_argument_suffix    = "required";
                std::string required_multiple_argument_indicator
                    = "val1 [val2 ...]";
                bool reset_indent_after_newline                = false;
                bool separate_required_from_optional_arguments = true;
                std::string suffix_end                         = "]";
                std::string suffix_start                       = "[";
                size_t title_indent                            = 0;
                std::string word_regex                         = "([^\\s]+)";

                virtual std::string build_suffix(
                    const std::string& body ) const {
                    return suffix_start + body + suffix_end;
                }

                std::vector<std::string> tokenize(
                    const std::string& input ) const {
                    // mostly this is to deal with strings already containing
                    // newlines
                    std::string replaced = std::regex_replace(
                        input, NCPA::strings::NEWLINE_REGEX,
                        HelpTextMarkers::marker( help_text_marker_t::NEWLINE,
                                                 true ) );
                    return NCPA::strings::split( replaced );
                }

                virtual std::vector<std::string> tokenize( const std::vector<std::string> sections ) const {
                    std::vector<std::string> tokens;
                    for (auto& str : sections) {
                        this->tokenize( str, tokens );
                    }
                    return tokens;
                }

                virtual void tokenize(
                    const std::string& input,
                    std::vector<std::string>& output ) const {
                    std::vector<std::string> tokens = this->tokenize( input );
                    // if (!output.empty()) {
                    //     output.push_back( HelpTextMarkers::marker(
                    //         help_text_marker_t::NEWLINE ) );
                    // }
                    output.insert( output.end(), tokens.begin(),
                                   tokens.end() );
                }
        };
    }  // namespace config
}  // namespace NCPA
