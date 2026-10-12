#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/strings.hpp"

#include <regex>
#include <string>
#include <unordered_map>

namespace NCPA {
    namespace config {

        class HelpTextMarkers {

            public:
                static const std::unordered_map<help_text_marker_t,
                                                std::string>
                    marker_map() {
                    return {
                        {              help_text_marker_t::BREAK, "br" },
                        {            help_text_marker_t::NEWLINE, "nl" },
                        { help_text_marker_t::SET_HANGING_INDENT, "hi" },
                        {              help_text_marker_t::SPACE, "sp" }
                    };
                }

                static help_text_marker_t identify_marker( std::string& m ) {
                    return parse_marker( m ).first;
                }

                static bool is_marker( std::string& m ) {
                    return ( parse_marker( m ).first
                             != help_text_marker_t::NO_MARKER );
                }

                static std::string marker( help_text_marker_t mt,
                                           bool spaced = false ) {
                    std::ostringstream oss;
                    if (spaced) {
                        oss << " ";
                    }
                    oss << options().marker_start << marker_tag( mt )
                        << options().marker_end;
                    if (spaced) {
                        oss << " ";
                    }
                    return oss.str();
                }

                static std::string marker( help_text_marker_t mt, int payload,
                                           bool spaced = false ) {
                    std::ostringstream oss;
                    if (spaced) {
                        oss << " ";
                    }
                    oss << options().marker_start << marker_tag( mt ) << ":"
                        << payload << options().marker_end;
                    if (spaced) {
                        oss << " ";
                    }
                    return oss.str();
                }

                static std::string marker_tag( help_text_marker_t mt ) {
                    return marker_map().at( mt );
                }

                static help_text_marker_t marker_type( std::string sm ) {
                    if (sm.empty()) {
                        return help_text_marker_t::NO_MARKER;
                    }
                    for (auto& markerpair : marker_map()) {
                        if (sm == markerpair.second) {
                            return markerpair.first;
                        }
                    }
                    return help_text_marker_t::NO_MARKER;
                }

                static std::regex marker_regex() {
                    std::ostringstream oss;
                    oss << "^"
                        << NCPA::strings::escape_regex( options().marker_start )
                        << options().marker_contents_regex
                        << NCPA::strings::escape_regex( options().marker_end )
                        << "$";
                    return std::regex( oss.str() );
                }

                static const help_text_marker_options options() {
                    return help_text_marker_options();
                } 

                static marker_contents_t parse_marker( std::string& token ) {
                    // std::string t                  = token;
                    int payload = 0;
                    std::smatch matches;
                    std::regex pattern = marker_regex();
                    if (std::regex_search( token, matches, pattern )) {
                        help_text_marker_t mt
                            = marker_type( matches[ 1 ].str() );
                        payload = matches[ 2 ].matched
                                    ? std::stoi( matches[ 2 ].str() )
                                    : payload;
                        return std::make_pair( mt, payload );
                    }
                    return std::make_pair( help_text_marker_t::NO_MARKER, 0 );
                    // throw BadMarker( token + " cannot be parsed as a marker"
                    // );
                }

                static std::string strip_marker( std::string& m ) {
                    std::string front_pattern  = "^";
                    front_pattern             += NCPA::strings::escape_regex(
                        options().marker_start );
                    std::string back_pattern
                        = NCPA::strings::escape_regex( options().marker_end )
                        + "$";
                    std::string stripped = std::regex_replace(
                        m, std::regex( front_pattern ), "" );
                    stripped = std::regex_replace(
                        stripped, std::regex( back_pattern ), "" );
                    return stripped;
                }
        };
    }  // namespace config
}  // namespace NCPA
