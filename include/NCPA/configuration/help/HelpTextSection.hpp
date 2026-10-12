#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/exceptions.hpp"
#include "NCPA/configuration/help/HelpTextFormatterOptions.hpp"
#include "NCPA/configuration/help/HelpTextMarkers.hpp"
#include "NCPA/strings.hpp"

#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        class HelpTextSection : public CloneBase<HelpTextSection> {
            public:
                // virtual HelpTextSection& add_tokens(
                //     const std::vector<std::string>& text ) = 0;
                // virtual std::vector<std::string> title_tokens(
                //     const HelpTextFormatterOptions& options ) const = 0;

                virtual HelpTextSection& add_text( const std::string& text )
                    = 0;
                virtual std::vector<std::string> tokens(
                    const HelpTextFormatterOptions& options ) const = 0;

                HelpTextSection() : HelpTextSection( "" ) {}

                HelpTextSection( const std::string& title ) :
                    _title { title } {}

                HelpTextSection( const HelpTextSection& other ) {
                    _title = other._title;
                    _depth = other._depth;
                    for (auto it = other._subsections.begin();
                         it != other._subsections.end(); ++it) {
                        _subsections.push_back( ( *it )->clone() );
                    }
                }

                HelpTextSection( HelpTextSection&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~HelpTextSection() {}

                friend void swap( HelpTextSection& a,
                                  HelpTextSection& b ) noexcept {
                    using std::swap;
                    swap( a._title, b._title );
                    swap( a._subsections, b._subsections );
                    swap( a._depth, b._depth );
                }

                virtual HelpTextSection& add_subsection(
                    const HelpTextSection& sub ) {
                    _subsections.push_back( sub.clone() );
                    _subsections.back()->increase_depth();
                    return *_subsections.back();
                }

                virtual size_t& depth() { return _depth; }

                virtual const size_t& depth() const { return _depth; }

                virtual void increase_depth() {
                    _depth++;
                    for (auto it = _subsections.begin();
                         it != _subsections.end(); ++it) {
                        ( *it )->increase_depth();
                    }
                }

                // virtual std::string make_marker(
                //     const std::string& contents ) const {
                //     std::ostringstream oss;
                //     oss << this->options().marker_start << contents
                //         << this->options().marker_end;
                //     return oss.str();
                // }

                // virtual HelpTextFormatterOptions& options() {
                //     return _options;
                // }

                // virtual const HelpTextFormatterOptions& options() const {
                //     return _options;
                // }

                // virtual HelpTextSection& set_text(
                //     const std::vector<std::string>& text ) {
                //     this->add_tokens( text );
                //     return *this;
                // }

                // virtual HelpTextSection& set_text( const std::string& text )
                // {
                //     return this->set_text( _split_text( text ) );
                // }

                virtual HelpTextSection& set_title(
                    const std::string& title ) {
                    _title = title;
                    return *this;
                }

                virtual HelpTextSection *subsection(
                    const std::string& title ) {
                    for (auto& ptr : _subsections) {
                        if (ptr->title() == title) {
                            return ptr.get();
                        }
                    }
                    throw SubsectionNotFound( "No subsection found for title "
                                              + title );
                }

                virtual std::vector<std::unique_ptr<HelpTextSection>>&
                    subsections() {
                    return _subsections;
                }

                virtual const std::vector<std::unique_ptr<HelpTextSection>>&
                    subsections() const {
                    return _subsections;
                }

                virtual std::string title() const { return _title; }

                // virtual std::vector<std::string> tokenize(
                //     const std::string& input ) const {
                //     // mostly this is to deal with strings already containing
                //     // newlines
                //     std::string replaced = std::regex_replace(
                //         input, NCPA::strings::NEWLINE_REGEX,
                //         HelpTextMarkers::marker( help_text_marker_t::NEWLINE,
                //                                  true ) );
                //     return NCPA::strings::split( replaced );
                // }

                // virtual std::vector<std::string> tokenize_title() const {
                //     // mostly this is to deal with strings already containing
                //     // newlines
                //     return this->tokenize( this->title() );
                // }

                

                // virtual std::vector<std::string> title_tokens() const {
                //     return this->tokenize( this->title() );
                // }

                // virtual std::vector<std::string> tokenize(
                //     const std::string& input ) const {
                //     // mostly this is to deal with strings already
                //     containing
                //     // newlines
                //     std::string replaced = std::regex_replace(
                //         input, NCPA::strings::NEWLINE_REGEX,
                //         _options.marker( help_text_marker_t::NEWLINE, true )
                //         );
                //     return NCPA::strings::split( replaced );
                // }

                // virtual void tokenize(
                //     const std::string& input,
                //     std::vector<std::string>& output ) const {
                //     std::vector<std::string> tokens = this->tokenize( input
                //     ); output.insert( output.end(), tokens.begin(),
                //                    tokens.end() );
                // }


            protected:
                std::string _title;
                // HelpTextFormatterOptions _options;

                // std::vector<std::string> _tokens;
                std::vector<std::unique_ptr<HelpTextSection>> _subsections;
                size_t _depth = 0;

                // std::vector<std::string> _split_text(
                //     const std::string& str ) const {
                //     return NCPA::strings::split( str );
                // }
        };
    }  // namespace config
}  // namespace NCPA
