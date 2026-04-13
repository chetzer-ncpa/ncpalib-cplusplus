#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/HelpTextSection.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        class HelpTextFormatter {
            public:
                HelpTextFormatter() : HelpTextFormatter( 2, 80 ) {}

                HelpTextFormatter( size_t indent, size_t maxwidth ) :
                    _indent { indent },
                    _maxwidth { maxwidth },
                    _linepos { 0 },
                    _linestart { true } {}

                HelpTextFormatter( const HelpTextFormatter& other ) :
                    _indent { other._indent },
                    _maxwidth { other._maxwidth },
                    _linepos { other._linepos },
                    _linestart { other._linestart },
                    _options { other._options } {}

                HelpTextFormatter( HelpTextFormatter&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~HelpTextFormatter() {}

                HelpTextFormatter& operator=( HelpTextFormatter other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( HelpTextFormatter& a,
                                  HelpTextFormatter& b ) noexcept;

                void format( const HelpTextSection& section,
                             std::ostream& buffer, size_t depth = 0 ) {
                    std::vector<std::string> words;
                    _linestart = true;
                    _linepos   = 0;
                    if (!section.title().empty()) {
                        // std::vector<std::string> titlewords =
                        // this->split(section.title());
                        std::vector<std::string> titlewords
                            = this->split_title( section );
                        for (auto it = titlewords.begin();
                             it != titlewords.end(); ++it) {
                            this->insert( buffer, *it, depth );
                        }
                        this->insert_newline( buffer );
                    }
                    if (!section.text().empty()) {
                        std::vector<std::string> textwords
                            = this->split( section.text() );
                        for (auto it = textwords.begin();
                             it != textwords.end(); ++it) {
                            this->insert( buffer, *it, depth );
                        }
                        this->insert_newline( buffer );
                    }
                    for (auto it = section.subsections().begin();
                         it != section.subsections().end(); ++it) {
                        this->format( **it, buffer, depth + 1 );
                    }
                }

                std::vector<std::string> split_title(
                    const HelpTextSection& section ) const {
                    std::vector<std::string> words
                        = this->split( section.title() );
                    if (this->options().newline_before_title) {
                        words.insert( words.begin(),
                                      this->options().newline_marker );
                    }
                    if (this->options().newline_after_title) {
                        words.push_back( this->options().newline_marker );
                    }
                    return words;
                }

                void insert( std::ostream& buffer, const std::string& word,
                             size_t depth ) {
                    if (word == this->options().newline_marker) {
                        this->insert_newline( buffer );
                    } else if (_linestart) {
                        if (word.size() + this->indent().size() > _maxwidth) {
                            std::ostringstream oss;
                            oss << "Word \"" << word
                                << "\" too long for maximum width of "
                                << _maxwidth;
                            throw std::out_of_range( oss.str() );
                        }
                        this->insert_indent( buffer, depth );
                        buffer << word;
                        _linepos   += word.size();
                        _linestart  = false;
                    } else {
                        if (_linepos + 1 + word.size() > _maxwidth) {
                            this->insert_newline( buffer );
                            this->insert( buffer, word, depth );
                        } else {
                            buffer << " " << word;
                            _linepos += word.size() + 1;
                        }
                    }
                }

                void insert_indent( std::ostream& buffer, size_t depth ) {
                    for (size_t i = 0; i < depth; ++i) {
                        buffer << indent();
                        _linepos += _indent;
                    }
                }

                void insert_newline( std::ostream& buffer ) {
                    buffer << "\n";
                    _linepos   = 0;
                    _linestart = true;
                }

                std::string indent() const {
                    return std::string( _indent, ' ' );
                }

                HelpTextFormatterOptions& options() { return _options; }

                const HelpTextFormatterOptions& options() const {
                    return _options;
                }

                virtual std::vector<std::string> split(
                    const std::string& contents ) const {
                    std::vector<std::string> words;
                    // std::string contents = this->text();

                    // first, split into words
                    std::regex word_regex( this->options().word_regex );
                    auto words_begin = std::sregex_iterator(
                        contents.begin(), contents.end(), word_regex );
                    auto words_end = std::sregex_iterator();

                    for (std::sregex_iterator i = words_begin; i != words_end;
                         ++i) {
                        std::smatch match     = *i;
                        std::string match_str = match.str();
                        size_t strbegin       = 0;
                        size_t nlpos
                            = match_str.find( this->options().newline_marker );
                        while (nlpos != std::string::npos) {
                            if (nlpos > 0) {
                                // std::cout << "Pushing \"" <<
                                // match_str.substr(0, nlpos)
                                //           << "\"" << std::endl;
                                words.push_back(
                                    match_str.substr( 0, nlpos ) );
                                match_str = match_str.substr( nlpos );
                                // std::cout << "match_str is now \"" <<
                                // match_str << "\""
                                //           << std::endl;
                                words.push_back(
                                    this->options().newline_marker );
                                match_str = match_str.substr(
                                    this->options().newline_marker.size() );
                                // std::cout << "match_str is now \"" <<
                                // match_str << "\""
                                //           << std::endl;
                            } else {
                                words.push_back(
                                    this->options().newline_marker );
                                match_str = match_str.substr(
                                    this->options().newline_marker.size() );
                            }
                            nlpos = match_str.find(
                                this->options().newline_marker );
                        }
                        if (!match_str.empty()) {
                            words.push_back( match_str );
                        }
                    }
                    return words;
                }

            protected:
                size_t _indent;
                size_t _maxwidth;

                size_t _linepos;
                bool _linestart;
                HelpTextFormatterOptions _options;
        };

        void swap( HelpTextFormatter& a, HelpTextFormatter& b ) noexcept {
            using std::swap;
            swap( a._indent, b._indent );
            swap( a._maxwidth, b._maxwidth );
            swap( a._linepos, b._linepos );
            swap( a._linestart, b._linestart );
            swap( a._options, b._options );
        }
    }  // namespace config
}  // namespace NCPA
