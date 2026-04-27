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
                HelpTextFormatter() {}

                HelpTextFormatter( std::ostream& stream ) {
                    _buffer = &stream;
                }

                HelpTextFormatter( std::ostream& stream, size_t indent,
                                   size_t maxwidth ) {
                    _buffer                 = &stream;
                    options().indent_spaces = indent;
                    options().max_width     = maxwidth;
                }

                HelpTextFormatter( const HelpTextFormatter& other ) :
                    _linepos { other._linepos },
                    _buffer { other._buffer },
                    _section { other._section },
                    _first_line { other._first_line },
                    _title_active { other._title_active },
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

                virtual void insert_indent( size_t n ) {
                    for (size_t i = 0; i < n; ++i) {
                        this->insert_space();
                    }
                }

                virtual void insert( const std::string& s ) {
                    *_buffer << s;
                    _linepos                += s.size();
                    _last_insert_whitespace  = false;
                }

                virtual void insert_newline() {
                    this->insert( "\n" );
                    _linepos                = 0;
                    _last_insert_whitespace = true;
                }

                virtual void insert_space() {
                    this->insert( " " );
                    _last_insert_whitespace = true;
                }

                virtual size_t remaining_in_line() const {
                    return ( _linepos <= options().max_width
                                 ? options().max_width - _linepos
                                 : 0 );
                }

                virtual void end_line() {
                    // *_buffer << "\n";
                    this->insert_newline();
                    if (_title_active
                        || _section->options().reset_indent_after_newline) {
                        _first_line = true;
                    } else {
                        _first_line = false;
                    }
                }

                virtual size_t base_indent() const {
                    return options().indent_spaces * _section->depth();
                }

                virtual void start_line() {
                    if (_title_active) {
                        this->insert_indent(
                            this->base_indent()
                            + _section->options().title_indent );
                    } else if (_first_line) {
                        this->insert_indent(
                            this->base_indent()
                            + _section->options().first_line_indent );
                    } else {
                        this->insert_indent(
                            this->base_indent()
                            + _section->options().hanging_indent );
                    }
                }

                virtual bool process_tag( const std::string& word ) {
                    if (word.empty()) {
                        return false;
                    }
                    if (word.front() != '<' || word.back() != '>') {
                        return false;
                    }
                    if (word == NEWLINE_MARKER) {
                        this->end_line();
                        return true;
                    }
                    if (word[ 1 ] == 'i') {
                        size_t relative_indent
                            = std::stoi( word.substr( 2, word.size() - 1 ) );
                        size_t total_indent
                            = relative_indent + this->base_indent();
                        while (_linepos < this->options().max_width
                               && _linepos < total_indent) {
                            this->insert_space();
                        }
                        return true;
                    }
                    return false;
                }

                virtual void insert_word( const std::string& word ) {
                    // if (inword == NEWLINE_MARKER) {
                    //     this->end_line();
                    //     return;
                    // }
                    if (this->process_tag( word )) {
                        return;
                    }

                    // std::string space = " ";
                    if (_linepos == 0) {
                        this->start_line();
                        // space = "";
                    }

                    if (word.size() > remaining_in_line()) {
                        this->end_line();
                        this->insert_word( word );
                    } else {
                        if (!_last_insert_whitespace) {
                            this->insert_space();
                        }
                        this->insert( word );
                    }
                }

                std::string _adjust_word( const std::string& inword ) {
                    std::string word = inword;
                    size_t colonpos  = word.find_first_of( ':' );
                    if (colonpos != word.npos) {
                        size_t width = std::stoi(
                            word.substr( colonpos + 1, word.npos ) );
                        word = word.substr( colonpos )
                             + std::string( width, ' ' );
                    }
                    return word;
                }

                // virtual std::string indent( const HelpTextSection& section )
                // {
                //     return std::string(
                //         options().indent_spaces * section.depth(), ' ' );
                // }

                // size_t insert( std::ostream& buffer,
                //                const std::string& item ) const {}

                // void stream( std::ostream& buffer, const std::string& text,
                //              const std::string& base_indent,
                //              const std::string first_line_indent = "",
                //              const std::string hanging_indent    = "" )
                //              const {
                //     std::vector<std::string> words = this->split( text );
                // }

                void stream( const HelpTextSection& section,
                             std::ostream& buffer ) {
                    _buffer = &buffer;
                    this->stream( section );
                }

                void stream( const HelpTextSection& section ) {
                    if (_buffer == nullptr) {
                        throw std::logic_error( "HelpTextFormatter.stream(): "
                                                "No stream has been set!" );
                    }
                    _section    = &section;
                    _first_line = true;
                    // _section_options = _section->options();

                    // title
                    if (!section.title().empty()) {
                        _title_active = true;
                        std::vector<std::string> titlewords
                            = this->split( section.title() );
                        if (_section->options().newline_before_title) {
                            titlewords.insert( titlewords.begin(),
                                               NEWLINE_MARKER );
                        }
                        titlewords.push_back( NEWLINE_MARKER );
                        if (_section->options().newline_after_title) {
                            titlewords.push_back( NEWLINE_MARKER );
                        }
                        for (auto word : titlewords) {
                            this->insert_word( word );
                        }
                        _title_active = false;
                    }

                    std::vector<std::string> words
                        = this->split( section.text() );
                    for (auto word : words) {
                        this->insert_word( word );
                    }
                    this->end_line();
                    _section = nullptr;

                    for (auto it = section.subsections().begin();
                         it != section.subsections().end(); ++it) {
                        this->stream( **it );
                    }
                }

                // std::vector<std::string> split_title(
                //     const HelpTextSection& section ) const {
                //     std::vector<std::string> words
                //         = this->split( section.title() );
                //     if (this->options().newline_before_title) {
                //         words.insert( words.begin(),
                //                       NEWLINE_MARKER );
                //     }
                //     if (this->options().newline_after_title) {
                //         words.push_back( NEWLINE_MARKER );
                //     }
                //     return words;
                // }

                // void insert( std::ostream& buffer, const std::string& word,
                //              size_t depth ) {
                //     if (word == NEWLINE_MARKER) {
                //         this->insert_newline( buffer );
                //     } else if (_linestart) {
                //         if (word.size() + this->indent().size() > _maxwidth)
                //         {
                //             std::ostringstream oss;
                //             oss << "Word \"" << word
                //                 << "\" too long for maximum width of "
                //                 << _maxwidth;
                //             throw std::out_of_range( oss.str() );
                //         }
                //         this->insert_indent( buffer, depth );
                //         buffer << word;
                //         _linepos   += word.size();
                //         _linestart  = false;
                //     } else {
                //         if (_linepos + 1 + word.size() > _maxwidth) {
                //             this->insert_newline( buffer );
                //             this->insert( buffer, word, depth );
                //         } else {
                //             buffer << " " << word;
                //             _linepos += word.size() + 1;
                //         }
                //     }
                // }

                // void insert_indent( std::ostream& buffer, size_t depth ) {
                //     for (size_t i = 0; i < depth; ++i) {
                //         buffer << indent();
                //         _linepos += _indent;
                //     }
                // }

                // void insert_newline( std::ostream& buffer ) {
                //     buffer << "\n";
                //     _linepos   = 0;
                //     _linestart = true;
                // }

                // std::string indent( bool firstline ) const {
                //     return std::string( _indent + , ' ' );
                // }

                help_text_formatter_options_t& options() { return _options; }

                const help_text_formatter_options_t& options() const {
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
                        size_t nlpos = match_str.find( NEWLINE_MARKER );
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
                                words.push_back( NEWLINE_MARKER );
                                match_str = match_str.substr(
                                    NEWLINE_MARKER.size() );
                                // std::cout << "match_str is now \"" <<
                                // match_str << "\""
                                //           << std::endl;
                            } else {
                                words.push_back( NEWLINE_MARKER );
                                match_str = match_str.substr(
                                    NEWLINE_MARKER.size() );
                            }
                            nlpos = match_str.find( NEWLINE_MARKER );
                        }
                        if (!match_str.empty()) {
                            words.push_back( match_str );
                        }
                    }
                    return words;
                }

            protected:
                // size_t _indent;
                // size_t _maxwidth;

                size_t _linepos                 = 0;
                std::ostream *_buffer           = nullptr;
                const HelpTextSection *_section = nullptr;
                bool _first_line                = true;
                bool _title_active              = false;
                bool _last_insert_whitespace    = true;

                help_text_formatter_options_t _options;
                // help_text_section_formatter_options_t _section_options;
        };

        void swap( HelpTextFormatter& a, HelpTextFormatter& b ) noexcept {
            using std::swap;
            swap( a._buffer, b._buffer );
            swap( a._section, b._section );
            swap( a._linepos, b._linepos );
            swap( a._first_line, b._first_line );
            swap( a._title_active, b._title_active );
            swap( a._options, b._options );
        }
    }  // namespace config
}  // namespace NCPA
