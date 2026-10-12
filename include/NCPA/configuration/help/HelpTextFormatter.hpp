#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/help/HelpTextFormatterOptions.hpp"
#include "NCPA/configuration/help/HelpTextMarkers.hpp"
#include "NCPA/configuration/help/HelpTextSection.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        class HelpTextFormatter {
            public:
                HelpTextFormatter() {
                    _buffer = &std::cout;
                }

                HelpTextFormatter( std::ostream& stream ) {
                    _buffer = &stream;
                }

                HelpTextFormatter( std::ostream& stream,
                                   HelpTextFormatterOptions options ) {
                    _buffer  = &stream;
                    _options = options;
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
                                  HelpTextFormatter& b ) noexcept {
                    using std::swap;
                    swap( a._buffer, b._buffer );
                    swap( a._section, b._section );
                    swap( a._linepos, b._linepos );
                    swap( a._first_line, b._first_line );
                    swap( a._title_active, b._title_active );
                    swap( a._options, b._options );
                }

                virtual bool process_marker( std::string word ) {
                    marker_contents_t contents
                        = HelpTextMarkers::parse_marker( word );
                    switch (contents.first) {
                        case help_text_marker_t::BREAK:
                            _end_line( false );
                            return true;
                        case help_text_marker_t::NO_MARKER:
                            return false;
                        case help_text_marker_t::NEWLINE:
                            _end_line( true );
                            return true;
                        case help_text_marker_t::SPACE:
                            _insert_space( contents.second >= 0
                                               ? (size_t)contents.second
                                               : 0 );
                            return true;
                        case help_text_marker_t::SET_HANGING_INDENT:
                            _hanging = contents.second >= 0
                                         ? (size_t)contents.second
                                         : 0;
                            return true;
                        default:
                            return false;
                    }
                }

                virtual void process_token( const std::string& word ) {
                    if (this->process_marker( word )) {
                        return;
                    }

                    // std::string space = " ";
                    if (_linepos == 0) {
                        _start_line();
                        // space = "";
                    }

                    if (word.size() > _remaining_in_line()) {
                        this->_end_line();
                        this->process_token( word );
                    } else {
                        if (!_last_insert_whitespace) {
                            _insert_space();
                        }
                        _insert( word );
                    }
                }

                void set_options( const HelpTextFormatterOptions& options ) {
                    _options = options;
                }

                // virtual std::vector<std::string> tokenize(
                //     const std::string& input ) const {
                //     // mostly this is to deal with strings already
                //     containing
                //     // newlines
                //     std::string replaced = std::regex_replace(
                //         input, std::regex( this->options().newline_regex ),
                //         HelpTextMarkers::marker(
                //         help_text_marker_t::NEWLINE,
                //                                  true ) );
                //     return NCPA::strings::split( replaced );
                // }

                // virtual void tokenize(
                //     const std::string& input,
                //     std::vector<std::string>& output ) const {
                //     std::vector<std::string> tokens = this->tokenize( input
                //     ); output.insert( output.end(), tokens.begin(),
                //                    tokens.end() );
                // }

                // void write( HelpTextSection& section, std::ostream& buffer )
                // {
                //     _buffer = &buffer;
                //     this->write( section );
                // }

                void format( const HelpTextSection& section,
                             const HelpTextFormatterOptions& options,
                             std::ostream& buffer = std::cout ) {
                    set_options( options );
                    _buffer = &buffer;
                    this->format( section );
                }

                void format( const HelpTextSection& section ) {
                    _set_section( section );
                    _first_line = true;

                    if (!_section->title().empty()) {
                        _title_active = true;
                        std::vector<std::string> titlewords = this->options().tokenize( _section->title() );
                            // = _section->tokenize_title();
                        // = _section->title_tokens( _options );
                        if (this->options().newline_before_title) {
                            titlewords.insert(
                                titlewords.begin(),
                                HelpTextMarkers::marker(
                                    help_text_marker_t::NEWLINE ) );
                        }
                        titlewords.push_back( HelpTextMarkers::marker(
                            help_text_marker_t::NEWLINE ) );
                        if (this->options().newline_after_title) {
                            titlewords.push_back( HelpTextMarkers::marker(
                                help_text_marker_t::NEWLINE ) );
                        }
                        for (auto word : titlewords) {
                            this->process_token( word );
                        }
                        _title_active = false;
                    }

                    for (auto& token : _section->tokens( options() )) {
                        this->process_token( token );
                    }
                    _end_line();
                    _section = nullptr;

                    for (auto it = section.subsections().begin();
                         it != section.subsections().end(); ++it) {
                        this->format( **it );
                    }
                }

                // void write( HelpTextSection& section ) {
                //     if (_buffer == nullptr) {
                //         throw std::logic_error( "HelpTextFormatter.stream():
                //         "
                //                                 "No stream has been set!" );
                //     }

                //     _set_section( section );
                //     _first_line = true;

                //     // title
                //     if (!_section->title().empty()) {
                //         _title_active = true;
                //         std::vector<std::string> titlewords
                //             = tokenize( _section->title() );
                //         if (this->options().newline_before_title) {
                //             titlewords.insert(
                //                 titlewords.begin(),
                //                 HelpTextMarkers::marker(
                //                     help_text_marker_t::NEWLINE ) );
                //         }
                //         titlewords.push_back( HelpTextMarkers::marker(
                //             help_text_marker_t::NEWLINE ) );
                //         if (this->options().newline_after_title) {
                //             titlewords.push_back( HelpTextMarkers::marker(
                //                 help_text_marker_t::NEWLINE ) );
                //         }
                //         for (auto word : titlewords) {
                //             this->process_token( word );
                //         }
                //         _title_active = false;
                //     }

                //     for (auto& token : _section->tokens(_options)) {
                //         this->process_token( token );
                //     }
                //     _end_line();
                //     _section = nullptr;

                //     for (auto it = section.subsections().begin();
                //          it != section.subsections().end(); ++it) {
                //         this->write( **it );
                //     }
                // }

                HelpTextFormatterOptions& options() { return _options; }

                const HelpTextFormatterOptions& options() const {
                    return _options;
                }

            protected:
                size_t _linepos              = 0;
                size_t _hanging              = 0;
                bool _first_line             = true;
                bool _title_active           = false;
                bool _last_insert_whitespace = true;

                std::ostream *_buffer           = nullptr;
                const HelpTextSection *_section = nullptr;
                HelpTextFormatterOptions _options;

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

                size_t _base_indent( size_t depth ) const {
                    return options().indent_spaces * depth;
                }

                void _end_line( bool hard = false ) {
                    _insert_newline();
                    if (hard || _title_active
                        || options().reset_indent_after_newline) {
                        _first_line = true;
                    } else {
                        _first_line = false;
                    }
                }

                size_t _hanging_indent() const {
                    if (options().hanging_indent < 0) {
                        return _hanging;
                    } else {
                        return options().hanging_indent;
                    }
                }

                void _insert_indent( size_t n ) {
                    for (size_t i = 0; i < n; ++i) {
                        _insert_space();
                    }
                }

                void _insert( const std::string& s ) {
                    *_buffer << s;
                    _linepos                += s.size();
                    _last_insert_whitespace  = false;
                }

                void _insert_newline() {
                    _insert( "\n" );
                    _linepos                = 0;
                    _last_insert_whitespace = true;
                }

                void _insert_space( size_t n = 1 ) {
                    n = std::min( n, _remaining_in_line() );
                    for (size_t i = 0; i < n; ++i) {
                        _insert( " " );
                    }
                    _last_insert_whitespace = true;
                }

                size_t _remaining_in_line() const {
                    return ( _linepos <= options().max_width
                                 ? options().max_width - _linepos
                                 : 0 );
                }

                void _set_section( const HelpTextSection& section ) {
                    _section = &section;
                }

                void _start_line() {
                    if (_title_active) {
                        _insert_indent( _base_indent( _section->depth() )
                                        + options().title_indent );
                    } else if (_first_line) {
                        _insert_indent( this->_base_indent( _section->depth() )
                                        + this->options().first_line_indent );
                    } else {
                        _insert_indent( this->_base_indent( _section->depth() )
                                        + _hanging_indent() );
                    }
                }
        };
    }  // namespace config
}  // namespace NCPA
