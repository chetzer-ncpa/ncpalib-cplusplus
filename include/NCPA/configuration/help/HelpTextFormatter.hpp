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
                HelpTextFormatter() : HelpTextFormatter( std::cout ) {}

                HelpTextFormatter( std::ostream& stream ) :
                    _stream { &stream } {}

                HelpTextFormatter( std::ostream& stream,
                                   HelpTextFormatterOptions options ) :
                    _stream { &stream }, _default_options { options } {}

                HelpTextFormatter( const HelpTextFormatter& other ) :
                    _linepos { other._linepos },
                    _hanging { other._hanging },
                    _first_line { other._first_line },
                    _title_on { other._title_on },
                    _last_insert_whitespace { other._last_insert_whitespace },
                    _box_on { other._box_on },
                    _stream { other._stream },
                    _section { other._section },
                    _default_options { other._default_options } {
                    _buffer.str( other._buffer.str() );
                }

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
                    swap( a._linepos, b._linepos );
                    swap( a._hanging, b._hanging );
                    swap( a._first_line, b._first_line );
                    swap( a._title_on, b._title_on );
                    swap( a._last_insert_whitespace,
                          b._last_insert_whitespace );
                    swap( a._box_on, b._box_on );
                    swap( a._buffer, b._buffer );
                    swap( a._stream, b._stream );
                    swap( a._section, b._section );
                    swap( a._default_options, b._default_options );
                }

                std::string align_buffer_contents( bool title = false ) {
                    help_text_align_t align
                        = title ? options().title_align : options().text_align;
                    size_t remaining = remaining_spaces_in_line();
                    size_t before = 0, after = 0;
                    std::ostringstream output;
                    switch (align) {
                        case help_text_align_t::LEFT:
                            before = this->hanging_indent();
                            after  = _box_on ? remaining : 0;
                            break;
                        case help_text_align_t::RIGHT:
                            before = this->hanging_indent() + remaining;
                            after  = 0;
                            break;
                        case help_text_align_t::CENTER:
                            before
                                = this->hanging_indent() + ( remaining / 2 );
                            after = _box_on ? remaining - before : 0;
                            break;
                    }
                    output << std::string( before, ' ' ) << _buffer.str()
                           << std::string( after, ' ' );
                    return output.str();
                }

                size_t base_indent() const {
                    return ( _section == nullptr ? 0
                                                 : options().indent_spaces
                                                       * _section->depth() );
                }

                std::string box_left() const {
                    std::ostringstream oss;
                    oss << options().box_side;
                    if (options().box_side_pad) {
                        oss << " ";
                    }
                    return oss.str();
                }

                std::string box_lid() {
                    size_t lid_size = this->useable_space( false );
                    size_t total    = 0;
                    std::ostringstream os;
                    while (total++ < lid_size) {
                        os << options().box_top;
                    }
                    return os.str();
                }

                std::string box_right() const {
                    std::ostringstream oss;
                    if (options().box_side_pad) {
                        oss << " ";
                    }
                    oss << options().box_side;
                    return oss.str();
                }

                void end_box() {
                    _box_on = false;
                    start_line();
                    stream( box_lid() );
                    end_line();
                }

                void end_section() {
                    if (options().box_style == help_text_box_t::ON) {
                        end_box();
                    }
                }

                size_t hanging_indent() const {
                    if (_title_on) {
                        return options().title_indent;
                    } else if (_first_line) {
                        return this->options().first_line_indent;
                    } else if (options().hanging_indent < 0) {
                        return _hanging;
                    } else {
                        return options().hanging_indent;
                    }
                }

                // this assumes that we've not inserted the left indent yet
                size_t remaining_spaces_in_line() const {
                    // size_t buffer_size = _buffer.str().size();
                    size_t useable = this->useable_space();
                    size_t linepos = _buffer.str().size();
                    return ( linepos < useable ? useable - linepos : 0 );
                }

                size_t useable_space( bool include_box_sides = true ) const {
                    size_t reserved = this->left_indent()
                                    + this->hanging_indent()
                                    + this->right_indent();
                    if (include_box_sides && _box_on) {
                        reserved += this->box_left().size()
                                  + this->box_right().size();
                        if (options().box_side_pad) {
                            reserved += 2;
                        }
                    }
                    if (reserved >= options().max_width) {
                        throw FormattingError( "Not enough formatting space "
                                               "available for content!" );
                    }
                    return options().max_width - reserved;
                }

                void start_line() {
                    if (_line_started) {
                        end_line();
                    }
                    indent_line();
                    _buffer.str( "" );
                    _line_started = true;
                }

                // end line:
                // if box_on and at start of line:
                //     insert box left
                //     insert optional pad
                // align buffer contents
                // insert aligned buffer contents
                // if box_on
                //     insert trailing spaces
                //     insert optional pad
                //     insert box right
                // insert newline

                void end_line( bool hard = false ) {
                    if (!_line_started) {
                        start_line();
                    }

                    if (_box_on) {
                        stream( box_left() );
                        if (options().box_side_pad) {
                            stream( " " );
                        }
                    }
                    stream( align_buffer_contents( _title_on ) );
                    _buffer.str( "" );
                    if (_box_on) {
                        if (options().box_side_pad) {
                            stream( " " );
                        }
                        stream( box_right() );
                    }
                    stream_newline();
                    if (hard || _title_on
                        || options().reset_indent_after_newline) {
                        _first_line = true;
                    } else {
                        _first_line = false;
                    }
                }

                void format( const HelpTextSection& section ) {
                    set_section( section );
                    start_section();
                    make_title();

                    _first_line = true;
                    for (auto& token : _section->tokens( options() )) {
                        this->process_token( token );
                    }
                    end_line();
                    _section = nullptr;

                    for (auto it = section.subsections().begin();
                         it != section.subsections().end(); ++it) {
                        this->format( **it, *_stream );
                    }
                    end_section();
                }

                void format( const HelpTextSection& section,
                             std::ostream& stream ) {
                    _stream = &stream;
                    this->format( section );
                }

                void buffer( const std::string& s ) { _buffer << s; }

                void stream( const std::string& s ) { *_stream << s; }

                // void insert( const std::string& s ) {
                //     this->insert( s, _buffer );
                // }

                // void insert( const std::string& s, std::ostream& oss ) {
                //     oss << s;
                //     _linepos                += s.size();
                //     _last_insert_whitespace  = false;
                // }

                void indent_line() {
                    std::string ind( this->left_indent(), ' ' );
                    stream( ind );
                }

                void stream_newline() {
                    stream( "\n" );
                    _linepos                = 0;
                    _last_insert_whitespace = true;
                    _line_started           = false;
                }

                void buffer_space( size_t n ) {
                    n = std::min( n, remaining_spaces_in_line() );
                    for (size_t i = 0; i < n; ++i) {
                        buffer( " " );
                    }
                    _last_insert_whitespace = true;
                }

                // void left_indent( size_t& before_box,
                //                   size_t& after_box ) const {
                //     before_box = base_indent();
                //     if (_title_on) {
                //         after_box = options().title_indent;
                //     } else if (_first_line) {
                //         after_box = this->options().first_line_indent;
                //     } else {
                //         after_box = hanging_indent();
                //     }
                // }

                size_t left_indent() const { return base_indent(); }

                void make_title() {
                    if (!_section->title().empty()) {
                        _title_on = true;
                        if (!_box_on && this->options().newline_before_title) {
                            end_line();
                        }
                        if (options().box_style
                            == help_text_box_t::TITLE_ONLY) {
                            start_box();
                        }
                        std::vector<std::string> titlewords
                            = this->options().tokenize( _section->title() );

                        titlewords.push_back( HelpTextMarkers::marker(
                            help_text_marker_t::NEWLINE ) );

                        for (auto word : titlewords) {
                            this->process_token( word );
                        }
                        if (options().box_style
                            == help_text_box_t::TITLE_ONLY) {
                            end_box();
                        }
                        if (this->options().newline_after_title) {
                            end_line();
                        }
                        _title_on = false;
                    }
                }

                const HelpTextFormatterOptions& options() const {
                    return ( _section == nullptr ? _default_options
                                                 : _section->options() );
                }

                virtual bool process_marker( std::string word ) {
                    marker_contents_t contents
                        = HelpTextMarkers::parse_marker( word );
                    switch (contents.first) {
                        case help_text_marker_t::BREAK:
                            end_line( false );
                            return true;
                        case help_text_marker_t::NO_MARKER:
                            return false;
                        case help_text_marker_t::NEWLINE:
                            end_line( true );
                            return true;
                        case help_text_marker_t::SPACE:
                            buffer_space( contents.second >= 0
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

                    if (!_line_started) {
                        start_line();
                    }

                    if (word.size() > remaining_spaces_in_line()) {
                        this->end_line();
                        this->process_token( word );
                    } else {
                        if (!_last_insert_whitespace) {
                            buffer_space( 1 );
                        }
                        buffer( word );
                        _last_insert_whitespace = false;
                    }
                }

                size_t right_indent() const {
                    return options().block_quote ? base_indent() : 0;
                }

                void set_section( const HelpTextSection& section ) {
                    _section = &section;
                }

                void start_box() {
                    _box_on = false;
                    start_line();
                    stream( box_lid() );
                    end_line();
                    _box_on = true;
                }

                void start_section() {
                    _buffer.str( "" );
                    if (options().box_style == help_text_box_t::ON) {
                        start_box();
                    }
                }


            protected:
                size_t _linepos              = 0;
                size_t _hanging              = 0;
                bool _first_line             = true;
                bool _title_on               = false;
                bool _last_insert_whitespace = true;
                bool _box_on                 = false;
                bool _line_started           = false;
                std::ostringstream _buffer;

                std::ostream *_stream = &std::cout;
                HelpTextFormatterOptions _default_options;
                const HelpTextSection *_section = nullptr;
        };
    }  // namespace config
}  // namespace NCPA
