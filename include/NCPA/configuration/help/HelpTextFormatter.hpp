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
                HelpTextFormatter() { _buffer = &std::cout; }

                HelpTextFormatter( std::ostream& stream ) {
                    _buffer = &stream;
                }

                HelpTextFormatter( std::ostream& stream,
                                   HelpTextFormatterOptions options ) {
                    _buffer = &stream;
                    // _options = options;
                }

                HelpTextFormatter( const HelpTextFormatter& other ) :
                    _linepos { other._linepos },
                    _buffer { other._buffer },
                    _section { other._section },
                    _first_line { other._first_line },
                    _title_on { other._title_on } {}

                // _options { other._options } {}

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
                    swap( a._title_on, b._title_on );
                    // swap( a._options, b._options );
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
                    _start_line();
                    *_buffer << box_lid();
                    _end_line();
                }

                void end_section() {
                    if (options().box_style == help_text_box_t::ON) {
                        end_box();
                    }
                }

                void format( const HelpTextSection& section ) {
                    _set_section( section );
                    start_section();
                    make_title();

                    _first_line = true;
                    for (auto& token : _section->tokens( options() )) {
                        this->process_token( token );
                    }
                    _end_line();
                    _section = nullptr;

                    for (auto it = section.subsections().begin();
                         it != section.subsections().end(); ++it) {
                        this->format( **it, *_buffer );
                    }
                    end_section();
                }

                void format( const HelpTextSection& section,
                             std::ostream& buffer ) {
                    _buffer = &buffer;
                    this->format( section );
                }

                void insert_left_indent() {
                    size_t b, a;
                    left_indent( b, a );
                    _insert_space( b );
                    if (_box_on) {
                        _insert( box_left() );
                    }
                    _insert_space( a );
                }

                void left_indent( size_t& before_box,
                                  size_t& after_box ) const {
                    before_box = _base_indent();
                    if (_title_on) {
                        after_box = options().title_indent;
                    } else if (_first_line) {
                        after_box = this->options().first_line_indent;
                    } else {
                        after_box = _hanging_indent();
                    }
                }

                size_t left_indent() const {
                    size_t b, a;
                    left_indent( b, a );
                    return b + a;
                }

                void make_title() {
                    if (!_section->title().empty()) {
                        _title_on = true;
                        if (this->options().newline_before_title) {
                            _insert_newline();
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
                            _insert_newline();
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

                    if (_linepos == 0) {
                        _start_line();
                    }

                    if (word.size() > _remaining_in_line()) {
                        this->_end_line();
                        this->process_token( word );
                    } else {
                        if (!_last_insert_whitespace) {
                            _insert_space( 1 );
                        }
                        _insert( word );
                    }
                }

                size_t right_indent() const {
                    return options().block_quote ? _base_indent() : 0;
                }

                void start_box() {
                    _box_on = false;
                    _start_line();
                    _insert( box_lid() );
                    _end_line();
                    _box_on = true;
                }

                void start_section() {
                    if (options().box_style == help_text_box_t::ON) {
                        start_box();
                    }
                }

                size_t useable_space( bool include_box_sides = true ) const {
                    size_t left_space  = this->left_indent();
                    size_t right_space = this->right_indent();
                    size_t box_space   = include_box_sides && _box_on
                                           ? this->box_left().size()
                                                 + this->box_right().size()
                                           : 0;
                    size_t reserved    = left_space + right_space + box_space;
                    return ( reserved >= options().max_width )
                             ? 0
                             : options().max_width - reserved;
                }


            protected:
                size_t _linepos              = 0;
                size_t _hanging              = 0;
                bool _first_line             = true;
                bool _title_on               = false;
                bool _last_insert_whitespace = true;
                bool _box_on                 = false;

                std::ostream *_buffer           = &std::cout;
                const HelpTextSection *_section = nullptr;
                HelpTextFormatterOptions _default_options;

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

                size_t _base_indent() const {
                    return ( _section == nullptr ? 0
                                                 : options().indent_spaces
                                                       * _section->depth() );
                }

                void _end_line( bool hard = false ) {
                    if (_box_on) {
                        if (_linepos == 0) {
                            _start_line();
                        }
                        _insert_space( _remaining_in_line() );
                        _insert( box_right() );
                    }
                    _insert_newline();
                    if (hard || _title_on
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
                    _insert_space( n );
                    // for (size_t i = 0; i < n; ++i) {
                    //     _insert_space();
                    // }
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

                void _insert_space( size_t n ) {
                    n = std::min( n, _remaining_in_line() );
                    for (size_t i = 0; i < n; ++i) {
                        _insert( " " );
                    }
                    _last_insert_whitespace = true;
                }

                size_t _remaining_in_line() const {
                    size_t remaining = options().max_width
                                     - ( _box_on ? box_right().size() : 0 )
                                     - right_indent();
                    return ( _linepos <= remaining ? remaining - _linepos
                                                   : 0 );
                }

                void _set_section( const HelpTextSection& section ) {
                    _section = &section;
                }

                void _start_line() {
                    // if (_title_on) {
                    //     _insert_indent( _base_indent( _section->depth() )
                    //                     + options().title_indent );
                    // } else if (_first_line) {
                    //     _insert_indent( this->_base_indent(
                    //     _section->depth() )
                    //                     + this->options().first_line_indent
                    //                     );
                    // } else {
                    //     _insert_indent( this->_base_indent(
                    //     _section->depth() )
                    //                     + _hanging_indent() );
                    // }

                    insert_left_indent();

                    // _insert_indent( left_indent() );
                    // if (_box_on) {
                    //     _insert( box_left() );
                    // }
                }
        };
    }  // namespace config
}  // namespace NCPA
