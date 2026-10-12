#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/arguments/Argument.hpp"
#include "NCPA/configuration/arguments/ArgumentSet.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/help/HelpTextSection.hpp"
#include "NCPA/logging.hpp"

#include <iostream>

namespace NCPA {
    namespace config {
        class HelpTextArgumentSection
            : public Cloneable<HelpTextArgumentSection, HelpTextSection> {
                using cloneable_t
                    = Cloneable<HelpTextArgumentSection, HelpTextSection>;

            public:
                HelpTextArgumentSection() : cloneable_t() {}

                HelpTextArgumentSection( const std::string& title ) :
                    cloneable_t( title ) {}

                HelpTextArgumentSection( const std::string& title,
                                         const Argument& arg ) :
                    cloneable_t( title ) {
                    this->add_argument( arg );
                }

                HelpTextArgumentSection( const std::string& title,
                                         const Argument *arg ) :
                    cloneable_t( title ) {
                    this->add_argument( arg );
                }

                HelpTextArgumentSection( const std::string& title,
                                         const argument_ptr_t& arg ) :
                    cloneable_t( title ) {
                    this->add_argument( arg );
                }

                HelpTextArgumentSection( const std::string& title,
                                         const ArgumentSet& argset ) :
                    cloneable_t( title ) {
                    this->add_arguments( argset );
                }

                HelpTextArgumentSection( const std::string& title,
                                         const ArgumentSet& argset,
                                         const std::string& tag ) :
                    cloneable_t( title ) {
                    this->add_argument( argset, tag );
                }

                HelpTextArgumentSection(
                    const std::string& title, const ArgumentSet& argset,
                    const std::vector<std::string>& tags ) :
                    cloneable_t( title ) {
                    this->add_arguments( argset, tags );
                }

                HelpTextArgumentSection(
                    const HelpTextArgumentSection& other ) :
                    cloneable_t( other ) {
                    _arguments.reserve( other._arguments.size() );
                    for (auto it = other._arguments.begin();
                         it != other._arguments.end(); ++it) {
                        _arguments.push_back( ( *it )->clone() );
                    }
                }

                HelpTextArgumentSection(
                    HelpTextArgumentSection&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~HelpTextArgumentSection() {}

                HelpTextArgumentSection& operator=(
                    HelpTextArgumentSection other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( HelpTextArgumentSection& a,
                                  HelpTextArgumentSection& b ) noexcept {
                    using std::swap;
                    swap( static_cast<cloneable_t&>( a ),
                          static_cast<cloneable_t&>( b ) );
                    swap( a._arguments, b._arguments );
                }

                virtual void add_argument( const Argument& arg ) {
                    _arguments.push_back( arg.clone() );
                    // std::string flagname  = arg.tag();
                    // std::string indicator = _indicator( arg );
                    // _hanging_indent       = std::max(
                    //     _hanging_indent,
                    //     std::min(
                    //         flagname.size()
                    //             + ( indicator.empty() ? 1
                    //                                   : indicator.size() + 2
                    //                                   ),
                    //         this->options().max_tag_length_without_newline )
                    //         );
                    // NCPA_DEBUG << "Setting hanging indent to "
                    //            << _hanging_indent << std::endl;
                }

                virtual void add_argument( const Argument *arg ) {
                    add_argument( *arg );
                }

                virtual void add_argument(
                    const std::unique_ptr<Argument>& arg ) {
                    add_argument( *arg );
                }

                virtual void add_argument( const ArgumentSet& argset,
                                           const std::string& tag ) {
                    this->add_argument( argset.argument( tag ) );
                }

                virtual void add_arguments( const ArgumentSet& argset ) {
                    for (auto& argpair : argset.arguments()) {
                        this->add_argument( argpair.second );
                    }
                }

                virtual void add_arguments(
                    const ArgumentSet& argset,
                    const std::vector<std::string>& tags ) {
                    for (auto& tag : tags) {
                        this->add_argument( argset.argument( tag ) );
                    }
                }

                // virtual HelpTextSection& add_tokens(
                //     const std::vector<std::string>& text ) override {
                //     _preface.insert( _preface.end(), text.begin(),
                //                      text.end() );
                //     return *this;
                // }

                // virtual HelpTextArgumentSection& set_text(
                //     const std::string& text ) override {
                //     _preface = tokenize( text );
                //     return *this;
                // }

                virtual HelpTextArgumentSection& add_text(
                    const std::string& text ) override {
                    _preface.push_back( text );
                    return *this;
                }

                virtual std::vector<std::string> tokens(
                    const HelpTextFormatterOptions& options ) const override {
                    std::vector<std::string> tok;

                    if (!_preface.empty()) {
                        tok = options.tokenize( _preface );
                    }

                    size_t hanging = _calculate_hanging_indent( options );
                    tok.push_back( HelpTextMarkers::marker(
                        help_text_marker_t::SET_HANGING_INDENT,
                        (int)hanging ) );
                    if (options.separate_required_from_optional_arguments) {
                        for (auto it = _arguments.begin();
                             it != _arguments.end(); ++it) {
                            if (( *it )->required()) {
                                _tokenize_argument_help_text(
                                    tok, **it, options, hanging );
                            }
                        }
                        for (auto it = _arguments.begin();
                             it != _arguments.end(); ++it) {
                            if (!( *it )->required()) {
                                _tokenize_argument_help_text(
                                    tok, **it, options, hanging );
                            }
                        }
                    } else {
                        for (auto it = _arguments.begin();
                             it != _arguments.end(); ++it) {
                            _tokenize_argument_help_text( tok, **it, options,
                                                          hanging );
                        }
                    }
                    return tok;
                }

            protected:
                std::vector<std::unique_ptr<Argument>> _arguments;
                std::vector<std::string> _preface;
                // size_t _hanging_indent = 0;

                size_t _calculate_hanging_indent(
                    const HelpTextFormatterOptions& options ) const {
                    size_t _hanging_indent = 0;
                    for (auto& arg : _arguments) {
                        std::string flagname  = arg->tag();
                        std::string indicator = _indicator( *arg, options );
                        _hanging_indent       = std::max(
                            _hanging_indent,
                            std::min(
                                flagname.size()
                                    + ( indicator.empty()
                                            ? 1
                                            : indicator.size() + 2 ),
                                options.max_tag_length_without_newline ) );
                    }
                    _hanging_indent
                        += ( options.print_argument_hyphens ? 2 : 0 );
                    return _hanging_indent;
                }

                void _tokenize_argument_help_text(
                    std::vector<std::string>& tok, Argument& arg,
                    const HelpTextFormatterOptions& options,
                    size_t hanging_indent ) const {
                    // std::string flagname = arg.tag();

                    std::string argtag = _argument_tag( arg, options );
                    tok.push_back( argtag );
                    size_t total_size = argtag.size() + 1;
                    // tok.push_back( this->options().marker(
                    //     help_text_marker_t::SET_HANGING_INDENT,
                    //     (int)_hanging_indent ) );

                    if (options.break_after_tag
                        || ( total_size
                             > options.max_tag_length_without_newline )) {
                        tok.push_back( HelpTextMarkers::marker(
                            help_text_marker_t::BREAK ) );
                    } else {
                        int space = (int)( hanging_indent - total_size );
                        if (space > 0) {
                            // add an extra space to account for the one
                            // automatically added after the indicator by the
                            // formatter
                            tok.push_back( HelpTextMarkers::marker(
                                help_text_marker_t::SPACE, space + 1 ) );
                        }
                    }

                    options.tokenize( arg.description(), tok );
                    _add_suffixes( tok, arg, options );
                    // if (arg.has_default_value()) {
                    //     tok.push_back( std::string(
                    //         "[" + this->options().default_value_suffix
                    //         + arg.default_as_string() + "]" ) );
                    // }
                    // if (arg.required()) {
                    //     tok.push_back(
                    //         this->options().required_argument_suffix );
                    // }
                    tok.push_back( HelpTextMarkers::marker(
                        help_text_marker_t::NEWLINE ) );
                }

                virtual void _add_suffixes(
                    std::vector<std::string>& tok, Argument& arg,
                    const HelpTextFormatterOptions& options ) const {
                    std::vector<std::string> suffixes;
                    if (arg.required()) {
                        suffixes.push_back( options.required_argument_suffix );
                    }
                    if (arg.has_default_value()) {
                        suffixes.push_back( options.default_value_suffix
                                            + arg.default_as_string() );
                    }
                    if (!suffixes.empty()) {
                        tok.push_back( options.build_suffix(
                            NCPA::strings::join( suffixes, "," ) ) );
                    }
                }

                virtual std::string _argument_tag(
                    const Argument& arg,
                    const HelpTextFormatterOptions& options ) const {
                    std::ostringstream oss;
                    if (options.print_argument_hyphens) {
                        oss << "--";
                    }
                    oss << arg.tag();
                    std::string indicator = _indicator( arg, options );
                    if (!indicator.empty()) {
                        oss << " " << indicator;
                    }
                    return oss.str();
                }

                std::string _indicator(
                    const Argument& arg,
                    const HelpTextFormatterOptions& options ) const {
                    switch (arg.value_count()) {
                        case value_count_t::ZERO:
                            return options.flag_argument_indicator;
                        case value_count_t::ONE:
                            return options.required_argument_indicator;
                        case value_count_t::ANY:
                            return options
                                .optional_multiple_argument_indicator;
                        case value_count_t::AT_LEAST_ONE:
                            return options
                                .required_multiple_argument_indicator;
                    }
                    return "";
                }

                // void _set_options() {
                //     this->options().first_line_indent          = 0;
                //     this->options().indent_subsections         = true;
                //     this->options().reset_indent_after_newline = true;
                // }

                // std::vector<std::string> _tokenize_preface(
                //     std::vector<std::string>& tok ) const {
                //     // std::vector<std::string> tok;
                //     if (this->options().newline_before_title) {
                //         tok.push_back( this->options().marker(
                //             help_text_marker_t::NEWLINE ) );
                //     }
                //     tok.insert( tok.begin(), _preface.begin(),
                //                 _preface.end() );
                //     tok.push_back( this->options().marker(
                //         help_text_marker_t::NEWLINE ) );
                //     if (this->options().newline_after_title) {
                //         tok.push_back( this->options().marker(
                //             help_text_marker_t::NEWLINE ) );
                //     }
                //     return tok;
                // }
        };

    }  // namespace config
}  // namespace NCPA
