#pragma once

#include "NCPA/configuration/Argument.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/HelpTextSection.hpp"

#include <iostream>

namespace NCPA {
    namespace config {
        class HelpTextArgumentSection : public HelpTextSection {
            public:
                HelpTextArgumentSection() : HelpTextSection() {
                    _set_options();
                }

                HelpTextArgumentSection( const std::string& title ) :
                    HelpTextSection( title ) {
                    _set_options();
                }

                HelpTextArgumentSection(
                    const HelpTextArgumentSection& other ) :
                    HelpTextSection( other ) {
                    _newline_after_tag = other._newline_after_tag;
                    _max_tag_length_without_newline
                        = other._max_tag_length_without_newline;
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
                                  HelpTextArgumentSection& b ) noexcept;

                virtual std::unique_ptr<HelpTextSection> clone() const {
                    return std::unique_ptr<HelpTextSection>(
                        new HelpTextArgumentSection( *this ) );
                }

                virtual HelpTextArgumentSection& set_text(
                    const std::string& text ) override {
                    return *this;
                }

                virtual std::string text() const override {
                    size_t arg_length = 0;
                    std::ostringstream oss;
                    if (_separate_required_from_optional) {
                        for (auto it = _arguments.begin();
                             it != _arguments.end(); ++it) {
                            if (( *it )->required()) {
                                _add_argument_help_text( oss, **it );
                            }
                        }
                        for (auto it = _arguments.begin();
                             it != _arguments.end(); ++it) {
                            if (!( *it )->required()) {
                                _add_argument_help_text( oss, **it );
                            }
                        }
                    } else {
                        for (auto it = _arguments.begin();
                             it != _arguments.end(); ++it) {
                            _add_argument_help_text( oss, **it );
                            // std::string flagname = ( *it )->tag();
                            // arg_length = std::max( arg_length,
                            // flagname.size() ); oss << "--" << flagname; if
                            // (_newline_after_tag
                            //     || ( flagname.size() + 3 )
                            //            > _max_tag_length_without_newline) {
                            //     oss << " " << NEWLINE_MARKER << " ";
                            // }
                            // oss << ( *it )->help_text() << NEWLINE_MARKER;
                        }
                    }
                    return oss.str();
                }

                virtual void add_argument( const Argument& arg ) {
                    _arguments.push_back( arg.clone() );
                    std::string flagname           = arg.tag();
                    size_t hanging = std::max(
                        this->options().hanging_indent,
                        std::min( flagname.size() + (arg.expects_value() ? 9 : 3),
                                  _max_tag_length_without_newline ) );
                    std::cout << "Setting hanging indent to " << hanging << std::endl;
                    this->options().hanging_indent = hanging;
                }

                virtual void add_argument( Argument *arg ) {
                    add_argument( *arg );
                }

                virtual void add_argument( std::unique_ptr<Argument>& arg ) {
                    add_argument( *arg );
                }

            protected:
                void _add_argument_help_text( std::ostream& os,
                                              Argument& arg ) const {
                    // std::string flagname = arg.tag();
                    os << "--" << arg.tag() << " ";
                    size_t total_size = arg.tag().size() + 3;
                    if (arg.expects_value()) {
                        os << "[val] ";
                        total_size += 6;
                    }
                    if (_newline_after_tag
                        || ( total_size ) > _max_tag_length_without_newline) {
                        os << NEWLINE_MARKER << " ";
                    }
                    os << "<i" << this->options().hanging_indent << "> " << arg.help_text();
                    if (arg.has_default()) { 
                        os << " [default=" << arg.default_string() << "]";
                    }
                    os << NEWLINE_MARKER;
                }

                std::vector<std::unique_ptr<Argument>> _arguments;

                bool _newline_after_tag                = false;
                bool _separate_required_from_optional  = true;
                size_t _max_tag_length_without_newline = 30;

                void _set_options() {
                    this->options().first_line_indent          = 0;
                    this->options().indent_subsections         = true;
                    this->options().reset_indent_after_newline = true;
                }
        };

        inline void swap( HelpTextArgumentSection& a,
                          HelpTextArgumentSection& b ) noexcept {
            using std::swap;
            swap( static_cast<HelpTextSection&>( a ),
                  static_cast<HelpTextSection&>( b ) );
            swap( a._arguments, b._arguments );
            swap( a._newline_after_tag, b._newline_after_tag );
            swap( a._max_tag_length_without_newline,
                  b._max_tag_length_without_newline );
        }
    }  // namespace config
}  // namespace NCPA
