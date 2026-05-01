#pragma once

#include "NCPA/configuration/declarations.hpp"

#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        class HelpTextSection {
            public:
                HelpTextSection() : HelpTextSection( "" ) {}

                HelpTextSection( const std::string& title ) :
                    _title { title } {}

                HelpTextSection( const HelpTextSection& other ) {
                    _title   = other._title;
                    _depth   = other._depth;
                    _options = other._options;
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
                                  HelpTextSection& b ) noexcept;

                virtual HelpTextSection& add_subsection(
                    const HelpTextSection& sub ) {
                    _subsections.push_back( sub.clone() );
                    if (this->options().indent_subsections) {
                        _subsections.back()->increase_depth();
                    }
                    return *this;
                }

                virtual HelpTextSection& set_title(
                    const std::string& title ) {
                    _title = title;
                    return *this;
                }

                virtual std::string title() const { return _title; }

                virtual std::vector<std::unique_ptr<HelpTextSection>>&
                    subsections() {
                    return _subsections;
                }

                virtual const std::vector<std::unique_ptr<HelpTextSection>>&
                    subsections() const {
                    return _subsections;
                }

                virtual size_t& depth() { return _depth; }

                virtual const size_t& depth() const { return _depth; }

                virtual help_text_section_formatter_options_t& options() {
                    return _options;
                }

                virtual const help_text_section_formatter_options_t& options()
                    const {
                    return _options;
                }

                // virtual size_t& first_line_indent() {
                //     return _first_line_indent;
                // }

                // virtual const size_t& first_line_indent() const {
                //     return _first_line_indent;
                // }

                // virtual size_t& hanging_indent() { return _hanging_indent; }

                // virtual const size_t& hanging_indent() const {
                //     return _hanging_indent;
                // }

                // virtual size_t& title_indent() { return _title_indent; }

                // virtual const size_t& title_indent() const {
                //     return _title_indent;
                // }

                // virtual bool& indent_subsections() {
                //     return _indent_subsections;
                // }

                // virtual const bool& indent_subsections() const {
                //     return _indent_subsections;
                // }

                virtual void increase_depth() {
                    _depth++;
                    for (auto it = _subsections.begin();
                         it != _subsections.end(); ++it) {
                        ( *it )->increase_depth();
                    }
                }

                virtual HelpTextSection& set_text( const std::string& text )
                    = 0;
                virtual std::string text() const                       = 0;
                virtual std::unique_ptr<HelpTextSection> clone() const = 0;

            protected:
                std::string _title;
                std::vector<std::unique_ptr<HelpTextSection>> _subsections;

                size_t _depth = 0;
                help_text_section_formatter_options_t _options;
                // size_t _first_line_indent        = 0;
                // size_t _hanging_indent           = 0;
                // size_t _title_indent             = 0;
                // bool _indent_subsections         = true;
                // bool _reset_indent_after_newline = false;
                // HelpTextSectionFormattingOptions _options;
        };

        inline void swap( HelpTextSection& a, HelpTextSection& b ) noexcept {
            using std::swap;
            swap( a._title, b._title );
            swap( a._subsections, b._subsections );
            swap( a._depth, b._depth );
            swap( a._options, b._options );
        }
    }  // namespace config
}  // namespace NCPA
