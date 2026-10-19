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
                virtual HelpTextSection& add_text( const std::string& text )
                    = 0;
                virtual std::vector<std::string> tokens(
                    const HelpTextFormatterOptions& options ) const = 0;

                HelpTextSection() : HelpTextSection( "" ) {}

                HelpTextSection( const std::string& title ) :
                    _title { title }, _subsections {} {}

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
                    const HelpTextSection& sub,
                    bool inherit_format_options = true ) {
                    _subsections.push_back( sub.clone() );
                    _subsections.back()->increase_depth();
                    if (inherit_format_options) {
                        _subsections.back()->options() = this->options();
                    }
                    _subsections.back()->override_options();
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

                HelpTextFormatterOptions& options() { return _options; }

                const HelpTextFormatterOptions& options() const {
                    return _options;
                }

                virtual void override_options() {}

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

            protected:
                std::string _title;
                std::vector<std::unique_ptr<HelpTextSection>> _subsections;
                size_t _depth = 0;

                HelpTextFormatterOptions _options;
        };
    }  // namespace config
}  // namespace NCPA
