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
                    _title = other._title;
                    // _text = other._text;
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

                virtual HelpTextSection& set_text( const std::string& text )
                    = 0;
                virtual std::string text() const                       = 0;
                virtual std::unique_ptr<HelpTextSection> clone() const = 0;

            protected:
                std::string _title;
                // std::string _text;
                std::vector<std::unique_ptr<HelpTextSection>> _subsections;
        };

        inline void swap( HelpTextSection& a, HelpTextSection& b ) noexcept {
            using std::swap;
            swap( a._title, b._title );
            // swap(a._text, b._text);
            swap( a._subsections, b._subsections );
        }
    }  // namespace config
}  // namespace NCPA
