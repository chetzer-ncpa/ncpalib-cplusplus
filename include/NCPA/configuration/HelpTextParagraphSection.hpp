#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/HelpTextSection.hpp"

#include <memory>
#include <string>

namespace NCPA {
    namespace config {
        class HelpTextParagraphSection : public HelpTextSection {
            public:
                HelpTextParagraphSection() : HelpTextSection() {
                    _set_options();
                }

                HelpTextParagraphSection( const std::string& title ) :
                    HelpTextSection( title ) {
                    _set_options();
                }

                HelpTextParagraphSection( const std::string& title, const std::string& text ) :
                    HelpTextSection( title ) {
                    _set_options();
                    _text = text;
                }

                HelpTextParagraphSection(
                    const HelpTextParagraphSection& other ) :
                    HelpTextSection( other ) {
                    _text = other._text;
                }

                HelpTextParagraphSection(
                    HelpTextParagraphSection&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~HelpTextParagraphSection() {}

                HelpTextParagraphSection& operator=(
                    HelpTextParagraphSection other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( HelpTextParagraphSection& a,
                                  HelpTextParagraphSection& b ) noexcept;

                virtual std::unique_ptr<HelpTextSection> clone() const {
                    return std::unique_ptr<HelpTextSection>(
                        new HelpTextParagraphSection( *this ) );
                }

                virtual HelpTextParagraphSection& set_text(
                    const std::string& text ) override {
                    _text = text;
                    return *this;
                }

                virtual std::string text() const override { return _text; }

            protected:
                void _set_options() {
                    this->options().first_line_indent = 0;
                    this->options().hanging_indent    = 0;
                    this->options().title_indent      = 0;
                }

            private:
                std::string _text;
        };

        inline void swap( HelpTextParagraphSection& a,
                          HelpTextParagraphSection& b ) noexcept {
            using std::swap;
            swap( static_cast<HelpTextSection&>( a ),
                  static_cast<HelpTextSection&>( b ) );
            swap( a._text, b._text );
        }
    }  // namespace config
}  // namespace NCPA
