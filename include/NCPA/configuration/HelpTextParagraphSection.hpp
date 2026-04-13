#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/HelpTextSection.hpp"

#include <memory>
#include <string>

namespace NCPA {
    namespace config {
        class HelpTextParagraphSection : public HelpTextSection {
            public:
                HelpTextParagraphSection() : HelpTextSection() {}

                HelpTextParagraphSection( const std::string& title ) :
                    HelpTextSection( title ) {}

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

                virtual std::string text() const override { return _text; }

                virtual HelpTextParagraphSection& set_text(
                    const std::string& text ) override {
                    _text = text;
                    return *this;
                }

            protected:
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
