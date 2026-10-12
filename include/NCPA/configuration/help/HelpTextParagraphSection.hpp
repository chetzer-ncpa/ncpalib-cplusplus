#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/help/HelpTextSection.hpp"

#include <memory>
#include <string>

namespace NCPA {
    namespace config {
        class HelpTextParagraphSection
            : public Cloneable<HelpTextParagraphSection, HelpTextSection> {
                using cloneable_t
                    = Cloneable<HelpTextParagraphSection, HelpTextSection>;

            public:
                HelpTextParagraphSection() = default;

                HelpTextParagraphSection( const std::string& title ) :
                    cloneable_t( title ) {
                    // _set_options();
                }

                HelpTextParagraphSection( const std::string& title,
                                          const std::string& text ) :
                    cloneable_t( title ) {
                    // _set_options();
                    this->add_text( text );
                }

                HelpTextParagraphSection(
                    const HelpTextParagraphSection& other ) = default;

                HelpTextParagraphSection(
                    HelpTextParagraphSection&& other ) noexcept = default;
                virtual ~HelpTextParagraphSection()             = default;

                HelpTextParagraphSection& operator=(
                    HelpTextParagraphSection other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( HelpTextParagraphSection& a,
                                  HelpTextParagraphSection& b ) noexcept {
                    using std::swap;
                    swap( static_cast<cloneable_t&>( a ),
                          static_cast<cloneable_t&>( b ) );
                    // swap( a._text, b._text );
                }

                // virtual HelpTextParagraphSection& add_tokens(
                //     const std::vector<std::string>& text ) override {
                //     _tokens.insert( _tokens.end(), text.begin(), text.end() );
                //     return *this;
                // }

                virtual HelpTextParagraphSection& add_text( const std::string& text ) override {
                    _paragraphs.push_back( text );
                    return *this;
                }

                // virtual void write( HelpTextFormatter& formatter,
                //                     std::ostream& os = std::cout ) const
                //                     override {
                //                         formatter.
                //                     }

                virtual std::vector<std::string> tokens(const HelpTextFormatterOptions& options) const override {
                    std::vector<std::string> tok;
                    for (auto& par : _paragraphs) {
                        options.tokenize( par, tok );
                    }
                    return tok;
                }

                // virtual std::string text() const override { return _text; }

            protected:
                std::vector<std::string> _paragraphs;
                std::vector<std::string> _tokens;

                // void _set_options() {
                //     this->options().first_line_indent = 0;
                //     this->options().hanging_indent    = 0;
                //     this->options().title_indent      = 0;
                // }
        };
    }  // namespace config
}  // namespace NCPA
