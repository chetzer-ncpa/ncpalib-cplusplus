#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/help/HelpTextFormatterOptions.hpp"
#include "NCPA/configuration/help/HelpTextParagraphSection.hpp"

#include <memory>
#include <string>

namespace NCPA {
    namespace config {
        class HelpTextHeaderSection
            : public Cloneable<HelpTextHeaderSection, HelpTextParagraphSection> {
                using cloneable_t
                    = Cloneable<HelpTextHeaderSection, HelpTextParagraphSection>;

            public:
                HelpTextHeaderSection() = default;

                HelpTextHeaderSection( const std::string& title ) :
                    cloneable_t( title ) {}

                HelpTextHeaderSection( const std::string& title,
                                       const std::string& text ) :
                    cloneable_t( title ) {
                    this->add_text( text );
                }

                HelpTextHeaderSection( const HelpTextHeaderSection& other )
                    = default;

                HelpTextHeaderSection( HelpTextHeaderSection&& other ) noexcept
                    = default;
                virtual ~HelpTextHeaderSection() = default;

                HelpTextHeaderSection& operator=(
                    HelpTextHeaderSection other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( HelpTextHeaderSection& a,
                                  HelpTextHeaderSection& b ) noexcept {
                    using std::swap;
                    swap( static_cast<cloneable_t&>( a ),
                          static_cast<cloneable_t&>( b ) );
                    swap( a._text, b._text );
                }

                virtual HelpTextHeaderSection& add_text(
                    const std::string& text ) override {
                    _text += " ";
                    _text += text;
                    return *this;
                }

                virtual void override_options() override {
                    _options.box_style = help_text_box_t::TITLE_ONLY;
                }

                virtual std::vector<std::string> tokens(
                    const HelpTextFormatterOptions& options ) const override {
                    return options.tokenize( _text );
                }

            protected:
                std::string _text;
        };
    }  // namespace config
}  // namespace NCPA
