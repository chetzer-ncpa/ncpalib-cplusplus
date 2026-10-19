#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/help/HelpTextHeaderSection.hpp"
#include "NCPA/configuration/help/HelpTextOrganizerSection.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace NCPA {
    namespace config {
        class HelpText : public Cloneable<HelpText, HelpTextOrganizerSection> {
                using cloneable_t
                    = Cloneable<HelpText, HelpTextOrganizerSection>;

            public:
                HelpText() = default;

                HelpText( const std::string& title ) : cloneable_t() {
                    HelpTextHeaderSection header( title );
                    this->add_subsection( header );
                }

                HelpText( const HelpText& other ) = default;

                HelpText( HelpText&& other ) noexcept = default;

                virtual ~HelpText() = default;

                HelpText& operator=( HelpText other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( HelpText& a, HelpText& b ) noexcept {
                    using std::swap;
                    swap( static_cast<cloneable_t&>( a ),
                          static_cast<cloneable_t&>( b ) );
                }

                virtual HelpText& add_text(
                    const std::string& text ) override {
                    return *this;
                }

                virtual HelpText& add_tokens(
                    const std::vector<std::string>& text ) {
                    return *this;
                }

                virtual std::vector<std::string> tokens(
                    const HelpTextFormatterOptions& options ) const override {
                    return std::vector<std::string>();
                }
        };

    }  // namespace config
}  // namespace NCPA
