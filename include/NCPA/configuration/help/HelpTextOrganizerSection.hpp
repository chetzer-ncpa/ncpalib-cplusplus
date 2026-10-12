#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/help/HelpTextSection.hpp"

#include <memory>
#include <string>

namespace NCPA {
    namespace config {
        class HelpTextOrganizerSection
            : public Cloneable<HelpTextOrganizerSection, HelpTextSection> {
                using cloneable_t
                    = Cloneable<HelpTextOrganizerSection, HelpTextSection>;

            public:
                HelpTextOrganizerSection() = default;

                HelpTextOrganizerSection( const std::string& title ) :
                    cloneable_t( title ) {}

                HelpTextOrganizerSection(
                    const HelpTextOrganizerSection& other ) = default;

                HelpTextOrganizerSection(
                    HelpTextOrganizerSection&& other ) noexcept = default;

                virtual ~HelpTextOrganizerSection() = default;

                HelpTextOrganizerSection& operator=(
                    HelpTextOrganizerSection other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( HelpTextOrganizerSection& a,
                                  HelpTextOrganizerSection& b ) noexcept {
                    using std::swap;
                    swap( static_cast<cloneable_t&>( a ),
                          static_cast<cloneable_t&>( b ) );
                }

                virtual HelpTextOrganizerSection& add_text(
                    const std::string& text ) override {
                    return *this;
                }

                virtual HelpTextOrganizerSection& add_tokens(
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
