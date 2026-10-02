#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/HelpTextSection.hpp"

#include <memory>
#include <string>

namespace NCPA {
    namespace config {
        class HelpTextOrganizerSection : public HelpTextSection {
            public:
                HelpTextOrganizerSection() = default;

                HelpTextOrganizerSection( const std::string& title ) :
                    HelpTextSection( title ) {}

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
                    swap( static_cast<HelpTextSection&>( a ),
                          static_cast<HelpTextSection&>( b ) );
                }

                virtual std::unique_ptr<HelpTextSection> clone() const {
                    return std::unique_ptr<HelpTextSection>(
                        new HelpTextOrganizerSection( *this ) );
                }

                virtual std::string text() const override { return ""; }

                virtual HelpTextOrganizerSection& set_text(
                    const std::string& text ) override {
                    return *this;
                }
        };
    }  // namespace config
}  // namespace NCPA
