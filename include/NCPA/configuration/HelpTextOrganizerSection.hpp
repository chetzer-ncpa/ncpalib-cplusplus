#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/HelpTextSection.hpp"

#include <memory>
#include <string>

namespace NCPA {
    namespace config {
        class HelpTextOrganizerSection : public HelpTextSection {
            public:
                HelpTextOrganizerSection() : HelpTextSection() {}

                HelpTextOrganizerSection( const std::string& title ) :
                    HelpTextSection( title ) {}

                HelpTextOrganizerSection(
                    const HelpTextOrganizerSection& other ) :
                    HelpTextSection( other ) {}

                HelpTextOrganizerSection(
                    HelpTextOrganizerSection&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~HelpTextOrganizerSection() {}

                HelpTextOrganizerSection& operator=(
                    HelpTextOrganizerSection other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( HelpTextOrganizerSection& a,
                                  HelpTextOrganizerSection& b ) noexcept;

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

        inline void swap( HelpTextOrganizerSection& a,
                          HelpTextOrganizerSection& b ) noexcept {
            using std::swap;
            swap( static_cast<HelpTextSection&>( a ),
                  static_cast<HelpTextSection&>( b ) );
        }
    }  // namespace config
}  // namespace NCPA
