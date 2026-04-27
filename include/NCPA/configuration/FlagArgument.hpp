#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/TypedArgument.hpp"

#include <stdexcept>
#include <string>

namespace NCPA {
    namespace config {
        class FlagArgument : public TypedArgument<bool> {
                using TypedArgument<bool>::set;

            public:
                FlagArgument() : TypedArgument<bool>() {}

                FlagArgument( const std::string& tag, const std::string& help_text ) :
                    TypedArgument<bool>( tag, help_text, false ) {}

                FlagArgument(
                    const std::string& tag, const std::string& help_text,
                    const mapping_ptr_t<bool, std::string>& mapping ) :
                    TypedArgument<bool>( tag, help_text, false, *mapping ) {}

                // FlagArgument(
                //     const std::string& tag, bool defaultval,
                //     const mapping_ptr_t<bool, std::string>& mapping ) :
                //     TypedArgument<bool>( tag, defaultval, *mapping ) {}

                FlagArgument( const std::string& tag,const std::string& help_text,
                              const Mapping<bool, std::string>& mapping ) :
                    TypedArgument<bool>( tag, help_text, false, mapping ) {}

                // FlagArgument( const std::string& tag, bool defaultval,
                //               const Mapping<bool, std::string>& mapping ) :
                //     TypedArgument<bool>( tag, defaultval, mapping ) {}

                FlagArgument( const FlagArgument& other ) :
                    TypedArgument<bool>( other ) {}

                FlagArgument( FlagArgument&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~FlagArgument() {}

                FlagArgument& operator=( FlagArgument other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( FlagArgument& a, FlagArgument& b ) noexcept;

                virtual std::unique_ptr<Argument> clone() const override {
                    return std::unique_ptr<Argument>(
                        new FlagArgument( *this ) );
                }

                virtual bool expects_value() const override { return false; }

                virtual TypedArgument<bool>& set() override {
                    return this->set( true );
                }

                virtual bool has_default() const override { return false; }

                virtual std::string default_string() const override {
                    return "false";
                }
        };

        void swap( FlagArgument& a, FlagArgument& b ) noexcept {
            using std::swap;
            swap( static_cast<TypedArgument<bool>&>( a ),
                  static_cast<TypedArgument<bool>&>( b ) );
        }
    }  // namespace config
}  // namespace NCPA
