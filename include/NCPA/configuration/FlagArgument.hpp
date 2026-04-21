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

                FlagArgument( const std::string& tag ) :
                    TypedArgument<bool>( tag ) {}

                FlagArgument( const std::string& tag, bool defaultval ) :
                    TypedArgument<bool>( tag, defaultval ) {}

                FlagArgument(
                    const std::string& tag,
                    const mapping_ptr_t<bool, std::string>& mapping ) :
                    TypedArgument<bool>( tag, *mapping ) {}

                FlagArgument(
                    const std::string& tag, bool defaultval,
                    const mapping_ptr_t<bool, std::string>& mapping ) :
                    TypedArgument<bool>( tag, defaultval, *mapping ) {}

                FlagArgument( const std::string& tag,
                              const Mapping<bool, std::string>& mapping ) :
                    TypedArgument<bool>( tag, mapping ) {}

                FlagArgument( const std::string& tag, bool defaultval,
                              const Mapping<bool, std::string>& mapping ) :
                    TypedArgument<bool>( tag, defaultval, mapping ) {}

                virtual ~FlagArgument() {}

                virtual bool expects_value() const override { return false; }

                virtual TypedArgument<bool>& set() override {
                    return this->set( true );
                }
        };
    }  // namespace config
}  // namespace NCPA
