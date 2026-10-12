#pragma once

#include "NCPA/configuration/arguments/Argument.hpp"
#include "NCPA/configuration/arguments/TypedArgument.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/exceptions.hpp"
#include "NCPA/configuration/parameters/VectorParameter.hpp"

#include <stdexcept>
#include <string>

namespace NCPA {
    namespace config {
        class FlagArgument
            : public Cloneable<FlagArgument, TypedArgument<bool>> {
            public:
                using parent_t = Cloneable<FlagArgument, TypedArgument<bool>>;
                using parent_t::set_value;

                FlagArgument() = default;

                template<typename... Args>
                explicit FlagArgument( const std::string& tag,
                                       const std::string& description,
                                       Args&&...args ) :
                    parent_t( tag, description, value_count_t::ZERO,
                              std::forward<Args>( args )... ) {}

                FlagArgument( const FlagArgument& other )     = default;
                FlagArgument( FlagArgument&& other ) noexcept = default;
                virtual ~FlagArgument()                       = default;

                FlagArgument& operator=( FlagArgument other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( FlagArgument& a, FlagArgument& b ) noexcept {
                    using std::swap;
                    swap( static_cast<parent_t&>( a ),
                          static_cast<parent_t&>( b ) );
                }

                virtual void set_value(
                    const std::vector<bool>& val ) override {
                    if (val.size() > 1) {
                        throw TooManyValues( "Attempted to set multiple "
                                             "values for flag argument "
                                             + this->tag() );
                    }
                    this->set_value( true );
                }
        };
    }  // namespace config
}  // namespace NCPA
