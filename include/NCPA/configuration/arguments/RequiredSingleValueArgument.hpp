#pragma once

#include "NCPA/configuration/arguments/TypedArgument.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace NCPA {
    namespace config {
        template<typename T>
        class RequiredSingleValueArgument
            : public Cloneable<RequiredSingleValueArgument<T>, TypedArgument<T>> {
            public:
                using parent_t
                    = Cloneable<RequiredSingleValueArgument<T>, TypedArgument<T>>;
                using parent_t::set_value;

                RequiredSingleValueArgument() = default;

                template<typename... Args>
                explicit RequiredSingleValueArgument( const std::string& tag,
                                              const std::string& description,
                                              Args&&...args ) :
                    parent_t( tag, description, value_count_t::ONE,
                              std::forward<Args>( args )... ) {}

                RequiredSingleValueArgument( const RequiredSingleValueArgument<T>& other )
                    = default;
                RequiredSingleValueArgument( RequiredSingleValueArgument<T>&& other ) noexcept
                    = default;
                virtual ~RequiredSingleValueArgument() = default;

                RequiredSingleValueArgument& operator=(
                    RequiredSingleValueArgument<T> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( RequiredSingleValueArgument<T>& a,
                                  RequiredSingleValueArgument<T>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<parent_t&>( a ),
                          static_cast<parent_t&>( b ) );
                }

                virtual void set_value( const std::vector<T>& val ) override {
                    if (val.size() > 1) {
                        throw TooManyValues(
                            "Attempted to set multiple values for "
                            "single-value argument "
                            + this->tag() );
                    }
                    this->set_value( val );
                }
        };
    }  // namespace config
}  // namespace NCPA
