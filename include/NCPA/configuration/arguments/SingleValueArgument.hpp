#pragma once

#include "NCPA/configuration/arguments/TypedArgument.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace NCPA {
    namespace config {
        template<typename T>
        class SingleValueArgument
            : public Cloneable<SingleValueArgument<T>, TypedArgument<T>> {
            public:
                using parent_t
                    = Cloneable<SingleValueArgument<T>, TypedArgument<T>>;
                using parent_t::set_value;

                SingleValueArgument() = default;

                template<typename... Args>
                explicit SingleValueArgument( const std::string& tag,
                                              const std::string& description,
                                              Args&&...args ) :
                    parent_t( tag, description, value_count_t::ZERO_OR_ONE,
                              std::forward<Args>( args )... ) {}

                SingleValueArgument( const SingleValueArgument<T>& other )
                    = default;
                SingleValueArgument( SingleValueArgument<T>&& other ) noexcept
                    = default;
                virtual ~SingleValueArgument() = default;

                SingleValueArgument& operator=(
                    SingleValueArgument<T> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( SingleValueArgument<T>& a,
                                  SingleValueArgument<T>& b ) noexcept {
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
