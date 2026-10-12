#pragma once

#include "NCPA/configuration/arguments/TypedArgument.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <stdexcept>
#include <string>

namespace NCPA {
    namespace config {
        template<typename T>
        class RequiredMultiValueArgument
            : public Cloneable<RequiredMultiValueArgument<T>, TypedArgument<T>> {
                using parent_t
                    = Cloneable<RequiredMultiValueArgument<T>, TypedArgument<T>>;

            public:
            RequiredMultiValueArgument() = default;

                template<typename... Args>
                explicit RequiredMultiValueArgument( const std::string& tag,
                                             const std::string& description,
                                             Args&&...args ) :
                    parent_t( tag, description, value_count_t::AT_LEAST_ONE,
                              std::forward<Args>( args )... ) {}

                RequiredMultiValueArgument( const RequiredMultiValueArgument<T>& other )
                    = default;
                RequiredMultiValueArgument( RequiredMultiValueArgument<T>&& other ) noexcept
                    = default;
                virtual ~RequiredMultiValueArgument() = default;

                RequiredMultiValueArgument& operator=( RequiredMultiValueArgument<T> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( RequiredMultiValueArgument<T>& a,
                                  RequiredMultiValueArgument<T>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<parent_t&>( a ),
                          static_cast<parent_t&>( b ) );
                }
        };
    }  // namespace config
}  // namespace NCPA
