#pragma once

#include "NCPA/configuration/arguments/TypedArgument.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <stdexcept>
#include <string>

namespace NCPA {
    namespace config {
        template<typename T>
        class MultiValueArgument
            : public Cloneable<MultiValueArgument<T>, TypedArgument<T>> {
                using parent_t = Cloneable<MultiValueArgument<T>, TypedArgument<T>> ;
            public:
                MultiValueArgument() = default;
                
                template<typename... Args>
                explicit MultiValueArgument( const std::string& tag,
                                       const std::string& description, Args&&...args ) :
                    parent_t( tag, description, value_count_t::ANY, std::forward<Args>( args )... ) {}

                MultiValueArgument( const MultiValueArgument<T>& other ) = default;
                MultiValueArgument( MultiValueArgument<T>&& other ) noexcept = default;
                virtual ~MultiValueArgument() = default;

                MultiValueArgument& operator=( MultiValueArgument<T> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( MultiValueArgument<T>& a,
                                  MultiValueArgument<T>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<parent_t&>( a ),
                          static_cast<parent_t&>( b ) );
                }
        };
    }  // namespace config
}  // namespace NCPA
