#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/Validation.hpp"

#include <functional>

namespace NCPA {
    namespace config {
        template<typename T>
        class TypedValidation : public Validation {
            public:

                TypedValidation() : Validation() {}

                virtual ~TypedValidation() {}

                TypedValidation( const TypedValidation<T>& other ) :
                    Validation( other ) {
                    validate = other.validate;
                }

                TypedValidation( TypedValidation<T>&& other ) noexcept :
                    TypedValidation<T>() {
                    swap( *this, other );
                }

                TypedValidation& operator=( TypedValidation<T> other ) noexcept {
                    swap( *this, other );
                    return *this;
                }

                TypedValidation( validation_function_t<T> v ) { validate = v; }

                friend void swap( TypedValidation<T>& a,
                                  TypedValidation<T>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<Validation&>( a ),
                          static_cast<Validation&>( b ) );
                    swap( a.validate, b.validate );
                }

                NCPA_CLONE_METHOD( TypedValidation<T>, Validation )

                validation_function_t<T> validate;
        };
    }  // namespace config
}  // namespace NCPA
