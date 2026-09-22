#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/Validation.hpp"

#include <functional>

namespace NCPA {
    namespace config {
        template<typename T>
        class TypedValidation : // public Validation,
        public Cloneable<TypedValidation<T>, Validation> {
            public:

                TypedValidation() {}

                virtual ~TypedValidation() {}

                TypedValidation( const TypedValidation<T>& other ) :
                    Cloneable<TypedValidation<T>, Validation> ( other ) {
                    validate = other.validate;
                }

                TypedValidation( TypedValidation<T>&& other ) noexcept {
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
                    swap( static_cast<Cloneable<TypedValidation<T>, Validation> &>( a ),
                          static_cast<Cloneable<TypedValidation<T>, Validation> &>( b ) );
                    swap( a.validate, b.validate );
                }

                validation_function_t<T> validate;
        };
    }  // namespace config
}  // namespace NCPA
