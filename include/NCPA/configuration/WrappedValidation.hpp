#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/TypedValidation.hpp"

#include <functional>

namespace NCPA {
    namespace config {
        template<typename T>
        class WrappedValidation : public TypedValidation<std::string> {
            public:
                WrappedValidation() : TypedValidation<std::string>() {}

                WrappedValidation( std::function<T( std::string )> conversion,
                                   validation_function_t<T> base ) :
                    TypedValidation<std::string>(),
                    convert { conversion },
                    _sub { TypedValidation<T>( base ) } {}

                virtual ~WrappedValidation() {}

                WrappedValidation( const WrappedValidation<T>& other ) :
                    TypedValidation<std::string>( other ) {
                    convert = other.convert;
                    _sub    = other._sub;
                    // validate = other.validate;
                }

                WrappedValidation( WrappedValidation<T>&& other ) noexcept :
                    TypedValidation<std::string>() {
                    swap( *this, other );
                }

                WrappedValidation& operator=(
                    WrappedValidation<T> other ) noexcept {
                    swap( *this, other );
                    return *this;
                }

                WrappedValidation( validation_function_t<T> v ) {
                    _sub._validate = v;
                }

                friend void swap( WrappedValidation<T>& a,
                                  WrappedValidation<T>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<TypedValidation<std::string>&>( a ),
                          static_cast<TypedValidation<std::string>&>( b ) );
                    swap( a._sub, b._sub );
                    swap( a.convert, b.convert );
                }

                virtual validation_status_t validate(
                    const std::string& val ) const override {
                    return _sub.validate( convert( val ) );
                }

                NCPA_CLONE_METHOD( WrappedValidation<std::string>, Validation )

                std::function<T( std::string )> convert;

            protected:
                TypedValidation<T> _sub;
        };
    }  // namespace config
}  // namespace NCPA
