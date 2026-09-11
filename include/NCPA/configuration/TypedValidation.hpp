#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/Validation.hpp"

#include <functional>

namespace NCPA {
    namespace config {
        template<typename T>
        class TypedValidation : public Validation {
            public:
                TypedValidation( const std::string& failmsg = "" ) :
                    Validation( failmsg ) {}

                virtual ~TypedValidation() {}

                TypedValidation( const TypedValidation<T>& other ) :
                    Validation( other ) {
                    _validate = other._validate;
                }

                TypedValidation( TypedValidation<T>&& other ) noexcept :
                    TypedValidation<T>() {
                    swap( *this, other );
                }

                TypedValidation& operator=(
                    TypedValidation<T> other ) noexcept {
                    swap( *this, other );
                    return *this;
                }

                TypedValidation( validation_function_t<T> v ) {
                    _validate = v;
                }

                friend void swap( TypedValidation<T>& a,
                                  TypedValidation<T>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<Validation&>( a ),
                          static_cast<Validation&>( b ) );
                    swap( a._validate, b._validate );
                }

                virtual validation_status_t validate( const T& val ) const {
                    validation_status_t status;
                    if (_validate( val )) {
                        status.result = test_result_t::PASSED;
                    } else {
                        status.result  = test_result_t::FAILED;
                        status.message = this->failure_message();
                    }
                    return status;
                }

                NCPA_CLONE_METHOD( TypedValidation<T>, Validation )

                validation_function_t<T> _validate;
        };
    }  // namespace config
}  // namespace NCPA
