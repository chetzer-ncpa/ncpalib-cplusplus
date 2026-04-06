#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/TypedParameter.hpp"
#include "NCPA/configuration/ValidationTest.hpp"

namespace NCPA {
    namespace config {
        template<typename T>
        class TypedValidationTest : public ValidationTest {
            public:
                TypedValidationTest() {}

                TypedValidationTest( const TypedValidationTest<T>& other ) {}

                virtual ~TypedValidationTest() {}

                friend void ::swap<>( TypedValidationTest<T>& a,
                                      TypedValidationTest<T>& b ) noexcept;

                virtual const T& parameter_value(
                    const Parameter *param ) const {
                    return dynamic_cast<const TypedParameter<T> *>( param )
                        ->value();
                }

                virtual T value( size_t n ) const = 0;
        };
    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap( NCPA::config::TypedValidationTest<T>& a,
           NCPA::config::TypedValidationTest<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::ValidationTest&>( a ),
            static_cast<NCPA::config::ValidationTest&>( b ) );
}