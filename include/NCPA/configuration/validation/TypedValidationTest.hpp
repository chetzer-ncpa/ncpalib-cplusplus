#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/declarations.hpp"
#include "NCPA/configuration/validation/ValidationTest.hpp"

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
                    return dynamic_cast<const Parameter<T> *>( param )
                        ->value();
                }

                virtual T value( size_t n ) const = 0;
        };
    }  // namespace config
}  // namespace NCPA
