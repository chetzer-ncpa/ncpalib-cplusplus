#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/declarations.hpp"
#include "NCPA/configuration/validation/UnaryValidationTest.hpp"

#include <memory>
#include <type_traits>

namespace NCPA {
    namespace config {
        template<typename T>
        class IsEqualToTest : public UnaryValidationTest<T> {
            public:
                IsEqualToTest( T value ) : UnaryValidationTest<T>( value ) {}

                virtual ~IsEqualToTest() {}

                friend void ::swap<>( IsEqualToTest<T>& a,
                                      IsEqualToTest<T>& b ) noexcept;

                virtual std::unique_ptr<ValidationTest> clone()
                    const override {
                    return std::unique_ptr<ValidationTest>(
                        new IsEqualToTest<T>( *this ) );
                }

                virtual std::string description() const override {
                    std::ostringstream oss;
                    oss << "is equal to " << this->value();
                    return oss.str();
                }

                virtual ValidationTest& test(
                    const Parameter *param ) override {
                    if (this->value() == this->parameter_value( param )) {
                        this->pass();
                    } else {
                        this->fail();
                    }
                    return *this;
                }
        };

    }  // namespace config
}  // namespace NCPA
