#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/UnaryValidationTest.hpp"

#include <memory>
#include <type_traits>

namespace NCPA {
    namespace config {
        template<typename T>
        class IsNotEqualToTest : public UnaryValidationTest<T> {
            public:
                IsNotEqualToTest( T value ) :
                    UnaryValidationTest<T>( value ) {}

                virtual ~IsNotEqualToTest() {}

                friend void ::swap<>( IsNotEqualToTest<T>& a,
                                      IsNotEqualToTest<T>& b ) noexcept;

                virtual std::unique_ptr<ValidationTest> clone()
                    const override {
                    return std::unique_ptr<ValidationTest>(
                        new IsNotEqualToTest<T>( *this ) );
                }

                virtual std::string description() const override {
                    std::ostringstream oss;
                    oss << "is not equal to " << this->value();
                    return oss.str();
                }

                virtual ValidationTest& test(
                    const BaseParameter *param ) override {
                    if (this->value() != this->parameter_value( param )) {
                        this->pass();
                    } else {
                        this->fail();
                    }
                    return *this;
                }
        };

    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap( NCPA::config::IsNotEqualToTest<T>& a,
           NCPA::config::IsNotEqualToTest<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::UnaryValidationTest<T>&>( a ),
            static_cast<NCPA::config::UnaryValidationTest<T>&>( b ) );
}
