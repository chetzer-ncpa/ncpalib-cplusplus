#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/TypedValidationTest.hpp"

#include <memory>
#include <type_traits>

namespace NCPA {
    namespace config {
        // specialization for integers
        template<typename T>
        class IsGreaterThanOrEqualToTest<
            T, typename std::enable_if<std::is_integral<T>::value>::type>
            : public UnaryValidationTest<T> {
            public:
                IsGreaterThanOrEqualToTest( T value ) :
                    UnaryValidationTest<T>( value ) {}

                virtual ~IsGreaterThanOrEqualToTest() {}

                friend void ::swap<>(
                    IsGreaterThanOrEqualToTest<T>& a,
                    IsGreaterThanOrEqualToTest<T>& b ) noexcept;

                virtual std::unique_ptr<ValidationTest> clone()
                    const override {
                    return std::unique_ptr<ValidationTest>(
                        new IsGreaterThanOrEqualToTest<T>( *this ) );
                }

                virtual std::string description() const override {
                    std::ostringstream oss;
                    oss << "is greater than or equal to " << this->value();
                    return oss.str();
                }

                virtual ValidationTest& test(
                    const BaseParameter *param ) override {
                    if (param->as_int() >= this->value()) {
                        this->pass();
                    } else {
                        this->fail();
                    }
                    return *this;
                }
        };

        // specialization for floating-point
        template<typename T>
        class IsGreaterThanOrEqualToTest<
            T, typename std::enable_if<std::is_floating_point<T>::value>::type>
            : public UnaryValidationTest<T> {
            public:
                IsGreaterThanOrEqualToTest( T value ) :
                    UnaryValidationTest<T>( value ) {}

                virtual ~IsGreaterThanOrEqualToTest() {}

                friend void ::swap<>(
                    IsGreaterThanOrEqualToTest<T>& a,
                    IsGreaterThanOrEqualToTest<T>& b ) noexcept;

                virtual std::unique_ptr<ValidationTest> clone()
                    const override {
                    return std::unique_ptr<ValidationTest>(
                        new IsGreaterThanOrEqualToTest<T>( *this ) );
                }

                virtual std::string description() const override {
                    std::ostringstream oss;
                    oss << "is greater than or equal to " << this->value();
                    return oss.str();
                }

                virtual ValidationTest& test(
                    const BaseParameter *param ) override {
                    if (param->as_double() >= this->value()) {
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
void swap( NCPA::config::IsGreaterThanOrEqualToTest<T>& a,
           NCPA::config::IsGreaterThanOrEqualToTest<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::UnaryValidationTest<T>&>( a ),
            static_cast<NCPA::config::UnaryValidationTest<T>&>( b ) );
}
