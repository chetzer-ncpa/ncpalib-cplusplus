#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/TypedValidationTest.hpp"

#include <memory>
#include <type_traits>

namespace NCPA {
    namespace config {
        // specialization for integers
        template<typename T>
        class IsGreaterThanTest<
            T, typename std::enable_if<std::is_integral<T>::value>::type>
            : public UnaryValidationTest<T> {
            public:
                IsGreaterThanTest( T value ) :
                    UnaryValidationTest<T>( value ) {}

                IsGreaterThanTest( const IsGreaterThanTest<T>& other ) {}

                virtual ~IsGreaterThanTest() {}

                friend void ::swap<>( IsGreaterThanTest<T>& a,
                                      IsGreaterThanTest<T>& b ) noexcept;

                virtual std::unique_ptr<ValidationTest> clone()
                    const override {
                    return std::unique_ptr<ValidationTest>(
                        new IsGreaterThanTest<T>( *this ) );
                }

                virtual std::string description() const override {
                    std::ostringstream oss;
                    oss << "is greater than " << this->value();
                    return oss.str();
                }

                virtual ValidationTest& test(
                    const Parameter *param ) override {
                    if (param->as_int() > this->value()) {
                        this->pass();
                    } else {
                        this->fail();
                    }
                    return *this;
                }
        };

        // specialization for floating-point
        template<typename T>
        class IsGreaterThanTest<
            T, typename std::enable_if<std::is_floating_point<T>::value>::type>
            : public UnaryValidationTest<T> {
            public:
                IsGreaterThanTest( T value ) :
                    UnaryValidationTest<T>( value ) {}

                IsGreaterThanTest( const IsGreaterThanTest<T>& other ) {}

                virtual ~IsGreaterThanTest() {}

                friend void ::swap<>( IsGreaterThanTest<T>& a,
                                      IsGreaterThanTest<T>& b ) noexcept;

                virtual std::unique_ptr<ValidationTest> clone()
                    const override {
                    return std::unique_ptr<ValidationTest>(
                        new IsGreaterThanTest<T>( *this ) );
                }

                virtual std::string description() const override {
                    std::ostringstream oss;
                    oss << "is greater than " << this->value();
                    return oss.str();
                }

                virtual ValidationTest& test(
                    const Parameter *param ) override {
                    if (param->as_double() > this->value()) {
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
void swap( NCPA::config::IsGreaterThanTest<T>& a,
           NCPA::config::IsGreaterThanTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::UnaryValidationTest<T>&>( a ),
        static_cast<NCPA::config::UnaryValidationTest<T>&>( b ) );
}