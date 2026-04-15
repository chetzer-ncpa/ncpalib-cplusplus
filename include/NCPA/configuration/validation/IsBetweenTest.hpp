#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/TypedValidationTest.hpp"

#include <memory>
#include <type_traits>

namespace NCPA {
    namespace config {
        // specialization for integers
            template<typename T>
            class IsBetweenTest<
                T, typename std::enable_if<std::is_integral<T>::value>::type>
                : public BinaryValidationTest<T> {
                public:
                    IsBetweenTest( T value1, T value2, bool inclusive ) :
                        BinaryValidationTest<T>( value1, value2 ),
                        _inclusive { inclusive } {}

                    IsBetweenTest( const IsBetweenTest<T>& other ) {
                        _inclusive = other._inclusive;
                    }

                    virtual ~IsBetweenTest() {}

                    friend void ::swap<>( IsBetweenTest<T>& a,
                                          IsBetweenTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsBetweenTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is between " << this->value( 0 ) << " and "
                            << this->value( 1 );
                        if (_inclusive) {
                            oss << " (inclusive)";
                        }
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const BaseParameter *param ) override {
                        int val = param->as_int();
                        if (_inclusive) {
                            if (val >= this->value( 0 )
                                && val <= this->value( 1 )) {
                                this->pass();
                            } else {
                                this->fail();
                            }
                        } else {
                            if (val > this->value( 0 )
                                && val < this->value( 1 )) {
                                this->pass();
                            } else {
                                this->fail();
                            }
                        }
                        return *this;
                    }

                private:
                    bool _inclusive;
            };

            // specialization for floating-point
            template<typename T>
            class IsBetweenTest<T, typename std::enable_if<
                                       std::is_floating_point<T>::value>::type>
                : public BinaryValidationTest<T> {
                public:
                    IsBetweenTest( T value1, T value2, bool inclusive ) :
                        BinaryValidationTest<T>( value1, value2 ),
                        _inclusive { inclusive } {}

                    IsBetweenTest( const IsBetweenTest<T>& other ) {
                        _inclusive = other._inclusive;
                    }

                    virtual ~IsBetweenTest() {}

                    friend void ::swap<>( IsBetweenTest<T>& a,
                                          IsBetweenTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsBetweenTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is between " << this->value( 0 ) << " and "
                            << this->value( 1 );
                        if (_inclusive) {
                            oss << " (inclusive)";
                        }
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const BaseParameter *param ) override {
                        double val = param->as_double();
                        if (_inclusive) {
                            if (val >= this->value( 0 )
                                && val <= this->value( 1 )) {
                                this->pass();
                            } else {
                                this->fail();
                            }
                        } else {
                            if (val > this->value( 0 )
                                && val < this->value( 1 )) {
                                this->pass();
                            } else {
                                this->fail();
                            }
                        }
                        return *this;
                    }

                private:
                    bool _inclusive;
            };
    }
}

template<typename T>
void swap( NCPA::config::IsBetweenTest<T>& a,
           NCPA::config::IsBetweenTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::BinaryValidationTest<T>&>( a ),
        static_cast<NCPA::config::BinaryValidationTest<T>&>( b ) );
    swap( a._inclusive, b._inclusive );
}