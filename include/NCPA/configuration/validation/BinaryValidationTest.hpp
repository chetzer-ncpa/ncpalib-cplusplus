#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/TypedValidationTest.hpp"

namespace NCPA {
    namespace config {
        template<typename T>
        class BinaryValidationTest : public TypedValidationTest<T> {
            public:
                BinaryValidationTest( T value1, T value2 ) :
                    TypedValidationTest<T>() {
                    _values.resize( 2 );
                    _values[ 0 ] = value1;
                    _values[ 1 ] = value2;
                }

                BinaryValidationTest( const BinaryValidationTest<T>& other ) {
                    _values = other._values;
                }

                virtual ~BinaryValidationTest() {}

                friend void ::swap<>( BinaryValidationTest<T>& a,
                                      BinaryValidationTest<T>& b ) noexcept;

                virtual T value( size_t n ) const override {
                    if (n > 1) {
                        throw std::out_of_range(
                            "Requested index out of range" );
                    }
                    return _values.at( n );
                }

                virtual const std::vector<T>& values() const {
                    return _values;
                }

            protected:
                std::vector<T> _values;
                // T _value1, _value2;
        };

    }  // namespace config
}  // namespace NCPA
