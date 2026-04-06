#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/Parameter.hpp"
#include "NCPA/configuration/validation/TypedValidationTest.hpp"

#include <vector>

namespace NCPA {
    namespace config {
        template<typename T>
        class UnaryValidationTest : public TypedValidationTest<T> {
            public:
                UnaryValidationTest( T value ) : TypedValidationTest<T>() {
                    _values.resize( 1 );
                    _values[ 0 ] = value;
                }

                UnaryValidationTest( const UnaryValidationTest<T>& other ) {
                    _values = other._values;
                }

                virtual ~UnaryValidationTest() {}

                friend void ::swap<>( UnaryValidationTest<T>& a,
                                      UnaryValidationTest<T>& b ) noexcept;

                virtual T value( size_t n = 0 ) const override {
                    if (n > 0) {
                        throw std::out_of_range(
                            "Requested index out of range" );
                    }
                    return _values.at( 0 );
                }

                virtual const std::vector<T>& values() const {
                    return _values;
                }

            protected:
                std::vector<T> _values;
        };

    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap( NCPA::config::UnaryValidationTest<T>& a,
           NCPA::config::UnaryValidationTest<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::TypedValidationTest<T>&>( a ),
            static_cast<NCPA::config::TypedValidationTest<T>&>( b ) );
    swap( a._values, b._values );
}
