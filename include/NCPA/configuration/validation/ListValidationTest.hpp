#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/declarations.hpp"
#include "NCPA/configuration/validation/TypedValidationTest.hpp"
#include <vector>

namespace NCPA {
    namespace config {
        template<typename T>
            class ListValidationTest : public TypedValidationTest<T> {
                public:
                    ListValidationTest( std::vector<T> values ) :
                        TypedValidationTest<T>(), _values { values } {}

                    ListValidationTest( std::initializer_list<T> values ) :
                        TypedValidationTest<T>(), _values { values } {}

                    ListValidationTest( const ListValidationTest<T>& other ) {
                        _values = other._values;
                    }

                    virtual ~ListValidationTest() {}

                    friend void ::swap<>( ListValidationTest<T>& a,
                                          ListValidationTest<T>& b ) noexcept;

                    virtual T value( size_t n ) const override {
                        if (n >= _values.size()) {
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
            };

    }  // namespace config
}  // namespace NCPA
