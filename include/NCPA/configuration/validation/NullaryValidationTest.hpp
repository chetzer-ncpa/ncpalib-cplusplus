#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/ValidationTest.hpp"

namespace NCPA {
    namespace config {
        class NullaryValidationTest : public ValidationTest {
                public:
                    NullaryValidationTest() : ValidationTest() {}

                    NullaryValidationTest(
                        const NullaryValidationTest& other ) {}

                    virtual ~NullaryValidationTest() {}

                    friend void ::swap( NullaryValidationTest& a,
                                        NullaryValidationTest& b ) noexcept;
            };

    }  // namespace config
}  // namespace NCPA

inline void swap(
    NCPA::config::NullaryValidationTest& a,
    NCPA::config::NullaryValidationTest& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::ValidationTest&>( a ),
            static_cast<NCPA::config::ValidationTest&>( b ) );
}