#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/declarations.hpp"
#include "NCPA/configuration/validation/ValidationTest.hpp"

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
