#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/Parameter.hpp"
#include "NCPA/configuration/validation/NullaryValidationTest.hpp"

#include <memory>
#include <type_traits>

namespace NCPA {
    namespace config {
        class WasSetTest : public NullaryValidationTest {
            public:
                WasSetTest() {}

                virtual ~WasSetTest() {}

                friend void ::swap( WasSetTest& a, WasSetTest& b ) noexcept;

                virtual std::unique_ptr<ValidationTest> clone()
                    const override {
                    return std::unique_ptr<ValidationTest>(
                        new WasSetTest( *this ) );
                }

                virtual std::string description() const override {
                    return "was set";
                }

                virtual ValidationTest& test(
                    const Parameter *param ) override {
                    if (param->was_set()) {
                        this->pass();
                    } else {
                        this->fail();
                    }
                    return *this;
                }
        };


    }  // namespace config
}  // namespace NCPA

inline void swap( NCPA::config::WasSetTest& a,
                  NCPA::config::WasSetTest& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::NullaryValidationTest&>( a ),
            static_cast<NCPA::config::NullaryValidationTest&>( b ) );
}
