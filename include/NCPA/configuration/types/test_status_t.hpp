#pragma once

#include <string>

namespace NCPA {
    namespace config {
        enum class test_status_t { NONE, PENDING, FAILED, PASSED };

        inline std::string to_string( test_status_t t ) {
            switch (t) {
                case test_status_t::NONE:
                    return "none";
                    break;
                case test_status_t::PENDING:
                    return "pending";
                    break;
                case test_status_t::FAILED:
                    return "failed";
                    break;
                case test_status_t::PASSED:
                    return "passed";
                    break;
                default:
                    return "undefined";
            }
        }
    }  // namespace config
}  // namespace NCPA
