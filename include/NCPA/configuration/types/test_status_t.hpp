#pragma once

#include <string>

namespace NCPA {
    namespace config {
        enum class test_status_t { UNDEF, NONE, PENDING, FAILED, PASSED };

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

inline void from_string( const std::string& s,
                         NCPA::config::test_status_t& t ) {
    std::string lcs = NCPA::strings::to_lower( s );
    if (lcs == "pending") {
        t = NCPA::config::test_status_t::PENDING;
    } else if (NCPA::strings::equals_up_to(lcs,"fail",4)) {
        t = NCPA::config::test_status_t::FAILED;
    } else if (NCPA::strings::equals_up_to(lcs,"pass",4)) {
        t = NCPA::config::test_status_t::PASSED;
    } else if (lcs == "none") {
        t = NCPA::config::test_status_t::NONE;
    } else {
        t = NCPA::config::test_status_t::UNDEF;
    }
}