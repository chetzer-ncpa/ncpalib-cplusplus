#pragma once

#include <stdexcept>
#include <string>

namespace NCPA {
    namespace config {
        enum class parse_status_t { UNDEF, PENDING, FAILURE, SUCCESS, NOT_FOUND };

        inline std::string to_string( parse_status_t status ) {
            switch (status) {
                case parse_status_t::PENDING:
                    return "pending";
                    break;
                case parse_status_t::FAILURE:
                    return "failure";
                    break;
                case parse_status_t::SUCCESS:
                    return "success";
                    break;
                case parse_status_t::NOT_FOUND:
                    return "not found";
                    break;
                default:
                    throw std::out_of_range( "Unrecognized or unsupported "
                                             "value of parse_status_t" );
            }
        }
    }  // namespace config
}  // namespace NCPA

inline void from_string( const std::string& s,
                         NCPA::config::parse_status_t& t ) {
    std::string lcs = NCPA::strings::to_lower( s );
    if (lcs == "pending") {
        t = NCPA::config::parse_status_t::PENDING;
    } else if (NCPA::strings::equals_up_to(lcs,"fail",4)) {
        t = NCPA::config::parse_status_t::FAILURE;
    } else if (NCPA::strings::equals_up_to(lcs,"succe",5)) {
        t = NCPA::config::parse_status_t::SUCCESS;
    } else if (lcs == "not found") {
        t = NCPA::config::parse_status_t::NOT_FOUND;
    } else {
        t = NCPA::config::parse_status_t::UNDEF;
    }
}