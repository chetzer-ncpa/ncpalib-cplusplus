#pragma once

#include "NCPA/strings.hpp"

#include <string>

namespace NCPA {
    namespace config {
        enum class parameter_form_t { UNDEF, SCALAR, VECTOR };

        inline std::string to_string( parameter_form_t t ) {
            switch (t) {
                case parameter_form_t::SCALAR:
                    return "scalar";
                    break;
                case parameter_form_t::VECTOR:
                    return "vector";
                    break;
                default:
                    return "undefined";
            }
        }
    }  // namespace config
}  // namespace NCPA

inline void from_string( const std::string& s,
                         NCPA::config::parameter_form_t& t ) {
    std::string lcs = NCPA::strings::to_lower( s );
    if (lcs == "scalar") {
        t = NCPA::config::parameter_form_t::SCALAR;
    } else if (lcs == "vector") {
        t = NCPA::config::parameter_form_t::VECTOR;
    } else {
        t = NCPA::config::parameter_form_t::UNDEF;
    }
}
