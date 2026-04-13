#pragma once

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
