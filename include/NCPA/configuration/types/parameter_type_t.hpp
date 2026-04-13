#pragma once

#include <string>

namespace NCPA {
    namespace config {
        enum class parameter_type_t {
            UNDEF,
            INTEGER,
            UNSIGNED_INTEGER,
            FLOAT,
            STRING,
            BOOLEAN,
            ENUM,
            COMPLEX,
            OTHER
        };

        inline std::string to_string( parameter_type_t t ) {
            switch (t) {
                case parameter_type_t::FLOAT:
                    return "float";
                    break;
                case parameter_type_t::INTEGER:
                    return "signed integer";
                    break;
                case parameter_type_t::UNSIGNED_INTEGER:
                    return "unsigned integer";
                    break;
                case parameter_type_t::STRING:
                    return "string";
                    break;
                case parameter_type_t::BOOLEAN:
                    return "boolean";
                    break;
                case parameter_type_t::COMPLEX:
                    return "complex";
                    break;
                case parameter_type_t::ENUM:
                    return "enum";
                    break;
                case parameter_type_t::OTHER:
                    return "other";
                    break;
                default:
                    return "undefined";
            }
        }
    }  // namespace config
}  // namespace NCPA
