#pragma once

#include "NCPA/strings.hpp"

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

inline void from_string( const std::string& s, NCPA::config::parameter_type_t& t ) {
            using NCPA::strings::equals_up_to;

            std::string lcs = NCPA::strings::to_lower( s );
            if (equals_up_to( lcs, "float", 5 )
                || equals_up_to( lcs, "double", 6 )) {
                t = NCPA::config::parameter_type_t::FLOAT;
            } else if (equals_up_to( lcs, "int", 3 )) {
                t = NCPA::config::parameter_type_t::INTEGER;
            } else if (lcs == "unsigned integer" || lcs == "size_t"
                       || equals_up_to( lcs, "uint", 4 )) {
                t = NCPA::config::parameter_type_t::UNSIGNED_INTEGER;
            } else if (lcs == "string" || equals_up_to( lcs, "char", 4 )
                       || lcs == "std::string" || lcs == "text") {
                t = NCPA::config::parameter_type_t::STRING;
            } else if (equals_up_to( lcs, "bool", 4 )) {
                t = NCPA::config::parameter_type_t::BOOLEAN;
            } else if (equals_up_to( lcs, "enum", 4 )) {
                t = NCPA::config::parameter_type_t::ENUM;
            } else if (equals_up_to( lcs, "complex", 7 )
                       || equals_up_to( lcs, "std::complex", 12 )) {
                t = NCPA::config::parameter_type_t::COMPLEX;
            } else {
                t = NCPA::config::parameter_type_t::OTHER;
            }
        }