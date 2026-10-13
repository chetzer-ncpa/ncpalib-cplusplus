#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/strings.hpp"
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        struct ArgumentParseResult {
                ArgumentParseResult( parse_result_t res,
                                     const std::string& msg ) :
                    result { res }, messages { 1, msg } {}

                ArgumentParseResult( parse_result_t res,
                                     const std::vector<std::string>& msg
                                     = {} ) :
                    result { res }, messages { msg } {}

                parse_result_t result;
                std::vector<std::string> messages;

                std::string message() const {
                    return strings::join( messages );
                }
        };
    }  // namespace config
}  // namespace NCPA
