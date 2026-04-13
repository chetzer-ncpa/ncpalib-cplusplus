#pragma once

#include <stdexcept>
#include <string>

namespace NCPA {
    namespace config {
        enum class parse_status_t { PENDING, FAILURE, SUCCESS };

        std::string to_string( parse_status_t status ) {
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
                default:
                    throw std::out_of_range( "Unrecognized or unsupported "
                                             "value of parse_status_t" );
            }
        }
    }  // namespace config
}  // namespace NCPA
