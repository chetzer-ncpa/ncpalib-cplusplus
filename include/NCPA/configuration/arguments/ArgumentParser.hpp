#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/arguments/Argument.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <regex>

// #define TAG_REGEX R""
// #define VALUE_WITH_UNITS_REGEX R"^(-?\d?\.?\d+)(Dw*)$"

namespace NCPA {
    namespace config {
        class ArgumentParser : public CloneBase<ArgumentParser> {
            public:
                

                enum class parser_expecting_t { TAG, VALUE, EITHER };
                enum class token_t { NONE, TAG, VALUE, VALUE_WITH_UNITS };

                struct Token {
                        token_t token_type;
                        std::string text;
                        std::string units_text = "";
                };

                typedef std::pair<std::string, std::vector<std::string>>
                    arg_value_set_t;

                ArgumentParser()                                  = default;
                virtual ~ArgumentParser()                         = default;
                ArgumentParser( const ArgumentParser& other )     = default;
                ArgumentParser( ArgumentParser&& other ) noexcept = default;

                friend void swap( ArgumentParser& a,
                                  ArgumentParser& b ) noexcept {}

                virtual ArgumentParseResult parse( ArgumentSet& argset ) const
                    = 0;
        };


    }  // namespace config
}  // namespace NCPA
