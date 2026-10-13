#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        // yes, this is basically a bool, but then it gets confused with
        // default int values
        enum class value_count_t { ZERO, ZERO_OR_ONE, ONE, ANY, AT_LEAST_ONE };
        enum class parse_result_t {
            NO_ARGUMENTS,
            SUCCESS,
            BARE_VALUES,
            PARSE_ERROR,
            HELP_REQUESTED
        };

        struct ArgumentParseResult;
        class Argument;

        template<typename T>
        class TypedArgument;
        template<typename T>
        class MultiValueArgument;
        template<typename T>
        class RequiredMultiValueArgument;
        template<typename T>
        class RequiredSingleValueArgument;
        template<typename T>
        class SingleValueArgument;
        class FlagArgument;

        typedef std::unique_ptr<Argument> argument_ptr_t;
        typedef std::map<std::string, argument_ptr_t> argument_map_t;

        class ArgumentSet;
    }  // namespace config
}  // namespace NCPA
