#pragma once

#include <memory>
#include <string>
#include <stdexcept>

#if __has_include( "nlohmann/json.hpp" )
#  define HAVE_NLOHMANN_JSON_HPP 1
#else
#  define HAVE_NLOHMANN_JSON_HPP 0
#endif

namespace NCPA {
    namespace config {
        enum class parameter_type_t { UNDEF, INTEGER, FLOAT, STRING, BOOLEAN, ENUM };
        inline std::string to_string( parameter_type_t form ) {
            switch (form) {
                case parameter_type_t::UNDEF:
                    return "undefined";
                    break;
                case parameter_type_t::INTEGER:
                    return "integer";
                    break;
                case parameter_type_t::FLOAT:
                    return "float";
                    break;
                case parameter_type_t::STRING:
                    return "string";
                    break;
                case parameter_type_t::BOOLEAN:
                    return "boolean";
                    break;
                case parameter_type_t::ENUM:
                    return "enum";
                    break;
                default:
                    throw std::out_of_range( "Unrecognized or unsupported parameter_type_t value!");
            }
        }

        enum class parameter_form_t { UNDEF, SCALAR, ARRAY };
        inline std::string to_string( parameter_form_t form ) {
            switch (form) {
                case parameter_form_t::UNDEF:
                    return "undefined";
                    break;
                case parameter_form_t::SCALAR:
                    return "scalar";
                    break;
                case parameter_form_t::ARRAY:
                    return "array";
                    break;
                default:
                    throw std::out_of_range( "Unrecognized or unsupported parameter_form_t value!");
            }
        }

        class Parameter;
        template<typename T>
        class ScalarParameter;
        template<typename T>
        class VectorParameter;

        class IntegerParameter;
        class DoubleParameter;
        class StringParameter;
        class BooleanParameter;
        class EnumParameter;
        class IntegerVectorParameter;
        class DoubleVectorParameter;
        class StringVectorParameter;
        class BooleanVectorParameter;
        class EnumVectorParameter;

        typedef std::unique_ptr<Parameter> parameter_ptr_t;
    }  // namespace config
}  // namespace NCPA
