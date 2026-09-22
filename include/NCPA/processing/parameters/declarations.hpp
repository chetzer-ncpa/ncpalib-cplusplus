#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace NCPA {
    namespace processing {
        enum class parameter_type_t { INTEGER, FLOAT, STRING, BOOLEAN, ENUM };

        enum class parameter_form_t { SCALAR, ARRAY };

        struct passthrough_parameter_t {
                passthrough_parameter_t( std::string t, std::string in,
                                         std::string out ) :
                    tag { t }, in_key { in }, out_key { out } {}

                std::string tag;
                std::string in_key;
                std::string out_key;
        };

        class Parameter;
        template<typename T>
        class ScalarParameter;
        template<typename T>
        class VectorParameter;

        class IntegerParameter;
        class DoubleParameter;
        class StringParameter;
        class BooleanParameter;
        class IntegerVectorParameter;
        class DoubleVectorParameter;
        class StringVectorParameter;
        class BooleanVectorParameter;

        class ParameterTree;

        typedef std::unique_ptr<Parameter> parameter_ptr_t;
        typedef std::unordered_map<std::string, std::vector<parameter_ptr_t>>
            parameter_tree_t;

    }  // namespace processing
}  // namespace NCPA
