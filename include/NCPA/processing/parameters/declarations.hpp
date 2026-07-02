#pragma once

#include <memory>
#include <string>
#include <unordered_map>

namespace NCPA {
    namespace processing {
        enum class parameter_type_t { INTEGER, FLOAT, STRING, BOOLEAN, ENUM };

        enum class parameter_form_t { SCALAR, ARRAY };

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

    }
}