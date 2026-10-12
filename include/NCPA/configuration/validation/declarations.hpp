#pragma once

#include <functional>
#include <stdexcept>
#include <vector>

namespace NCPA {
    namespace config {
        enum class test_result_t { NONE = 0, PENDING, FAILED, PASSED };

        class ValidationStatus;

        template<typename T>
        using validation_function_t = std::function<ValidationStatus( T )>;

        class AbstractValidation;

        template<typename T>
        class TypedValidation;

        template<typename T>
        using validation_ptr_t = std::unique_ptr<TypedValidation<T>>;

        class Validation;

    }  // namespace config
}  // namespace NCPA
