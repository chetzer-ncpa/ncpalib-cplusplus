#pragma once

namespace NCPA {
    namespace params {
        template<typename TESTTYPE>
        class Parameter;

        template<typename TESTTYPE>
        class Validator;

        template<typename TESTTYPE>
        class RangeValidator;

        namespace details {
            template<typename TESTTYPE>
            class ValidationTest;

            template<typename TESTTYPE>
            class GreaterTest;

            template<typename TESTTYPE>
            class LessTest;

            template<typename TESTTYPE>
            class EqualTest;

            template<typename TESTTYPE>
            class OneOfTest;

        }  // namespace details

    }  // namespace params
}  // namespace NCPA
