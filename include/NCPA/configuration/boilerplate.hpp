#pragma once

#ifdef NCPA_CONFIGURATION_SCALARPARAMETER_PUBLIC_BOILERPLATE
#  undef NCPA_CONFIGURATION_SCALARPARAMETER_PUBLIC_BOILERPLATE
#endif

#ifdef NCPA_CONFIGURATION_VECTORPARAMETER_PUBLIC_BOILERPLATE
#  undef NCPA_CONFIGURATION_VECTORPARAMETER_PUBLIC_BOILERPLATE
#endif

// #define NCPA_CONFIGURATION_SCALARPARAMETER_PUBLIC_BOILERPLATE                 \
//     ScalarParameter() : hidden::_base_scalar_parameter<PARAMTYPE>() {}        \
//     ScalarParameter(PARAMTYPE defaultval)                                     \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(defaultval) {}            \
//     ScalarParameter(const std::vector<PARAMTYPE>& defaultval)                 \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(defaultval) {}            \
//     ScalarParameter(const ValidationTest& newtest)                            \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(newtest) {}               \
//     ScalarParameter(const test_ptr_t& newtest)                                \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(newtest) {}               \
//     ScalarParameter(std::initializer_list<test_ptr_t> new_tests)              \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(new_tests) {}             \
//     ScalarParameter(PARAMTYPE defaultval, const ValidationTest& newtest)      \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(defaultval, newtest) {}   \
//     ScalarParameter(PARAMTYPE defaultval, const test_ptr_t& newtest)          \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(defaultval, newtest) {}   \
//     ScalarParameter(PARAMTYPE defaultval,                                     \
//                     std::initializer_list<test_ptr_t> new_tests)              \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(defaultval, new_tests) {} \
//     ScalarParameter(const std::vector<PARAMTYPE>& defaultval,                 \
//                     const ValidationTest& newtest)                            \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(defaultval, newtest) {}   \
//     ScalarParameter(const std::vector<PARAMTYPE>& defaultval,                 \
//                     const test_ptr_t& newtest)                                \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(defaultval, newtest) {}   \
//     ScalarParameter(const std::vector<PARAMTYPE>& defaultval,                 \
//                     std::initializer_list<test_ptr_t> new_tests)              \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(defaultval, new_tests) {} \
//     ScalarParameter(const ScalarParameter<PARAMTYPE>& other)                  \
//         : hidden::_base_scalar_parameter<PARAMTYPE>(other) {}                 \
//     ScalarParameter(ScalarParameter<PARAMTYPE>&& other) noexcept              \
//         : hidden::_base_scalar_parameter<PARAMTYPE>() {                       \
//         ::swap(*this, other);                                                 \
//     }                                                                         \
//     virtual ~ScalarParameter() {}                                             \
//     ScalarParameter<PARAMTYPE>& operator=(ScalarParameter<PARAMTYPE> other) { \
//         ::swap(*this, other);                                                 \
//         return *this;                                                         \
//     }                                                                         \
//     friend void ::swap<>(ScalarParameter<PARAMTYPE> & a,                      \
//                          ScalarParameter<PARAMTYPE> & b) noexcept;            \
//     virtual param_ptr_t clone() const override {                              \
//         return param_ptr_t(new ScalarParameter<PARAMTYPE>(*this));            \
//     }


#define NCPA_CONFIGURATION_VECTORPARAMETER_PUBLIC_BOILERPLATE                 \
    VectorParameter() : hidden::_base_vector_parameter<PARAMTYPE>() {}        \
    VectorParameter(PARAMTYPE defaultval)                                     \
        : hidden::_base_vector_parameter<PARAMTYPE>(defaultval) {}            \
    VectorParameter(const std::vector<PARAMTYPE>& defaultval)                 \
        : hidden::_base_vector_parameter<PARAMTYPE>(defaultval) {}            \
    VectorParameter(const ValidationTest& newtest)                            \
        : hidden::_base_vector_parameter<PARAMTYPE>(newtest) {}               \
    VectorParameter(const test_ptr_t& newtest)                                \
        : hidden::_base_vector_parameter<PARAMTYPE>(newtest) {}               \
    VectorParameter(std::initializer_list<test_ptr_t> new_tests)              \
        : hidden::_base_vector_parameter<PARAMTYPE>(new_tests) {}             \
    VectorParameter(PARAMTYPE defaultval, const ValidationTest& newtest)      \
        : hidden::_base_vector_parameter<PARAMTYPE>(defaultval, newtest) {}   \
    VectorParameter(PARAMTYPE defaultval, const test_ptr_t& newtest)          \
        : hidden::_base_vector_parameter<PARAMTYPE>(defaultval, newtest) {}   \
    VectorParameter(PARAMTYPE defaultval,                                     \
                    std::initializer_list<test_ptr_t> new_tests)              \
        : hidden::_base_vector_parameter<PARAMTYPE>(defaultval, new_tests) {} \
    VectorParameter(const std::vector<PARAMTYPE>& defaultval,                 \
                    const ValidationTest& newtest)                            \
        : hidden::_base_vector_parameter<PARAMTYPE>(defaultval, newtest) {}   \
    VectorParameter(const std::vector<PARAMTYPE>& defaultval,                 \
                    const test_ptr_t& newtest)                                \
        : hidden::_base_vector_parameter<PARAMTYPE>(defaultval, newtest) {}   \
    VectorParameter(const std::vector<PARAMTYPE>& defaultval,                 \
                    std::initializer_list<test_ptr_t> new_tests)              \
        : hidden::_base_vector_parameter<PARAMTYPE>(defaultval, new_tests) {} \
    VectorParameter(const VectorParameter<PARAMTYPE>& other)                  \
        : hidden::_base_vector_parameter<PARAMTYPE>(other) {}                 \
    VectorParameter(VectorParameter<PARAMTYPE>&& other) noexcept              \
        : hidden::_base_vector_parameter<PARAMTYPE>() {                       \
        ::swap(*this, other);                                                 \
    }                                                                         \
    virtual ~VectorParameter() {}                                             \
    VectorParameter<PARAMTYPE>& operator=(VectorParameter<PARAMTYPE> other) { \
        ::swap(*this, other);                                                 \
        return *this;                                                         \
    }                                                                         \
    friend void ::swap<>(VectorParameter<PARAMTYPE> & a,                      \
                         VectorParameter<PARAMTYPE> & b) noexcept;            \
    virtual param_ptr_t clone() const override {                              \
        return param_ptr_t(new VectorParameter<PARAMTYPE>(*this));            \
    }


    