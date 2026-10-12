#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/parameters/BaseParameter.hpp"
#include "NCPA/configuration/parameters/ScalarParameter.hpp"
#include "NCPA/configuration/parameters/ScalarParameterWithUnits.hpp"
#include "NCPA/configuration/parameters/TypedParameter.hpp"
#include "NCPA/configuration/parameters/VectorParameter.hpp"
#include "NCPA/configuration/validation/Validation.hpp"
#include "NCPA/units.hpp"

#include <memory>

namespace NCPA {
    namespace config {
        class Parameter {
            public:
                template<typename PARAMTYPE>
                static param_ptr_t scalar( const std::string& tag ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( tag ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( const std::string& tag,
                                           NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t(
                            new ScalarParameter<PARAMTYPE>( tag ) );
                    } else {
                        return param_ptr_t(
                            new ScalarParameterWithUnits<PARAMTYPE>( tag, 0.0,
                                                                     ptr ) );
                    }
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar( const std::string& tag,
                                           const PARAMTYPE& defaultval ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( tag, defaultval ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( const std::string& tag,
                                           const PARAMTYPE& defaultval,
                                           NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                            tag, defaultval ) );
                    } else {
                        return param_ptr_t(
                            new ScalarParameterWithUnits<PARAMTYPE>(
                                tag, defaultval, ptr ) );
                    }
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( tag, defaultval ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                            tag, defaultval ) );
                    } else {
                        return param_ptr_t(
                            new ScalarParameterWithUnits<PARAMTYPE>(
                                tag, defaultval, ptr ) );
                    }
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const TypedValidation<PARAMTYPE>& newtest ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( tag, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( tag, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::string& tag,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( tag, new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar( const std::string& tag,
                                           PARAMTYPE defaultval,
                                           const Validation& newtest ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        tag, defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::string& tag, PARAMTYPE defaultval,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        tag, defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::string& tag, PARAMTYPE defaultval,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        tag, defaultval, new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    const Validation& newtest ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        tag, defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        tag, defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        tag, defaultval, new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( const std::string& tag,
                                           NCPA::units::units_ptr_t ptr,
                                           const Validation& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>( tag, ptr,
                                                                 newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::string& tag, NCPA::units::units_ptr_t ptr,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>( tag, ptr,
                                                                 newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::string& tag, NCPA::units::units_ptr_t ptr,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>( tag, ptr,
                                                                 new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( const std::string& tag,
                                           PARAMTYPE defaultval,
                                           NCPA::units::units_ptr_t ptr,
                                           const Validation& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::string& tag, PARAMTYPE defaultval,
                    NCPA::units::units_ptr_t ptr,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::string& tag, PARAMTYPE defaultval,
                    NCPA::units::units_ptr_t ptr,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr, const Validation& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const ScalarParameter<PARAMTYPE>& other ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( other ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::string& tag,
                    const ScalarParameter<PARAMTYPE>& other,
                    NCPA::units::units_ptr_t u ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            tag, other.get(), u ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const ScalarParameterWithUnits<PARAMTYPE>& other ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>( other ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector( const std::string& tag ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( tag ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector( const std::string& tag,
                                           const PARAMTYPE& defaultval ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( tag, defaultval ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( const std::string& tag,
                                           NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t(
                            new VectorParameter<PARAMTYPE>( tag ) );
                    } else {
                        return param_ptr_t(
                            new VectorParameterWithUnits<PARAMTYPE>( tag,
                                                                     ptr ) );
                    }
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( const std::string& tag,
                                           const PARAMTYPE& defaultval,
                                           NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t( new VectorParameter<PARAMTYPE>(
                            tag, defaultval ) );
                    } else {
                        return param_ptr_t(
                            new VectorParameterWithUnits<PARAMTYPE>(
                                tag, defaultval, ptr ) );
                    }
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( tag, defaultval ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag,
                    const NCPA::units::VectorWithUnits<PARAMTYPE>&
                        defaultval ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            tag, defaultval ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector( const std::string& tag,
                                           const Validation& newtest ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( tag, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::string& tag,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( tag, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::string& tag,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( tag, new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector( const std::string& tag,
                                           PARAMTYPE defaultval,
                                           const Validation& newtest ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        tag, defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::string& tag, PARAMTYPE defaultval,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        tag, defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::string& tag, PARAMTYPE defaultval,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        tag, defaultval, new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    const Validation& newtest ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        tag, defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        tag, defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        tag, defaultval, new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( const std::string& tag,
                                           NCPA::units::units_ptr_t ptr,
                                           const Validation& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( tag, ptr,
                                                                 newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag, NCPA::units::units_ptr_t ptr,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( tag, ptr,
                                                                 newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag, NCPA::units::units_ptr_t ptr,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( tag, ptr,
                                                                 new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( const std::string& tag,
                                           PARAMTYPE defaultval,
                                           NCPA::units::units_ptr_t ptr,
                                           const Validation& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag, PARAMTYPE defaultval,
                    NCPA::units::units_ptr_t ptr,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag, PARAMTYPE defaultval,
                    NCPA::units::units_ptr_t ptr,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr, const Validation& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr,
                    const validation_ptr_t<PARAMTYPE>& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr,
                    std::initializer_list<validation_ptr_t<PARAMTYPE>>
                        new_tests ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            tag, defaultval, ptr, new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const VectorParameter<PARAMTYPE>& other ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( other ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::string& tag,
                    const VectorParameter<PARAMTYPE>& other,
                    NCPA::units::units_ptr_t u ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            tag, other.get(), u ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const VectorParameterWithUnits<PARAMTYPE>& other ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( other ) );
                }
        };
    }  // namespace config
}  // namespace NCPA
