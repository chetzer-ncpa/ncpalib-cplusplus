#pragma once

#include "NCPA/configuration/BaseParameter.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/ScalarParameter.hpp"
#include "NCPA/configuration/ScalarParameterWithUnits.hpp"
#include "NCPA/configuration/TypedParameter.hpp"
#include "NCPA/configuration/VectorParameter.hpp"
#include "NCPA/units.hpp"

#include <memory>

namespace NCPA {
    namespace config {
        class Parameter {
            public:
                template<typename PARAMTYPE>
                static param_ptr_t scalar() {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>() );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t( new ScalarParameter<PARAMTYPE>() );
                    } else {
                        return param_ptr_t(
                            new ScalarParameterWithUnits<PARAMTYPE>( 0.0,
                                                                     ptr ) );
                    }
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar( const PARAMTYPE& defaultval ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( defaultval ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( const PARAMTYPE& defaultval,
                                           NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t(
                            new ScalarParameter<PARAMTYPE>( defaultval ) );
                    } else {
                        return param_ptr_t(
                            new ScalarParameterWithUnits<PARAMTYPE>(
                                defaultval, ptr ) );
                    }
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::vector<PARAMTYPE>& defaultval ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( defaultval ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t(
                            new ScalarParameter<PARAMTYPE>( defaultval ) );
                    } else {
                        return param_ptr_t(
                            new ScalarParameterWithUnits<PARAMTYPE>(
                                defaultval, ptr ) );
                    }
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar( const ValidationTest& newtest ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar( const test_ptr_t& newtest ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar( PARAMTYPE defaultval,
                                           const ValidationTest& newtest ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar( PARAMTYPE defaultval,
                                           const test_ptr_t& newtest ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        defaultval, new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::vector<PARAMTYPE>& defaultval,
                    const ValidationTest& newtest ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::vector<PARAMTYPE>& defaultval,
                    const test_ptr_t& newtest ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t scalar(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t( new ScalarParameter<PARAMTYPE>(
                        defaultval, new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( NCPA::units::units_ptr_t ptr,
                                           const ValidationTest& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>( ptr,
                                                                 newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( NCPA::units::units_ptr_t ptr,
                                           const test_ptr_t& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>( ptr,
                                                                 newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    NCPA::units::units_ptr_t ptr,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>( ptr,
                                                                 new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( PARAMTYPE defaultval,
                                           NCPA::units::units_ptr_t ptr,
                                           const ValidationTest& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar( PARAMTYPE defaultval,
                                           NCPA::units::units_ptr_t ptr,
                                           const test_ptr_t& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    PARAMTYPE defaultval, NCPA::units::units_ptr_t ptr,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr,
                    const ValidationTest& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr, const test_ptr_t& newtest ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t scalar(
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, new_tests ) );
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
                    const ScalarParameter<PARAMTYPE>& other,
                    NCPA::units::units_ptr_t u ) {
                    return param_ptr_t(
                        new ScalarParameterWithUnits<PARAMTYPE>( other.get(),
                                                                 u ) );
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
                static param_ptr_t vector() {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>() );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector( const PARAMTYPE& defaultval ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( defaultval ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t( new VectorParameter<PARAMTYPE>() );
                    } else {
                        return param_ptr_t(
                            new VectorParameterWithUnits<PARAMTYPE>( ptr ) );
                    }
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( const PARAMTYPE& defaultval,
                                           NCPA::units::units_ptr_t ptr ) {
                    if (ptr == nullptr) {
                        return param_ptr_t(
                            new VectorParameter<PARAMTYPE>( defaultval ) );
                    } else {
                        return param_ptr_t(
                            new VectorParameterWithUnits<PARAMTYPE>(
                                defaultval, ptr ) );
                    }
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::vector<PARAMTYPE>& defaultval ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( defaultval ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( defaultval,
                                                                 ptr ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const NCPA::units::VectorWithUnits<PARAMTYPE>& defaultval ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( defaultval ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector( const ValidationTest& newtest ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector( const test_ptr_t& newtest ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector( PARAMTYPE defaultval,
                                           const ValidationTest& newtest ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector( PARAMTYPE defaultval,
                                           const test_ptr_t& newtest ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        defaultval, new_tests ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::vector<PARAMTYPE>& defaultval,
                    const ValidationTest& newtest ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::vector<PARAMTYPE>& defaultval,
                    const test_ptr_t& newtest ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        defaultval, newtest ) );
                }

                template<typename PARAMTYPE>
                static param_ptr_t vector(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t( new VectorParameter<PARAMTYPE>(
                        defaultval, new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( NCPA::units::units_ptr_t ptr,
                                           const ValidationTest& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( ptr,
                                                                 newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( NCPA::units::units_ptr_t ptr,
                                           const test_ptr_t& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( ptr,
                                                                 newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    NCPA::units::units_ptr_t ptr,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( ptr,
                                                                 new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( PARAMTYPE defaultval,
                                           NCPA::units::units_ptr_t ptr,
                                           const ValidationTest& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector( PARAMTYPE defaultval,
                                           NCPA::units::units_ptr_t ptr,
                                           const test_ptr_t& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    PARAMTYPE defaultval, NCPA::units::units_ptr_t ptr,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, new_tests ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr,
                    const ValidationTest& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr, const test_ptr_t& newtest ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, newtest ) );
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                static param_ptr_t vector(
                    const std::vector<PARAMTYPE>& defaultval,
                    NCPA::units::units_ptr_t ptr,
                    std::initializer_list<test_ptr_t> new_tests ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>(
                            defaultval, ptr, new_tests ) );
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
                    const VectorParameter<PARAMTYPE>& other,
                    NCPA::units::units_ptr_t u ) {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( other.get(),
                                                                 u ) );
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
