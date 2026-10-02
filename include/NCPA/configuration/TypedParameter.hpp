#pragma once

#include "NCPA/configuration/BaseParameter.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/TypedValidation.hpp"
#include "NCPA/types.hpp"
#include "NCPA/units.hpp"

#include <limits>
#include <type_traits>

namespace NCPA {
    namespace config {

        template<typename PARAMTYPE>
        class TypedParameter : public BaseParameter {
            public:
                using value_type = PARAMTYPE;

                explicit TypedParameter() : BaseParameter() {}

                TypedParameter( const std::string& tag         = "",
                                const std::string& description = "" ) :
                    BaseParameter( tag, description ) {}

                TypedParameter( const TypedValidation<PARAMTYPE>& v ) :
                    TypedParameter<PARAMTYPE>() {
                    _validations.push_back( v );
                }

                TypedParameter( const TypedValidation<PARAMTYPE> *v ) :
                    TypedParameter<PARAMTYPE>( *v ) {}

                TypedParameter(
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    TypedParameter<PARAMTYPE>() {
                    for (auto& test : tests) {
                        _validations.push_back( test );
                    }
                }

                virtual ~TypedParameter() {}

                TypedParameter( const TypedParameter<PARAMTYPE>& other ) :
                    BaseParameter( other ) {
                    _validations = other._validations;
                }

                TypedParameter( TypedParameter<PARAMTYPE>&& other ) noexcept :
                    TypedParameter() {
                    swap( *this, other );
                }

                TypedParameter<PARAMTYPE>& operator=(
                    TypedParameter<PARAMTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( TypedParameter<PARAMTYPE>& a,
                                  TypedParameter<PARAMTYPE>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<NCPA::config::BaseParameter&>( a ),
                          static_cast<NCPA::config::BaseParameter&>( b ) );
                    swap( a._validations, b._validations );
                }

                virtual parameter_type_t type() const override {
                    return parameter_type<PARAMTYPE>();
                }

                virtual PARAMTYPE get( size_t n = 0 ) const = 0;

                virtual std::vector<PARAMTYPE> get_vector() const = 0;

                template<typename ASTYPE,
                         typename std::enable_if<
                             std::is_convertible<PARAMTYPE, ASTYPE>::value,
                             int>::type = 0>
                ASTYPE as() const {
                    return static_cast<ASTYPE>( this->get() );
                }

                template<typename ASTYPE,
                         typename std::enable_if<
                             std::is_convertible<PARAMTYPE, ASTYPE>::value,
                             int>::type = 0>
                ScalarWithUnits<ASTYPE> as_with_units() const {
                    return ScalarWithUnits<ASTYPE>(
                        static_cast<ASTYPE>( this->get() ),
                        this->get_units() );
                }

                template<
                    typename ASTYPE,
                    typename std::enable_if<
                        ( !( std::is_convertible<PARAMTYPE, ASTYPE>::value ) ),
                        int>::type = 0>
                ScalarWithUnits<ASTYPE> as_with_units() const {
                    return ScalarWithUnits<ASTYPE>(
                        static_cast<ASTYPE>( this->as_double() ),
                        this->get_units() );
                }

                virtual TypedParameter<PARAMTYPE>& add_test(
                    Validation&& v ) override {
                    try {
                        auto& typed_ref
                            = dynamic_cast<TypedValidation<PARAMTYPE>&>( v );
                        _validations.push_back( std::move( typed_ref ) );
                    } catch (const std::bad_cast&) {
                        throw std::logic_error(
                            "Error in Parameter: Can't cast to TypedParameter "
                            "of proper type" );
                    }
                    return *this;
                }

                virtual TypedParameter<PARAMTYPE>& add_test(
                    Validation& v ) override {
                    try {
                        auto& typed_ref
                            = dynamic_cast<TypedValidation<PARAMTYPE>&>( v );
                        _validations.push_back(
                            TypedValidation<PARAMTYPE>( typed_ref ) );
                    } catch (const std::bad_cast&) {
                        throw std::logic_error(
                            "Error in Parameter: Can't cast to TypedParameter "
                            "of proper type" );
                    }
                    return *this;
                }

                virtual validation_status_t validate(
                    bool short_circuit = false ) const override {
                    validation_status_t status;
                    if (_validations.empty()) {
                        status.result = test_result_t::NONE;
                    } else {
                        status.result = test_result_t::PENDING;
                        for (const TypedValidation<PARAMTYPE>& test :
                             _validations) {
                            validation_status_t teststatus
                                = test.validate( this->get() );
                            switch (teststatus.result) {
                                case test_result_t::FAILED:
                                    status.result = test_result_t::FAILED;
                                    status.message
                                        += "\n" + teststatus.message;
                                    if (short_circuit) {
                                        return status;
                                    }
                                    break;
                                default:
                                    break;
                            }
                        }
                        if (status.result == test_result_t::PENDING) {
                            status.result = test_result_t::PASSED;
                        }
                    }
                    return status;
                }

            protected:
                std::vector<TypedValidation<PARAMTYPE>> _validations;

                template<typename T = PARAMTYPE,
                         typename std::enable_if<
                             NCPA::types::has_to_string<T>::value, int>::type
                         = 0>
                std::string _as_string( size_t n            = 0,
                                        bool throw_on_error = false ) const {
                    return to_string( this->get( n ) );
                }

                template<typename T = PARAMTYPE,
                         typename std::enable_if<
                             !( NCPA::types::has_to_string<T>::value ),
                             int>::type = 0>
                std::string _as_string( size_t n            = 0,
                                        bool throw_on_error = false ) const {
                    if (throw_on_error) {
                        throw std::out_of_range(
                            "No as_string() conversion defined!" );
                    } else {
                        return "<no string conversion defined>";
                    }
                }
        };
    }  // namespace config
}  // namespace NCPA
