#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/parameters/BaseParameter.hpp"
#include "NCPA/configuration/validation/TypedValidation.hpp"
#include "NCPA/configuration/validation/ValidationStatus.hpp"
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

                virtual void add_value( const PARAMTYPE& v )      = 0;
                virtual PARAMTYPE get( size_t n = 0 ) const       = 0;
                virtual std::vector<PARAMTYPE> get_vector() const = 0;
                virtual void set_value( const PARAMTYPE& val )    = 0;
                virtual void set_value( const std::vector<PARAMTYPE>& vals )
                    = 0;

                explicit TypedParameter( const std::string& tag = "" ) :
                    BaseParameter( tag ) {}

                TypedParameter( const std::string& tag,
                                const TypedValidation<PARAMTYPE>& v ) :
                    TypedParameter<PARAMTYPE>( tag ) {
                    _validations.push_back( v );
                }

                TypedParameter( const std::string& tag,
                                const TypedValidation<PARAMTYPE> *v ) :
                    TypedParameter<PARAMTYPE>( tag, *v ) {}

                TypedParameter(
                    const std::string& tag,
                    const std::vector<TypedValidation<PARAMTYPE>>& tests ) :
                    TypedParameter<PARAMTYPE>( tag ) {
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
                    const Validation& v ) override {
                    try {
                        auto& typed_ref
                            = dynamic_cast<const TypedValidation<PARAMTYPE>&>(
                                v );
                        _validations.push_back(
                            TypedValidation<PARAMTYPE>( typed_ref ) );
                    } catch (const std::bad_cast&) {
                        throw std::logic_error(
                            "Error in Parameter: Can't cast to TypedParameter "
                            "of proper type" );
                    }
                    return *this;
                }

                virtual ValidationStatus validate() const override {
                    ValidationStatus status;
                    for (PARAMTYPE val : this->get_vector()) {
                        status
                            = status && _run_validations( _validations, val );
                    }
                    if (status.result == test_result_t::PENDING) {
                        status.result = test_result_t::PASSED;
                    }
                    return status;
                }

                template<typename T>
                ValidationStatus _run_validations(
                    const std::vector<TypedValidation<T>>& validations,
                    const T& value ) const {
                    ValidationStatus status;
                    if (validations.empty()) {
                        status.result = test_result_t::NONE;
                    } else {
                        status.result = test_result_t::PENDING;
                        for (const TypedValidation<PARAMTYPE>& test :
                             validations) {
                            status = status && test.validate( value );
                        }
                        if (status.result == test_result_t::PENDING) {
                            status.result = test_result_t::PASSED;
                        }
                    }
                    for (auto& msg : status.messages) {
                        msg = this->tag() + ": " + msg;
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
                    using std::to_string;
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
