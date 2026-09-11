#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/TypedValidation.hpp"
#include "NCPA/configuration/Validation.hpp"
#include "NCPA/configuration/functions.hpp"

#include <memory>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        template<typename BASETYPE>
        class Validated {
            public:
                virtual validation_status_t validate() const = 0;

                Validated() {}

                virtual ~Validated() {}

                virtual BASETYPE& add_test( const Validation& v ) {
                    _validations.push_back( v.clone() );
                    return static_cast<BASETYPE&>( *this );
                }

                virtual bool failed() const {
                    validation_status_t res = this->validate();
                    return ( res.result == test_result_t::FAILED );
                }

                virtual bool passed() const {
                    validation_status_t res = this->validate();
                    return ( res.result == test_result_t::NONE
                             || res.result == test_result_t::PASSED );
                }

                template<typename PARAMTYPE>
                validation_status_t validate_value( const PARAMTYPE& value ) const {
                    validation_status_t status;
                    if (_validations.empty()) {
                        status.result = test_result_t::NONE;
                    } else {
                        status.result = test_result_t::PENDING;
                        for (const std::unique_ptr<Validation>& gtest :
                             _validations) {
                            try {
                                auto& test = dynamic_cast<
                                    const TypedValidation<PARAMTYPE>&>(
                                    *gtest );
                                validation_status_t teststatus
                                    = test.validate( value );
                                if (teststatus.result
                                    == test_result_t::FAILED) {
                                    status.result = test_result_t::FAILED;
                                    status.message
                                        += "\n" + teststatus.message;
                                }
                            } catch (const std::bad_cast&) {
                                throw std::logic_error(
                                    "Error in validation: Can't cast test to "
                                    "proper type" );
                            }
                        }
                        if (status.result == test_result_t::PENDING) {
                            status.result = test_result_t::PASSED;
                        }
                    }
                    return status;
                }

                size_t validation_count() const noexcept {
                    return _validations.size();
                }


            private:
                std::vector<std::unique_ptr<Validation>> _validations;
        };
    }  // namespace config
}  // namespace NCPA
