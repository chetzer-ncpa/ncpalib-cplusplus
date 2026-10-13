#pragma once

#include "NCPA/configuration/declarations.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        class ValidationStatus {
            public:
                test_result_t result;
                std::vector<std::string> messages;

                ValidationStatus() :
                    ValidationStatus( test_result_t::PENDING,
                                      std::vector<std::string> {} ) {}

                ValidationStatus( test_result_t res ) : result { res } {}

                ValidationStatus( test_result_t res,
                                  const std::vector<std::string>& msg ) :
                    result { res }, messages { msg } {}

                ValidationStatus( test_result_t res, const std::string& msg ) :
                    result { res }, messages { 1, msg } {}

                ValidationStatus( const ValidationStatus& other ) = default;
                ValidationStatus( ValidationStatus&& other ) noexcept
                    = default;
                ValidationStatus& operator=( ValidationStatus status ) {
                    swap(*this,status);
                    return *this;
                }

                friend void swap( ValidationStatus& a,
                                  ValidationStatus& b ) noexcept {
                    using std::swap;
                    swap( a.result, b.result );
                    swap( a.messages, b.messages );
                }

                ValidationStatus& add_message( const std::string& msg ) {
                    messages.push_back( msg );
                    return *this;
                }

                void fail( const std::string& msg ) {
                    result = test_result_t::FAILED;
                    this->add_message( msg );
                }

                const std::string message() const {
                    std::ostringstream oss;
                    bool first = true;
                    for (auto& s : messages) {
                        if (!first) {
                            oss << "\n";
                            first = false;
                        }
                        oss << s;
                    }
                    return oss.str();
                }

                ValidationStatus& operator&=( const ValidationStatus& other ) {
                    ValidationStatus status = *this && other;
                    swap( *this, status );
                    return *this;
                }

                ValidationStatus operator&&(
                    const ValidationStatus& other ) const {
                    ValidationStatus v( *this );
                    switch (other.result) {
                        case test_result_t::NONE:
                            break;
                        case test_result_t::FAILED:
                            v.result = test_result_t::FAILED;
                            break;
                        case test_result_t::PENDING:
                            v.result = ( v.result == test_result_t::FAILED
                                             ? test_result_t::FAILED
                                             : test_result_t::PENDING );
                            break;
                        case test_result_t::PASSED:
                            if (v.result == test_result_t::NONE
                                || v.result == test_result_t::PASSED) {
                                v.result = test_result_t::PASSED;
                            }
                            break;
                        default:
                            throw std::out_of_range(
                                "Unsupported or unrecognized test status "
                                "received" );
                    }
                    v.messages.insert( v.messages.end(),
                                       other.messages.cbegin(),
                                       other.messages.cend() );
                    return v;
                }
        };
    }  // namespace config
}  // namespace NCPA
