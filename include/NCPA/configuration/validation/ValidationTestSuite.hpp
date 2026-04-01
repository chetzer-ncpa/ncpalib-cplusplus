#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation/declarations.hpp"
#include "NCPA/configuration/validation/ValidationTest.hpp"

#include <vector>

namespace NCPA {
    namespace config {
        class ValidationTestSuite {
            public:
                ValidationTestSuite() {}

                ValidationTestSuite( const ValidationTestSuite& other ) :
                    ValidationTestSuite() {
                    for (auto it = other._tests.begin();
                         it != other._tests.end(); ++it) {
                        this->append( **it );
                    }
                }

                ValidationTestSuite( ValidationTestSuite&& other ) noexcept {
                    ::swap( *this, other );
                }

                virtual ~ValidationTestSuite() {}

                ValidationTestSuite& operator=( ValidationTestSuite other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( ValidationTestSuite& a,
                                    ValidationTestSuite& b ) noexcept;

                virtual ValidationTestSuite& append(
                    const ValidationTest& newtest ) {
                    _tests.push_back( newtest.clone() );
                    return *this;
                }

                virtual ValidationTestSuite& append(
                    const ValidationTest *newtest ) {
                    _tests.push_back( newtest->clone() );
                    return *this;
                }

                virtual ValidationTestSuite& append(
                    const test_ptr_t& newtest ) {
                    _tests.push_back( newtest->clone() );
                    return *this;
                }

                virtual bool failed() const {
                    return ( this->status() == test_status_t::FAILED );
                }

                virtual std::vector<const ValidationTest *> failed_tests()
                    const {
                    std::vector<const ValidationTest *> failed_tests;
                    for (auto it = _tests.cbegin(); it != _tests.cend();
                         ++it) {
                        if (( *it )->failed()) {
                            failed_tests.push_back( it->get() );
                        }
                    }
                    return failed_tests;
                }

                virtual bool passed() const {
                    return ( this->status() == test_status_t::PASSED );
                }

                virtual bool pending() const {
                    return ( this->status() == test_status_t::PENDING );
                }

                virtual ValidationTestSuite& prepend(
                    const ValidationTest& newtest ) {
                    _tests.emplace( _tests.begin(), newtest.clone() );
                    return *this;
                }

                virtual ValidationTestSuite& prepend(
                    const test_ptr_t& newtest ) {
                    _tests.insert( _tests.begin(), newtest->clone() );
                    return *this;
                }

                virtual ValidationTestSuite& run_tests(
                    const Parameter *param,
                    bool short_circuit = false ) {
                    // _status = test_status_t::PENDING;
                    bool pass = true;
                    for (auto it = _tests.begin(); it != _tests.end(); ++it) {
                        pass = ( ( *it )->test( param ).passed() ) && pass;
                        if (short_circuit && !pass) {
                            return *this;
                        }
                    }
                    return *this;
                }

                virtual size_t size() const { return _tests.size(); }

                // if any tests have failed, return FAILED.  Otherwise, if any
                // tests are pending, return PENDING.  Otherwise return PASSED.
                virtual test_status_t status() const {
                    if (_tests.empty()) {
                        return test_status_t::PASSED;
                    }
                    bool pending = false;
                    for (auto it = _tests.cbegin(); it != _tests.cend();
                         ++it) {
                        switch (( *it )->status()) {
                            case test_status_t::FAILED:
                                return test_status_t::FAILED;
                                break;
                            case test_status_t::PENDING:
                                pending = true;
                                break;
                            case test_status_t::PASSED:
                                break;
                            default:
                                throw std::range_error(
                                    "ValidationTestSuite.status(): "
                                    "Unrecognized status type!" );
                        }
                    }
                    if (pending) {
                        return test_status_t::PENDING;
                    } else {
                        return test_status_t::PASSED;
                    }
                }

                virtual std::vector<std::string> status_messages(
                    bool failed_only = false ) {
                    std::vector<std::string> messages;
                    for (auto it = _tests.cbegin(); it != _tests.cend();
                         ++it) {
                        if (( !failed_only ) || ( *it )->failed()) {
                            messages.push_back( ( *it )->status_message() );
                        }
                    }
                    return messages;
                }

            private:
                // test_status_t _status;
                std::vector<std::unique_ptr<ValidationTest>> _tests;
        };
    }
}