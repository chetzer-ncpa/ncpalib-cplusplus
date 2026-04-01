#pragma once

#include "NCPA/configuration/declarations.hpp"

#include <sstream>

namespace NCPA {
    namespace config {

        class ValidationTest {
            public:
                ValidationTest() : _status { test_status_t::PENDING } {}

                ValidationTest( const ValidationTest& other ) :
                    ValidationTest() {
                    _status = other._status;
                }

                ValidationTest( ValidationTest&& other ) noexcept {
                    ::swap( *this, other );
                }

                virtual ~ValidationTest() {}

                friend void ::swap( ValidationTest& a,
                                    ValidationTest& b ) noexcept;

                virtual ValidationTest& clear() {
                    _status = test_status_t::PENDING;
                    return *this;
                }

                virtual ValidationTest& fail() {
                    _status = test_status_t::FAILED;
                    return *this;
                }

                virtual bool failed() const {
                    return ( _status == test_status_t::FAILED );
                }

                virtual ValidationTest& pass() {
                    _status = test_status_t::PASSED;
                    return *this;
                }

                virtual bool passed() const {
                    return ( _status == test_status_t::PASSED );
                }

                virtual bool pending() const {
                    return ( _status == test_status_t::PENDING );
                }

                virtual test_status_t status() const { return _status; }

                virtual std::string status_message() const {
                    std::ostringstream oss;
                    oss << this->description() << ": ";
                    switch (this->status()) {
                        case test_status_t::PENDING:
                            oss << "PENDING";
                            break;
                        case test_status_t::PASSED:
                            oss << "PASSED";
                            break;
                        case test_status_t::FAILED:
                            oss << "FAILED";
                    }
                    return oss.str();
                }

                virtual ValidationTest& test(
                    const std::unique_ptr<Parameter>& param ) {
                    return this->test( param.get() );
                }

                virtual std::string description() const = 0;
                virtual ValidationTest& test(
                    const Parameter *param )           = 0;
                virtual std::unique_ptr<ValidationTest> clone() const = 0;

            protected:
                test_status_t _status;
        };
    }  // namespace config
}  // namespace NCPA
