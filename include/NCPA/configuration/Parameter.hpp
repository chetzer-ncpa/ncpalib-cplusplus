#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/ValidationTest.hpp"
#include "NCPA/configuration/ValidationTestSuite.hpp"
#include "NCPA/units.hpp"

#include <string>
#include <type_traits>
#include <vector>

namespace NCPA {
    namespace config {
        class Parameter {
            public:
                Parameter() {}

                Parameter( const ValidationTest& newtest ) : Parameter() {
                    this->append_test( newtest );
                }

                Parameter( const test_ptr_t& newtest ) : Parameter() {
                    this->append_test( newtest );
                }

                Parameter( std::initializer_list<test_ptr_t> new_tests ) :
                    Parameter() {
                    this->append_tests( new_tests );
                }

                Parameter( const Parameter& other ) { 
                    _tests = other._tests;
                 }

                Parameter( Parameter&& other ) noexcept : Parameter() {
                    ::swap( *this, other );
                }

                virtual ~Parameter() {}

                friend void ::swap( Parameter& a, Parameter& b ) noexcept;

                virtual Parameter& append_test(
                    const ValidationTest& newtest ) {
                    _tests.append( newtest );
                    return *this;
                }

                virtual Parameter& append_test( const test_ptr_t& newtest ) {
                    _tests.append( newtest );
                    return *this;
                }

                virtual Parameter& append_tests(
                    std::initializer_list<test_ptr_t> new_tests ) {
                    for (auto it = new_tests.begin(); it != new_tests.end();
                         ++it) {
                        this->append_test( *it );
                    }
                    return *this;
                }

                virtual const ValidationTestSuite& tests() const {
                    return _tests;
                }

                virtual std::vector<const ValidationTest *> failed_tests()
                    const {
                    return _tests.failed_tests();
                }

                virtual Parameter& prepend_test(
                    const ValidationTest& newtest ) {
                    _tests.prepend( newtest );
                    return *this;
                }

                virtual std::ostream& validation_report(
                    std::ostream& os, bool newline = true,
                    const std::string& prepend = "" ) const {
                    if (this->passed()) {
                        os << prepend << "All tests passed";
                    } else {
                        auto f = this->failed_tests();
                        for (auto it = f.begin(); it != f.end(); ++it) {
                            if (it != f.begin()) {
                                os << std::endl;
                            }
                            os << prepend
                               << "Failed: " << ( *it )->description();
                        }
                    }
                    if (newline) {
                        os << std::endl;
                    }
                    return os;
                }

                virtual Parameter& prepend_tests(
                    std::initializer_list<ValidationTest> new_tests ) {
                    for (auto it = new_tests.begin(); it != new_tests.end();
                         ++it) {
                        this->prepend_test( *it );
                    }
                    return *this;
                }

                virtual Parameter& validate( bool short_circuit = false ) {
                    _tests.run_tests( this, short_circuit );
                    return *this;
                }

                virtual bool failed() const { return _tests.failed(); }

                virtual bool passed() const { return _tests.passed(); }

                virtual bool pending() const { return _tests.pending(); }

                virtual test_status_t validation_status() const {
                    return _tests.status();
                }

                virtual bool was_set() const { return true; }

                virtual parameter_form_t form() const                   = 0;
                virtual parameter_type_t type() const                   = 0;
                virtual size_t size() const                             = 0;
                // virtual bool was_set() const                        = 0;
                virtual param_ptr_t clone() const                       = 0;
                virtual bool as_bool( size_t n = 0 ) const              = 0;
                virtual std::string as_string( size_t n = 0 ) const     = 0;
                virtual long long as_int( size_t n = 0 ) const          = 0;
                virtual unsigned long long as_unsigned_int( size_t n
                                                            = 0 ) const = 0;
                virtual double as_double( size_t n = 0 ) const          = 0;
                virtual std::complex<double> as_complex( size_t n = 0 ) const
                    = 0;

                virtual void from_bool( bool b )                       = 0;
                virtual void from_string( const std::string& s )       = 0;
                virtual void from_int( long long i )                   = 0;
                virtual void from_unsigned_int( unsigned long long n ) = 0;
                virtual void from_double( double d )                   = 0;
                virtual void from_complex( std::complex<double> c )    = 0;
                virtual void from_bool( const std::vector<bool>& b )   = 0;
                virtual void from_string( const std::vector<std::string>& s )
                    = 0;
                virtual void from_int( const std::vector<long long>& i ) = 0;
                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& n )           = 0;
                virtual void from_double( const std::vector<double>& d ) = 0;
                virtual void from_complex(
                    const std::vector<std::complex<double>>& c ) = 0;

                bool is_scalar() const {
                    return ( this->form() == parameter_form_t::SCALAR );
                }

                bool is_vector() const {
                    return ( this->form() == parameter_form_t::VECTOR );
                }

                virtual std::vector<double> as_double_vector() const {
                    std::vector<double> v( this->size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        v.at( i ) = this->as_double( i );
                    }
                    return v;
                }

            protected:
                ValidationTestSuite _tests;

        };
    }  // namespace config
}  // namespace NCPA

void swap( NCPA::config::Parameter& a, NCPA::config::Parameter& b ) noexcept {
    using std::swap;
    swap( a._tests, b._tests );
}
