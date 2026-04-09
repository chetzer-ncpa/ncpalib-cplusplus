#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/ValidationTest.hpp"
#include "NCPA/configuration/ValidationTestSuite.hpp"
#include "NCPA/units.hpp"

#include <string>
#include <type_traits>
#include <vector>

#ifndef NCPA_CONFIG_USE_STRICT_UNITS
#  define NCPA_CONFIG_USE_STRICT_UNITS false
#endif

namespace NCPA {
    namespace config {
        using namespace NCPA::units;

        class BaseParameter {
            public:
                BaseParameter( units_ptr_t u = nullptr ) {}

                BaseParameter( const ValidationTest& newtest ) :
                    BaseParameter() {
                    this->append_test( newtest );
                }

                BaseParameter( const test_ptr_t& newtest ) : BaseParameter() {
                    this->append_test( newtest );
                }

                BaseParameter( std::initializer_list<test_ptr_t> new_tests ) :
                    BaseParameter() {
                    this->append_tests( new_tests );
                }

                BaseParameter( const BaseParameter& other ) {
                    _tests = other._tests;
                }

                BaseParameter( BaseParameter&& other ) noexcept :
                    BaseParameter() {
                    ::swap( *this, other );
                }

                virtual ~BaseParameter() {}

                friend void ::swap( BaseParameter& a,
                                    BaseParameter& b ) noexcept;

                virtual BaseParameter& append_test(
                    const ValidationTest& newtest ) {
                    _tests.append( newtest );
                    return *this;
                }

                virtual BaseParameter& append_test(
                    const test_ptr_t& newtest ) {
                    _tests.append( newtest );
                    return *this;
                }

                virtual BaseParameter& append_tests(
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

                virtual BaseParameter& prepend_test(
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

                virtual BaseParameter& prepend_tests(
                    std::initializer_list<ValidationTest> new_tests ) {
                    for (auto it = new_tests.begin(); it != new_tests.end();
                         ++it) {
                        this->prepend_test( *it );
                    }
                    return *this;
                }

                virtual BaseParameter& validate( bool short_circuit = false ) {
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

                // virtual void add_units( units_ptr_t u ) = 0;

                virtual bool has_units() const { return false; }

                virtual void set_units( units_ptr_t u ) {
                    this->_check_has_units();
                }

                virtual units_ptr_t get_units() const {
                    this->_check_has_units();
                    return nullptr;
                }

                virtual BaseParameter& convert_units( units_ptr_t u ) {
                    this->_check_has_units();
                    return *this;
                }

                virtual BaseParameter& convert_units( const std::string& s ) {
                    return this->convert_units(
                        NCPA::units::Units::from_string( s ) );
                }

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

                // virtual VectorWithUnits<double> as_double_vector_with_units()
                //     const {
                //     _check_has_units();
                //     return VectorWithUnits( this->as_double_vector(),
                //                             this->get_units() );
                // }

                // virtual ScalarWithUnits<double> as_double_with_units(
                //     size_t n = 0 ) const {
                //     return ScalarWithUnits<double>( this->as_double( n ),
                //                                     this->get_units() );
                // }

                // virtual ScalarWithUnits<long long> as_int_with_units(
                //     size_t n = 0 ) const {
                //     return ScalarWithUnits<long long>( this->as_int( n ),
                //                                        this->get_units() );
                // }

                // virtual ScalarWithUnits<unsigned long long>
                //     as_unsigned_int_with_units( size_t n = 0 ) const {
                //     return ScalarWithUnits<unsigned long long>(
                //         this->as_unsigned_int( n ), this->get_units() );
                // }

                // abstract API
                virtual bool as_bool( size_t n = 0 ) const = 0;
                virtual std::complex<double> as_complex( size_t n = 0 ) const
                    = 0;
                virtual double as_double( size_t n = 0 ) const = 0;
                virtual long long as_int( size_t n = 0 ) const = 0;
                virtual std::string as_string( size_t n = 0 ) const     = 0;
                virtual unsigned long long as_unsigned_int( size_t n
                                                            = 0 ) const = 0;


                virtual param_ptr_t clone() const     = 0;
                virtual parameter_form_t form() const = 0;

                virtual void from_bool( bool b )                     = 0;
                virtual void from_bool( const std::vector<bool>& b ) = 0;

                virtual void from_complex( std::complex<double> c ) = 0;
                virtual void from_complex(
                    const std::vector<std::complex<double>>& c ) = 0;

                virtual void from_double( double d )                     = 0;
                // virtual void from_double( double d, units_ptr_t u )      =
                // 0;
                virtual void from_double( const std::vector<double>& d ) = 0;
                // virtual void from_double( const std::vector<double>& d,
                //                           units_ptr_t u )                =
                //                           0;


                virtual void from_int( long long i ) = 0;
                // virtual void from_int( long long i, units_ptr_t u ) = 0;

                virtual void from_int( const std::vector<long long>& i ) = 0;
                // virtual void from_int( const std::vector<long long>& i,
                //                        units_ptr_t u )                   =
                //                        0;


                virtual void from_string( const std::string& s ) = 0;
                virtual void from_string( const std::vector<std::string>& s )
                    = 0;

                virtual void from_unsigned_int( unsigned long long n ) = 0;
                // virtual void from_unsigned_int( unsigned long long n,
                //                                 units_ptr_t u )        = 0;
                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& n ) = 0;
                // virtual void from_unsigned_int(
                //     const std::vector<unsigned long long>& n, units_ptr_t u
                //     ) = 0;

                virtual size_t size() const           = 0;
                virtual parameter_type_t type() const = 0;

            protected:
                ValidationTestSuite _tests;

                virtual void _check_has_units() const {
                    if (NCPA_CONFIG_USE_STRICT_UNITS && !this->has_units()) {
                        throw std::logic_error(
                            "Parameter does not have units!" );
                    }
                }
        };

    }  // namespace config
}  // namespace NCPA

void swap( NCPA::config::BaseParameter& a,
           NCPA::config::BaseParameter& b ) noexcept {
    using std::swap;
    swap( a._tests, b._tests );
}
