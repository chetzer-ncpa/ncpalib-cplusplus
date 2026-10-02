#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/ValidationTest.hpp"
#include "NCPA/configuration/ValidationTestSuite.hpp"
// #include "NCPA/json.hpp"
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

        class BaseParameter : public virtual CloneBase<BaseParameter> {
            public:
                BaseParameter() {}

                BaseParameter( const BaseParameter& other ) {}

                BaseParameter( BaseParameter&& other ) noexcept :
                    BaseParameter() {
                    swap( *this, other );
                }

                virtual ~BaseParameter() {}

                friend void swap( BaseParameter& a,
                                  BaseParameter& b ) noexcept {
                    using std::swap;
                }

                virtual bool as_bool() const { return this->as_bool( 0 ); }

                virtual std::complex<double> as_complex() const {
                    return this->as_complex( 0 );
                }

                virtual double as_double() const {
                    return this->as_double( 0 );
                }

                virtual std::vector<double> as_double_vector() const {
                    std::vector<double> v( this->size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        v.at( i ) = this->as_double( i );
                    }
                    return v;
                }

                virtual long long as_int() const { return this->as_int( 0 ); }

                virtual std::string as_string(bool throw_on_error = false) const {
                    return this->as_string( 0, throw_on_error );
                }

                virtual unsigned long long as_unsigned_int() const {
                    return this->as_unsigned_int( 0 );
                }

                virtual BaseParameter& convert_units( units_ptr_t u ) {
                    this->_check_has_units();
                    return *this;
                }

                virtual BaseParameter& convert_units( const std::string& s ) {
                    return this->convert_units(
                        NCPA::units::Units::from_string( s ) );
                }

                virtual double diff() const {
                    return ( this->size() > 1
                                 ? this->as_double( 1 ) - this->as_double( 0 )
                                 : 0.0 );
                }

                virtual bool failed() const {
                    validation_status_t res = this->validate();
                    return ( res.result == test_result_t::FAILED );
                }

                virtual units_ptr_t get_units() const {
                    this->_check_has_units();
                    return nullptr;
                }

                virtual bool has_units() const { return false; }

                virtual bool is_scalar() const {
                    return ( this->form() == parameter_form_t::SCALAR );
                }

                virtual bool is_vector() const {
                    return ( this->form() == parameter_form_t::VECTOR );
                }

                virtual bool passed() const {
                    validation_status_t res = this->validate();
                    return ( res.result == test_result_t::NONE
                             || res.result == test_result_t::PASSED );
                }

                virtual void set_units( units_ptr_t u ) {
                    this->_check_has_units();
                }

                virtual bool was_set() const { return true; }

// #if HAVE_NCPA_JSON
//                 virtual NCPA::json::JSONItem as_json( const std::string& key, const std::string& comment = "" ) const {
//                     if (this->is_scalar()) {
//                     switch (this->json_type()) {
//                         case NCPA::json::json_type_t::BOOLEAN:
//                             return NCPA::json::JSONItem( key, this->as_bool(), this->json_form(), this->json_type(), comment );
//                     }
//                 }
//                 }

//                 virtual NCPA::json::json_form_t json_form() const {
//                     switch (this->form()) {
//                         case parameter_form_t::SCALAR:
//                             return NCPA::json::json_form_t::SCALAR;
//                         case parameter_form_t::VECTOR:
//                             return NCPA::json::json_form_t::VECTOR;
//                         default:
//                             throw std::out_of_range(
//                                 "Parameter form not set" );
//                     }
//                 }

//                 virtual NCPA::json::json_type_t json_type() const {
//                     switch (this->type()) {
//                         case parameter_type_t::BOOLEAN:
//                             return NCPA::json::json_type_t::BOOLEAN;
//                         case parameter_type_t::FLOAT:
//                             return NCPA::json::json_type_t::FLOAT;
//                         case parameter_type_t::INTEGER:
//                             return NCPA::json::json_type_t::INTEGER;
//                         case parameter_type_t::STRING:
//                             return NCPA::json::json_type_t::STRING;
//                         default:
//                             // this will throw if it can't be done
//                             std::string s = this->as_string( true );
//                             return NCPA::json::json_type_t::STRING;
//                     }
//                 }
// #endif

                // abstract API
                virtual bool as_bool( size_t n ) const                    = 0;
                virtual std::complex<double> as_complex( size_t n ) const = 0;
                virtual double as_double( size_t n ) const                = 0;
                virtual long long as_int( size_t n ) const                = 0;
                virtual std::string as_string( size_t n, bool throw_on_error
                                                         = false ) const  = 0;
                virtual unsigned long long as_unsigned_int( size_t n ) const
                    = 0;


                virtual BaseParameter& add_test( Validation& v )      = 0;
                virtual BaseParameter& add_test( Validation&& v )     = 0;
                virtual validation_status_t validate( bool short_circuit
                                                      = false ) const = 0;

                // virtual param_ptr_t clone() const     = 0;
                virtual parameter_form_t form() const = 0;

                virtual void from_bool( bool b )                     = 0;
                virtual void from_bool( const std::vector<bool>& b ) = 0;

                virtual void from_complex( std::complex<double> c ) = 0;
                virtual void from_complex(
                    const std::vector<std::complex<double>>& c ) = 0;

                virtual void from_double( double d )                     = 0;
                virtual void from_double( const std::vector<double>& d ) = 0;

                virtual void from_int( long long i ) = 0;

                virtual void from_int( const std::vector<long long>& i ) = 0;


                virtual void from_string( const std::string& s ) = 0;
                virtual void from_string( const std::vector<std::string>& s )
                    = 0;

                virtual void from_unsigned_int( unsigned long long n ) = 0;
                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& n ) = 0;

                virtual size_t size() const           = 0;
                virtual parameter_type_t type() const = 0;

            protected:
                virtual void _check_has_units() const {
                    if (NCPA_CONFIG_USE_STRICT_UNITS && !this->has_units()) {
                        throw std::logic_error(
                            "Parameter does not have units!" );
                    }
                }
        };

    }  // namespace config
}  // namespace NCPA
