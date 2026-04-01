#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/validation.hpp"
#include "NCPA/units.hpp"

#include <string>
#include <type_traits>
#include <vector>

namespace NCPA {
    namespace config {
        class Parameter {
            public:
                Parameter() {}

                Parameter( const Parameter& other ) { _tests = other._tests; }

                Parameter( Parameter&& other ) noexcept : Parameter() {
                    ::swap( *this, other );
                }

                Parameter& operator=( Parameter other ) {
                    ::swap( *this, other );
                    return *this;
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

                // template specializations, scalar form
                // floating point
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                PARAMTYPE as() const {
                    return static_cast<PARAMTYPE>( this->as_double() );
                }

                // signed integer
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_signed<PARAMTYPE>::value ),
                             int>::type = 0>
                PARAMTYPE as() const {
                    return static_cast<PARAMTYPE>( this->as_int() );
                }

                // unsigned integer
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_unsigned<PARAMTYPE>::value ),
                             int>::type = 0>
                PARAMTYPE as() const {
                    return static_cast<PARAMTYPE>( this->as_size_t() );
                }

                // boolean
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_same<PARAMTYPE, bool>::value, int>::type
                         = 0>
                PARAMTYPE as() const {
                    return this->as_bool();
                }

                // string
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<PARAMTYPE>::value )
                               && std::is_convertible<PARAMTYPE,
                                                      std::string>::value ),
                             int>::type = 0>
                PARAMTYPE as() const {
                    std::string s = this->as_string();
                    return s;
                }

                // complex
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_scalar<PARAMTYPE>::value )
                               && NCPA::types::is_complex<PARAMTYPE>::value ),
                             int>::type = 0>
                PARAMTYPE as() const {
                    return PARAMTYPE( this->as_complex().real(),
                                      this->as_complex().imag() );
                }

                // everything else
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<PARAMTYPE>::value
                                  || std::is_convertible<PARAMTYPE,
                                                         std::string>::value
                                  || ( !( std::is_scalar<PARAMTYPE>::value )
                                       && NCPA::types::is_complex<
                                           PARAMTYPE>::value ) ) ),
                             int>::type = 0>
                PARAMTYPE as() const {
                    return this->typed<PARAMTYPE>().get();
                }

                // template specializations, scalar form
                // floating point
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                std::vector<PARAMTYPE> as_vector() const {
                    std::vector<PARAMTYPE> newvec( this->size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        newvec[ i ]
                            = static_cast<PARAMTYPE>( this->as_double( i ) );
                    }
                    return newvec;
                }

                // signed integer
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_signed<PARAMTYPE>::value ),
                             int>::type = 0>
                std::vector<PARAMTYPE> as_vector() const {
                    std::vector<PARAMTYPE> newvec( this->size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        newvec[ i ]
                            = static_cast<PARAMTYPE>( this->as_int( i ) );
                    }
                    return newvec;
                }

                // unsigned integer
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_unsigned<PARAMTYPE>::value ),
                             int>::type = 0>
                std::vector<PARAMTYPE> as_vector() const {
                    std::vector<PARAMTYPE> newvec( this->size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        newvec[ i ]
                            = static_cast<PARAMTYPE>( this->as_size_t( i ) );
                    }
                    return newvec;
                }

                // boolean
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_same<PARAMTYPE, bool>::value, int>::type
                         = 0>
                std::vector<PARAMTYPE> as_vector() const {
                    std::vector<PARAMTYPE> newvec( this->size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        newvec[ i ] = this->as_bool( i );
                    }
                    return newvec;
                }

                // string
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<PARAMTYPE>::value )
                               && std::is_convertible<PARAMTYPE,
                                                      std::string>::value ),
                             int>::type = 0>
                std::vector<PARAMTYPE> as_vector() const {
                    std::vector<PARAMTYPE> newvec( this->size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        newvec[ i ] = this->as_string( i );
                    }
                    return newvec;
                }

                // complex
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_scalar<PARAMTYPE>::value )
                               && NCPA::types::is_complex<PARAMTYPE>::value ),
                             int>::type = 0>
                std::vector<PARAMTYPE> as_vector() const {
                    std::vector<PARAMTYPE> newvec( this->size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        newvec[ i ] = this->as_complex( i );
                    }
                    return newvec;
                }

                virtual size_t as_size_t( size_t n = 0 ) const {
                    return static_cast<size_t>( this->as_long_int( n ) );
                }

                // template specializations, scalar form
                // floating point
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_scalar<PARAMTYPE>::value
                               && std::is_floating_point<PARAMTYPE>::value ),
                             int>::type = 0>
                void from( PARAMTYPE input ) const {
                    this->from_double( static_cast<double>( input ) );
                }

                // signed integer
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_scalar<PARAMTYPE>::value
                               && std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_signed<PARAMTYPE>::value ),
                             int>::type = 0>
                void from( PARAMTYPE input ) const {
                    this->from_int( static_cast<long long>( input ) );
                }

                // unsigned integer
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_scalar<PARAMTYPE>::value
                               && std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_unsigned<PARAMTYPE>::value ),
                             int>::type = 0>
                void from( PARAMTYPE input ) const {
                    this->from_size_t( static_cast<size_t>( input ) );
                }

                // boolean
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_same<PARAMTYPE, bool>::value, int>::type
                         = 0>
                void from( PARAMTYPE input ) const {
                    this->from_bool( input );
                }

                // string
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<PARAMTYPE>::value )
                               && std::is_convertible<PARAMTYPE,
                                                      std::string>::value ),
                             int>::type = 0>
                void from( PARAMTYPE input ) const {
                    this->from_string( input );
                }

                // complex
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_scalar<PARAMTYPE>::value )
                               && NCPA::types::is_complex<PARAMTYPE>::value ),
                             int>::type = 0>
                void from( PARAMTYPE input ) const {
                    this->from_complex( std::complex<double>(
                        static_cast<double>( input.real() ),
                        static_cast<double>( input.imag() ) ) );
                }

                // everything else scalar
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( ( std::is_scalar<PARAMTYPE>::value
                                 && !( std::is_arithmetic<PARAMTYPE>::value
                                       || std::is_convertible<
                                           PARAMTYPE, std::string>::value ) )
                               || ( !( std::is_scalar<PARAMTYPE>::value )
                                    && !( NCPA::types::is_complex<
                                          PARAMTYPE>::value ) ) ),
                             int>::type = 0>
                void from( PARAMTYPE input ) const {
                    this->scalar<PARAMTYPE>().set( input );
                }

                // template specializations, vector form
                // floating point
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_scalar<PARAMTYPE>::value
                               && std::is_floating_point<PARAMTYPE>::value ),
                             int>::type = 0>
                void from_vector( std::vector<PARAMTYPE> input ) const {
                    this->from_double( NCPA::arrays::cast_vector<PARAMTYPE,double>( input ) );
                }

                // signed integer
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_scalar<PARAMTYPE>::value
                               && std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_signed<PARAMTYPE>::value ),
                             int>::type = 0>
                void from_vector( std::vector<PARAMTYPE> input ) const {
                    this->from_int( NCPA::arrays::cast_vector<PARAMTYPE,long long>( input ) );
                }

                // unsigned integer
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_scalar<PARAMTYPE>::value
                               && std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_unsigned<PARAMTYPE>::value ),
                             int>::type = 0>
                void from_vector( std::vector<PARAMTYPE> input ) const {
                    this->from_size_t( NCPA::arrays::cast_vector<PARAMTYPE,size_t>( input ) );
                }

                // boolean
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_same<PARAMTYPE, bool>::value, int>::type
                         = 0>
                void from_vector( std::vector<PARAMTYPE> input ) const {
                    this->from_bool( input);
                }

                // string
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<PARAMTYPE>::value )
                               && std::is_convertible<PARAMTYPE,
                                                      std::string>::value ),
                             int>::type = 0>
                void from_vector( std::vector<PARAMTYPE> input ) const {
                    std::vector<std::string> newvec(input.size());
                    for (size_t i = 0; i < input.size(); ++i) {
                        newvec[i] = input.at(i);
                    }
                    this->from_string( newvec );
                }

                // complex
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_scalar<PARAMTYPE>::value )
                               && NCPA::types::is_complex<PARAMTYPE>::value ),
                             int>::type = 0>
                void from_vector( std::vector<PARAMTYPE> input ) const {
                    std::vector<std::complex<double>> newvec(input.size());
                    for (size_t i = 0; i < input.size(); ++i) {
                        newvec[i].real( static_cast<double>( input.at(i).real() ) );
                        newvec[i].imag( static_cast<double>( input.at(i).imag() ) );
                    }
                    this->from_complex( newvec );
                }

                // everything else vector
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( ( std::is_scalar<PARAMTYPE>::value
                                 && !( std::is_arithmetic<PARAMTYPE>::value
                                       || std::is_convertible<
                                           PARAMTYPE, std::string>::value ) )
                               || ( !( std::is_scalar<PARAMTYPE>::value )
                                    && !( NCPA::types::is_complex<
                                          PARAMTYPE>::value ) ) ),
                             int>::type = 0>
                void from_vector( const std::vector<PARAMTYPE>& input ) const {
                    this->vector<PARAMTYPE>().set( input );
                }


                virtual parameter_form_t form() const               = 0;
                virtual parameter_type_t type() const               = 0;
                virtual size_t size() const                         = 0;
                // virtual bool was_set() const                        = 0;
                virtual param_ptr_t clone() const                   = 0;
                virtual bool as_bool( size_t n = 0 ) const          = 0;
                virtual std::string as_string( size_t n = 0 ) const = 0;
                virtual int as_int( size_t n = 0 ) const            = 0;
                virtual long long as_long_int( size_t n = 0 ) const = 0;
                virtual double as_double( size_t n = 0 ) const      = 0;
                virtual std::complex<double> as_complex( size_t n = 0 ) const
                    = 0;

                virtual void from_bool( bool b )                     = 0;
                virtual void from_string( const std::string& s )     = 0;
                virtual void from_int( long long i )                 = 0;
                virtual void from_unsigned_int( size_t n )           = 0;
                virtual void from_double( double d )                 = 0;
                virtual void from_complex( std::complex<double> c )  = 0;
                virtual void from_bool( const std::vector<bool>& b ) = 0;
                virtual void from_string( const std::vector<std::string>& s )
                    = 0;
                virtual void from_int( const std::vector<long long>& i ) = 0;
                virtual void from_unsigned_int( const std::vector<size_t>& n )
                    = 0;
                virtual void from_double( const std::vector<double>& d ) = 0;
                virtual void from_complex(
                    const std::vector<std::complex<double>>& c ) = 0;

                template<typename T>
                TypedParameter<T>& typed() {
                    return *dynamic_cast<TypedParameter<T> *>( this );
                }

                template<typename T>
                ScalarParameter<T>& scalar() {
                    return *dynamic_cast<ScalarParameter<T> *>( this );
                }

                template<typename T>
                VectorParameter<T>& vector() {
                    return *dynamic_cast<VectorParameter<T> *>( this );
                }

                bool is_scalar() const {
                    return ( this->form() == parameter_form_t::SCALAR );
                }

                bool is_vector() const {
                    return ( this->form() == parameter_form_t::VECTOR );
                }


            protected:
                ValidationTestSuite _tests;

                // const parameter_form_t _form;
                // const parameter_type_t _type;
        };
    }  // namespace config
}  // namespace NCPA

void swap( NCPA::config::Parameter& a, NCPA::config::Parameter& b ) noexcept {
    using std::swap;
    swap( a._tests, b._tests );
}
