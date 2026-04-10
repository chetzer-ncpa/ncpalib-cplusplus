#pragma once

#include "NCPA/configuration/BaseParameter.hpp"
#include "NCPA/configuration/boilerplate.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/TypedParameter.hpp"

#include <regex>
#include <vector>

namespace NCPA {
    namespace config {
        namespace hidden {
            template<typename T>
            class _base_scalar_parameter;
        }
    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap( NCPA::config::hidden::_base_scalar_parameter<T>& a,
           NCPA::config::hidden::_base_scalar_parameter<T>& b ) noexcept;

namespace NCPA {
    namespace config {

        namespace hidden {
            template<typename PARAMTYPE>
            class _base_scalar_parameter : public TypedParameter<PARAMTYPE> {
                public:
                    _base_scalar_parameter() :
                        _base_scalar_parameter<PARAMTYPE>() {}

                    _base_scalar_parameter( PARAMTYPE defaultval ) :
                        TypedParameter<PARAMTYPE>(), _value { defaultval } {}

                    _base_scalar_parameter(
                        const std::vector<PARAMTYPE>& defaultval ) :
                        TypedParameter<PARAMTYPE>() {
                        _get_first_value_from_vector( defaultval );
                    }

                    _base_scalar_parameter( const ValidationTest& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ) {}

                    _base_scalar_parameter( const test_ptr_t& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ) {}

                    _base_scalar_parameter(
                        std::initializer_list<test_ptr_t> new_tests ) :
                        TypedParameter<PARAMTYPE>( new_tests ) {}

                    _base_scalar_parameter( PARAMTYPE defaultval,
                                            const ValidationTest& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ),
                        _value { defaultval } {}

                    _base_scalar_parameter( PARAMTYPE defaultval,
                                            const test_ptr_t& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ),
                        _value { defaultval } {}

                    _base_scalar_parameter(
                        PARAMTYPE defaultval,
                        std::initializer_list<test_ptr_t> new_tests ) :
                        TypedParameter<PARAMTYPE>( new_tests ),
                        _value { defaultval } {}

                    _base_scalar_parameter(
                        const std::vector<PARAMTYPE>& defaultval,
                        const ValidationTest& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ) {
                        _get_first_value_from_vector( defaultval );
                    }

                    _base_scalar_parameter(
                        const std::vector<PARAMTYPE>& defaultval,
                        const test_ptr_t& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ) {
                        _get_first_value_from_vector( defaultval );
                    }

                    _base_scalar_parameter(
                        const std::vector<PARAMTYPE>& defaultval,
                        std::initializer_list<test_ptr_t> new_tests ) :
                        TypedParameter<PARAMTYPE>( new_tests ) {
                        _get_first_value_from_vector( defaultval );
                    }

                    virtual ~_base_scalar_parameter() {}

                    _base_scalar_parameter(
                        const _base_scalar_parameter<PARAMTYPE>& other ) :
                        TypedParameter<PARAMTYPE>( other ) {
                        _value = other._value;
                    }

                    _base_scalar_parameter(
                        _base_scalar_parameter<PARAMTYPE>&& other ) noexcept :
                        TypedParameter<PARAMTYPE>() {
                        ::swap( *this, other );
                    }

                    _base_scalar_parameter<PARAMTYPE>& operator=(
                        _base_scalar_parameter<PARAMTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    friend void ::swap<>(
                        _base_scalar_parameter<PARAMTYPE>& a,
                        _base_scalar_parameter<PARAMTYPE>& b ) noexcept;

                    virtual parameter_form_t form() const override {
                        return parameter_form_t::SCALAR;
                    }

                    virtual size_t size() const override { return 1; }

                    virtual PARAMTYPE get( size_t n = 0 ) const override {
                        return this->_value;
                    }

                    virtual std::vector<PARAMTYPE> get_vector()
                        const override {
                        return std::vector<PARAMTYPE> { this->_value };
                    }

                protected:
                    PARAMTYPE _value;

                    void _get_first_value_from_vector(
                        const std::vector<PARAMTYPE>& defaultval ) {
                        if (defaultval.size() > 0) {
                            _value = defaultval.at( 0 );
                        }
                    }

                    template<typename T = PARAMTYPE>
                    typename std::enable_if<
                        NCPA::types::has_to_string<T>::value,
                        std::string>::type
                        _as_string( size_t n = 0 ) const {
                        return this->get( n ).to_string();
                    }

                    template<typename T = PARAMTYPE>
                    typename std::enable_if<
                        ( !( NCPA::types::has_to_string<T>::value )
                          && NCPA::types::can_use_std_to_string<T>::value ),
                        std::string>::type
                        _as_string( size_t n = 0 ) const {
                        return std::to_string( this->get( n ) );
                    }

                    template<typename T = PARAMTYPE>
                    typename std::enable_if<
                        ( !( NCPA::types::can_use_std_to_string<T>::value
                             || NCPA::types::has_to_string<T>::value )
                          && NCPA::types::can_use_to_string<T>::value ),
                        std::string>::type
                        _as_string( size_t n = 0 ) const {
                        return to_string( this->get( n ) );
                    }

                    template<typename T = PARAMTYPE>
                    typename std::enable_if<
                        !( NCPA::types::can_use_to_string<T>::value
                           || NCPA::types::has_to_string<T>::value
                           || NCPA::types::can_use_std_to_string<T>::value ),
                        std::string>::type
                        _as_string( size_t n = 0 ) const {
                            return "<no string conversion defined>";
                        // throw std::out_of_range(
                        //     "No as_string() conversion defined!" );
                    }
            };
        }  // namespace hidden

        // floating point
        template<typename PARAMTYPE>
        class ScalarParameter<
            PARAMTYPE, typename std::enable_if<
                           std::is_floating_point<PARAMTYPE>::value>::type>
            : public hidden::_base_scalar_parameter<PARAMTYPE> {
            public:
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_bool;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_complex;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_double;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_int;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_string;
                using hidden::_base_scalar_parameter<
                    PARAMTYPE>::as_unsigned_int;

                ScalarParameter() :
                    hidden::_base_scalar_parameter<PARAMTYPE>( 0.0 ) {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter( const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter(
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( new_tests ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const ScalarParameter<PARAMTYPE>& other ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( other ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    hidden::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                virtual ~ScalarParameter() {}

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ScalarParameter<PARAMTYPE>& a,
                                      ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual long long as_int( size_t n ) const override {
                    return static_cast<long long>(
                        this->get( n )
                        + std::numeric_limits<PARAMTYPE>::epsilon() );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n ) const override {
                    return static_cast<unsigned long long>(
                        this->get( n )
                        + std::numeric_limits<PARAMTYPE>::epsilon() );
                }

                virtual bool as_bool( size_t n ) const override {
                    return ( this->get( n ) != 0 );
                }

                virtual std::string as_string( size_t n ) const override {
                    return std::to_string( this->get( n ) );
                }

                virtual double as_double( size_t n ) const override {
                    return static_cast<double>( this->get( n ) );
                }

                virtual std::complex<double> as_complex(
                    size_t n ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_string( const std::string& x ) override {
                    this->_value = static_cast<PARAMTYPE>( std::stod( x ) );
                }

                virtual void from_bool( bool x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_int( long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_double( double x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value = static_cast<PARAMTYPE>( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_bool( vec.at( 0 ) );
                    }
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_int( vec.at( 0 ) );
                    }
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_unsigned_int( vec.at( 0 ) );
                    }
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_double( vec.at( 0 ) );
                    }
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_complex( vec.at( 0 ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    if (vec.size() > 0) {
                        std::string s = vec.at( 0 );
                        this->from_string( s );
                    }
                }
        };

        // signed integer
        template<typename PARAMTYPE>
        class ScalarParameter<PARAMTYPE,
                              typename std::enable_if<(
                                  std::is_integral<PARAMTYPE>::value
                                  && !( std::is_same<PARAMTYPE, bool>::value )
                                  && std::is_signed<PARAMTYPE>::value )>::type>
            : public hidden::_base_scalar_parameter<PARAMTYPE> {
            public:
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_bool;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_complex;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_double;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_int;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_string;
                using hidden::_base_scalar_parameter<
                    PARAMTYPE>::as_unsigned_int;

                ScalarParameter() :
                    hidden::_base_scalar_parameter<PARAMTYPE>( 0 ) {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter( const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter(
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( new_tests ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const ScalarParameter<PARAMTYPE>& other ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( other ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    hidden::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                virtual ~ScalarParameter() {}

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ScalarParameter<PARAMTYPE>& a,
                                      ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n ) const override {
                    return ( this->get( n ) != 0 );
                }

                virtual std::string as_string( size_t n ) const override {
                    return std::to_string( this->get( n ) );
                }

                virtual long long as_int( size_t n ) const override {
                    return static_cast<long long>( this->get( n ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n ) const override {
                    return static_cast<unsigned long long>( this->get( n ) );
                }

                virtual double as_double( size_t n ) const override {
                    return static_cast<double>( this->get( n ) );
                }

                virtual std::complex<double> as_complex(
                    size_t n ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_string( const std::string& x ) override {
                    this->_value = static_cast<PARAMTYPE>( std::stol( x ) );
                }

                virtual void from_bool( bool x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_int( long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_double( double x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value = static_cast<PARAMTYPE>( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_bool( vec.at( 0 ) );
                    }
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_int( vec.at( 0 ) );
                    }
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_unsigned_int( vec.at( 0 ) );
                    }
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_double( vec.at( 0 ) );
                    }
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_complex( vec.at( 0 ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    if (vec.size() > 0) {
                        std::string s = vec.at( 0 );
                        this->from_string( s );
                    }
                }
        };

        // unsigned integer
        template<typename PARAMTYPE>
        class ScalarParameter<
            PARAMTYPE, typename std::enable_if<(
                           std::is_integral<PARAMTYPE>::value
                           && !( std::is_same<PARAMTYPE, bool>::value )
                           && std::is_unsigned<PARAMTYPE>::value )>::type>
            : public hidden::_base_scalar_parameter<PARAMTYPE> {
            public:
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_bool;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_complex;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_double;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_int;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_string;
                using hidden::_base_scalar_parameter<
                    PARAMTYPE>::as_unsigned_int;

                ScalarParameter() :
                    hidden::_base_scalar_parameter<PARAMTYPE>( 0 ) {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter( const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter(
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( new_tests ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const ScalarParameter<PARAMTYPE>& other ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( other ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    hidden::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                virtual ~ScalarParameter() {}

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ScalarParameter<PARAMTYPE>& a,
                                      ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n ) const override {
                    return ( this->get( n ) != 0 );
                }

                virtual std::string as_string( size_t n ) const override {
                    return std::to_string( this->get( n ) );
                }

                virtual long long as_int( size_t n ) const override {
                    return static_cast<long long>( this->get( n ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n ) const override {
                    return static_cast<unsigned long long>( this->get( n ) );
                }

                virtual double as_double( size_t n ) const override {
                    return static_cast<double>( this->get( n ) );
                }

                virtual std::complex<double> as_complex(
                    size_t n ) const override {
                    return std::complex<double>( this->as_double(), 0.0 );
                }

                virtual void from_string( const std::string& x ) override {
                    this->_value = static_cast<PARAMTYPE>( std::stoull( x ) );
                }

                virtual void from_bool( bool x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_int( long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_double( double x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value = static_cast<PARAMTYPE>( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_bool( vec.at( 0 ) );
                    }
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_int( vec.at( 0 ) );
                    }
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_unsigned_int( vec.at( 0 ) );
                    }
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_double( vec.at( 0 ) );
                    }
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_complex( vec.at( 0 ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    if (vec.size() > 0) {
                        std::string s = vec.at( 0 );
                        this->from_string( s );
                    }
                }
        };

        // // boolean
        template<typename PARAMTYPE>
        class ScalarParameter<PARAMTYPE, typename std::enable_if<std::is_same<
                                             PARAMTYPE, bool>::value>::type>
            : public hidden::_base_scalar_parameter<PARAMTYPE> {
            public:
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_bool;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_complex;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_double;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_int;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_string;
                using hidden::_base_scalar_parameter<
                    PARAMTYPE>::as_unsigned_int;

                ScalarParameter() :
                    hidden::_base_scalar_parameter<PARAMTYPE>( false ) {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter( const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter(
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( new_tests ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const ScalarParameter<PARAMTYPE>& other ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( other ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    hidden::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                virtual ~ScalarParameter() {}

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ScalarParameter<PARAMTYPE>& a,
                                      ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n ) const override {
                    return this->get( n );
                }

                virtual std::string as_string( size_t n ) const override {
                    return ( this->get( n ) ? "true" : "false" );
                }

                virtual long long as_int( size_t n ) const override {
                    return static_cast<long long>( this->get( n ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n ) const override {
                    return static_cast<unsigned long long>( this->get( n ) );
                }

                virtual double as_double( size_t n ) const override {
                    return static_cast<double>( this->get( n ) );
                }

                virtual std::complex<double> as_complex(
                    size_t n ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_bool( bool x ) override { this->_value = x; }

                virtual void from_string( const std::string& x ) override {
                    this->_value = ( x == "true" );
                }

                virtual void from_int( long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_double( double x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value = static_cast<PARAMTYPE>( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_bool( vec.at( 0 ) );
                    }
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_int( vec.at( 0 ) );
                    }
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_unsigned_int( vec.at( 0 ) );
                    }
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_double( vec.at( 0 ) );
                    }
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_complex( vec.at( 0 ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    if (vec.size() > 0) {
                        std::string s = vec.at( 0 );
                        this->from_string( s );
                    }
                }

                //    private:
                //     PARAMTYPE _value;
        };

        // string
        template<typename PARAMTYPE>
        class ScalarParameter<
            PARAMTYPE,
            typename std::enable_if<(
                !( std::is_arithmetic<PARAMTYPE>::value )
                && std::is_convertible<PARAMTYPE, std::string>::value )>::type>
            : public hidden::_base_scalar_parameter<PARAMTYPE> {
            public:
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_bool;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_complex;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_double;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_int;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_string;
                using hidden::_base_scalar_parameter<
                    PARAMTYPE>::as_unsigned_int;

                ScalarParameter() :
                    hidden::_base_scalar_parameter<PARAMTYPE>( "" ) {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter( const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter(
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( new_tests ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const ScalarParameter<PARAMTYPE>& other ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( other ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    hidden::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                virtual ~ScalarParameter() {}

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ScalarParameter<PARAMTYPE>& a,
                                      ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n ) const override {
                    std::string s = this->get( n );
                    std::transform( s.begin(), s.end(), s.begin(),
                                    []( unsigned char c ) {
                                        return std::tolower( c );
                                    }  // correct
                    );
                    return ( s == "true" );
                }

                virtual std::string as_string( size_t n ) const override {
                    std::string s = this->get( n );
                    return s;
                }

                virtual long long as_int( size_t n ) const override {
                    return std::stoll( this->get( n ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n ) const override {
                    return std::stoull( this->get( n ) );
                }

                virtual double as_double( size_t n ) const override {
                    return std::stod( this->get( n ) );
                }

                virtual std::complex<double> as_complex(
                    size_t n ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_bool( bool x ) override {
                    this->_value = ( x ? "true" : "false" );
                }

                virtual void from_string( const std::string& x ) override {
                    this->_value = x;
                }

                virtual void from_int( long long x ) override {
                    this->_value = std::to_string( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value = std::to_string( x );
                }

                virtual void from_double( double x ) override {
                    this->_value = std::to_string( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value = std::to_string( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_bool( vec.at( 0 ) );
                    }
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_int( vec.at( 0 ) );
                    }
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_unsigned_int( vec.at( 0 ) );
                    }
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_double( vec.at( 0 ) );
                    }
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_complex( vec.at( 0 ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    if (vec.size() > 0) {
                        std::string s = vec.at( 0 );
                        this->from_string( s );
                    }
                }
        };

        // complex
        template<typename PARAMTYPE>
        class ScalarParameter<
            PARAMTYPE,
            typename std::enable_if<(
                !( std::is_scalar<PARAMTYPE>::value )
                && NCPA::types::is_complex<PARAMTYPE>::value )>::type>
            : public hidden::_base_scalar_parameter<PARAMTYPE> {
            public:
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_bool;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_complex;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_double;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_int;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_string;
                using hidden::_base_scalar_parameter<
                    PARAMTYPE>::as_unsigned_int;

                ScalarParameter() :
                    hidden::_base_scalar_parameter<PARAMTYPE>(
                        PARAMTYPE { 0, 0 } ) {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter( const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter(
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( new_tests ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const ScalarParameter<PARAMTYPE>& other ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( other ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    hidden::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                virtual ~ScalarParameter() {}

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ScalarParameter<PARAMTYPE>& a,
                                      ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n ) const override {
                    return ( std::abs( this->get( n ) ) != 0.0 );
                }

                virtual std::string as_string( size_t n ) const override {
                    std::ostringstream oss;
                    oss << std::to_string( this->get( n ).real() ) << " "
                        << ( this->get( n ).imag() < 0.0 ? "- " : "+ " )
                        << std::to_string( this->get( n ).imag() );
                    return oss.str();
                }

                virtual long long as_int( size_t n ) const override {
                    return static_cast<long long>(
                        std::abs( this->get( n ) ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n ) const override {
                    return static_cast<unsigned long long>(
                        std::abs( this->get( n ) ) );
                }

                virtual double as_double( size_t n ) const override {
                    return static_cast<double>( std::abs( this->get( n ) ) );
                }

                virtual std::complex<double> as_complex(
                    size_t n ) const override {
                    return std::complex<double>( this->get( n ).real(),
                                                 this->get( n ).imag() );
                }

                virtual void from_string( const std::string& x ) override {
                    throw std::out_of_range(
                        "from_string() not yet implemented for complex "
                        "numbers" );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->get().real( x.real() );
                    this->get().imag( x.imag() );
                }

                virtual void from_bool( bool x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_int( long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_double( double x ) override {
                    this->_value = static_cast<PARAMTYPE>( x );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_bool( vec.at( 0 ) );
                    }
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_int( vec.at( 0 ) );
                    }
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_unsigned_int( vec.at( 0 ) );
                    }
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_double( vec.at( 0 ) );
                    }
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    if (vec.size() > 0) {
                        this->from_complex( vec.at( 0 ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    if (vec.size() > 0) {
                        std::string s = vec.at( 0 );
                        this->from_string( s );
                    }
                }
        };

        // everything else
        template<typename PARAMTYPE>
        class ScalarParameter<
            PARAMTYPE,
            typename std::enable_if<( !(
                std::is_arithmetic<PARAMTYPE>::value
                || std::is_convertible<PARAMTYPE, std::string>::value
                || ( !( std::is_scalar<PARAMTYPE>::value )
                     && NCPA::types::is_complex<PARAMTYPE>::value ) ) )>::type>
            : public hidden::_base_scalar_parameter<PARAMTYPE> {
            public:
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_bool;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_complex;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_double;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_int;
                using hidden::_base_scalar_parameter<PARAMTYPE>::as_string;
                using hidden::_base_scalar_parameter<
                    PARAMTYPE>::as_unsigned_int;

                ScalarParameter() :
                    hidden::_base_scalar_parameter<PARAMTYPE>() {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval ) {}

                ScalarParameter( const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter( const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( newtest ) {}

                ScalarParameter(
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( new_tests ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( PARAMTYPE defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                ScalarParameter(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                ScalarParameter( const ScalarParameter<PARAMTYPE>& other ) :
                    hidden::_base_scalar_parameter<PARAMTYPE>( other ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    hidden::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                virtual ~ScalarParameter() {}

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ScalarParameter<PARAMTYPE>& a,
                                      ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n ) const override {
                    throw std::out_of_range(
                        "No as_bool() conversion defined!" );
                }

                virtual std::string as_string( size_t n ) const override {
                    // throw std::out_of_range("No as_string() conversion
                    // defined!");
                    return this->_as_string( n );
                }

                virtual long long as_int( size_t n ) const override {
                    throw std::out_of_range(
                        "No as_int() conversion defined!" );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n ) const override {
                    throw std::out_of_range(
                        "No as_unsigned_int() conversion defined!" );
                }

                virtual double as_double( size_t n ) const override {
                    throw std::out_of_range(
                        "No as_double() conversion defined!" );
                }

                virtual std::complex<double> as_complex(
                    size_t n ) const override {
                    throw std::out_of_range(
                        "No as_complex() conversion defined!" );
                }

                virtual void from_bool( bool x ) override {
                    throw std::out_of_range(
                        "No from_bool() conversion defined!" );
                }

                virtual void from_string( const std::string& x ) override {
                    throw std::out_of_range(
                        "No from_string() conversion defined!" );
                }

                virtual void from_int( long long x ) override {
                    throw std::out_of_range(
                        "No from_int() conversion defined!" );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    throw std::out_of_range(
                        "No from_unsigned_int() conversion defined!" );
                }

                virtual void from_double( double x ) override {
                    throw std::out_of_range(
                        "No from_double() conversion defined!" );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    throw std::out_of_range(
                        "No from_complex() conversion defined!" );
                }

                virtual void from_bool( const std::vector<bool>& x ) override {
                    throw std::out_of_range(
                        "No from_bool() conversion defined!" );
                }

                virtual void from_string(
                    const std::vector<std::string>& x ) override {
                    throw std::out_of_range(
                        "No from_string() conversion defined!" );
                }

                virtual void from_int(
                    const std::vector<long long>& x ) override {
                    throw std::out_of_range(
                        "No from_int() conversion defined!" );
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& x ) override {
                    throw std::out_of_range(
                        "No from_unsigned_int() conversion defined!" );
                }

                virtual void from_double(
                    const std::vector<double>& x ) override {
                    throw std::out_of_range(
                        "No from_double() conversion defined!" );
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& x ) override {
                    throw std::out_of_range(
                        "No from_complex() conversion defined!" );
                }

            protected:
                using hidden::_base_scalar_parameter<PARAMTYPE>::_as_string;
        };

    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap( NCPA::config::hidden::_base_scalar_parameter<T>& a,
           NCPA::config::hidden::_base_scalar_parameter<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::TypedParameter<T>&>( a ),
            static_cast<NCPA::config::TypedParameter<T>&>( b ) );
    swap( a._value, b._value );
}

template<typename T>
void swap( NCPA::config::ScalarParameter<T>& a,
           NCPA::config::ScalarParameter<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::hidden::_base_scalar_parameter<T>&>( a ),
        static_cast<NCPA::config::hidden::_base_scalar_parameter<T>&>( b ) );
}
