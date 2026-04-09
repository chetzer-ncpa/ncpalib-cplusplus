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
            class _base_vector_parameter;
        }
    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap( NCPA::config::hidden::_base_vector_parameter<T>& a,
           NCPA::config::hidden::_base_vector_parameter<T>& b ) noexcept;

namespace NCPA {
    namespace config {

        namespace hidden {
            template<typename PARAMTYPE>
            class _base_vector_parameter : public TypedParameter<PARAMTYPE> {
                public:
                    _base_vector_parameter() : TypedParameter<PARAMTYPE>() {}

                    _base_vector_parameter( PARAMTYPE defaultval ) :
                        TypedParameter<PARAMTYPE>(),
                        _value { std::vector<PARAMTYPE> { defaultval } } {}

                    _base_vector_parameter(
                        const std::vector<PARAMTYPE>& defaultval ) :
                        TypedParameter<PARAMTYPE>(), _value { defaultval } {}

                    _base_vector_parameter( const ValidationTest& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ) {}

                    _base_vector_parameter( const test_ptr_t& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ) {}

                    _base_vector_parameter(
                        std::initializer_list<test_ptr_t> new_tests ) :
                        TypedParameter<PARAMTYPE>( new_tests ) {}

                    _base_vector_parameter( PARAMTYPE defaultval,
                                            const ValidationTest& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ),
                        _value { std::vector<PARAMTYPE> { defaultval } } {}

                    _base_vector_parameter( PARAMTYPE defaultval,
                                            const test_ptr_t& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ),
                        _value { std::vector<PARAMTYPE> { defaultval } } {}

                    _base_vector_parameter(
                        PARAMTYPE defaultval,
                        std::initializer_list<test_ptr_t> new_tests ) :
                        TypedParameter<PARAMTYPE>( new_tests ),
                        _value { std::vector<PARAMTYPE> { defaultval } } {}

                    _base_vector_parameter(
                        const std::vector<PARAMTYPE>& defaultval,
                        const ValidationTest& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ),
                        _value { defaultval } {}

                    _base_vector_parameter(
                        const std::vector<PARAMTYPE>& defaultval,
                        const test_ptr_t& newtest ) :
                        TypedParameter<PARAMTYPE>( newtest ),
                        _value { defaultval } {}

                    _base_vector_parameter(
                        const std::vector<PARAMTYPE>& defaultval,
                        std::initializer_list<test_ptr_t> new_tests ) :
                        TypedParameter<PARAMTYPE>( new_tests ),
                        _value { defaultval } {}

                    virtual ~_base_vector_parameter() {}

                    _base_vector_parameter(
                        const _base_vector_parameter<PARAMTYPE>& other ) :
                        TypedParameter<PARAMTYPE>( other ) {
                        _value = other._value;
                    }

                    _base_vector_parameter(
                        _base_vector_parameter<PARAMTYPE>&& other ) noexcept :
                        TypedParameter<PARAMTYPE>() {
                        ::swap( *this, other );
                    }

                    _base_vector_parameter<PARAMTYPE>& operator=(
                        _base_vector_parameter<PARAMTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    friend void ::swap<>(
                        _base_vector_parameter<PARAMTYPE>& a,
                        _base_vector_parameter<PARAMTYPE>& b ) noexcept;

                    virtual parameter_form_t form() const override {
                        return parameter_form_t::VECTOR;
                    }

                    virtual size_t size() const override {
                        return _value.size();
                    }

                    virtual PARAMTYPE get( size_t n = 0 ) const override {
                        return ( _value.empty() ? this->_nullval()
                                                : _value.at( n ) );
                    }

                    virtual std::vector<PARAMTYPE> get_vector()
                        const override {
                        return _value;
                    }

                protected:
                    std::vector<PARAMTYPE> _value;

                    virtual PARAMTYPE _nullval() const = 0;

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
                        throw std::out_of_range(
                            "No as_string() conversion defined!" );
                    }
            };
        }  // namespace hidden

        // floating point
        template<typename PARAMTYPE>
        class VectorParameter<
            PARAMTYPE, typename std::enable_if<
                           std::is_floating_point<PARAMTYPE>::value>::type>
            : public hidden::_base_vector_parameter<PARAMTYPE> {
            public:
                VectorParameter() :
                    hidden::_base_vector_parameter<PARAMTYPE>() {}

                VectorParameter( PARAMTYPE defaultval ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( defaultval ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( defaultval ) {}

                VectorParameter( const ValidationTest& newtest ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( newtest ) {}

                VectorParameter( const test_ptr_t& newtest ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( newtest ) {}

                VectorParameter(
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( new_tests ) {}

                VectorParameter( PARAMTYPE defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                VectorParameter( PARAMTYPE defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                VectorParameter(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const ValidationTest& newtest ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval,
                                 const test_ptr_t& newtest ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( defaultval,
                                                               newtest ) {}

                VectorParameter(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( defaultval,
                                                               new_tests ) {}

                VectorParameter( const VectorParameter<PARAMTYPE>& other ) :
                    hidden::_base_vector_parameter<PARAMTYPE>( other ) {}

                VectorParameter( VectorParameter<PARAMTYPE>&& other ) noexcept
                    :
                    hidden::_base_vector_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                virtual ~VectorParameter() {}

                VectorParameter<PARAMTYPE>& operator=(
                    VectorParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( VectorParameter<PARAMTYPE>& a,
                                      VectorParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( *this ) );
                }

                virtual long long as_int( size_t n = 0 ) const override {
                    return static_cast<long long>(
                        this->get( n )
                        + std::numeric_limits<PARAMTYPE>::epsilon() );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n = 0 ) const override {
                    return static_cast<unsigned long long>(
                        this->get( n )
                        + std::numeric_limits<PARAMTYPE>::epsilon() );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    return ( this->get( n ) != 0 );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    return std::to_string( this->get( n ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return static_cast<double>( this->get( n ) );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_string( const std::string& x ) override {
                    this->_value[ 0 ]
                        = static_cast<PARAMTYPE>( std::stod( x ) );
                }

                virtual void from_bool( bool x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_int( long long x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_double( double x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value[ 0 ]
                        = static_cast<PARAMTYPE>( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, bool>( vec );
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, long long>(
                            vec );
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE,
                                                    unsigned long long>( vec );
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, double>( vec );
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = static_cast<PARAMTYPE>(
                            std::abs( vec.at( i ) ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = static_cast<PARAMTYPE>(
                            std::stod( vec.at( i ) ) );
                    }
                }

            protected:
                virtual PARAMTYPE _nullval() const override { return 0.0; }
        };

        // signed integer
        template<typename PARAMTYPE>
        class VectorParameter<PARAMTYPE,
                              typename std::enable_if<(
                                  std::is_integral<PARAMTYPE>::value
                                  && !( std::is_same<PARAMTYPE, bool>::value )
                                  && std::is_signed<PARAMTYPE>::value )>::type>
            : public hidden::_base_vector_parameter<PARAMTYPE> {
            public:
                NCPA_CONFIGURATION_VECTORPARAMETER_PUBLIC_BOILERPLATE

                virtual long long as_int( size_t n = 0 ) const override {
                    return static_cast<long long>( this->get( n ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n = 0 ) const override {
                    return static_cast<unsigned long long>( this->get( n ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    return ( this->get( n ) != 0 );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    return std::to_string( this->get( n ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return static_cast<double>( this->get( n ) );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_string( const std::string& x ) override {
                    this->_value[ 0 ]
                        = static_cast<PARAMTYPE>( std::stoll( x ) );
                }

                virtual void from_bool( bool x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_int( long long x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_double( double x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value[ 0 ]
                        = static_cast<PARAMTYPE>( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, bool>( vec );
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, long long>(
                            vec );
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE,
                                                    unsigned long long>( vec );
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, double>( vec );
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = static_cast<PARAMTYPE>(
                            std::abs( vec.at( i ) ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = static_cast<PARAMTYPE>(
                            std::stoll( vec.at( i ) ) );
                    }
                }

            protected:
                virtual PARAMTYPE _nullval() const override { return 0; }
        };

        // unsigned integer
        template<typename PARAMTYPE>
        class VectorParameter<
            PARAMTYPE, typename std::enable_if<(
                           std::is_integral<PARAMTYPE>::value
                           && !( std::is_same<PARAMTYPE, bool>::value )
                           && std::is_unsigned<PARAMTYPE>::value )>::type>
            : public hidden::_base_vector_parameter<PARAMTYPE> {
            public:
                NCPA_CONFIGURATION_VECTORPARAMETER_PUBLIC_BOILERPLATE

                virtual long long as_int( size_t n = 0 ) const override {
                    return static_cast<long long>( this->get( n ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n = 0 ) const override {
                    return static_cast<unsigned long long>( this->get( n ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    return ( this->get( n ) != 0 );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    return std::to_string( this->get( n ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return static_cast<double>( this->get( n ) );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_string( const std::string& x ) override {
                    this->get( 0 )
                        = static_cast<PARAMTYPE>( std::stoull( x ) );
                }

                virtual void from_bool( bool x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_int( long long x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_double( double x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value[ 0 ]
                        = static_cast<PARAMTYPE>( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, bool>( vec );
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, long long>(
                            vec );
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE,
                                                    unsigned long long>( vec );
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, double>( vec );
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = static_cast<PARAMTYPE>(
                            std::abs( vec.at( i ) ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = static_cast<PARAMTYPE>(
                            std::stoull( vec.at( i ) ) );
                    }
                }

            protected:
                virtual PARAMTYPE _nullval() const override { return 0; }
        };

        // // boolean
        template<typename PARAMTYPE>
        class VectorParameter<PARAMTYPE, typename std::enable_if<std::is_same<
                                             PARAMTYPE, bool>::value>::type>
            : public hidden::_base_vector_parameter<PARAMTYPE> {
            public:
                NCPA_CONFIGURATION_VECTORPARAMETER_PUBLIC_BOILERPLATE

                virtual bool as_bool( size_t n = 0 ) const override {
                    return this->get( n );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    return ( this->get( n ) ? "true" : "false" );
                }

                virtual long long as_int( size_t n = 0 ) const override {
                    return static_cast<long long>( this->get( n ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n = 0 ) const override {
                    return static_cast<unsigned long long>( this->get( n ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return static_cast<double>( this->get( n ) );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_bool( bool x ) override {
                    this->_value[ 0 ] = x;
                }

                virtual void from_string( const std::string& x ) override {
                    this->_value[ 0 ] = ( x == "true" );
                }

                virtual void from_int( long long x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_double( double x ) override {
                    this->_value[ 0 ] = static_cast<PARAMTYPE>( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value[ 0 ]
                        = static_cast<PARAMTYPE>( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    this->_value = vec;
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, long long>(
                            vec );
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE,
                                                    unsigned long long>( vec );
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    this->_value
                        = NCPA::arrays::cast_vector<PARAMTYPE, double>( vec );
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = static_cast<PARAMTYPE>(
                            std::abs( vec.at( i ) ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i )
                            = ( NCPA::strings::to_lower( vec.at( i ) )
                                == "true" );
                    }
                }

            protected:
                virtual PARAMTYPE _nullval() const override { return false; }

                //    private:
                //     PARAMTYPE _value;
        };

        // string
        template<typename PARAMTYPE>
        class VectorParameter<
            PARAMTYPE,
            typename std::enable_if<(
                !( std::is_arithmetic<PARAMTYPE>::value )
                && std::is_convertible<PARAMTYPE, std::string>::value )>::type>
            : public hidden::_base_vector_parameter<PARAMTYPE> {
            public:
                NCPA_CONFIGURATION_VECTORPARAMETER_PUBLIC_BOILERPLATE

                virtual bool as_bool( size_t n = 0 ) const override {
                    std::string s = this->get( n );
                    std::transform( s.begin(), s.end(), s.begin(),
                                    []( unsigned char c ) {
                                        return std::tolower( c );
                                    }  // correct
                    );
                    return ( s == "true" );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    std::string s = this->get( n );
                    return s;
                }

                virtual long long as_int( size_t n = 0 ) const override {
                    return std::stoll( this->get( n ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n = 0 ) const override {
                    return std::stoull( this->get( n ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return std::stod( this->get( n ) );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_bool( bool x ) override {
                    this->_value[ 0 ] = ( x ? "true" : "false" );
                }

                virtual void from_string( const std::string& x ) override {
                    this->_value[ 0 ] = x;
                }

                virtual void from_int( long long x ) override {
                    this->_value[ 0 ] = std::to_string( x );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->_value[ 0 ] = std::to_string( x );
                }

                virtual void from_double( double x ) override {
                    this->_value[ 0 ] = std::to_string( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->_value[ 0 ] = std::to_string( std::abs( x ) );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i )
                            = ( vec.at( i ) ? "true" : "false" );
                    }
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = std::to_string( vec.at( i ) );
                    }
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = std::to_string( vec.at( i ) );
                    }
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = std::to_string( vec.at( i ) );
                    }
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i )
                            = std::to_string( std::abs( vec.at( i ) ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ) = vec.at( i );
                    }
                }

            protected:
                virtual PARAMTYPE _nullval() const override { return ""; }
        };

        // complex
        template<typename PARAMTYPE>
        class VectorParameter<
            PARAMTYPE,
            typename std::enable_if<(
                !( std::is_scalar<PARAMTYPE>::value )
                && NCPA::types::is_complex<PARAMTYPE>::value )>::type>
            : public hidden::_base_vector_parameter<PARAMTYPE> {
            public:
                NCPA_CONFIGURATION_VECTORPARAMETER_PUBLIC_BOILERPLATE

                virtual bool as_bool( size_t n = 0 ) const override {
                    return ( std::abs( this->get( n ) ) != 0.0 );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    std::ostringstream oss;
                    oss << std::to_string( this->get( n ).real() ) << " "
                        << ( this->get( n ).imag() < 0.0 ? "- " : "+ " )
                        << std::to_string( this->get( n ).imag() );
                    return oss.str();
                }

                virtual long long as_int( size_t n = 0 ) const override {
                    return static_cast<long long>(
                        std::abs( this->get( n ) ) );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n = 0 ) const override {
                    return static_cast<unsigned long long>(
                        std::abs( this->get( n ) ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return static_cast<double>( std::abs( this->get( n ) ) );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->get( n ).real(),
                                                 this->get( n ).imag() );
                }

                virtual void from_string( const std::string& x ) override {
                    throw std::out_of_range(
                        "from_string() not yet implemented for complex "
                        "numbers" );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->get( 0 ).real( x.real() );
                    this->get( 0 ).imag( x.imag() );
                }

                virtual void from_bool( bool x ) override {
                    this->get( 0 ).real( static_cast<double>( x ) );
                }

                virtual void from_int( long long x ) override {
                    this->get( 0 ).real( static_cast<double>( x ) );
                }

                virtual void from_unsigned_int(
                    unsigned long long x ) override {
                    this->get( 0 ).real( static_cast<double>( x ) );
                }

                virtual void from_double( double x ) override {
                    this->get( 0 ).real( x );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ).real(
                            static_cast<double>( vec.at( i ) ) );
                    }
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ).real(
                            static_cast<double>( vec.at( i ) ) );
                    }
                }

                virtual void from_unsigned_int(
                    const std::vector<unsigned long long>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ).real(
                            static_cast<double>( vec.at( i ) ) );
                    }
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ).real( vec.at( i ) );
                    }
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    this->_value.resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->_value.at( i ).real(
                            static_cast<double>( vec.at( i ).real() ) );
                        this->_value.at( i ).imag(
                            static_cast<double>( vec.at( i ).imag() ) );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    throw std::out_of_range(
                        "from_string() not yet implemented for complex "
                        "numbers" );
                }

            protected:
                virtual PARAMTYPE _nullval() const override {
                    return PARAMTYPE { 0, 0 };
                }
        };

        // everything else
        template<typename PARAMTYPE>
        class VectorParameter<
            PARAMTYPE,
            typename std::enable_if<( !(
                std::is_arithmetic<PARAMTYPE>::value
                || std::is_convertible<PARAMTYPE, std::string>::value
                || ( !( std::is_scalar<PARAMTYPE>::value )
                     && NCPA::types::is_complex<PARAMTYPE>::value ) ) )>::type>
            : public hidden::_base_vector_parameter<PARAMTYPE> {
            public:
                NCPA_CONFIGURATION_VECTORPARAMETER_PUBLIC_BOILERPLATE

                virtual bool as_bool( size_t n = 0 ) const override {
                    throw std::out_of_range(
                        "No as_bool() conversion defined!" );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    // throw std::out_of_range("No as_string() conversion
                    // defined!");
                    return this->_as_string( n );
                }

                virtual long long as_int( size_t n = 0 ) const override {
                    throw std::out_of_range(
                        "No as_int() conversion defined!" );
                }

                virtual unsigned long long as_unsigned_int(
                    size_t n = 0 ) const override {
                    throw std::out_of_range(
                        "No as_unsigned_int() conversion defined!" );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    throw std::out_of_range(
                        "No as_double() conversion defined!" );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
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
                using hidden::_base_vector_parameter<PARAMTYPE>::_as_string;

                virtual PARAMTYPE _nullval() const override {
                    throw std::out_of_range( "Vector is empty and no null "
                                             "value has been defined!" );
                }
        };

    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap( NCPA::config::hidden::_base_vector_parameter<T>& a,
           NCPA::config::hidden::_base_vector_parameter<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::config::TypedParameter<T>&>( a ),
          static_cast<NCPA::config::TypedParameter<T>&>( b ) );
    swap( a._value, b._value );
}

template<typename T>
void swap( NCPA::config::VectorParameter<T>& a,
           NCPA::config::VectorParameter<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::config::hidden::_base_vector_parameter<T>&>( a ),
          static_cast<NCPA::config::hidden::_base_vector_parameter<T>&>( b ) );
}
