#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/Parameter.hpp"
#include "NCPA/configuration/TypedParameter.hpp"

#include <vector>

namespace NCPA {
    namespace config {
        namespace implementation {
            template<typename T>
            class _base_vector_parameter;
        }
    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap(
    NCPA::config::implementation::_base_vector_parameter<T>& a,
    NCPA::config::implementation::_base_vector_parameter<T>& b ) noexcept;

namespace NCPA {
    namespace config {
        namespace implementation {
            template<typename PARAMTYPE>
            class _base_vector_parameter : public TypedParameter<PARAMTYPE> {
                public:
                    _base_vector_parameter() :
                        TypedParameter<PARAMTYPE>(),
                        _value { std::vector<PARAMTYPE>( 1 ) } {}

                    _base_vector_parameter( PARAMTYPE defaultval ) :
                        TypedParameter<PARAMTYPE>(),
                        _value { std::vector<PARAMTYPE>( 1, defaultval ) } {}

                    _base_vector_parameter(
                        const std::vector<PARAMTYPE>& defaultval ) :
                        TypedParameter<PARAMTYPE>(), _value { defaultval } {}

                    virtual ~_base_vector_parameter() {}

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

                    virtual ~_base_vector_parameter() {}

                    friend void ::swap(
                        _base_vector_parameter<PARAMTYPE>& a,
                        _base_vector_parameter<PARAMTYPE>& b ) noexcept;

                    virtual parameter_form_t form() const override {
                        return parameter_form_t::VECTOR;
                    }

                    virtual size_t size() const override {
                        return _value.size();
                    }

                    virtual std::vector<PARAMTYPE>& value() { return _value; }

                    virtual const std::vector<PARAMTYPE>& value() const {
                        return _value;
                    }

                    virtual PARAMTYPE value( size_t n ) const {
                        return _value.at( n );
                    }

                    virtual bool as_bool( size_t n = 0 ) const override {
                        return ( this->value( n ) != 0 );
                    }

                    virtual std::string as_string( size_t n
                                                   = 0 ) const override {
                        return std::to_string( this->value( n ) );
                    }

                    virtual int as_int( size_t n = 0 ) const override {
                        return static_cast<int>( this->value( n ) );
                    }

                    virtual long long as_long_int( size_t n
                                                   = 0 ) const override {
                        return static_cast<long long>( this->value( n ) );
                    }

                    virtual size_t as_size_t( size_t n = 0 ) const override {
                        return static_cast<size_t>( this->as_long_int( n ) );
                    }

                    virtual double as_double( size_t n = 0 ) const override {
                        return static_cast<double>( this->value( n ) );
                    }

                    virtual std::complex<double> as_complex(
                        size_t n = 0 ) const override {
                        return std::complex<double>( this->as_double( n ),
                                                     0.0 );
                    }

                    virtual void from_bool( bool x ) override {
                        this->from_bool( std::vector<bool>( 1, x ) );
                    }

                    virtual void from_string( const std::string& x ) override {
                        this->from_string( std::vector<std::string>( 1, x ) );
                    }

                    virtual void from_int( long long x ) override {
                        this->from_int( std::vector<long long>( 1, x ) );
                    }

                    virtual void from_unsigned_int( size_t x ) override {
                        this->from_unsigned_int( std::vector<size_t>( 1, x ) );
                    }

                    virtual void from_double( double x ) override {
                        this->from_double( std::vector<double>( 1, x ) );
                    }

                    virtual void from_complex(
                        std::complex<double> x ) override {
                        this->from_complex(
                            std::vector<std::complex<double>>( 1, x ) );
                    }

                    virtual void from_bool(
                        const std::vector<bool>& vec ) override {
                        this->value().resize( vec.size() );
                        for (size_t i = 0; i < vec.size(); ++i) {
                            this->value().at( i )
                                = static_cast<PARAMTYPE>( vec.at( i ) );
                        }
                    }

                    virtual void from_int(
                        const std::vector<long long>& vec ) override {
                        this->value().resize( vec.size() );
                        for (size_t i = 0; i < vec.size(); ++i) {
                            this->value().at( i )
                                = static_cast<PARAMTYPE>( vec.at( i ) );
                        }
                    }

                    virtual void from_unsigned_int(
                        const std::vector<size_t>& vec ) override {
                        this->value().resize( vec.size() );
                        for (size_t i = 0; i < vec.size(); ++i) {
                            this->value().at( i )
                                = static_cast<PARAMTYPE>( vec.at( i ) );
                        }
                    }

                    virtual void from_double(
                        const std::vector<double>& vec ) override {
                        this->value().resize( vec.size() );
                        for (size_t i = 0; i < vec.size(); ++i) {
                            this->value().at( i )
                                = static_cast<PARAMTYPE>( vec.at( i ) );
                        }
                    }

                    virtual void from_complex(
                        const std::vector<std::complex<double>>& vec )
                        override {
                        this->value().resize( vec.size() );
                        for (size_t i = 0; i < vec.size(); ++i) {
                            this->value().at( i ) = static_cast<PARAMTYPE>(
                                std::abs( vec.at( i ) ) );
                        }
                    }

                protected:
                    std::vector<PARAMTYPE> _value;
            };
        }  // namespace implementation

        // floating point vector
        template<typename PARAMTYPE>
        class VectorParameter<
            PARAMTYPE, typename std::enable_if<
                           std::is_floating_point<PARAMTYPE>::value>::type>
            : public implementation::_base_vector_parameter {
            public:
                VectorParameter() :
                    implementation::_base_vector_parameter<PARAMTYPE>() {}

                VectorParameter( PARAMTYPE defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( VectorParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_vector_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                VectorParameter<PARAMTYPE>& operator=(
                    VectorParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_vector_parameter() {}

                friend void ::swap( VectorParameter<PARAMTYPE>& a,
                                    VectorParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( *this ) );
                }

                virtual int as_int( size_t n = 0 ) const override {
                    return static_cast<int>(
                        this->value( n )
                        + std::numeric_limits<PARAMTYPE>::epsilon() );
                }

                virtual long long as_long_int( size_t n = 0 ) const override {
                    return static_cast<long long>(
                        this->value( n )
                        + std::numeric_limits<PARAMTYPE>::epsilon() );
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i ) = static_cast<PARAMTYPE>(
                            std::stod( vec.at( i ) ) );
                    }
                }
        };

        // signed integer vector
        template<typename PARAMTYPE>
        class VectorParameter<PARAMTYPE,
                              typename std::enable_if<(
                                  std::is_integral<PARAMTYPE>::value
                                  && !( std::is_same<PARAMTYPE, bool>::value )
                                  && std::is_signed<PARAMTYPE>::value )>::type>
            : public implementation::_base_vector_parameter {
            public:
                VectorParameter() :
                    implementation::_base_vector_parameter<PARAMTYPE>() {}

                VectorParameter( PARAMTYPE defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( VectorParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_vector_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                VectorParameter<PARAMTYPE>& operator=(
                    VectorParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_vector_parameter() {}

                friend void ::swap( VectorParameter<PARAMTYPE>& a,
                                    VectorParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( *this ) );
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i ) = static_cast<PARAMTYPE>(
                            std::stoll( vec.at( i ) ) );
                    }
                }

                //    private:
                //     PARAMTYPE _value;
        };

        // unsigned integer
        template<typename PARAMTYPE>
        class VectorParameter<
            PARAMTYPE, typename std::enable_if<(
                           std::is_integral<PARAMTYPE>::value
                           && !( std::is_same<PARAMTYPE, bool>::value )
                           && std::is_unsigned<PARAMTYPE>::value )>::type>
            : public implementation::_base_vector_parameter {
            public:
                VectorParameter() :
                    implementation::_base_vector_parameter<PARAMTYPE>() {}

                VectorParameter( PARAMTYPE defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( VectorParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_vector_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                VectorParameter<PARAMTYPE>& operator=(
                    VectorParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_vector_parameter() {}

                friend void ::swap( VectorParameter<PARAMTYPE>& a,
                                    VectorParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( *this ) );
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i ) = static_cast<PARAMTYPE>(
                            std::stoull( vec.at( i ) ) );
                    }
                }
        };

        // boolean
        template<typename PARAMTYPE>
        class VectorParameter<PARAMTYPE, typename std::enable_if<std::is_same<
                                             PARAMTYPE, bool>::value>::type>
            : public implementation::_base_vector_parameter {
            public:
                VectorParameter() :
                    implementation::_base_vector_parameter<PARAMTYPE>() {}

                VectorParameter( PARAMTYPE defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( VectorParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_vector_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                VectorParameter<PARAMTYPE>& operator=(
                    VectorParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_vector_parameter() {}

                friend void ::swap( VectorParameter<PARAMTYPE>& a,
                                    VectorParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    return this->value( n );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    return ( this->value( n ) ? "true" : "false" );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    this->value() = vec;
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i )
                            = ( NCPA::strings::to_lower( vec.at( i ) )
                                == "true" );
                    }
                }

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
            : public implementation::_base_vector_parameter {
            public:
                VectorParameter() :
                    implementation::_base_vector_parameter<PARAMTYPE>() {}

                VectorParameter( PARAMTYPE defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( VectorParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_vector_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                VectorParameter<PARAMTYPE>& operator=(
                    VectorParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_vector_parameter() {}

                friend void ::swap( VectorParameter<PARAMTYPE>& a,
                                    VectorParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    std::string s = this->value( n );
                    std::transform( s.begin(), s.end(), s.begin(),
                                    []( unsigned char c ) {
                                        return std::tolower( c );
                                    }  // correct
                    );
                    return s == "true";
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    std::string s = this->value( n );
                    return s;
                }

                virtual int as_int( size_t n = 0 ) const override {
                    return std::stoi( this->value( n ) );
                }

                virtual long long as_long_int( size_t n = 0 ) const override {
                    return std::stol( this->value( n ) );
                }

                virtual size_t as_size_t( size_t n = 0 ) const override {
                    return static_cast<size_t>( this->as_long_int( n ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return std::stod( this->value( n ) );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->as_double( n ), 0.0 );
                }

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i )
                            = ( vec.at( i ) ? "true" : "false" );
                    }
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i ) = vec.at( i );
                    }
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i ) = std::to_string( vec.at( i ) );
                    }
                }

                virtual void from_unsigned_int(
                    const std::vector<size_t>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i ) = std::to_string( vec.at( i ) );
                    }
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i ) = std::to_string( vec.at( i ) );
                    }
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i )
                            = std::to_string( vec.at( i ).real() )
                            + ( vec.at( i ).imag() < 0.0 ? " - " : " + " )
                            + std::to_string( std::abs( vec.at( i ).imag() ) );
                    }
                }

                //    private:
                //     PARAMTYPE _value;
        };

        // complex
        template<typename PARAMTYPE>
        class VectorParameter<
            PARAMTYPE,
            typename std::enable_if<(
                !( std::is_scalar<PARAMTYPE>::value )
                && NCPA::types::is_complex<PARAMTYPE>::value )>::type>
            : public implementation::_base_vector_parameter {
            public:
                VectorParameter() :
                    implementation::_base_vector_parameter<PARAMTYPE>() {}

                VectorParameter( PARAMTYPE defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( VectorParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_vector_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                VectorParameter<PARAMTYPE>& operator=(
                    VectorParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_vector_parameter() {}

                friend void ::swap( VectorParameter<PARAMTYPE>& a,
                                    VectorParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    return ( std::abs( this->value( n ) ) != 0.0 );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    std::ostringstream oss;
                    oss << std::to_string( this->value( n ).real() ) << " "
                        << ( this->value( n ).imag() < 0.0 ? "- " : "+ " )
                        << std::to_string( this->value( n ).imag() );
                    return oss.str();
                }

                virtual int as_int( size_t n = 0 ) const override {
                    return static_cast<int>( std::abs( this->value( n ) ) );
                }

                virtual long long as_long_int( size_t n = 0 ) const override {
                    return static_cast<long long>(
                        std::abs( this->value( n ) ) );
                }

                virtual size_t as_size_t( size_t n = 0 ) const override {
                    return static_cast<size_t>( this->as_long_int( n ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return static_cast<double>( std::abs( this->value( n ) ) );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->value( n ).real(),
                                                 this->value( n ).imag() );
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    throw std::out_of_range(
                        "No from_string() conversion implemented for complex "
                        "numbers" );
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    this->value().resize( vec.size() );
                    for (size_t i = 0; i < vec.size(); ++i) {
                        this->value().at( i ) = vec.at( i );
                    }
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
            : public implementation::_base_vector_parameter {
            public:
                VectorParameter() :
                    implementation::_base_vector_parameter<PARAMTYPE>() {}

                VectorParameter( PARAMTYPE defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_vector_parameter<PARAMTYPE>(
                        defaultval ) {}

                VectorParameter( VectorParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_vector_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                VectorParameter<PARAMTYPE>& operator=(
                    VectorParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_vector_parameter() {}

                friend void ::swap( VectorParameter<PARAMTYPE>& a,
                                    VectorParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new VectorParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    throw std::out_of_range(
                        "No as_bool() conversion defined!" );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    return this->_as_string( n );
                }

                virtual int as_int( size_t n = 0 ) const override {
                    throw std::out_of_range(
                        "No as_int() conversion defined!" );
                }

                virtual long long as_long_int( size_t n = 0 ) const override {
                    throw std::out_of_range(
                        "No as_long_int() conversion defined!" );
                }

                virtual size_t as_size_t( size_t n = 0 ) const override {
                    return static_cast<size_t>( this->as_long_int( n ) );
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

                virtual void from_bool(
                    const std::vector<bool>& vec ) override {
                    throw std::out_of_range(
                        "No from_bool() conversion defined!" );
                }

                virtual void from_string(
                    const std::vector<std::string>& vec ) override {
                    throw std::out_of_range(
                        "No from_string() conversion defined!" );
                }

                virtual void from_int(
                    const std::vector<long long>& vec ) override {
                    throw std::out_of_range(
                        "No from_int() conversion defined!" );
                }

                virtual void from_unsigned_int(
                    const std::vector<size_t>& vec ) override {
                    throw std::out_of_range(
                        "No from_unsigned_int() conversion defined!" );
                }

                virtual void from_double(
                    const std::vector<double>& vec ) override {
                    throw std::out_of_range(
                        "No from_double() conversion defined!" );
                }

                virtual void from_complex(
                    const std::vector<std::complex<double>>& vec ) override {
                    throw std::out_of_range(
                        "No from_complex() conversion defined!" );
                }
        };

        using DoubleVectorParameter = VectorParameter<double>;
    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap(
    NCPA::config::implementation::_base_vector_parameter<T>& a,
    NCPA::config::implementation::_base_vector_parameter<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::config::TypedParameter<T>&>( a ),
          static_cast<NCPE::config::TypedParameter<T>&>( b ) );
    swap( a._value, b._value );
}

template<typename T>
void swap( NCPA::config::VectorParameter<T>& a,
           NCPA::config::VectorParameter<T>& b ) noexcept {
    using std::swap;
    swap(
        static_cast<NCPA::config::implementation::_base_vector_parameter<T>&>(
            a ),
        static_cast<NCPE::config::implementation::_base_vector_parameter<T>&>(
            b ) );
}
