#pragma once
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/Parameter.hpp"
#include "NCPA/configuration/TypedParameter.hpp"

#include <regex>
#include <vector>

namespace NCPA {
    namespace config {
        namespace implementation {
            template<typename T>
            class _base_scalar_parameter;
        }
    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap(
    NCPA::config::implementation::_base_scalar_parameter<T>& a,
    NCPA::config::implementation::_base_scalar_parameter<T>& b ) noexcept;

namespace NCPA {
    namespace config {
        namespace implementation {
            template<typename PARAMTYPE>
            class _base_scalar_parameter : public TypedParameter<PARAMTYPE> {
                public:
                    _base_scalar_parameter() : TypedParameter<PARAMTYPE>() {}

                    _base_scalar_parameter( PARAMTYPE defaultval ) :
                        TypedParameter<PARAMTYPE>(), _value { defaultval } {}

                    _base_scalar_parameter(
                        const std::vector<PARAMTYPE>& defaultval ) :
                        TypedParameter<PARAMTYPE>() {
                        if (defaultval.size() > 0) {
                            _value = defaultval.at( 0 );
                        }
                    }

                    virtual ~_base_scalar_parameter() {}

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

                    virtual ~_base_scalar_parameter() {}

                    friend void ::swap(
                        _base_scalar_parameter<PARAMTYPE>& a,
                        _base_scalar_parameter<PARAMTYPE>& b ) noexcept;

                    virtual parameter_form_t form() const override {
                        return parameter_form_t::SCALAR;
                    }

                    virtual size_t size() const override { return 1; }

                    virtual PARAMTYPE& value( size_t n = 0 ) { return _value; }

                    virtual const PARAMTYPE& value( size_t n = 0 ) const {
                        return _value;
                    }

                    virtual bool as_bool( size_t n = 0 ) const override {
                        return ( this->value() != 0 );
                    }

                    virtual std::string as_string( size_t n
                                                   = 0 ) const override {
                        return std::to_string( this->value() );
                    }

                    virtual int as_int( size_t n = 0 ) const override {
                        return static_cast<int>( this->value() );
                    }

                    virtual long long as_long_int( size_t n
                                                   = 0 ) const override {
                        return static_cast<long long>( this->value() );
                    }

                    virtual size_t as_size_t( size_t n = 0 ) const override {
                        return static_cast<size_t>( this->as_long_int( n ) );
                    }

                    virtual double as_double( size_t n = 0 ) const override {
                        return static_cast<double>( this->value() );
                    }

                    virtual std::complex<double> as_complex(
                        size_t n = 0 ) const override {
                        return std::complex<double>( this->as_double(), 0.0 );
                    }

                    virtual void from_bool( bool x ) override {
                        this->value() = static_cast<PARAMTYPE>( x );
                    }

                    // virtual void from_string( const std::string& x )
                    // override {
                    //     this->value() = static_cast<PARAMTYPE>( std::stod( x
                    //     ) );
                    // }

                    virtual void from_int( long long x ) override {
                        this->value() = static_cast<PARAMTYPE>( x );
                    }

                    virtual void from_unsigned_int( size_t x ) override {
                        this->value() = static_cast<PARAMTYPE>( x );
                    }

                    virtual void from_double( double x ) override {
                        this->value() = static_cast<PARAMTYPE>( x );
                    }

                    virtual void from_complex(
                        std::complex<double> x ) override {
                        this->value()
                            = static_cast<PARAMTYPE>( std::abs( x ) );
                    }

                    virtual void from_bool(
                        const std::vector<bool>& vec ) override {
                        if (vec.size() > 0) {
                            this->from_bool( vec.at( 0 ) );
                        }
                    }

                    virtual void from_string(
                        const std::vector<std::string>& vec ) override {
                        if (vec.size() > 0) {
                            this->from_string( vec.at( 0 ) );
                        }
                    }

                    virtual void from_int(
                        const std::vector<long long>& vec ) override {
                        if (vec.size() > 0) {
                            this->from_int( vec.at( 0 ) );
                        }
                    }

                    virtual void from_unsigned_int(
                        const std::vector<size_t>& vec ) override {
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
                        const std::vector<std::complex<double>>& vec )
                        override {
                        if (vec.size() > 0) {
                            this->from_complex( vec.at( 0 ) );
                        }
                    }

                protected:
                    PARAMTYPE _value;

                    template<typename T = PARAMTYPE,
                             typename std::enable_if<
                                 NCPA::types::has_to_string<T>::value,
                                 int>::type = 0>
                    std::string _as_string( size_t n = 0 ) const {
                        return to_string( this->value() );
                    }

                    template<typename T = PARAMTYPE,
                             typename std::enable_if<
                                 !( NCPA::types::has_to_string<T>::value ),
                                 int>::type = 0>
                    std::string _as_string( size_t n = 0 ) const {
                        throw std::out_of_range(
                            "No as_string() conversion defined!" );
                    }
            };
        }  // namespace implementation

        // floating point
        template<typename PARAMTYPE>
        class ScalarParameter<
            PARAMTYPE, typename std::enable_if<
                           std::is_floating_point<PARAMTYPE>::value>::type>
            : public implementation::_base_scalar_parameter<PARAMTYPE> {
            public:
                ScalarParameter() :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_scalar_parameter() {}

                friend void ::swap( ScalarParameter<PARAMTYPE>& a,
                                    ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                // virtual bool as_bool( size_t n = 0 ) const override {
                //     return ( this->value() != 0.0 );
                // }

                // virtual std::string as_string( size_t n = 0 ) const override
                // {
                //     return std::to_string( this->value() );
                // }

                virtual int as_int( size_t n = 0 ) const override {
                    return static_cast<int>(
                        this->value()
                        + std::numeric_limits<PARAMTYPE>::epsilon() );
                }

                virtual long long as_long_int( size_t n = 0 ) const override {
                    return static_cast<long long>(
                        this->value()
                        + std::numeric_limits<PARAMTYPE>::epsilon() );
                }

                // virtual size_t as_size_t( size_t n = 0 ) const override {
                //     return static_cast<size_t>( this->as_long_int( n ) );
                // }

                // virtual double as_double( size_t n = 0 ) const override {
                //     return static_cast<double>( this->value() );
                // }

                // virtual std::complex<double> as_complex( size_t n
                //                                          = 0 ) const
                //                                          override {
                //     return std::complex<double>( this->as_double(), 0.0 );
                // }

                virtual void from_string( const std::string& x ) override {
                    this->value() = static_cast<PARAMTYPE>( std::stod( x ) );
                }
        };

        // signed integer
        template<typename PARAMTYPE>
        class ScalarParameter<PARAMTYPE,
                              typename std::enable_if<(
                                  std::is_integral<PARAMTYPE>::value
                                  && !( std::is_same<PARAMTYPE, bool>::value )
                                  && std::is_signed<PARAMTYPE>::value )>::type>
            : public implementation::_base_scalar_parameter<PARAMTYPE> {
            public:
                ScalarParameter() :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_scalar_parameter() {}

                friend void ::swap( ScalarParameter<PARAMTYPE>& a,
                                    ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                // virtual bool as_bool( size_t n = 0 ) const override {
                //     return ( this->value() != 0 );
                // }

                // virtual std::string as_string( size_t n = 0 ) const override
                // {
                //     return std::to_string( this->value() );
                // }

                // virtual int as_int( size_t n = 0 ) const override {
                //     return static_cast<int>( this->value() );
                // }

                // virtual long long as_long_int( size_t n = 0 ) const override
                // {
                //     return static_cast<long long>( this->value() );
                // }

                // virtual size_t as_size_t( size_t n = 0 ) const override {
                //     return static_cast<size_t>( this->as_long_int( n ) );
                // }

                // virtual double as_double( size_t n = 0 ) const override {
                //     return static_cast<double>( this->value() );
                // }

                // virtual std::complex<double> as_complex( size_t n
                //                                          = 0 ) const
                //                                          override {
                //     return std::complex<double>( this->as_double(), 0.0 );
                // }

                // virtual void from_bool( bool x ) override {
                //     this->value() = static_cast<PARAMTYPE>( x );
                // }

                virtual void from_string( const std::string& x ) override {
                    this->value() = static_cast<PARAMTYPE>( std::stol( x ) );
                }

                // virtual void from_int( long long x ) override {
                //     this->value() = static_cast<PARAMTYPE>( x );
                // }

                // virtual void from_unsigned_int( size_t x ) override {
                //     this->value() = static_cast<PARAMTYPE>( x );
                // }

                // virtual void from_double( double x ) override {
                //     this->value() = static_cast<PARAMTYPE>( x );
                // }

                // virtual void from_complex( std::complex<double> x ) override
                // {
                //     this->value() = static_cast<PARAMTYPE>( std::abs( x ) );
                // }
        };

        // unsigned integer
        template<typename PARAMTYPE>
        class ScalarParameter<
            PARAMTYPE, typename std::enable_if<(
                           std::is_integral<PARAMTYPE>::value
                           && !( std::is_same<PARAMTYPE, bool>::value )
                           && std::is_unsigned<PARAMTYPE>::value )>::type>
            : public implementation::_base_scalar_parameter<PARAMTYPE> {
            public:
                ScalarParameter() :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_scalar_parameter() {}

                friend void ::swap( ScalarParameter<PARAMTYPE>& a,
                                    ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                // virtual bool as_bool( size_t n = 0 ) const override {
                //     return ( this->value() != 0 );
                // }

                // virtual std::string as_string( size_t n = 0 ) const override
                // {
                //     return std::to_string( this->value() );
                // }

                // virtual int as_int( size_t n = 0 ) const override {
                //     return static_cast<int>( this->value() );
                // }

                // virtual long long as_long_int( size_t n = 0 ) const override
                // {
                //     return static_cast<long long>( this->value() );
                // }

                // virtual size_t as_size_t( size_t n = 0 ) const override {
                //     return static_cast<size_t>( this->as_long_int( n ) );
                // }

                // virtual double as_double( size_t n = 0 ) const override {
                //     return static_cast<double>( this->value() );
                // }

                // virtual std::complex<double> as_complex( size_t n
                //                                          = 0 ) const
                //                                          override {
                //     return std::complex<double>( this->as_double(), 0.0 );
                // }

                // virtual void from_bool( bool x ) override {
                //     this->value() = static_cast<PARAMTYPE>( x );
                // }

                virtual void from_string( const std::string& x ) override {
                    this->value() = static_cast<PARAMTYPE>( std::stoull( x ) );
                }

                // virtual void from_int( long long x ) override {
                //     this->value() = static_cast<PARAMTYPE>( x );
                // }

                // virtual void from_unsigned_int( size_t x ) override {
                //     this->value() = static_cast<PARAMTYPE>( x );
                // }

                // virtual void from_double( double x ) override {
                //     this->value() = static_cast<PARAMTYPE>( x );
                // }

                // virtual void from_complex( std::complex<double> x ) override
                // {
                //     this->value() = static_cast<PARAMTYPE>( std::abs( x ) );
                // }
        };

        // boolean
        template<typename PARAMTYPE>
        class ScalarParameter<PARAMTYPE, typename std::enable_if<std::is_same<
                                             PARAMTYPE, bool>::value>::type>
            : public implementation::_base_scalar_parameter<PARAMTYPE> {
            public:
                ScalarParameter() :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_scalar_parameter() {}

                friend void ::swap( ScalarParameter<PARAMTYPE>& a,
                                    ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    return this->value();
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    return ( this->value() ? "true" : "false" );
                }

                // virtual int as_int( size_t n = 0 ) const override {
                //     return static_cast<int>( this->value() );
                // }

                // virtual long long as_long_int( size_t n = 0 ) const override
                // {
                //     return static_cast<long long>( this->value() );
                // }

                // virtual size_t as_size_t( size_t n = 0 ) const override {
                //     return static_cast<size_t>( this->value() );
                // }

                // virtual double as_double( size_t n = 0 ) const override {
                //     return static_cast<double>( this->value() );
                // }

                // virtual std::complex<double> as_complex( size_t n
                //                                          = 0 ) const
                //                                          override {
                //     return std::complex<double>( this->as_double(), 0.0 );
                // }

                virtual void from_bool( bool x ) override {
                    this->value() = x;
                }

                virtual void from_string( const std::string& x ) override {
                    this->value() = ( NCPA::strings::to_lower( x ) == "true" );
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
            : public implementation::_base_scalar_parameter<PARAMTYPE> {
            public:
                ScalarParameter() :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_scalar_parameter() {}

                friend void ::swap( ScalarParameter<PARAMTYPE>& a,
                                    ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    std::string s = this->value();
                    std::transform( s.begin(), s.end(), s.begin(),
                                    []( unsigned char c ) {
                                        return std::tolower( c );
                                    }  // correct
                    );
                    return s == "true";
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    std::string s = this->value();
                    return s;
                }

                virtual int as_int( size_t n = 0 ) const override {
                    return std::stoi( this->value() );
                }

                virtual long long as_long_int( size_t n = 0 ) const override {
                    return std::stoll( this->value() );
                }

                virtual size_t as_size_t( size_t n = 0 ) const override {
                    return static_cast<size_t>( this->as_long_int( n ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return std::stod( this->value() );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->as_double(), 0.0 );
                }

                virtual void from_bool( bool x ) override {
                    this->value() = ( x ? "true" : "false" );
                }

                virtual void from_string( const std::string& x ) override {
                    this->value() = x;
                }

                virtual void from_int( long long x ) override {
                    this->value() = std::to_string( x );
                }

                virtual void from_unsigned_int( size_t x ) override {
                    this->value() = std::to_string( x );
                }

                virtual void from_double( double x ) override {
                    this->value() = std::to_string( x );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->value() = std::to_string( std::abs( x ) );
                }

                //    private:
                //     PARAMTYPE _value;
        };

        // complex
        template<typename PARAMTYPE>
        class ScalarParameter<
            PARAMTYPE,
            typename std::enable_if<(
                !( std::is_scalar<PARAMTYPE>::value )
                && NCPA::types::is_complex<PARAMTYPE>::value )>::type>
            : public implementation::_base_scalar_parameter<PARAMTYPE> {
            public:
                ScalarParameter() :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_scalar_parameter() {}

                friend void ::swap( ScalarParameter<PARAMTYPE>& a,
                                    ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    return ( std::abs( this->value() ) != 0.0 );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    std::ostringstream oss;
                    oss << std::to_string( this->value().real() ) << " "
                        << ( this->value().imag() < 0.0 ? "- " : "+ " )
                        << std::to_string( this->value().imag() );
                    return oss.str();
                }

                virtual int as_int( size_t n = 0 ) const override {
                    return static_cast<int>( std::abs( this->value() ) );
                }

                virtual long long as_long_int( size_t n = 0 ) const override {
                    return static_cast<long long>( std::abs( this->value() ) );
                }

                virtual size_t as_size_t( size_t n = 0 ) const override {
                    return static_cast<size_t>( this->as_long_int( n ) );
                }

                virtual double as_double( size_t n = 0 ) const override {
                    return static_cast<double>( std::abs( this->value() ) );
                }

                virtual std::complex<double> as_complex( size_t n
                                                         = 0 ) const override {
                    return std::complex<double>( this->value().real(),
                                                 this->value().imag() );
                }

                virtual void from_string( const std::string& x ) override {
                    throw std::out_of_range(
                        "from_string() not yet implemented for complex "
                        "numbers" );
                }

                virtual void from_complex( std::complex<double> x ) override {
                    this->value()
                        = std::complex<PARAMTYPE>( x.real(), x.imag() );
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
            : public implementation::_base_scalar_parameter<PARAMTYPE> {
            public:
                ScalarParameter() :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {}

                ScalarParameter( PARAMTYPE defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( const std::vector<PARAMTYPE>& defaultval ) :
                    implementation::_base_scalar_parameter<PARAMTYPE>(
                        defaultval ) {}

                ScalarParameter( ScalarParameter<PARAMTYPE>&& other ) noexcept
                    :
                    implementation::_base_scalar_parameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                ScalarParameter<PARAMTYPE>& operator=(
                    ScalarParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~_base_scalar_parameter() {}

                friend void ::swap( ScalarParameter<PARAMTYPE>& a,
                                    ScalarParameter<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new ScalarParameter<PARAMTYPE>( *this ) );
                }

                virtual bool as_bool( size_t n = 0 ) const override {
                    throw std::out_of_range(
                        "No as_bool() conversion defined!" );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    // throw std::out_of_range("No as_string() conversion
                    // defined!");
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

                virtual void from_unsigned_int( size_t x ) override {
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

                //    private:
                //     PARAMTYPE _value;
        };

        using DoubleParameter  = ScalarParameter<double>;
        using IntegerParameter = ScalarParameter<int>;
        using StringParameter  = ScalarParameter<std::string>;
        using BooleanParameter = ScalarParameter<bool>;
        // using UnitsParameter        =
        // ScalarParameter<NCPA::units::units_ptr_t>;

    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap(
    NCPA::config::implementation::_base_scalar_parameter<T>& a,
    NCPA::config::implementation::_base_scalar_parameter<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::config::TypedParameter<T>&>( a ),
          static_cast<NCPE::config::TypedParameter<T>&>( b ) );
    swap( a._value, b._value );
}

template<typename T>
void swap( NCPA::config::ScalarParameter<T>& a,
           NCPA::config::ScalarParameter<T>& b ) noexcept {
    using std::swap;
    swap(
        static_cast<NCPA::config::implementation::_base_scalar_parameter<T>&>(
            a ),
        static_cast<NCPE::config::implementation::_base_scalar_parameter<T>&>(
            b ) );
}
