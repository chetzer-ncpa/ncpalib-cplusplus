#pragma once

#include "NCPA/configuration/declarations.hpp"

#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#if HAVE_NLOHMANN_JSON_HPP
#  include "nlohmann/json.hpp"
#endif

void swap( NCPA::config::Parameter& a, NCPA::config::Parameter& b ) noexcept;
template<typename T>
void swap( NCPA::config::ScalarParameter<T>& a,
           NCPA::config::ScalarParameter<T>& b ) noexcept;
template<typename T>
void swap( NCPA::config::VectorParameter<T>& a,
           NCPA::config::VectorParameter<T>& b ) noexcept;
void swap( NCPA::config::DoubleParameter& a,
           NCPA::config::DoubleParameter& b ) noexcept;
void swap( NCPA::config::IntegerParameter& a,
           NCPA::config::IntegerParameter& b ) noexcept;
void swap( NCPA::config::StringParameter& a,
           NCPA::config::StringParameter& b ) noexcept;
void swap( NCPA::config::BooleanParameter& a,
           NCPA::config::BooleanParameter& b ) noexcept;
void swap( NCPA::config::DoubleVectorParameter& a,
           NCPA::config::DoubleVectorParameter& b ) noexcept;
void swap( NCPA::config::IntegerVectorParameter& a,
           NCPA::config::IntegerVectorParameter& b ) noexcept;
void swap( NCPA::config::StringVectorParameter& a,
           NCPA::config::StringVectorParameter& b ) noexcept;
void swap( NCPA::config::BooleanVectorParameter& a,
           NCPA::config::BooleanVectorParameter& b ) noexcept;

namespace NCPA {
    namespace config {
        class Parameter {
            public:
                Parameter() :
                    Parameter( parameter_form_t::UNDEF,
                               parameter_type_t::UNDEF, "", "" ) {}

                // Parameter( const std::string& key ) :
                //     _key { key }, _description { "" } {}

                // Parameter( const std::string& key,
                //            const std::string& description ) :

                //     _key { key }, _description { description } {}

                Parameter( parameter_form_t paramform,
                           parameter_type_t paramtype, const std::string& key,
                           const std::string& description ) :
                    _form { paramform },
                    _type { paramtype },
                    _key { key },
                    _description { description } {}

                Parameter( const Parameter& other ) :
                    _form { other._form },
                    _type { other._type },
                    _key { other._key },
                    _description { other._description } {}

                Parameter( Parameter&& other ) noexcept {
                    ::swap( *this, other );
                }

                virtual ~Parameter() {}

                friend void ::swap( Parameter& a, Parameter& b ) noexcept;

                virtual const std::string& key() const { return _key; }

                virtual const std::string& description() const {
                    return _description;
                }

                virtual size_t asSizeT() const {
                    int val = this->asInt();
                    if (val < 0) {
                        throw std::out_of_range(
                            "Value is negative, can't retrieve as size_t!" );
                    } else {
                        return static_cast<size_t>( val );
                    }
                }

                virtual std::vector<size_t> asSizeTVector() const {
                    std::vector<int> vi = this->asIntVector();
                    std::vector<size_t> vst( vi.size() );
                    for (size_t i = 0; i < vi.size(); ++i) {
                        if (vi[ i ] < 0) {
                            throw std::out_of_range(
                                "Value is negative, can't retrieve as "
                                "size_t!" );
                        } else {
                            vst[ i ] = static_cast<size_t>( vi[ i ] );
                        }
                    }
                    return vst;
                }

                virtual bool is_scalar() const {
                    return this->data_form() == parameter_form_t::SCALAR;
                }

                virtual bool is_vector() const {
                    return this->data_form() == parameter_form_t::ARRAY;
                }

                virtual std::string data_form_str() const {
                    return to_string( this->data_form() );
                }

                virtual std::string data_type_str() const {
                    return to_string( this->data_type() );
                }

                virtual parameter_form_t data_form() const { return _form; }

                virtual parameter_type_t data_type() const { return _type; }

                virtual Parameter& clear()                              = 0;
                virtual std::unique_ptr<Parameter> clone() const        = 0;
                virtual std::string json_value_str() const              = 0;
                virtual bool asBool() const                             = 0;
                virtual int asInt() const                               = 0;
                virtual double asDouble() const                         = 0;
                virtual std::string asString() const                    = 0;
                virtual std::vector<bool> asBoolVector() const          = 0;
                virtual std::vector<int> asIntVector() const            = 0;
                virtual std::vector<double> asDoubleVector() const      = 0;
                virtual std::vector<std::string> asStringVector() const = 0;

                // virtual Parameter& set( bool b )                        = 0;
                // virtual Parameter& set( int i )                         = 0;
                // virtual Parameter& set( double d )                      = 0;
                // virtual Parameter& set( const std::string& s )          = 0;
                // virtual Parameter& set( std::vector<bool>& b )          = 0;
                // virtual Parameter& set( std::vector<int>& i )           = 0;
                // virtual Parameter& set( std::vector<double>& d )        = 0;
                // virtual Parameter& set( std::vector<std::string>& s )
                //     = 0;

                template<typename T>
                ScalarParameter<T>& as_scalar() {
                    return dynamic_cast<ScalarParameter<T>&>( *this );
                }

                template<typename T>
                const ScalarParameter<T>& as_scalar() const {
                    return dynamic_cast<const ScalarParameter<T>&>( *this );
                }

                template<typename T>
                VectorParameter<T>& as_vector() {
                    return dynamic_cast<VectorParameter<T>&>( *this );
                }

                template<typename T>
                const VectorParameter<T>& as_vector() const {
                    return dynamic_cast<const VectorParameter<T>&>( *this );
                }

                virtual std::string as_json_str( bool pretty = false,
                                                 int indent  = 1,
                                                 char space  = ' ' ) const {
                    if (pretty) {
                        return this->as_json_str_pretty( indent, space );
                    } else {
                        return this->as_json_str_plain();
                    }
                }

#if HAVE_NLOHMANN_JSON_HPP
                virtual nlohmann::json as_json() const {
                    nlohmann::json json;
                    json[ "comment" ] = this->description();
                    json[ "form" ]    = this->data_form_str();
                    json[ "type" ]    = this->data_type_str();
                    this->add_value_to_json( json );
                    return json;
                }

                virtual std::string as_json_str_plain() const {
                    return this->as_json().dump( -1, 0, true );
                }

                virtual std::string as_json_str_pretty( int indent = 1,
                                                        char space
                                                        = ' ' ) const {
                    return this->as_json().dump( indent, space, true );
                }

                virtual void add_value_to_json( nlohmann::json& json ) const
                    = 0;
#else
                virtual std::string as_json( bool pretty = false,
                                             int indent  = 1,
                                             char space  = ' ' ) const {
                    if (pretty) {
                        return this->as_json_str_pretty( indent, space );
                    } else {
                        return this->as_json_str_plain();
                    }
                }

                virtual std::string as_json_str_plain() const {
                    std::ostringstream json;
                    json << "\"" << this->key() << "\": {"
                         << "\"comment\": \"" << this->description() << "\", "
                         << "\"form\": \"" << this->data_form_str() << "\", "
                         << "\"type\": \"" << this->data_type_str() << "\", "
                         << "\"value\": " << this->json_value_str() << " }";
                    return json.str();
                }

                virtual std::string as_json_str_pretty( int n_indent = 1,
                                                        char space
                                                        = ' ' ) const {
                    std::ostringstream json;
                    std::string indent = "";
                    if (n_indent > 0) {
                        std::ostringstream indentstr;
                        for (int i = 0; i < n_indent; ++i) {
                            indentstr << space;
                        }
                        indent = indentstr.str();
                    }

                    std::string newline = " ";
                    newline[ 0 ]        = '\n';
                    json << indent << "\"" << this->key() << "\": {"
                         << ( n_indent >= 0 ? newline : "" ) << indent
                         << indent << "\"comment\": \"" << this->description()
                         << "\"," << ( n_indent >= 0 ? newline : "" ) << indent
                         << indent << "\"form\": \"" << this->data_form_str()
                         << "\"," << ( n_indent >= 0 ? newline : "" ) << indent
                         << indent << "\"type\": \"" << this->data_type_str()
                         << "\"," << ( n_indent >= 0 ? newline : "" ) << indent
                         << indent << "\"value\": " << this->json_value_str()
                         << ( n_indent >= 0 ? newline : "" ) << "}";
                    return json.str();
                }
#endif

            private:
                parameter_form_t _form;
                parameter_type_t _type;
                std::string _key;
                std::string _description;
        };

        template<typename T>
        class ScalarParameter : public Parameter {
            public:
                ScalarParameter( parameter_type_t paramtype,
                                 const std::string& key, const T& value,
                                 const std::string& description = "" ) :
                    Parameter( parameter_form_t::SCALAR, paramtype, key,
                               description ),
                    _value { value } {}

                virtual ~ScalarParameter() {}

                virtual const T& value() const { return _value; }

                virtual ScalarParameter& set( const T& val ) {
                    _value = val;
                    return *this;
                }

#if HAVE_NLOHMANN_JSON_HPP
                virtual void add_value_to_json(
                    nlohmann::json& json ) const override {
                    json[ "value" ] = this->value();
                }
#endif

                // virtual bool is_scalar() const override { return true; }

                virtual std::vector<bool> asBoolVector() const override {
                    return std::vector<bool> { this->asBool() };
                }

                virtual std::vector<int> asIntVector() const override {
                    return std::vector<int> { this->asInt() };
                }

                virtual std::vector<size_t> asSizeTVector() const override {
                    return std::vector<size_t> { this->asSizeT() };
                }

                virtual std::vector<double> asDoubleVector() const override {
                    return std::vector<double> { this->asDouble() };
                }

                virtual std::vector<std::string> asStringVector()
                    const override {
                    return std::vector<std::string> { this->asString() };
                }

                // virtual parameter_form_t data_form() const override {
                //     return parameter_form_t::SCALAR;
                // }

                // virtual std::string data_form_str() const override {
                //     return "scalar";
                // }

                virtual std::string json_value_str() const override {
                    return this->asString();
                }

            private:
                T _value;
        };

        class IntegerParameter : public ScalarParameter<int> {
            public:
                IntegerParameter( const std::string& key, int value,
                                  const std::string& description = "" ) :
                    ScalarParameter<int>( parameter_type_t::INTEGER, key, value, description ) {}

                virtual ~IntegerParameter() {}

                virtual Parameter& clear() override {
                    this->set( 0 );
                    return *this;
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new IntegerParameter( this->key(), this->value() ) );
                }

                // virtual parameter_type_t data_type() const override {
                //     return parameter_type_t::INTEGER;
                // }

                virtual bool asBool() const override {
                    return ( this->value() != 0 );
                }

                virtual int asInt() const override { return this->value(); }

                virtual double asDouble() const override {
                    return static_cast<int>( this->value() );
                }

                virtual std::string asString() const override {
                    std::ostringstream oss;
                    oss << this->value();
                    return oss.str();
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, int value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, int value ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerParameter( key, value ) );
                }

                // virtual std::string data_type_str() const override {
                //     return "integer";
                // }
        };

        class DoubleParameter : public ScalarParameter<double> {
            public:
                DoubleParameter( const std::string& key, double value,
                                 const std::string& description = "" ) :
                    ScalarParameter<double>( parameter_type_t::FLOAT, key, value, description ) {}

                virtual ~DoubleParameter() {}

                virtual Parameter& clear() override {
                    this->set( 0.0 );
                    return *this;
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new DoubleParameter( this->key(), this->value() ) );
                }

                // virtual parameter_type_t data_type() const override {
                //     return parameter_type_t::FLOAT;
                // }

                virtual bool asBool() const override {
                    return ( this->value() != 0.0 );
                }

                virtual int asInt() const override {
                    return static_cast<int>(
                        this->value()
                        + std::numeric_limits<double>::epsilon() );
                }

                virtual double asDouble() const override {
                    return this->value();
                }

                virtual std::string asString() const override {
                    std::ostringstream oss;
                    oss << this->value();
                    return oss.str();
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, double value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, double value ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleParameter( key, value ) );
                }

                // virtual std::string data_type_str() const override {
                //     return "float";
                // }
        };

        class StringParameter : public ScalarParameter<std::string> {
            public:
                StringParameter( const std::string& key, std::string value,
                                 const std::string& description = "" ) :
                    ScalarParameter<std::string>( parameter_type_t::STRING, key, value, description ) {}

                virtual ~StringParameter() {}

                virtual Parameter& clear() override {
                    this->set( "" );
                    return *this;
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new StringParameter( this->key(), this->value() ) );
                }

                // virtual parameter_type_t data_type() const override {
                //     return parameter_type_t::STRING;
                // }

                virtual bool asBool() const override {
                    return ( this->value() != "" );
                }

                virtual int asInt() const override {
                    return std::stoi( this->value() );
                }

                virtual double asDouble() const override {
                    return std::stod( this->value() );
                }

                virtual std::string asString() const override {
                    return this->value();
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::string& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new StringParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::string& value ) {
                    return std::unique_ptr<Parameter>(
                        new StringParameter( key, value ) );
                }

                // virtual std::string data_type_str() const override {
                //     return "string";
                // }

                virtual std::string json_value_str() const override {
                    return "\"" + this->asString() + "\"";
                }
        };

        class BooleanParameter : public ScalarParameter<bool> {
            public:
                BooleanParameter( const std::string& key, bool value,
                                  const std::string& description = "" ) :
                    ScalarParameter<bool>( parameter_type_t::BOOLEAN, key, value, description ) {}

                virtual ~BooleanParameter() {}

                virtual Parameter& clear() override {
                    this->set( false );
                    return *this;
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new BooleanParameter( this->key(), this->value() ) );
                }

                // virtual parameter_type_t data_type() const override {
                //     return parameter_type_t::BOOLEAN;
                // }

                virtual bool asBool() const override { return this->value(); }

                virtual int asInt() const override {
                    return ( this->value() ? 1 : 0 );
                }

                virtual double asDouble() const override {
                    return ( this->value() ? 1.0 : 0.0 );
                }

                virtual std::string asString() const override {
                    return ( this->value() ? "true" : "false" );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, bool value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, bool value ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanParameter( key, value ) );
                }

                // virtual std::string data_type_str() const override {
                //     return "boolean";
                // }
        };

        template<typename ENUMTYPE>
        class EnumParameter : public ScalarParameter<ENUMTYPE> {
            public:
                EnumParameter( const std::string& key, ENUMTYPE value,
                                 const std::string& description = "" ) :
                    ScalarParameter<ENUMTYPE>( parameter_type_t::ENUM, key, value, description ), _default{value} {}

                EnumParameter( const std::string& key, ENUMTYPE value, ENUMTYPE defaultval,
                                 const std::string& description = "" ) :
                    ScalarParameter<ENUMTYPE>( parameter_type_t::ENUM, key, value, description ), _default{defaultval} {}

                virtual ~EnumParameter() {}

                virtual Parameter& clear() override {
                    this->set( _default );
                    return *this;
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new EnumParameter( this->key(), this->value(), _default, this->description() ) );
                }

                // virtual parameter_type_t data_type() const override {
                //     return parameter_type_t::STRING;
                // }

                virtual bool asBool() const override {
                    return ( this->value() != "" );
                }

                virtual int asInt() const override {
                    return std::stoi( this->value() );
                }

                virtual double asDouble() const override {
                    return std::stod( this->value() );
                }

                virtual std::string asString() const override {
                    return this->value();
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::string& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new StringParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::string& value ) {
                    return std::unique_ptr<Parameter>(
                        new StringParameter( key, value ) );
                }

                // virtual std::string data_type_str() const override {
                //     return "string";
                // }

                virtual std::string json_value_str() const override {
                    return "\"" + this->asString() + "\"";
                }

            private:
                E _default;
        };

        template<typename T>
        class VectorParameter : public Parameter {
            public:
                VectorParameter( parameter_type_t paramtype,
                    const std::string& key,
                                 const std::vector<T>& value,
                                 const std::string& description = "" ) :
                    Parameter( parameter_form_t::ARRAY, paramtype, key, description ), _value { value } {}

                VectorParameter( parameter_type_t paramtype, const std::string& key,
                                 const std::initializer_list<T>& value,
                                 const std::string& description = "" ) :
                    Parameter( parameter_form_t::ARRAY, paramtype, key, description ),
                    _value { std::vector<T> { value } } {}

                virtual ~VectorParameter() {}

                virtual Parameter& clear() override {
                    this->set( {} );
                    return *this;
                }

                // virtual bool is_scalar() const override { return false; }

                virtual const std::vector<T>& value() const { return _value; }

                virtual VectorParameter& set( const std::vector<T>& val ) {
                    _value = val;
                    return *this;
                }

                virtual VectorParameter& set(
                    const std::initializer_list<T>& val ) {
                    _value = std::vector<T> { val };
                    return *this;
                }

#if HAVE_NLOHMANN_JSON_HPP
                virtual void add_value_to_json(
                    nlohmann::json& json ) const override {
                    json[ "value" ] = nlohmann::json( this->value() );
                }
#endif

                virtual bool asBool() const override {
                    return this->value().size() > 0
                             ? this->asBoolVector().at( 0 )
                             : false;
                }

                virtual int asInt() const override {
                    return this->value().size() > 0
                             ? this->asIntVector().at( 0 )
                             : 0;
                }

                virtual double asDouble() const override {
                    return this->value().size() > 0
                             ? this->asDoubleVector().at( 0 )
                             : 0.0;
                }

                virtual std::string asString() const override {
                    return this->value().size() > 0
                             ? this->asStringVector().at( 0 )
                             : "";
                }

                // virtual parameter_form_t data_form() const override {
                //     return parameter_form_t::ARRAY;
                // }

                // virtual std::string data_form_str() const override {
                //     return "array";
                // }

                virtual std::string json_value_str() const override {
                    std::ostringstream oss;
                    oss << "[ ";
                    for (auto it = _value.begin(); it != _value.end(); ++it) {
                        if (it != _value.begin()) {
                            oss << ", ";
                        }
                        oss << *it;
                    }
                    oss << " ]";
                    return oss.str();
                }

            private:
                std::vector<T> _value;
        };

        class DoubleVectorParameter : public VectorParameter<double> {
            public:
                DoubleVectorParameter( const std::string& key,
                                       const std::vector<double>& value,
                                       const std::string& description = "" ) :
                    VectorParameter<double>( parameter_type_t::FLOAT, key, value, description ) {}

                DoubleVectorParameter( const std::string& key,
                                       std::initializer_list<double>& value,
                                       const std::string& description = "" ) :
                    VectorParameter<double>( parameter_type_t::FLOAT, key, value, description ) {}

                virtual ~DoubleVectorParameter() {}

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( this->key(), this->value(),
                                                   this->description() ) );
                }

                // virtual parameter_type_t data_type() const override {
                //     return parameter_type_t::FLOAT;
                // }

                virtual std::vector<bool> asBoolVector() const override {
                    std::vector<bool> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = ( this->value().at( i ) != 0.0 );
                    }
                    return vec;
                }

                virtual std::vector<int> asIntVector() const override {
                    std::vector<int> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = static_cast<int>(
                            this->value().at( i )
                            + std::numeric_limits<double>::epsilon() );
                    }
                    return vec;
                }

                virtual std::vector<double> asDoubleVector() const override {
                    return this->value();
                }

                virtual std::vector<std::string> asStringVector()
                    const override {
                    std::vector<std::string> vec;
                    std::ostringstream oss;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        oss << this->value().at( i );
                        vec[ i ] = oss.str();
                        oss.str( "" );
                    }
                    return vec;
                }

                // virtual std::string data_type_str() const override {
                //     return "float";
                // }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<double>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::vector<double>& value ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( key, value ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<double>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<double>& value ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( key, value ) );
                }
        };

        class IntegerVectorParameter : public VectorParameter<int> {
            public:
                IntegerVectorParameter( const std::string& key,
                                        const std::vector<int>& value,
                                        const std::string& description = "" ) :
                    VectorParameter<int>( parameter_type_t::INTEGER, key, value, description ) {}

                IntegerVectorParameter( const std::string& key,
                                        std::initializer_list<int>& value,
                                        const std::string& description = "" ) :
                    VectorParameter<int>( parameter_type_t::INTEGER, key, value, description ) {}

                virtual ~IntegerVectorParameter() {}

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( this->key(), this->value(),
                                                    this->description() ) );
                }

                // virtual parameter_type_t data_type() const override {
                //     return parameter_type_t::INTEGER;
                // }

                virtual std::vector<bool> asBoolVector() const override {
                    std::vector<bool> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = ( this->value().at( i ) != 0 );
                    }
                    return vec;
                }

                virtual std::vector<int> asIntVector() const override {
                    return this->value();
                }

                virtual std::vector<double> asDoubleVector() const override {
                    std::vector<double> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ]
                            = static_cast<double>( this->value().at( i ) );
                    }
                    return vec;
                }

                virtual std::vector<std::string> asStringVector()
                    const override {
                    std::vector<std::string> vec;
                    std::ostringstream oss;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        oss << this->value().at( i );
                        vec[ i ] = oss.str();
                        oss.str( "" );
                    }
                    return vec;
                }

                // virtual std::string data_type_str() const override {
                //     return "integer";
                // }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<int>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( key, value,
                                                    description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<int>& value ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( key, value ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<int>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( key, value,
                                                    description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<int>& value ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( key, value ) );
                }
        };

        class StringVectorParameter : public VectorParameter<std::string> {
            public:
                StringVectorParameter( const std::string& key,
                                       const std::vector<std::string>& value,
                                       const std::string& description = "" ) :
                    VectorParameter<std::string>( parameter_type_t::STRING, key, value, description ) {}

                StringVectorParameter(
                    const std::string& key,
                    std::initializer_list<std::string>& value,
                    const std::string& description = "" ) :
                    VectorParameter<std::string>( parameter_type_t::STRING, key, value, description ) {}

                virtual ~StringVectorParameter() {}

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( this->key(), this->value(),
                                                   this->description() ) );
                }

                // virtual parameter_type_t data_type() const override {
                //     return parameter_type_t::STRING;
                // }

                virtual std::vector<bool> asBoolVector() const override {
                    std::vector<bool> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = ( this->value().at( i ) != "" );
                    }
                    return vec;
                }

                virtual std::vector<int> asIntVector() const override {
                    std::vector<int> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = std::stoi( this->value().at( i ) );
                    }
                    return vec;
                }

                virtual std::vector<double> asDoubleVector() const override {
                    std::vector<double> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = std::stod( this->value().at( i ) );
                    }
                    return vec;
                }

                virtual std::vector<std::string> asStringVector()
                    const override {
                    return this->value();
                }

                // virtual std::string data_type_str() const override {
                //     return "string";
                // }

                virtual std::string json_value_str() const override {
                    std::ostringstream oss;
                    oss << "[ ";
                    for (auto it = this->value().begin();
                         it != this->value().end(); ++it) {
                        if (it != this->value().begin()) {
                            oss << ", ";
                        }
                        oss << "\"" << *it << "\"";
                    }
                    oss << " ]";
                    return oss.str();
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::vector<std::string>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::vector<std::string>& value ) {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( key, value ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<std::string>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<std::string>& value ) {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( key, value ) );
                }
        };

        class BooleanVectorParameter : public VectorParameter<bool> {
            public:
                BooleanVectorParameter( const std::string& key,
                                        const std::vector<bool>& value,
                                        const std::string& description = "" ) :
                    VectorParameter<bool>( parameter_type_t::BOOLEAN, key, value, description ) {}

                BooleanVectorParameter( const std::string& key,
                                        std::initializer_list<bool>& value,
                                        const std::string& description = "" ) :
                    VectorParameter<bool>(parameter_type_t::BOOLEAN,  key, value, description ) {}

                virtual ~BooleanVectorParameter() {}

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( this->key(), this->value(),
                                                    this->description() ) );
                }

                // virtual parameter_type_t data_type() const override {
                //     return parameter_type_t::BOOLEAN;
                // }

                virtual std::vector<bool> asBoolVector() const override {
                    return this->value();
                }

                virtual std::vector<int> asIntVector() const override {
                    std::vector<int> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = this->value().at( i ) ? 1 : 0;
                    }
                    return vec;
                }

                virtual std::vector<double> asDoubleVector() const override {
                    std::vector<double> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = this->value().at( i ) ? 1.0 : 0.0;
                    }
                    return vec;
                }

                virtual std::vector<std::string> asStringVector()
                    const override {
                    std::vector<std::string> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = this->value().at( i ) ? "true" : "false";
                    }
                    return vec;
                }

                // virtual std::string data_type_str() const override {
                //     return "boolean";
                // }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<bool>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( key, value,
                                                    description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<bool>& value ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( key, value ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<bool>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( key, value,
                                                    description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<bool>& value ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( key, value ) );
                }
        };

        // template<typename T,
        //          typename std::enable_if_t<std::is_scalar<T>::value>>
        // static std::unique_ptr<Parameter> make_parameter(
        //     const std::string& key, const T& val ) {
        //     return std::unique_ptr<Parameter>(
        //         new ScalarParameter<T>( key, val ) );
        // }

        // template<typename T>
        // static std::unique_ptr<Parameter> make_parameter(
        //     const std::string& key, const std::vector<T>& val ) {
        //     return std::unique_ptr<Parameter>(
        //         new VectorParameter<T>( key, val ) );
        // }

        // template<typename T>
        // static std::unique_ptr<Parameter> make_parameter(
        //     const std::string& key, const std::initializer_list<T>& val ) {
        //     return std::unique_ptr<Parameter>(
        //         new VectorParameter<T>( key, val ) );
        // }
    }  // namespace config
}  // namespace NCPA

void swap( NCPA::config::Parameter& a, NCPA::config::Parameter& b ) noexcept {
    using std::swap;
    swap( a._key, b._key );
    swap( a._description, b._description );
}

template<typename T>
void swap( NCPA::config::ScalarParameter<T>& a,
           NCPA::config::ScalarParameter<T>& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::Parameter&>( a ),
            dynamic_cast<NCPA::config::Parameter&>( b ) );
}

template<typename T>
void swap( NCPA::config::VectorParameter<T>& a,
           NCPA::config::VectorParameter<T>& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::Parameter&>( a ),
            dynamic_cast<NCPA::config::Parameter&>( b ) );
}

void swap( NCPA::config::IntegerParameter& a,
           NCPA::config::IntegerParameter& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::ScalarParameter<int>&>( a ),
            dynamic_cast<NCPA::config::ScalarParameter<int>&>( b ) );
}

void swap( NCPA::config::DoubleParameter& a,
           NCPA::config::DoubleParameter& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::ScalarParameter<double>&>( a ),
            dynamic_cast<NCPA::config::ScalarParameter<double>&>( b ) );
}

void swap( NCPA::config::StringParameter& a,
           NCPA::config::StringParameter& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::ScalarParameter<std::string>&>( a ),
            dynamic_cast<NCPA::config::ScalarParameter<std::string>&>( b ) );
}

void swap( NCPA::config::BooleanParameter& a,
           NCPA::config::BooleanParameter& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::ScalarParameter<bool>&>( a ),
            dynamic_cast<NCPA::config::ScalarParameter<bool>&>( b ) );
}

void swap( NCPA::config::IntegerVectorParameter& a,
           NCPA::config::IntegerVectorParameter& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::VectorParameter<int>&>( a ),
            dynamic_cast<NCPA::config::VectorParameter<int>&>( b ) );
}

void swap( NCPA::config::DoubleVectorParameter& a,
           NCPA::config::DoubleVectorParameter& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::VectorParameter<double>&>( a ),
            dynamic_cast<NCPA::config::VectorParameter<double>&>( b ) );
}

void swap( NCPA::config::StringVectorParameter& a,
           NCPA::config::StringVectorParameter& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::VectorParameter<std::string>&>( a ),
            dynamic_cast<NCPA::config::VectorParameter<std::string>&>( b ) );
}

void swap( NCPA::config::BooleanVectorParameter& a,
           NCPA::config::BooleanVectorParameter& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::config::VectorParameter<bool>&>( a ),
            dynamic_cast<NCPA::config::VectorParameter<bool>&>( b ) );
}
