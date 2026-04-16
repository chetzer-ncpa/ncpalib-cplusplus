#pragma once

#include "NCPA/configuration/BaseParameter.hpp"
#include "NCPA/configuration/ConfigurationMap.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/Parameter.hpp"
#include "NCPA/configuration/ScalarParameter.hpp"
#include "NCPA/configuration/ScalarParameterWithUnits.hpp"
#include "NCPA/configuration/VectorParameter.hpp"
#include "NCPA/configuration/VectorParameterWithUnits.hpp"
#include "NCPA/units.hpp"

#include <algorithm>
#include <string>
#include <type_traits>
#include <vector>

namespace NCPA {
    namespace config {
        using NCPA::units::units_ptr_t;

        // using UnitsParameter = ScalarParameter<units_ptr_t>;

        template<typename KEYTYPE>
        class Configurable {
            public:
                Configurable() {}

                Configurable( const Configurable<KEYTYPE>& other ) :
                    Configurable<KEYTYPE>() {
                    _parameters = other._parameters;
                }

                virtual ~Configurable() {}

                friend void swap<>( Configurable<KEYTYPE>& a,
                                    Configurable<KEYTYPE>& b ) noexcept;

                virtual bool parameter_has_units( KEYTYPE key ) const {
                    return this->parameter( key ).has_units();
                }

                // template<typename K = KEYTYPE,
                //          typename std::enable_if<
                //              std::is_convertible<K, std::string>::value,
                //              int>::type = 0>
                virtual double get_as(
                    KEYTYPE key, const NCPA::units::units_ptr_t u ) const {
                    return this->parameter( key ).get_units()->convert_to(
                        this->parameter( key ).as_double(), *u );
                    // return this->parameter(keymain).get_units()->
                    // return this->get<NCPA::units::units_ptr_t>( keyunits )
                    //     ->convert_to( this->parameter( key ).as_double(), *u
                    //     );
                }

                // virtual Configurable<KEYTYPE>& convert_parameter(
                //     KEYTYPE key, const NCPA::units::units_ptr_t ufrom,
                //     const NCPA::units::units_ptr_t uto ) {
                //     if (this->parameter(key).has_units() && ) {

                //     }
                //     if (this->parameter( key ).is_scalar()) {
                //         this->parameter( key ).from_double(
                //         ufrom->convert_to(
                //             this->parameter( key ).as_double(), uto ) );
                //     } else {
                //         this->parameter( key ).from_double(
                //         ufrom->convert_to(
                //             this->parameter( key ).as_double_vector(), uto )
                //             );
                //     }
                //     return *this;
                // }

                // convert assuming units are stored as "<key>_units"
                // template<typename K = KEYTYPE,
                //          typename std::enable_if<
                //              std::is_convertible<K, std::string>::value,
                //              int>::type = 0>
                BaseParameter& convert_parameter(
                    KEYTYPE key, const NCPA::units::units_ptr_t uto ) {
                    this->parameter( key ).convert_units( uto );
                    return this->parameter( key );
                    // std::string ukey = _make_units_key( key );
                    // return this->convert_parameter(
                    //     key, this->get<NCPA::units::units_ptr_t>( ukey ),
                    //     uto );
                }

                template<typename PARAMTYPE>
                BaseParameter *add_empty_parameter( KEYTYPE key ) {
                    _parameters[ key ] = Parameter::scalar<PARAMTYPE>();
                    // = param_ptr_t( new ScalarParameter<PARAMTYPE>() );
                    return _parameters[ key ].get();
                }

                template<typename PARAMTYPE>
                // , typename K = KEYTYPE,
                //          typename std::enable_if<
                //              std::is_convertible<K, std::string>::value,
                //              int>::type = 0>
                BaseParameter *add_empty_parameter(
                    KEYTYPE key, const NCPA::units::units_ptr_t units ) {
                    // _parameters[ _make_units_key( key ) ]
                    //     = param_ptr_t( new UnitsParameter( units ) );
                    // return this->add_empty_parameter<PARAMTYPE>( key );
                    _parameters[ key ] = Parameter::scalar<PARAMTYPE>( units );
                    return _parameters[ key ].get();
                }

                template<typename PARAMTYPE>
                BaseParameter *add_empty_vector_parameter( KEYTYPE key ) {
                    // _parameters[ key ]
                    //     = param_ptr_t( new VectorParameter<PARAMTYPE>() );
                    // return _parameters[ key ].get();
                    _parameters[ key ] = Parameter::vector<PARAMTYPE>();
                    return _parameters[ key ].get();
                }

                template<typename PARAMTYPE>
                // , typename K = KEYTYPE,
                //          typename std::enable_if<
                //              std::is_convertible<K, std::string>::value,
                //              int>::type = 0>
                BaseParameter *add_empty_vector_parameter(
                    KEYTYPE key, const NCPA::units::units_ptr_t units ) {
                    // _parameters[ _make_units_key( key ) ]
                    //     = param_ptr_t( new UnitsParameter( units ) );
                    // return this->add_empty_vector_parameter<PARAMTYPE>( key
                    // );
                    _parameters[ key ] = Parameter::vector<PARAMTYPE>( units );
                    return _parameters[ key ].get();
                }

                BaseParameter *add_parameter( KEYTYPE key,
                                              const BaseParameter *param ) {
                    _parameters[ key ] = param->clone();
                    return _parameters[ key ].get();
                }

                BaseParameter *add_parameter( KEYTYPE key,
                                              const param_ptr_t param ) {
                    // NCPA_DEBUG << "Parameter passed to " << key << ": " << param->as_string() << std::endl;
                    _parameters[ key ] = param->clone();
                    // NCPA_DEBUG << "Cloned parameter " << key << ": " << _parameters[ key ]->as_string() << std::endl;
                    return _parameters[ key ].get();
                }

                BaseParameter *add_parameter( KEYTYPE key,
                                              const BaseParameter& param ) {
                    _parameters[ key ] = param.clone();
                    return _parameters[ key ].get();
                }

                // BaseParameter *add_parameter(
                //     KEYTYPE key, const BaseParameter *param,
                //     const NCPA::units::units_ptr_t units ) {
                //     if (param->is_scalar)


                //     std::string ukey = _make_units_key( key );
                //     this->add_parameter( ukey, UnitsParameter( units ) );
                //     return this->add_parameter( key, param );
                // }

                // BaseParameter *add_parameter(
                //     KEYTYPE key, const param_ptr_t param,
                //     const NCPA::units::units_ptr_t units ) {
                //     std::string ukey = _make_units_key( key );
                //     this->add_parameter( ukey, UnitsParameter( units ) );
                //     return this->add_parameter( key, param );
                // }

                // BaseParameter *add_parameter(
                //     KEYTYPE key, const BaseParameter& param,
                //     const NCPA::units::units_ptr_t units ) {
                //     std::string ukey = _make_units_key( key );
                //     this->add_parameter( ukey, UnitsParameter( units ) );
                //     return this->add_parameter( key, param );
                // }

                BaseParameter& parameter( KEYTYPE key ) {
                    // this->init();
                    return *( _parameters.at( key ).get() );
                }

                const BaseParameter& parameter( KEYTYPE key ) const {
                    return *( _parameters.at( key ) );
                }

                Configurable<KEYTYPE>& copy_parameter(
                    const KEYTYPE& key, const param_ptr_t& ptr ) {
                    return this->copy_parameter( key, *ptr );
                }

                Configurable<KEYTYPE>& copy_parameter(
                    const KEYTYPE& key, const BaseParameter& ptr ) {
                    _parameters[ key ] = ptr.clone();
                    return *this;
                }

                Configurable<KEYTYPE>& validate_parameters() {
                    // this->init();
                    for (auto it = _parameters.cbegin();
                         it != _parameters.cend(); ++it) {
                        it->second->validate();
                    }
                    return *this;
                }

                std::string validation_report() const {
                    std::ostringstream oss;
                    validation_report( oss, false );
                    return oss.str();
                }

                std::ostream& validation_report( std::ostream& os,
                                                 bool newline = true ) const {
                    bool firsttime = true;
                    for (auto it = _parameters.begin();
                         it != _parameters.end(); ++it) {
                        if (firsttime) {
                            firsttime = false;
                        } else {
                            os << std::endl;
                        }
                        os << it->first << ":" << std::endl;
                        it->second->validation_report( os, false, "  " );
                    }
                    if (newline) {
                        os << std::endl;
                    }
                    return os;
                }

                bool passed() const {
                    bool pass = true;
                    for (auto it = _parameters.cbegin();
                         it != _parameters.cend(); ++it) {
                        pass = pass && it->second->passed();
                    }
                    return pass;
                }

                bool failed() const { return !( this->passed() ); }

                std::vector<Parameter *> invalid() const {
                    std::vector<Parameter *> inv;
                    for (auto it = _parameters.cbegin();
                         it != _parameters.cend(); ++it) {
                        if (it->second->failed()) {
                            inv.push_back( it->second.get() );
                        }
                    }
                    return inv;
                }

                // get() methods: specialize template by requested output type
                // case 1: floating point
                // criteria: is_floating_point == true
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                PARAMTYPE get( KEYTYPE key, size_t n = 0 ) const {
                    return static_cast<PARAMTYPE>(
                        this->parameter( key ).as_double( n ) );
                }

                // case 2: signed integer
                // criteria: is_integral == true,
                //           is_same<bool> == false,
                //           is_signed == true
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_signed<PARAMTYPE>::value ),
                             int>::type = 0>
                PARAMTYPE get( KEYTYPE key, size_t n = 0 ) const {
                    return static_cast<PARAMTYPE>(
                        this->parameter( key ).as_int( n ) );
                }

                // case 3: unsigned integer
                // criteria: is_integral == true,
                //           is_same<bool> == false,
                //           is_unsigned == true
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_integral<PARAMTYPE>::value
                               && !( std::is_same<PARAMTYPE, bool>::value )
                               && std::is_unsigned<PARAMTYPE>::value ),
                             int>::type = 0>
                PARAMTYPE get( KEYTYPE key, size_t n = 0 ) const {
                    return static_cast<PARAMTYPE>(
                        this->parameter( key ).as_unsigned_int( n ) );
                }

                // case 4: boolean
                // criteria: is_same<bool> == true
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_same<PARAMTYPE, bool>::value, int>::type
                         = 0>
                PARAMTYPE get( KEYTYPE key, size_t n = 0 ) const {
                    return this->parameter( key ).as_bool( n );
                }

                // case 5: string
                // criteria: is_arithmetic == false,
                //           is_convertible<string> == true
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<PARAMTYPE>::value )
                               && std::is_convertible<PARAMTYPE,
                                                      std::string>::value ),
                             int>::type = 0>
                PARAMTYPE get( KEYTYPE key, size_t n = 0 ) const {
                    std::string s = this->parameter( key ).as_string( n );
                    return s;
                }

                // case 6: complex
                // criteria: is_scalar == false,
                //           is_complex (custom) == true
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_scalar<PARAMTYPE>::value )
                               && NCPA::types::is_complex<PARAMTYPE>::value ),
                             int>::type = 0>
                PARAMTYPE get( KEYTYPE key, size_t n = 0 ) const {
                    PARAMTYPE output;
                    output.real( static_cast<typename PARAMTYPE::value_type>(
                        this->parameter( key ).as_complex( n ).real() ) );
                    output.imag( static_cast<typename PARAMTYPE::value_type>(
                        this->parameter( key ).as_complex( n ).imag() ) );
                    return output;
                }

                // case 7: everything else
                // we cast it to the appropriate TypedParameter and use get()
                // criteria: is_arithmetic == false,
                //           is_convertible<string> == false,
                //           ( is_scalar == false && is_complex == true)
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<PARAMTYPE>::value
                                  || std::is_convertible<PARAMTYPE,
                                                         std::string>::value
                                  || ( !( std::is_scalar<PARAMTYPE>::value )
                                       && NCPA::types::is_complex<
                                           PARAMTYPE>::value ) ) ),
                             int>::type = 0>
                PARAMTYPE get( KEYTYPE key, size_t n = 0 ) const {
                    return dynamic_cast<const TypedParameter<PARAMTYPE> *>(
                               &this->parameter( key ) )
                        ->get( n );
                }

                template<typename PARAMTYPE>
                std::vector<PARAMTYPE> get_vector( KEYTYPE key ) const {
                    size_t n = this->parameter( key ).size();
                    // NCPA_DEBUG << "Vector size is " << n << std::endl;
                    std::vector<PARAMTYPE> vec( n );
                    for (size_t i = 0; i < n; ++i) {
                        // NCPA_DEBUG << "Getting " << key << "[" << n << "]"
                        // << std::endl;
                        vec[ i ] = this->get<PARAMTYPE>( key, i );
                    }
                    return vec;
                }

                template<typename PARAMTYPE = double,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                NCPA::units::ScalarWithUnits<PARAMTYPE> get_with_units(
                    KEYTYPE key, size_t n = 0 ) const {
                    return ScalarWithUnits<PARAMTYPE>(
                        static_cast<PARAMTYPE>(
                            this->parameter( key ).as_double( n ) ),
                        this->parameter( key ).get_units() );
                }

                // template<typename PARAMTYPE = double,
                //          typename std::enable_if<
                //              std::is_convertible<KEYTYPE,
                //              std::string>::value, int>::type = 0>
                // NCPA::units::ScalarWithUnits<PARAMTYPE> get_with_units(
                //     KEYTYPE key, size_t n = 0 ) const {
                //     return NCPA::units::ScalarWithUnits<PARAMTYPE>(
                //         this->get<PARAMTYPE>( key, n ),
                //         this->get<NCPA::units::units_ptr_t>( key
                //                                              + "_units" ) );
                // }

                template<typename PARAMTYPE = double,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             int>::type = 0>
                NCPA::units::VectorWithUnits<PARAMTYPE> get_vector_with_units(
                    KEYTYPE key ) const {
                    return NCPA::units::VectorWithUnits<PARAMTYPE>(
                        this->get_vector<PARAMTYPE>( key ),
                        this->parameter( key ).get_units() );
                }

                // case 1: scalar or complex, not stringifyable
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( std::is_scalar<PARAMTYPE>::value
                               || NCPA::types::is_complex<PARAMTYPE>::value ),
                             int>::type = 0>
                Configurable<KEYTYPE>& set( KEYTYPE key, PARAMTYPE val ) {
                    this->add_parameter( key,
                                         Parameter::scalar<PARAMTYPE>( val ) );
                    // param_ptr_t( new ScalarParameter<PARAMTYPE>( val ) ) );
                    return *this;
                }

                // case 2: not scalar, is a string
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_scalar<PARAMTYPE>::value )
                               && std::is_convertible<PARAMTYPE,
                                                      std::string>::value ),
                             int>::type = 0>
                Configurable<KEYTYPE>& set( KEYTYPE key, PARAMTYPE val ) {
                    this->add_parameter(
                        key, Parameter::scalar<std::string>( val ) );
                    // key, param_ptr_t(
                    //          new ScalarParameter<std::string>( val ) ) );
                    return *this;
                }

                // case 3 is neither scalar, complex, nor a string, therefore
                // vector
                template<typename PARAMTYPE,
                         typename std::enable_if<
                             ( !( std::is_scalar<PARAMTYPE>::value
                                  || NCPA::types::is_complex<PARAMTYPE>::value
                                  || std::is_convertible<
                                      PARAMTYPE, std::string>::value ) ),
                             int>::type = 0>
                Configurable<KEYTYPE>& set( KEYTYPE key, PARAMTYPE val ) {
                    this->add_parameter(
                        key, Parameter::vector<typename PARAMTYPE::value_type>(
                                 val ) );
                    // key, param_ptr_t( new VectorParameter<
                    //                   typename PARAMTYPE::value_type>(
                    //          val ) ) );
                    return *this;
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<PARAMTYPE>::value,
                             //  std::is_convertible<KEYTYPE,
                             //  std::string>::value,
                             int>::type = 0>
                Configurable<KEYTYPE>& set( KEYTYPE key, PARAMTYPE val,
                                            NCPA::units::units_ptr_t units ) {
                    this->add_parameter(
                        key, Parameter::scalar<PARAMTYPE>( val, units ) );
                    return *this;
                }

                template<typename PARAMTYPE,
                         typename std::enable_if<

                             ( !( std::is_scalar<PARAMTYPE>::value
                                  || NCPA::types::is_complex<PARAMTYPE>::value
                                  || std::is_convertible<PARAMTYPE,
                                                         std::string>::value )
                               && std::is_floating_point<
                                   typename PARAMTYPE::value_type>::value ),
                             int>::type = 0>
                Configurable<KEYTYPE>& set( KEYTYPE key, PARAMTYPE val,
                                            NCPA::units::units_ptr_t units ) {
                    // NCPA_DEBUG << "Setting new " << key << " vector parameter.  Passed vector is size " << val.size() << std::endl;
                    this->add_parameter(
                        key, Parameter::vector<typename PARAMTYPE::value_type>(
                                 val, units ) );
                    return *this;
                }

                bool has_parameter( KEYTYPE key ) const {
                    return ( _parameters.find( key ) != _parameters.cend() );
                }

                Configurable<KEYTYPE>& copy_parameters_from(
                    const Configurable<KEYTYPE>& other,
                    bool create_if_missing = true ) {
                    // this->init();
                    for (auto it = other.parameters().cbegin();
                         it != other.parameters().cend(); ++it) {
                        if (create_if_missing
                            || this->has_parameter( it->first )) {
                            std::cout << "Setting parameter " << it->first
                                      << std::endl;
                            this->copy_parameter( it->first,
                                                  it->second->clone() );
                        }
                    }
                    return *this;
                }

                const Configurable<KEYTYPE>& copy_parameters_to(
                    Configurable<KEYTYPE>& other,
                    bool create_if_missing = true ) const {
                    other.copy_parameters_from( *this, create_if_missing );
                    return *this;
                }

                virtual ConfigurationMap<KEYTYPE>& parameters() {
                    // this->init();
                    return _parameters;
                }

                virtual const ConfigurationMap<KEYTYPE>& parameters() const {
                    return _parameters;
                }

                virtual std::vector<std::string> parameter_keys(
                    bool sorted = false ) const {
                    std::vector<std::string> keys;
                    for (auto it = _parameters.cbegin();
                         it != _parameters.cend(); ++it) {
                        keys.push_back( it->first );
                    }
                    if (sorted) {
                        std::sort( keys.begin(), keys.end() );
                    }
                    return keys;
                }

                virtual void print_configuration( std::ostream& os ) const {
                    size_t longest = 0;
                    std::vector<std::string> keys
                        = this->parameter_keys( true );
                    for (std::string key : keys) {
                        longest = std::max( key.size(), longest );
                    }
                    for (std::string key : keys) {
                        os.width( longest );
                        os << key << " : ";
                        os.width( 0 );
                        try {
                            os << parameter( key ).as_string();
                        } catch (std::out_of_range& oor) {
                            os << "<no string>";
                        }
                        os << std::endl;
                    }
                }

            private:
                ConfigurationMap<KEYTYPE> _parameters;

                // template<typename K = KEYTYPE,
                //          typename std::enable_if<
                //              std::is_convertible<K, std::string>::value,
                //              int>::type = 0>
                // std::string _make_units_key( const std::string& key ) const
                // {
                //     std::string ukey = key + "_units";
                //     return ukey;
                // }
        };
    }  // namespace config
}  // namespace NCPA

template<typename KEYTYPE>
void swap( NCPA::config::Configurable<KEYTYPE>& a,
           NCPA::config::Configurable<KEYTYPE>& b ) noexcept {
    using std::swap;
    swap( a._parameters, b._parameters );
}
