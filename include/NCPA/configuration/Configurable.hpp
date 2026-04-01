#pragma once

#include "NCPA/configuration/ConfigurationMap.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/Parameter.hpp"
#include "NCPA/configuration/ScalarParameter.hpp"
#include "NCPA/configuration/VectorParameter.hpp"
#include "NCPA/units.hpp"

#include <string>
#include <type_traits>
#include <vector>

namespace NCPA {
    namespace config {
        template<typename KEYTYPE>
        class Configurable {
            public:
                Configurable() {}

                Configurable( const Configurable<KEYTYPE>& other ) :
                    Configurable<KEYTYPE>() {
                    _parameters = other._parameters;
                }

                virtual ~Configurable() {}

                template<typename K = KEYTYPE,
                         typename std::enable_if<
                             std::is_convertible<K, std::string>::value,
                             int>::type = 0>
                double get_as( KEYTYPE key,
                               const NCPA::units::units_ptr_t u ) const {
                    std::string keymain  = key;
                    std::string keyunits = key + "_units";
                    return this->get<NCPA::units::units_ptr_t>( keyunits )
                        ->convert_to( this->parameter( key ).as_double(), *u );
                }

                void convert_parameter( KEYTYPE key,
                                        const NCPA::units::units_ptr_t ufrom,
                                        const NCPA::units::units_ptr_t uto ) {
                    if (this->parameter( key ).is_scalar()) {
                        this->parameter( key ).from<double>( ufrom->convert_to(
                            this->parameter( key ).as_double(), uto ) );
                    } else {
                        this->parameter( key ).from_vector<double>( ufrom->convert_to(
                            this->parameter( key ).as_vector<double>(), uto ) );
                    }
                }

                // converting to scalar
                template<typename TOTYPE                    = double,
                         typename std::enable_if<std::is_scalar<TOTYPE>::value,
                                                 int>::type = 0>
                void convert_parameter( KEYTYPE key,
                                        const NCPA::units::units_ptr_t ufrom,
                                        const NCPA::units::units_ptr_t uto ) {
                    this->add_parameter(
                        key, Parameter<TOTYPE>(
                                 static_cast<TOTYPE>( ufrom->convert_to(
                                     this->parameter( key ).as_double(),
                                     uto ) ) ) );
                }

                // converting to vector
                template<typename TOTYPE = double,
                         typename std::enable_if<
                             !( std::is_scalar<TOTYPE>::value
                                || NCPA::types::is_complex<TOTYPE>::value ),
                             int>::type = 0>
                void convert_parameter( KEYTYPE key,
                                        const NCPA::units::units_ptr_t ufrom,
                                        const NCPA::units::units_ptr_t uto ) {
                    TOTYPE newvec;
                    std::vector<double> as_double
                        = this->parameter( key ).as_double_vector();
                    for (auto it = as_double.begin(); it != as_double.end();
                         ++it) {
                        newvec.push_back( ufrom->convert_to( *it, uto ) );
                    }
                    this->add_parameter<TOTYPE>( key, newvec );
                }

                template<typename TOTYPE = double, typename K = KEYTYPE,
                         typename std::enable_if<
                             std::is_convertible<K, std::string>::value,
                             int>::type = 0>
                void convert_parameter( KEYTYPE key,
                                        const NCPA::units::units_ptr_t uto ) {
                    std::string ukey = _make_units_key( key );
                    this->convert_parameter<TOTYPE>(
                        key, this->get<NCPA::units::units_ptr_t>( ukey ),
                        uto );
                }

                template<typename PARAMTYPE>
                Parameter *add_empty_parameter( KEYTYPE key ) {
                    _parameters[ key ]
                        = param_ptr_t( new Parameter<PARAMTYPE>() );
                    return _parameters[ key ].get();
                }

                template<typename PARAMTYPE, typename K = KEYTYPE,
                         typename std::enable_if<
                             std::is_convertible<K, std::string>::value,
                             int>::type = 0>
                Parameter *add_empty_parameter(
                    KEYTYPE key, const NCPA::units::units_ptr_t units ) {
                    _parameters[ _make_units_key( key ) ]
                        = param_ptr_t( new UnitsParameter( units ) );
                    return this->add_parameter<PARAMTYPE>( key );
                }

                // template<typename PARAMTYPE>
                // Parameter *add_parameter( KEYTYPE key,
                //                                          PARAMTYPE value ) {
                //     _parameters[ key ]
                //         = param_ptr_t( new Parameter<PARAMTYPE>( value ) );
                //     return _parameters[ key ].get();
                // }

                // template<typename PARAMTYPE, typename K = KEYTYPE,
                //          typename std::enable_if<
                //              std::is_convertible<K, std::string>::value,
                //              int>::type = 0>
                // Parameter *add_parameter(
                //     KEYTYPE key, PARAMTYPE value,
                //     const NCPA::units::units_ptr_t units ) {
                //     _parameters[ _make_units_key( key ) ]
                //         = param_ptr_t( new UnitsParameter( units ) );
                //     return this->add_parameter<PARAMTYPE>( key, value );
                // }

                Parameter *add_parameter( KEYTYPE key,
                                          const Parameter *param ) {
                    _parameters[ key ] = param->clone();
                    return _parameters[ key ].get();
                }

                Parameter *add_parameter( KEYTYPE key,
                                          const param_ptr_t param ) {
                    _parameters[ key ] = param->clone();
                    return _parameters[ key ].get();
                }

                Parameter *add_parameter( KEYTYPE key,
                                          const Parameter& param ) {
                    _parameters[ key ] = param.clone();
                    return _parameters[ key ].get();
                }

                template<typename K = KEYTYPE,
                         typename std::enable_if<
                             std::is_convertible<K, std::string>::value,
                             int>::type = 0>
                Parameter *add_parameter(
                    KEYTYPE key, const Parameter *param,
                    const NCPA::units::units_ptr_t units ) {
                    std::string ukey = _make_units_key( key );
                    this->add_parameter( ukey, UnitsParameter( units ) );
                    return this->add_parameter( key, param );
                }

                template<typename K = KEYTYPE,
                         typename std::enable_if<
                             std::is_convertible<K, std::string>::value,
                             int>::type = 0>
                Parameter *add_parameter(
                    KEYTYPE key, const param_ptr_t param,
                    const NCPA::units::units_ptr_t units ) {
                    std::string ukey = _make_units_key( key );
                    this->add_parameter( ukey, UnitsParameter( units ) );
                    return this->add_parameter( key, param );
                }

                template<typename K = KEYTYPE,
                         typename std::enable_if<
                             std::is_convertible<K, std::string>::value,
                             int>::type = 0>
                Parameter *add_parameter(
                    KEYTYPE key, const Parameter& param,
                    const NCPA::units::units_ptr_t units ) {
                    std::string ukey = _make_units_key( key );
                    this->add_parameter( ukey, UnitsParameter( units ) );
                    return this->add_parameter( key, param );
                }

                Parameter& parameter( KEYTYPE key ) {
                    // this->init();
                    return *( _parameters.at( key ).get() );
                }

                const Parameter& parameter( KEYTYPE key ) const {
                    return *( _parameters.at( key ).get() );
                }

                void copy_parameter( const KEYTYPE& key,
                                     const param_ptr_t& ptr ) {
                    this->copy_parameter( key, *ptr );
                }

                void copy_parameter( const KEYTYPE& key,
                                     const Parameter& ptr ) {
                    _parameters[ key ] = ptr.clone();
                }

                // virtual void define_parameters() = 0;

                // virtual void init() {
                //     if (_parameters.size() == 0) {
                //         this->define_parameters();
                //     }
                // }

                void validate_parameters() {
                    // this->init();
                    for (auto it = _parameters.cbegin();
                         it != _parameters.cend(); ++it) {
                        it->second->validate();
                    }
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

                template<typename PARAMTYPE, typename K = KEYTYPE,
                         typename std::enable_if<
                             std::is_convertible<K, std::string>::value,
                             int>::type = 0>
                Parameter *set( KEYTYPE key, PARAMTYPE value ) {
                    // this->init();
                    if (!has_parameter( key )) {
                        return this->add_parameter(
                            key,
                            param_ptr_t( new Parameter<PARAMTYPE>( value ) ) );
                    } else {
                        if (auto sub = dynamic_cast<Parameter<PARAMTYPE> *>(
                                &this->parameter( key ) )) {
                            return sub->set( value );
                        } else {
                            throw std::logic_error( "Can't cast parameter "
                                                    + key
                                                    + " to requested type!" );
                        }
                    }
                }

                template<typename PARAMTYPE, typename K = KEYTYPE,
                         typename std::enable_if<
                             !( std::is_convertible<K, std::string>::value ),
                             int>::type = 0>
                Parameter *set( KEYTYPE key, PARAMTYPE value ) {
                    // this->init();
                    if (!has_parameter( key )) {
                        return this->add_parameter(
                            key,
                            param_ptr_t( new Parameter<PARAMTYPE>( value ) ) );
                    } else {
                        if (auto sub = dynamic_cast<Parameter<PARAMTYPE> *>(
                                &this->parameter( key ) )) {
                            return sub->set( value );
                        } else {
                            throw std::logic_error(
                                "Can't cast parameter to requested type!" );
                        }
                    }
                }

                template<typename PARAMTYPE, typename K = KEYTYPE,
                         typename std::enable_if<
                             std::is_convertible<K, std::string>::value,
                             int>::type = 0>
                const PARAMTYPE& get( KEYTYPE key ) const {
                    if (!has_parameter( key )) {
                        throw std::out_of_range( "Parameter " + key
                                                 + " not found!" );
                    }
                    if (auto sub = dynamic_cast<const Parameter<PARAMTYPE> *>(
                            &this->parameter( key ) )) {
                        return sub->value();
                    } else {
                        throw std::logic_error( "Can't cast parameter " + key
                                                + " to requested type!" );
                    }
                }

                template<typename PARAMTYPE, typename K = KEYTYPE,
                         typename std::enable_if<
                             !( std::is_convertible<K, std::string>::value ),
                             int>::type = 0>
                const PARAMTYPE& get( KEYTYPE key ) const {
                    if (!has_parameter( key )) {
                        throw std::out_of_range( "Parameter not found!" );
                    }
                    if (auto sub = dynamic_cast<const Parameter<PARAMTYPE> *>(
                            &this->parameter( key ) )) {
                        return sub->value();
                    } else {
                        throw std::logic_error(
                            "Can't cast parameter to requested type!" );
                    }
                }

                bool has_parameter( KEYTYPE key ) const {
                    return ( _parameters.find( key ) != _parameters.cend() );
                }

                void copy_parameters_from( const Configurable<KEYTYPE>& other,
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
                }

                void copy_parameters_to( Configurable<KEYTYPE>& other,
                                         bool create_if_missing
                                         = true ) const {
                    other.copy_parameters_from( *this );
                }

                virtual ConfigurationMap<KEYTYPE>& parameters() {
                    // this->init();
                    return _parameters;
                }

                virtual const ConfigurationMap<KEYTYPE>& parameters() const {
                    return _parameters;
                }

            private:
                // std::unordered_map<KEYTYPE,
                //                    std::unique_ptr<Parameter>>
                //     _parameters;
                ConfigurationMap<KEYTYPE> _parameters;

                template<typename K = KEYTYPE,
                         typename std::enable_if<
                             std::is_convertible<K, std::string>::value,
                             int>::type = 0>
                std::string _make_units_key( const std::string& key ) const {
                    std::string ukey = key + "_units";
                    return ukey;
                }
        };
    }  // namespace config
}  // namespace NCPA
