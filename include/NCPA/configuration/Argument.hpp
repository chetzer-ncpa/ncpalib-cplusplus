#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/Validated.hpp"
#include "NCPA/logging.hpp"

#include <string>

namespace NCPA {
    namespace config {
        class Argument : public Cloneable<Argument>,
                         public Validated<Argument> {
            protected:
                std::string _tag;
                bool _required;
                argument_multiplicity_t _multiplicity;
                std::string _default;
                std::vector<std::string> _values;
                bool _set = false;
                // std::vector<TypedValidation<std::string>> _validations;

            public:
                Argument() :
                    Argument( "", false, argument_multiplicity_t::ZERO ) {}

                Argument( const std::string& tag, bool req,
                          const std::string& defaultvalue
                          = NCPA_ARGUMENT_DEFAULT_STRING_VALUE ) :
                    Argument( tag, req,
                              ( req ? argument_multiplicity_t::SINGLE
                                    : argument_multiplicity_t::ZERO ),
                              defaultvalue ) {}

                Argument( const std::string& tag, bool req,
                          argument_multiplicity_t mult,
                          const std::string& defaultvalue
                          = NCPA_ARGUMENT_DEFAULT_STRING_VALUE ) :
                    _tag { tag },
                    _required { req },
                    _multiplicity { mult },
                    _default { defaultvalue },
                    _values {} {}

                Argument( const Argument& other ) : Argument() {
                    _tag          = other._tag;
                    _required     = other._required;
                    _multiplicity = other._multiplicity;
                    _default      = other._default;
                    _values       = other._values;
                    _set          = other._set;
                    // _validations  = other._validations;
                }

                virtual ~Argument() {}

                friend void swap( Argument& a, Argument& b ) noexcept {
                    using std::swap;
                    swap( a._tag, b._tag );
                    swap( a._required, b._required );
                    swap( a._multiplicity, b._multiplicity );
                    swap( a._default, b._default );
                    swap( a._values, b._values );
                    swap( a._set, b._set );
                    // swap( a._validations, b._validations );
                }

                NCPA_CLONE_METHOD( Argument, Argument )

                virtual Argument& add_value_string(
                    const std::string& newval ) {
                    if (this->can_accept_value()) {
                        NCPA_DEBUG << "Adding value " << newval
                                   << " to Argument " << this->tag()
                                   << std::endl;
                        _values.push_back( newval );
                        _set = true;
                    } else {
                        throw TooManyValuesAssigned(
                            "Argument cannot accept more values" );
                    }
                    return *this;
                }

                virtual Argument& add_value_strings(
                    const std::vector<std::string>& newvals ) {
                    for (auto val : newvals) {
                        this->add_value_string( val );
                    }
                    return *this;
                }

                bool can_accept_value() noexcept {
                    switch (this->multiplicity()) {
                        case argument_multiplicity_t::ZERO:
                            return false;
                            break;
                        case argument_multiplicity_t::SINGLE:
                            return ( this->value_count() == 0 );
                            break;
                        default:
                            return true;
                    }
                }

                virtual Argument& clip_values() {
                    switch (this->multiplicity()) {
                        case argument_multiplicity_t::ZERO:
                            _values.clear();
                            break;
                        case argument_multiplicity_t::SINGLE:
                            if (this->value_count() > 1) {
                                _values.resize( 1 );
                            }
                            break;
                    }
                    return *this;
                }

                virtual bool expects_value() const {
                    return ( _multiplicity != argument_multiplicity_t::ZERO );
                }

                virtual bool has_default_value() const noexcept {
                    return ( !_default.empty() );
                }

                virtual Argument& mark_set() {
                    _set = true;
                    return *this;
                }

                // virtual int max_values() const {
                //     switch (this->multiplicity()) {
                //         case argument_multiplicity_t::ZERO:
                //             return 0;
                //         case argument_multiplicity_t::SINGLE:
                //             return 1;
                //         default:
                //             return 1e12;
                //     }
                // }

                virtual argument_multiplicity_t multiplicity() const {
                    return _multiplicity;
                }

                virtual bool required() const { return _required; }

                virtual Argument& set_multiplicity(
                    argument_multiplicity_t mult ) {
                    _multiplicity = mult;
                    if (mult == argument_multiplicity_t::ZERO) {
                        _values.clear();
                        _required = false;
                    }
                    this->clip_values();
                    return *this;
                }

                virtual Argument& set_required( bool req ) {
                    _required = req;
                    return *this;
                }

                virtual Argument& set_tag( const std::string& newtag ) {
                    _tag = newtag;
                    return *this;
                }

                virtual Argument& set_value_string(
                    const std::string& strval ) {
                    if (this->expects_value()) {
                        _values.resize( 1, strval );
                        _set = true;
                    } else {
                        throw TooManyValuesAssigned(
                            "Cannot set a value on a flag argument" );
                    }
                    return *this;
                }

                virtual std::string tag() const noexcept { return _tag; }

                virtual validation_status_t validate() const override {
                    validation_status_t status;
                    if (this->required() && !this->was_set()) {
                        status.result = test_result_t::FAILED;
                        status.message
                            = this->tag() + ": required but not set";
                    } else {
                        status = merge( status, this->validate_value(
                                                    this->value_string() ) );
                    }
                    return status;
                }

                virtual size_t value_count() const noexcept {
                    return _values.size();
                }

                virtual std::string value_string( size_t n = 0 ) const {
                    if (this->multiplicity()
                        == argument_multiplicity_t::ZERO) {
                        throw ValueRequestedFromFlag(
                            "Value requested from flag option "
                            + this->tag() );
                    } else if (_values.empty()) {
                        if (this->has_default_value()) {
                            return _default;
                        } else {
                            throw NoDefaultValue(
                                "Value requested from option " + this->tag()
                                + " but it has not been set and no default "
                                  "value has been defined" );
                        }
                    } else if (this->value_count() <= n) {
                        throw ValueIndexOutOfRange(
                            "Requested value index " + std::to_string( n )
                            + " out of range for option " + this->tag() );
                    }
                    return _values.at( n );
                }

                virtual const std::vector<std::string>& value_strings() const {
                    return _values;
                }

                virtual bool was_set() const noexcept { return _set; }
        };
    }  // namespace config
}  // namespace NCPA
