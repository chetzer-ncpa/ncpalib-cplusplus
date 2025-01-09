#pragma once

#include "NCPA/defines.hpp"
#include "NCPA/parameters/declarations.hpp"

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace NCPA {
    namespace params {
        template<typename TESTTYPE>
        class Parameter {
            public:
                Parameter() : _required { true } {}

                Parameter( TESTTYPE initval ) :
                    _required { false }, _value { initval } {}

                virtual ~Parameter() {}

                virtual bool required() const { return _required; }

                virtual bool overridden() const { return _overridden; }

                virtual Parameter<TESTTYPE>& add_validator(
                    const Validator<TESTTYPE>& val ) {
                    _validators.push_back(val.clone() );
                    return *this;
                }

                virtual bool validate() {
                    bool ok = true;
                    for ( auto it = _validators.begin();
                          it != _validators.end(); ++it ) {
                        ok = ( *it )->validate( *this ) && ok;
                    }
                    return ok;
                }

                virtual Parameter<TESTTYPE>& required( bool tf ) {
                    _required = tf;
                    return *this;
                }

                virtual TESTTYPE value() const { return _value; }

                virtual Parameter<TESTTYPE>& set( TESTTYPE val ) {
                    _value = val;
                    return *this;
                }

                DECLARE_METHOD_ENABLED_IF_REAL( as_int, int, TESTTYPE ) const {
                    int val = (int)_value;
                    if ( _value - (double)val > 0.5 ) {
                        val++;
                    }
                    return val;
                }

                DECLARE_METHOD_ENABLED_IF_STRING( as_int, int, TESTTYPE )
                const {
                    return std::stoi( _value );
                }

                DECLARE_METHOD_ENABLED_IF_INTEGRAL( as_int, int, TESTTYPE )
                const {
                    return (int)_value;
                }

                DECLARE_METHOD_ENABLED_IF_REAL( as_double, double, TESTTYPE )
                const {
                    return (double)_value;
                }

                DECLARE_METHOD_ENABLED_IF_INTEGRAL( as_double, double,
                                                    TESTTYPE )
                const {
                    return (double)_value;
                }

                DECLARE_METHOD_ENABLED_IF_STRING( as_double, double, TESTTYPE )
                const {
                    return std::stod( _value );
                }

                std::string as_string() const {
                    return std::to_string( _value );
                }

            protected:
                bool _required = true;
                TESTTYPE _value;
                bool _overridden = false;

                std::vector<std::unique_ptr<Validator<TESTTYPE>>> _validators;
        };
    }  // namespace params
}  // namespace NCPA
