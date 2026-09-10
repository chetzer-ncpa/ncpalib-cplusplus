#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/Validated.hpp"

#include <string>

namespace NCPA {
    namespace config {
        class Argument : public Cloneable<Argument>,
                         public Validated<Argument> {
            public:
                Argument() {}

                virtual ~Argument() {}

                virtual bool expects_value() const { return _expects_value; }

                virtual bool required() const { return _required; }

                virtual Argument& set_expects_value( bool expects ) {
                    _expects_value = expects;
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

                virtual std::string tag() const { return _tag; }

                virtual std::string value_string() const {
                    return ( this->expects_value() ? _value : "" );
                }

            protected:
                std::string _tag;
                std::string _value;
                bool _required;
                bool _expects_value;
                bool _set;

                // std::vector<TypedValidation<std::string>> _validations;
        };
    }  // namespace config
}  // namespace NCPA
