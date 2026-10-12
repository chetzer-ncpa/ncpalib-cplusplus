#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/parameters/BaseParameter.hpp"

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace NCPA {
    namespace config {
        class Argument : public CloneBase<Argument> {
            public:
                virtual std::vector<bool> as_bools() const          = 0;
                virtual std::vector<double> as_doubles() const      = 0;
                virtual std::vector<int> as_ints() const            = 0;
                virtual std::vector<std::string> as_strings() const = 0;
                virtual bool as_bool() const                        = 0;
                virtual double as_double() const                    = 0;
                virtual int as_int() const                          = 0;
                virtual std::string as_string() const               = 0;
                virtual std::string default_as_string() const       = 0;
                virtual void from_string( const std::string& val )  = 0;
                virtual void from_string( const std::vector<std::string>& val )
                    = 0;
                virtual bool has_default_value() const              = 0;
                virtual ValidationStatus validate() const           = 0;

                Argument() = default;

                Argument( const std::string& tag,
                          const std::string& description,
                          value_count_t value_count
                          = value_count_t::ZERO_OR_ONE ) :
                    _tag { tag },
                    _description { description },
                    _value_count { value_count } {}

                Argument( const Argument& other ) :
                    _tag { other._tag },
                    _description { other._description },
                    _value_count { other._value_count },
                    _was_set { other._was_set } {}

                Argument( Argument&& other ) noexcept { swap( *this, other ); }

                virtual ~Argument() = default;

                friend void swap( Argument& a, Argument& b ) noexcept {
                    using std::swap;
                    swap( a._tag, b._tag );
                    swap( a._description, b._description );
                    swap( a._value_count, b._value_count );
                    swap( a._was_set, b._was_set );
                }

                virtual std::string description() const {
                    return _description;
                }

                virtual bool expects_values() const {
                    return this->value_count()
                        != value_count_t::ZERO;
                }

                virtual bool required() const {
                    switch (_value_count) {
                        case value_count_t::ONE:
                        case value_count_t::AT_LEAST_ONE:
                            return true;
                        default:
                            return false;
                    }
                }

                virtual void set_value_count( value_count_t count ) {
                    _value_count = count;
                }

                virtual std::string tag() const { return _tag; }

                virtual value_count_t value_count() const { return _value_count; }

                virtual bool was_set() const { return _was_set; }

            protected:
                std::string _tag         = "";
                std::string _description = "";
                value_count_t _value_count;
                // argument_required_t _required =
                // argument_required_t::OPTIONAL;

                bool _was_set = false;
        };
    }  // namespace config
}  // namespace NCPA
