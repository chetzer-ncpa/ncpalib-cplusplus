#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/arguments/Argument.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/parameters/BaseParameter.hpp"

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace NCPA {
    namespace config {
        template<typename BASETYPE>
        class TypedArgument
            : public Cloneable<TypedArgument<BASETYPE>, Argument> {
                using cloneable_t
                    = Cloneable<TypedArgument<BASETYPE>, Argument>;

            public:
                TypedArgument() = default;

                TypedArgument( const std::string& tag,
                               const std::string& description ) :
                    cloneable_t( tag, description ),
                    _underlying { tag },
                    _default {} {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               value_count_t count ) :
                    cloneable_t( tag, description, count ),
                    _underlying { tag },
                    _default {} {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               value_count_t count,
                               const BASETYPE& defaultval ) :
                    cloneable_t( tag, description, count ),
                    _underlying { tag, defaultval },
                    _default { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               value_count_t count,
                               const std::vector<BASETYPE>& defaultval ) :
                    cloneable_t( tag, description, count ),
                    _underlying { tag, defaultval },
                    _default { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               const BASETYPE& defaultval ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval },
                    _default { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               const std::vector<BASETYPE>& defaultval ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval },
                    _default { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               const TypedValidation<BASETYPE>& v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, v },
                    _default {} {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               const BASETYPE& defaultval,
                               const TypedValidation<BASETYPE>& v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval, v },
                    _default { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               const std::vector<BASETYPE>& defaultval,
                               const TypedValidation<BASETYPE>& v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval, v },
                    _default { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               const TypedValidation<BASETYPE> *v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, v },
                    _default {} {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               value_count_t required,
                               const TypedValidation<BASETYPE> *v ) :
                    cloneable_t( tag, description, required ),
                    _underlying { tag, v },
                    _default {} {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               const BASETYPE& defaultval,
                               const TypedValidation<BASETYPE> *v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval, v },
                    _default { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               const std::vector<BASETYPE>& defaultval,
                               const TypedValidation<BASETYPE> *v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval, v },
                    _default { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               const BASETYPE& defaultval,
                               value_count_t required,
                               const TypedValidation<BASETYPE> *v ) :
                    cloneable_t( tag, description, required ),
                    _underlying { tag, defaultval, v },
                    _default { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& description,
                               value_count_t required,
                               const std::vector<BASETYPE>& defaultval,
                               const TypedValidation<BASETYPE> *v ) :
                    cloneable_t( tag, description, required ),
                    _underlying { tag, defaultval, v },
                    _default { defaultval } {}

                TypedArgument(
                    const std::string& tag, const std::string& description,
                    const std::unique_ptr<TypedValidation<BASETYPE>>& v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, v },
                    _default {} {}

                TypedArgument(
                    const std::string& tag, const std::string& description,
                    value_count_t required,
                    const std::unique_ptr<TypedValidation<BASETYPE>>& v ) :
                    cloneable_t( tag, description, required ),
                    _underlying { tag, v },
                    _default {} {}

                TypedArgument(
                    const std::string& tag, const std::string& description,
                    const BASETYPE& defaultval,
                    const std::unique_ptr<TypedValidation<BASETYPE>>& v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval, v },
                    _default { defaultval } {}

                TypedArgument(
                    const std::string& tag, const std::string& description,
                    const std::vector<BASETYPE>& defaultval,
                    const std::unique_ptr<TypedValidation<BASETYPE>>& v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval, v },
                    _default { defaultval } {}

                TypedArgument(
                    const std::string& tag, const std::string& description,
                    std::initializer_list<TypedValidation<BASETYPE>> v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, v },
                    _default {} {}

                TypedArgument(
                    const std::string& tag, const std::string& description,
                    value_count_t required,
                    std::initializer_list<TypedValidation<BASETYPE>> v ) :
                    cloneable_t( tag, description, required ),
                    _underlying { tag, v },
                    _default {} {}

                TypedArgument(
                    const std::string& tag, const std::string& description,
                    const BASETYPE& defaultval,
                    std::initializer_list<TypedValidation<BASETYPE>> v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval, v },
                    _default { 1, defaultval } {}

                TypedArgument(
                    const std::string& tag, const std::string& description,
                    const std::vector<BASETYPE>& defaultval,
                    std::initializer_list<TypedValidation<BASETYPE>> v ) :
                    cloneable_t( tag, description ),
                    _underlying { tag, defaultval, v },
                    _default { defaultval } {}

                TypedArgument( const TypedArgument& other ) :
                    cloneable_t( other ),
                    _underlying { other._underlying },
                    _default { other._default } {}

                TypedArgument( TypedArgument&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~TypedArgument() = default;

                friend void swap( TypedArgument<BASETYPE>& a,
                                  TypedArgument<BASETYPE>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<cloneable_t&>( a ),
                          static_cast<cloneable_t&>( b ) );
                    swap( a._underlying, b._underlying );
                    swap( a._default, b._default );
                }

                virtual std::vector<bool> as_bools() const override {
                    std::vector<bool> out( _underlying.size() );
                    for (size_t i = 0; i < _underlying.size(); ++i) {
                        out[ i ] = _underlying.as_bool( i );
                    }
                    return out;
                }

                virtual std::vector<double> as_doubles() const override {
                    std::vector<double> out( _underlying.size() );
                    for (size_t i = 0; i < _underlying.size(); ++i) {
                        out[ i ] = _underlying.as_double( i );
                    }
                    return out;
                }

                virtual std::vector<int> as_ints() const override {
                    std::vector<int> out( _underlying.size() );
                    for (size_t i = 0; i < _underlying.size(); ++i) {
                        out[ i ] = _underlying.as_int( i );
                    }
                    return out;
                }

                virtual std::vector<std::string> as_strings() const override {
                    std::vector<std::string> out( _underlying.size() );
                    for (size_t i = 0; i < _underlying.size(); ++i) {
                        out[ i ] = _underlying.as_string( i );
                    }
                    return out;
                }

                virtual bool as_bool() const override {
                    return _underlying.as_bool();
                }

                virtual double as_double() const override {
                    return _underlying.as_double();
                }

                virtual int as_int() const override {
                    return _underlying.as_int();
                }

                virtual std::string as_string() const override {
                    return _underlying.as_string();
                }

                virtual std::string default_as_string() const override {
                    VectorParameter<BASETYPE> defparam( _underlying );
                    defparam.set_value( _default );
                    return defparam.as_string();
                }

                virtual void from_string( const std::string& val ) override {
                    _underlying.from_string( val );
                }

                virtual void from_string(
                    const std::vector<std::string>& val ) override {
                    _underlying.from_string( val );
                }

                virtual bool has_default_value() const override {
                    return !_default.empty();
                }

                virtual void set_value( const BASETYPE& val ) {
                    _underlying.set_value( val );
                    _was_set = true;
                }

                virtual void set_value( const std::vector<BASETYPE>& val ) {
                    _underlying.set_value( val );
                    _was_set = true;
                }

                virtual TypedArgument<BASETYPE>& add_test(
                    TypedValidation<BASETYPE>& v ) {
                    _underlying.add_test( v );
                    return *this;
                }

                virtual TypedArgument<BASETYPE>& add_test(
                    TypedValidation<BASETYPE>&& v ) {
                    _underlying.add_test( std::move( v ) );
                    return *this;
                }

                virtual ValidationStatus validate() const override {
                    ValidationStatus status;
                    if (this->required() && !this->was_set()) {
                        status.result = test_result_t::FAILED;
                        status.messages.push_back( "Required parameter " + _tag
                                                   + " was not set." );
                        return status;
                    } else {
                        return _underlying.validate();
                    }
                }

            protected:
                VectorParameter<BASETYPE> _underlying;
                std::vector<BASETYPE> _default;

                using Argument::_tag;
                using Argument::_was_set;
        };
    }  // namespace config
}  // namespace NCPA
