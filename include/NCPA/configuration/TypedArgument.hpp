#pragma once

#include "NCPA/configuration/Argument.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        template<typename INTYPE>
        void swap( TypedArgument<INTYPE>& a,
                   TypedArgument<INTYPE>& b ) noexcept;

        template<typename INTYPE>
        class TypedArgument : public Argument {
            public:
                TypedArgument() : Argument() {}

                // TypedArgument( const std::string& tag ) : Argument( tag ) {}

                TypedArgument( const std::string& tag,
                               const std::string& help_text ) :
                    Argument( tag, help_text ) {}

                // TypedArgument( const std::string& tag,
                //                const INTYPE& defaultval ) :
                //     Argument( tag ),
                //     _default { defaultval },
                //     _value { defaultval } {}

                TypedArgument( const std::string& tag,
                               const std::string& help_text,
                               const INTYPE& defaultval ) :
                    Argument( tag, help_text ),
                    _default { defaultval },
                    _value { defaultval }, _has_default{true} {}

                TypedArgument(
                    const std::string& tag, const std::string& help_text,
                    const mapping_ptr_t<INTYPE, std::string>& mapping ) :
                    Argument( tag, help_text ) {
                    _mappings.push_back( mapping.clone() );
                }

                TypedArgument(
                    const std::string& tag, const std::string& help_text,
                    const INTYPE& defaultval,
                    const mapping_ptr_t<INTYPE, std::string>& mapping ) :
                    Argument( tag, help_text ),
                    _default { defaultval },
                    _value { defaultval }, _has_default{true} {
                    _mappings.push_back( mapping.clone() );
                }

                TypedArgument( const std::string& tag,const std::string& help_text,
                               const Mapping<INTYPE, std::string>& mapping )
                               :
                    Argument( tag, help_text ) {
                    _mappings.push_back( mapping.clone() );
                }

                TypedArgument( const std::string& tag,const std::string& help_text,
                               const INTYPE& defaultval,
                               const Mapping<INTYPE, std::string>& mapping )
                               :
                    Argument( tag, help_text ),
                    _default { defaultval },
                    _value { defaultval }, _has_default{true} {
                    _mappings.push_back( mapping.clone() );
                }

                TypedArgument( const TypedArgument<INTYPE>& other ) :
                    Argument( other ),
                    _default { other._default },
                    _value { other._value }, _has_default{true},
                    _was_set { other._was_set } {
                    _mappings.reserve( other._mappings.size() );
                    for (auto it = other._mappings.begin();
                         it != other._mappings.end(); ++it) {
                        _mappings.push_back( ( *it )->clone() );
                    }
                }

                TypedArgument( TypedArgument<INTYPE>&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~TypedArgument() {}

                TypedArgument<INTYPE>& operator=(
                    TypedArgument<INTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap<>( TypedArgument<INTYPE>& a,
                                    TypedArgument<INTYPE>& b ) noexcept;

                virtual void apply() override {
                    for (auto it = _mappings.begin(); it != _mappings.end();
                         ++it) {
                        ( *it )->apply( _value );
                    }
                }

                virtual std::unique_ptr<Argument> clone() const override {
                    return std::unique_ptr<Argument>(
                        new TypedArgument<INTYPE>( *this ) );
                }

                virtual parse_status_t parse(
                    std::vector<std::string>& args ) override {
                    for (auto it = args.begin(); it != args.end(); ++it) {
                        std::string item
                            = it->substr( it->find_first_not_of( "-" ) );
                        if (item.size() > 0 && item == this->tag()) {
                            if (this->expects_value()) {
                                if (( it + 1 ) == args.end()) {
                                    throw std::out_of_range(
                                        "Argument " + item
                                        + " expects a value but none was "
                                          "found (reached end of argument "
                                          "list)!" );
                                } else if (( it + 1 )->at( 0 ) == '-') {
                                    throw std::out_of_range(
                                        "Argument " + item
                                        + " expects a value but none was "
                                          "found (next argument "
                                        + *( it + 1 )
                                        + " appears to be a flag)!" );
                                }
                                this->set( parse_value( *( it + 1 ) ) );
                                args.erase( it, it + 2 );
                            } else {
                                this->set();
                                args.erase( it );
                            }
                            return parse_status_t::SUCCESS;
                        }
                    }
                    if (this->required()) {
                        return parse_status_t::FAILURE;
                    } else {
                        return parse_status_t::NOT_FOUND;
                    }
                }

                static INTYPE parse_value( const std::string& input ) {
                    INTYPE val;
                    parse_string( input, val );
                    return val;
                }

                virtual TypedArgument& set() {
                    throw std::out_of_range( "Cannot set argument "
                                             + this->tag()
                                             + " to an undefined value!" );
                }

                virtual TypedArgument& set( const INTYPE& val ) {
                    _value   = val;
                    _was_set = true;
                    return *this;
                }

                virtual bool was_set() const override { return _was_set; }

                virtual bool has_default() const override { return _has_default; }

                virtual bool expects_value() const override { return true; }

                virtual std::string default_string() const override {
                    using std::to_string;
                    return to_string( _default );
                }

            protected:
                std::vector<mapping_ptr_t<INTYPE, std::string>> _mappings;
                INTYPE _default;
                INTYPE _value;
                bool _has_default = false;

                bool _was_set = false;
                
        };

        template<typename INTYPE>
        void swap( TypedArgument<INTYPE>& a,
                   TypedArgument<INTYPE>& b ) noexcept {
            using std::swap;
            swap( static_cast<Argument&>( a ), static_cast<Argument&>( b ) );
            swap( a._mappings, b._mappings );
            swap( a._default, b._default );
            swap( a._value, b._value );
            swap( a._was_set, b._was_set );
            swap( a._has_default, b._has_default );
        }
    }  // namespace config
}  // namespace NCPA
