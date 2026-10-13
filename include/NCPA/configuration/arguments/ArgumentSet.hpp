#pragma once

#include "NCPA/configuration/arguments/Argument.hpp"
#include "NCPA/configuration/arguments/ArgumentParser.hpp"
#include "NCPA/configuration/arguments/declarations.hpp"
#include "NCPA/configuration/arguments/FlagArgument.hpp"
#include "NCPA/configuration/arguments/MultiValueArgument.hpp"
#include "NCPA/configuration/arguments/SingleValueArgument.hpp"
#include "NCPA/configuration/arguments/TypedArgument.hpp"
#include "NCPA/configuration/exceptions.hpp"

#include <map>
#include <type_traits>

namespace NCPA {
    namespace config {
        class ArgumentSet {
            public:
                ArgumentSet()          = default;
                virtual ~ArgumentSet() = default;

                ArgumentSet( const ArgumentSet& other ) {
                    for (auto& argpair : other._args) {
                        _args.insert( std::make_pair(
                            argpair.first, argpair.second->clone() ) );
                    }
                }

                ArgumentSet( ArgumentSet&& other ) noexcept {
                    swap( *this, other );
                }

                ArgumentSet& operator=( ArgumentSet other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( ArgumentSet& a, ArgumentSet& b ) noexcept {
                    using std::swap;
                    swap( a._args, b._args );
                }

                Argument& add_argument( const Argument *arg ) {
                    _args.insert( std::make_pair( arg->tag(), arg->clone() ) );
                    return *_args.at( arg->tag() );
                }

                Argument& add_argument( const Argument& arg ) {
                    _args.insert( std::make_pair( arg.tag(), arg.clone() ) );
                    return *_args.at( arg.tag() );
                }

                Argument& add_argument(
                    const std::unique_ptr<Argument>& arg ) {
                    return this->add_argument( arg.get() );
                }

                template<typename ARGTYPE = int, typename... Args>
                Argument& add_argument( value_count_t valcount,
                                        const std::string& tag,
                                        const std::string& description,
                                        Args&&...args ) {
                    _args.insert( std::make_pair(
                        tag,
                        std::unique_ptr<Argument>( new TypedArgument<ARGTYPE>(
                            tag, description, valcount,
                            std::forward<Args>( args )... ) ) ) );

                    return this->argument( tag );
                }

                template<typename ARGTYPE, typename... Args>
                Argument& add_argument( const std::string& tag,
                                        const std::string& description,
                                        Args&&...args ) {
                    return this->add_argument<ARGTYPE>(
                        value_count_t::ZERO_OR_ONE, tag, description,
                        std::forward<Args>( args )... );
                }

                template<typename... Args>
                Argument& add_flag( const std::string& tag,
                                    const std::string& description,
                                    Args&&...args ) {
                    return this->add_argument<bool>(
                        value_count_t::ZERO, tag, description,
                        std::forward<Args>( args )... );
                }

                Argument& add_help_flag( const std::string& tag = "help" ) {
                    return this->add_flag( tag, "Display help text and exit" );
                }

                template<typename ARGTYPE, typename... Args>
                Argument& add_multivalue_argument(
                    const std::string& tag, const std::string& description,
                    Args&&...args ) {
                    return this->add_argument<ARGTYPE>(
                        value_count_t::ANY, tag, description,
                        std::forward<Args>( args )... );
                }

                template<typename ARGTYPE, typename... Args>
                Argument& add_optional_argument(
                    const std::string& tag, const std::string& description,
                    Args&&...args ) {
                    return this->add_argument<ARGTYPE>(
                        value_count_t::ZERO_OR_ONE, tag, description,
                        std::forward<Args>( args )... );
                }

                template<typename ARGTYPE, typename... Args>
                Argument& add_required_multivalue_argument(
                    const std::string& tag, const std::string& description,
                    Args&&...args ) {
                    return this->add_argument<ARGTYPE>(
                        value_count_t::AT_LEAST_ONE, tag, description,
                        std::forward<Args>( args )... );
                }

                template<typename ARGTYPE, typename... Args>
                Argument& add_required_argument(
                    const std::string& tag, const std::string& description,
                    Args&&...args ) {
                    return this->add_argument<ARGTYPE>(
                        value_count_t::ONE, tag, description,
                        std::forward<Args>( args )... );
                }

                Argument& argument( const std::string& key ) {
                    try {
                        return *_args.at( key );
                    } catch (std::out_of_range& e) {
                        throw ArgumentNotDefined( "Argument " + key
                                                  + " not defined" );
                    }
                }

                const Argument& argument( const std::string& key ) const {
                    try {
                        return *_args.at( key );
                    } catch (std::out_of_range& e) {
                        throw ArgumentNotDefined( "Argument " + key
                                                  + " not defined" );
                    }
                }

                const argument_map_t& arguments() const { return _args; }

                size_t count() const { return _args.size(); }

                ValidationStatus validate() const {
                    ValidationStatus status;
                    for (auto& arg : _args) {
                        if (arg.second->required() && !arg.second->was_set()) {
                            status.fail( arg.first + ": Required argument not found");
                        }
                        status &= arg.second->validate();
                    }
                    return status;
                }

            protected:
                argument_map_t _args;
        };
    }  // namespace config
}  // namespace NCPA
