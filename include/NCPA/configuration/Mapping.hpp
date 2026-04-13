/*
Mapping.hpp: Defines classes that take a value (presumably from an Argument),
turn it into something else, and then send the converted value to the set()
method of a Configurable.  You can make your own Mapping subclasses in a few
different ways.

Easiest:
Use the DECLARE_MAPPING_TYPE macro.  This is used as:

DECLARE_MAPPING_TYPE( SUBCLASS, INTYPE, OUTTYPE, KEYTYPE, LAMBDA )

where:
    * SUBCLASS is the name of the new class derived from mapping,
    * INTYPE is the input type (e.g. provided by the Argument),
    * OUTTYPE is the type expected by the Configurable,
    * KEYTYPE is the type that the Configurable uses as a key (usually
std::string), and
    * LAMBDA is the conversion function, which must be provided inside
parentheses.

For example, to declare a Mapping that takes a boolean value and
outputs its integer equivalent, you would use

DECLARE_MAPPING_TYPE( BoolToIntMapping, bool, int, std::string,
    ( [](bool b){ return (int)b; } ) )

This can also be used to declare other templated mappings, as in

template<typename T, typename KEYTYPE = std::string>
DECLARE_MAPPING_TYPE( DirectMapping, T, T, KEYTYPE,
                            ( []( T d ) { return d; } ) )

which declares a Mapping template that just passes its value directly, and can
be instantiated as

DirectMapping<double> direct;


More flexible:
Derive from Mapping directly.  This will allow you to use conversions that may
not be easily defined as a lambda.  You need to provide a default constructor,
a full constructor, and a virtual destructor.  To define the above
DirectMapping in this way, use:

template<typename T, typename KEYTYPE = std::string>
class DirectMapping : public Mapping<T, T, KEYTYPE> {
    public:
        DirectMapping() : Mapping<T,T,KEYTYPE>() {}
        DirectMapping( KEYTYPE key ) :
            Mapping<T, T, KEYTYPE>( key, []( T t ) { return t; } )
                    {}
        virtual ~DirectMapping() {}
};

*/
#pragma once

#include "NCPA/configuration/Configurable.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <functional>
#include <stdexcept>
#include <string>

#define DECLARE_MAPPING_TYPE( _SUBCLASS_, _INTYPE_, _OUTTYPE_, _KEYTYPE_,   \
                              _LAMBDA_ )                                    \
    class _SUBCLASS_ : public Mapping<_INTYPE_, _OUTTYPE_, _KEYTYPE_> {     \
        public:                                                             \
            _SUBCLASS_() : Mapping<_INTYPE_, _OUTTYPE_, _KEYTYPE_>() {}     \
            _SUBCLASS_( _KEYTYPE_ key ) :                                   \
                Mapping<_INTYPE_, _OUTTYPE_, _KEYTYPE_>( key, _LAMBDA_ ) {} \
            virtual ~_SUBCLASS_() {}                                        \
    };

namespace NCPA {
    namespace config {

        template<typename INTYPE, typename KEYTYPE>
        class Mapping {
            public:
                Mapping() {}

                Mapping( const KEYTYPE& key ) : _target_key { key } {}

                Mapping( const Mapping<INTYPE, KEYTYPE>& other ) :
                    _target_key { other._target_key } {}

                virtual ~Mapping() {}

                virtual void apply( const INTYPE& in ) = 0;

                virtual const KEYTYPE& key() const { return _target_key; }

            protected:
                KEYTYPE _target_key;
        };

        template<typename INTYPE, typename OUTTYPE,
                 typename KEYTYPE = std::string>
        class ConfigurationMapping : public Mapping<INTYPE, KEYTYPE> {
            public:
                ConfigurationMapping() : Mapping<INTYPE, KEYTYPE>() {}

                ConfigurationMapping(
                    const KEYTYPE& key,
                    std::function<OUTTYPE( INTYPE )> converter, Configurable<KEYTYPE> &target ) :
                    Mapping<INTYPE, KEYTYPE>( key ),
                    _converter { converter }, _target{ &target } {}

                virtual ~ConfigurationMapping() {}

                virtual void apply( const INTYPE& in ) override {
                    std::cout << "Applying value of " << this->convert( in )
                              << " to key " << this->key() << std::endl;
                    _target->set( this->key(), this->convert( in ) );
                }

                virtual OUTTYPE convert( INTYPE in ) const {
                    if (_converter) {
                        return _converter( in );
                    } else {
                        throw std::logic_error(
                            "No conversion function provided!" );
                    }
                }

            protected:
                std::function<OUTTYPE( INTYPE )> _converter;
                Configurable<KEYTYPE> *_target;
        };

        /*
        template<typename INTYPE, typename OUTTYPE, typename KEYTYPE>
        class Mapping {
            public:
                Mapping() {}

                Mapping( const std::string& key,
                         std::function<OUTTYPE( INTYPE )> converter ) :
                    _target_key { key }, _converter { converter } {}

                virtual ~Mapping() {}

                virtual void apply( const INTYPE& in,
                                    Configurable<KEYTYPE>& target ) {
                    target.set<OUTTYPE>( _target_key, this->convert( in ) );
                }

                virtual OUTTYPE convert( INTYPE in ) const {
                    if (_converter) {
                        return _converter( in );
                    } else {
                        throw std::logic_error(
                            "No conversion function provided!" );
                    }
                }

            protected:
                KEYTYPE _target_key;
                std::function<OUTTYPE( INTYPE )> _converter;
        };

        // template<typename T, typename KEYTYPE = std::string>
        // class DirectMapping : public Mapping<T, T, KEYTYPE> {
        //     public:
        //         DirectMapping( KEYTYPE key ) :
        //             Mapping<T, T, KEYTYPE>( key, []( T t ) { return t; } )
        //             {}
        // };
        template<typename T, typename KEYTYPE = std::string>
        DECLARE_MAPPING_TYPE( DirectMapping, T, T, KEYTYPE,
                              ( []( T d ) { return d; } ) )
                                */


    }  // namespace config
}  // namespace NCPA
