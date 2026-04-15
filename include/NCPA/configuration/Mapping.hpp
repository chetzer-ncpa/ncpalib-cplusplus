/*
Mapping.hpp: Defines classes that take a value (presumably from an Argument),
turn it into something else, and then send the converted value to the set()
method of a Configurable.
*/
#pragma once

#include "NCPA/configuration/Configurable.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/pointers.hpp"

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

// #define DECLARE_MAPPING_TYPE( _SUBCLASS_, _INTYPE_, _OUTTYPE_, _KEYTYPE_,   \
//                               _LAMBDA_ )                                    \
//     class _SUBCLASS_ : public Mapping<_INTYPE_, _KEYTYPE_> {                \
//         public:                                                             \
//             typedef Mapping<_INTYPE_, _KEYTYPE_> parent_t;                  \
//                                                                             \
//             _SUBCLASS_() : Mapping<_INTYPE_, _OUTTYPE_, _KEYTYPE_>() {}     \
//             _SUBCLASS_( _KEYTYPE_ key ) :                                   \
//                 Mapping<_INTYPE_, _OUTTYPE_, _KEYTYPE_>( key, _LAMBDA_ ) {} \
//             virtual ~_SUBCLASS_() {}                                        \
//     };

namespace NCPA {
    namespace config {

        using namespace NCPA::pointers;

        template<typename INTYPE, typename KEYTYPE>
        class Mapping : public Cloneable<Mapping<INTYPE, KEYTYPE>> {
            public:
                Mapping() {}

                Mapping( const KEYTYPE& key ) : _target_key { key } {}

                Mapping( const Mapping<INTYPE, KEYTYPE>& other ) :
                    _target_key { other._target_key } {}

                virtual ~Mapping() {}

                virtual void apply( const INTYPE& in ) = 0;

                virtual const KEYTYPE& key() const { return _target_key; }

                virtual void set_key( const std::string& newkey ) {
                    _target_key = newkey;
                }

            protected:
                KEYTYPE _target_key;
        };

        template<typename INTYPE, typename OUTTYPE,
                 typename KEYTYPE = std::string>
        class ConfigurationMapping : public Mapping<INTYPE, KEYTYPE> {
            public:
                typedef Mapping<INTYPE, KEYTYPE> parent_t;
                typedef ConfigurationMapping<INTYPE, OUTTYPE, KEYTYPE> this_t;

                ConfigurationMapping() : parent_t() {}

                ConfigurationMapping(
                    const KEYTYPE& key,
                    std::function<OUTTYPE( INTYPE )> converter,
                    Configurable<KEYTYPE>& target ) :
                    parent_t( key ),
                    _converter { converter },
                    _target { &target } {}

                ConfigurationMapping( const this_t& other ) :
                    parent_t( other ),
                    _converter { other._converter },
                    _target { other._target } {}

                virtual ~ConfigurationMapping() {}

                NCPA_CLONE_METHOD( this_t, parent_t )

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

                virtual void set_target( Configurable<KEYTYPE>& newtarget ) {
                    _target = &newtarget;
                }

                virtual void set_converter( std::function<OUTTYPE( INTYPE )> converter ) {
                    _converter = converter;
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
