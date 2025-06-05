#pragma once
#include "NCPA/defines.hpp"
#include "NCPA/math.hpp"
#include "NCPA/strings.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

// namespace NCPA {
//     namespace units {
//         class Unit;
//     }
// }

// std::ostream& operator<<( std::ostream& output,
//                           const NCPA::units::Unit& U );

namespace NCPA {
    namespace units {
        class Unit;

        template<typename T = NCPA::units::Unit>
        class invalid_conversion : public std::out_of_range {
            public:
                invalid_conversion() :
                    std::out_of_range( "Invalid conversion" ) {}

                invalid_conversion( const std::string& msg ) :
                    std::out_of_range( msg ) {}

                invalid_conversion( const T& from, const T& to ) :
                    std::out_of_range( "Invalid conversion from " + from.name()
                                       + " to " + to.name() ) {}

                invalid_conversion( const T *from, const T *to ) :
                    invalid_conversion( *from, *to ) {}
        };

        class Unit {
            public:
                Unit( const std::string& name                 = "",
                      const std::vector<std::string>& aliases = {},
                      const Unit *reference = nullptr, double scale = 1.0,
                      double prescale = 0.0, double postscale = 0.0 ) :
                    _name { name },
                    _aliases { aliases },
                    _scale { scale },
                    _prescale_offset { prescale },
                    _postscale_offset { postscale },
                    _reference { reference } {}

                Unit( const char *name,
                      const std::vector<std::string>& aliases = {},
                      const Unit *reference = nullptr, double scale = 1.0,
                      double prescale = 0.0, double postscale = 0.0 ) :
                    _name { name },
                    _aliases { aliases },
                    _scale { scale },
                    _prescale_offset { prescale },
                    _postscale_offset { postscale },
                    _reference { reference } {}

                virtual ~Unit() {}

                virtual const std::string& name() const { return _name; }

                virtual const std::vector<std::string>& aliases() const {
                    return _aliases;
                }

                virtual const Unit *reference() const {
                    return ( _reference == nullptr ? this : _reference );
                    // return _reference;
                }

                virtual const Unit *base_reference() const {
                    if (is_base_reference()) {
                        return this;
                    } else {
                        return _reference->base_reference();
                    }
                }

                virtual bool equals( const Unit *other ) const {
                    bool isbase1 = this->is_base_reference();
                    bool isbase2 = other->is_base_reference();

                    if (isbase1 != isbase2) {
                        return false;
                    } else if (isbase1 && isbase2) {
                        return ( this->name() == other->name() );

                    } else {
                        return ( this->reference()->equals(
                                   other->reference() ) )
                            && ( this->_scale == other->_scale )
                            && ( this->_prescale_offset
                                 == other->_prescale_offset )
                            && ( this->_postscale_offset
                                 == other->_postscale_offset );
                    }
                }

                virtual bool equals( const Unit& other ) const {
                    return this->equals( &other );
                }

                virtual bool is_convertible_to( const Unit& other ) const {
                    return base_reference()->equals(
                        *( other.base_reference() ) );
                }

                virtual bool is_base_reference() const {
                    return ( _reference == nullptr );
                }

                template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
                //  ENABLE_FUNCTION_IF_REAL(T)>
                T convert_to( T value, const Unit& target ) const {
                    if (!this->is_convertible_to( target )) {
                        throw invalid_conversion<>( *this, target );
                    }
                    if (target.equals( this )) {
                        return value;
                    } else if (target.equals( reference() )) {
                        return this->to_reference<T>( value );
                    } else {
                        return target.from_base_reference<T>(
                            this->to_base_reference<T>( value ) );
                    }
                }

                template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
                //  ENABLE_FUNCTION_IF_REAL(T)>
                T convert_to( T value, const Unit *target ) const {
                    return this->convert_to( value, *target );
                }

                // template<typename T, ENABLE_IF( NCPA::types::is_iterable<T>
                // )>
                template<typename T, ENABLE_FUNCTION_IF_ITERABLE( T )>
                T convert_to( const T& values, const Unit& target ) const {
                    T converted;
                    converted.reserve( values.size() );
                    for (auto it = values.cbegin(); it != values.cend();
                         ++it) {
                        converted.push_back( this->convert_to( *it, target ) );
                    }
                    return converted;
                }

                // template<typename T, ENABLE_IF( NCPA::types::is_iterable<T>
                // )>
                template<typename T, ENABLE_FUNCTION_IF_ITERABLE( T )>
                T convert_to( const T& values, const Unit *target ) const {
                    return this->convert_to( values, *target );
                }

                template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
                void convert_to( size_t N, T *values, const Unit& target,
                                 T *& converted ) const {
                    for (size_t i = 0; i < N; i++) {
                        converted[ i ]
                            = this->convert_to<T>( values[ i ], target );
                    }
                }

                template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
                void convert_to( size_t N, T *values, const Unit *target,
                                 T *& converted ) const {
                    this->convert_to( N, values, *target, converted );
                }

                template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
                T to_reference( T value ) const {
                    return ( value + (T)_prescale_offset ) * (T)_scale
                         + (T)_postscale_offset;
                }

                template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
                T to_base_reference( T value ) const {
                    if (is_base_reference()) {
                        return value;
                    } else {
                        return _reference->to_base_reference<T>(
                            this->to_reference<T>( value ) );
                    }
                }

                template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
                T from_reference( T value ) const {
                    return ( ( value - (T)_postscale_offset ) / (T)_scale )
                         - (T)_prescale_offset;
                }

                template<typename T = double, ENABLE_FUNCTION_IF_REAL( T )>
                T from_base_reference( T value ) const {
                    if (is_base_reference()) {
                        return value;
                    } else {
                        return from_reference(
                            _reference->from_base_reference( value ) );
                    }
                }

                friend std::ostream& operator<<( std::ostream& output,
                                                 const Unit& D ) {
                    output << D._name;
                    return output;
                }

            private:
                std::string _name;
                std::vector<std::string> _aliases;
                double _scale;
                double _prescale_offset;
                double _postscale_offset;

                const Unit *_reference = nullptr;
        };

        // class NullUnit : public Unit {
        //     public:
        //         NullUnit() {}
        //         virtual ~NullUnit() {}
        //         virtual const std::string& name() const override { return
        //         "N/A"; } virtual const std::vector<std::string>& aliases()
        //         const override { return std::vector<std::string>(); }
        //         virtual const Unit *reference() const override { return
        //         this; } virtual bool equals( const Unit *other ) const
        //         override { return false; } virtual bool is_convertible_to(
        //         const Unit& other ) const override { return false; } virtual
        //         bool is_base_reference() const override { return true; }

        // };

        // Standard units
        // const static NullUnit UNITLESS();
        const static Unit METERS( "m", { "meters" } );
        const static Unit KELVIN( "K", { "Kelvin", "degK" } );
        const static Unit METERS_PER_SECOND( "m/s",
                                             { "meters per second", "mps",
                                               "meters/second", "m/sec" } );
        const static Unit PASCALS( "Pa", { "Pascals" } );
        const static Unit KILOGRAMS_PER_CUBIC_METER(
            "kg/m3", { "kilograms per cubic meter", "kgpm3" } );
        const static Unit DECIBELS_PER_METER( "dB/m", { "decibels per meter",
                                                        "dB/meter" } );
        const static Unit KILOGRAMS( "kg", { "kilograms" } );
        const static Unit SECONDS( "s", { "seconds", "sec" } );
        const static Unit RADIANS( "rad", { "radians" } );
        // const static Unit HERTZ( "Hz", {"Hertz"} );

        // derived units
        const static Unit KILOMETERS( "km", { "kilometers" }, &METERS,
                                      1000.0 );
        const static Unit MILLIMETERS( "mm", { "millimeters" }, &METERS,
                                       0.001 );
        const static Unit CELSIUS( "C", { "Celsius", "degrees C", "degC" },
                                   &KELVIN, 1.0, 273.15 );
        const static Unit FAHRENHEIT( "F",
                                      { "degrees F", "Fahrenheit", "degF" },
                                      &CELSIUS, 5.0 / 9.0, -32.0 );
        const static Unit KILOMETERS_PER_SECOND( "km/s",
                                                 { "kilometers per second",
                                                   "kmps", "km/sec",
                                                   "kilometers/second" },
                                                 &METERS_PER_SECOND, 1000.0 );
        const static Unit HECTOPASCALS( "hPa", { "hectopascals" }, &PASCALS,
                                        100.0 );
        const static Unit MILLIBARS( "mbar", { "millibars" }, &HECTOPASCALS,
                                     1.0 );
        const static Unit ATMOSPHERES( "atm", { "atmospheres" }, &PASCALS,
                                       101325.0 );
        const static Unit GRAMS_PER_CUBIC_CENTIMETER(
            "g/cm3", { "grams per cubic centimeter", "gpcm3" },
            &KILOGRAMS_PER_CUBIC_METER, 1000.0 );
        const static Unit NEPERS_PER_METER( "np/m",
                                            { "nepers/meter",
                                              "nepers per meter" },
                                            &DECIBELS_PER_METER, 8.685889638 );
        const static Unit DECIBELS_PER_KILOMETER(
            "dB/km", { "dB/kilometer", "dB per kilometer" },
            &DECIBELS_PER_METER, 1000.0 );
        const static Unit GRAMS( "g", { "grams" }, &KILOGRAMS, 0.001 );
        const static Unit DAYS( "day", {}, &SECONDS, 86400.0 );
        const static Unit MINUTES( "min", { "minutes" }, &SECONDS, 60.0 );
        const static Unit HOURS( "hour", { "hr" }, &SECONDS, 3600.0 );
        const static Unit DEGREES( "deg", { "degrees" }, &RADIANS,
                                   NCPA::math::PI / 180.0 );

        // const static Unit ANGULAR_FREQUENCY( "omega", {}, &HERTZ, 2.0 *
        // NCPA::math::PI );

        namespace details {
            class units_map_t
                : public std::unordered_map<std::string, const Unit *> {
                public:
                    units_map_t()  {
                        this->init_units();
                    }

                    void init_units() {
                        if (this->size() == 0) {
                            // register_unit( UNITLESS );
                            register_unit( METERS );
                            register_unit( KELVIN );
                            register_unit( METERS_PER_SECOND );
                            register_unit( PASCALS );
                            register_unit( KILOGRAMS_PER_CUBIC_METER );
                            register_unit( DECIBELS_PER_METER );
                            register_unit( KILOGRAMS );
                            register_unit( SECONDS );
                            register_unit( RADIANS );
                            // register_unit( HERTZ);

                            register_unit( KILOMETERS );
                            register_unit( MILLIMETERS );
                            register_unit( CELSIUS );
                            register_unit( FAHRENHEIT );
                            register_unit( KILOMETERS_PER_SECOND );
                            register_unit( HECTOPASCALS );
                            register_unit( MILLIBARS );
                            register_unit( ATMOSPHERES );
                            register_unit( GRAMS_PER_CUBIC_CENTIMETER );
                            register_unit( NEPERS_PER_METER );
                            register_unit( DECIBELS_PER_KILOMETER );
                            register_unit( GRAMS );
                            register_unit( DAYS );
                            register_unit( MINUTES );
                            register_unit( HOURS );
                            register_unit( DEGREES );
                            // register_unit( AZIMUTH );
                            // register_unit( AZIMUTH_RADIANS );
                            // register_unit( ANGULAR_FREQUENCY );
                        }
                    }

                    void register_unit(const Unit& u) {
                        this->_add_unit_to_map( u.name(), &u );
                        for (auto it = u.aliases().cbegin();
                             it != u.aliases().cend(); ++it) {
                            this->_add_unit_to_map( *it, &u );
                        }
                    }

                    void _add_unit_to_map( const std::string& key,
                                           const Unit *u ) {
                        std::string capskey = NCPA::strings::to_upper( key );
                        if (this->find( capskey ) == this->end()) {
                            this->insert( { capskey, u } );
                            // NCPA::units::details::_units_map[ capskey ] = u;
                        } else {
                            std::ostringstream oss;
                            oss << "Unit with name or alias " << key
                                << " is already registered";
                            throw std::invalid_argument( oss.str() );
                        }
                    }
            };

            // static std::unordered_map<std::string, const Unit *> _units_map;
            static units_map_t _units_map;
        }  // namespace details

        class Units {
            public:
                // iterable version, string lookup
                // template<typename T, ENABLE_IF( NCPA::types::is_iterable<T>
                // )>
                template<typename T, ENABLE_FUNCTION_IF_ITERABLE( T )>
                static T convert( T values, const std::string& from,
                                  const std::string& to ) {
                    return Units::from_string( from )->convert_to<T>(
                        values, *( Units::from_string( to ) ) );
                }

                // scalar version, string lookup
                template<typename T, ENABLE_FUNCTION_IF_REAL( T )>
                static T convert( T value, const std::string& from,
                                  const std::string& to ) {
                    return Units::from_string( from )->convert_to<T>(
                        value, *( Units::from_string( to ) ) );
                }

                // array version, string lookup
                template<typename T, ENABLE_FUNCTION_IF_REAL( T )>
                static void convert( size_t N, T *values,
                                     const std::string& from,
                                     const std::string& to, T *& converted ) {
                    Units::from_string( from )->convert_to<T>(
                        N, values, *( Units::from_string( to ) ), converted );
                }

                // static void register_unit( const Unit& u ) {
                //     _add_unit_to_map( u.name(), &u );
                //     for (auto it = u.aliases().cbegin();
                //          it != u.aliases().cend(); ++it) {
                //         _add_unit_to_map( *it, &u );
                //     }
                // }

                static const Unit *from_string( const std::string& key ) {
                    // init_map();
                    return NCPA::units::details::_units_map.at(
                        NCPA::strings::to_upper( key ) );
                }

                static bool can_convert( const std::string& u1,
                                         const std::string& u2 ) {
                    return from_string( u1 )->is_convertible_to(
                        *from_string( u2 ) );
                }

                // static void init_map() {
                //     if (NCPA::units::details::_units_map.size() == 0) {
                //         // register_unit( UNITLESS );
                //         register_unit( METERS );
                //         register_unit( KELVIN );
                //         register_unit( METERS_PER_SECOND );
                //         register_unit( PASCALS );
                //         register_unit( KILOGRAMS_PER_CUBIC_METER );
                //         register_unit( DECIBELS_PER_METER );
                //         register_unit( KILOGRAMS );
                //         register_unit( SECONDS );
                //         register_unit( RADIANS );
                //         // register_unit( HERTZ);

                //         register_unit( KILOMETERS );
                //         register_unit( MILLIMETERS );
                //         register_unit( CELSIUS );
                //         register_unit( FAHRENHEIT );
                //         register_unit( KILOMETERS_PER_SECOND );
                //         register_unit( HECTOPASCALS );
                //         register_unit( MILLIBARS );
                //         register_unit( ATMOSPHERES );
                //         register_unit( GRAMS_PER_CUBIC_CENTIMETER );
                //         register_unit( NEPERS_PER_METER );
                //         register_unit( DECIBELS_PER_KILOMETER );
                //         register_unit( GRAMS );
                //         register_unit( DAYS );
                //         register_unit( MINUTES );
                //         register_unit( HOURS );
                //         register_unit( DEGREES );
                //         // register_unit( AZIMUTH );
                //         // register_unit( AZIMUTH_RADIANS );
                //         // register_unit( ANGULAR_FREQUENCY );
                //     }
                // }

            private:
                static void _add_unit_to_map( const std::string& key,
                                              const Unit *u ) {
                    std::string capskey = NCPA::strings::to_upper( key );
                    if (NCPA::units::details::_units_map.find( capskey )
                        == NCPA::units::details::_units_map.end()) {
                        NCPA::units::details::_units_map[ capskey ] = u;
                    } else {
                        std::ostringstream oss;
                        oss << "Unit with name or alias " << key
                            << " is already registered";
                        throw std::invalid_argument( oss.str() );
                    }
                }
        };

        typedef const NCPA::units::Unit *units_ptr_t;
    }  // namespace units
}  // namespace NCPA
