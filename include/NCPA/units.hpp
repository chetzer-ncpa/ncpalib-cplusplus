#pragma once
#include "NCPA/strings.hpp"
#include "NCPA/types.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

// forward declarations
namespace NCPA {
    namespace units {
        class Unit;

        template<typename T, ENABLE_IF( std::is_floating_point<T> )>
        class ScalarWithUnits;

        template<typename T, ENABLE_IF( std::is_floating_point<T> )>
        class VectorWithUnits;
    }  // namespace units
}  // namespace NCPA

// template<typename T, ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
void swap( NCPA::units::ScalarWithUnits<T>&,
           NCPA::units::ScalarWithUnits<T>& ) noexcept;

// template<typename T, ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator==( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator!=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator>=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator<=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator>( const NCPA::units::ScalarWithUnits<T>& a,
                const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator<( const NCPA::units::ScalarWithUnits<T>& a,
                const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
std::ostream& operator<<( std::ostream& output,
                          const NCPA::units::ScalarWithUnits<T>& D );

template<typename T>
void swap( NCPA::units::VectorWithUnits<T>&,
           NCPA::units::VectorWithUnits<T>& ) noexcept;

namespace NCPA {
    namespace units {
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

                // std::out_of_range( "Invalid conversion from " + from->name()
                //                    + " to " + to->name() ) {}
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

                const std::string& name() const { return _name; }

                const std::vector<std::string>& aliases() const {
                    return _aliases;
                }

                const Unit *reference() const {
                    return ( _reference == nullptr ? this : _reference );
                    // return _reference;
                }

                const Unit *base_reference() const {
                    if ( is_base_reference() ) {
                        return this;
                    } else {
                        return _reference->base_reference();
                    }
                }

                bool equals( const Unit *other ) const {
                    bool isbase1 = this->is_base_reference();
                    bool isbase2 = other->is_base_reference();

                    if ( isbase1 != isbase2 ) {
                        return false;
                    } else if ( isbase1 && isbase2 ) {
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

                bool equals( const Unit& other ) const {
                    return this->equals( &other );
                }

                bool is_convertible_to( const Unit& other ) const {
                    return base_reference()->equals(
                        *( other.base_reference() ) );
                }

                bool is_base_reference() const {
                    return ( _reference == nullptr );
                }

                template<typename T = double,
                         ENABLE_IF( std::is_floating_point<T> )>
                T convert_to( T value, const Unit& target ) const {
                    if ( !this->is_convertible_to( target ) ) {
                        throw invalid_conversion<>( *this, target );
                    }
                    if ( target.equals( this ) ) {
                        return value;
                    } else if ( target.equals( reference() ) ) {
                        return this->to_reference<T>( value );
                    } else {
                        return target.from_base_reference<T>(
                            this->to_base_reference<T>( value ) );
                    }
                }

                template<typename T = double,
                         ENABLE_IF( std::is_floating_point<T> )>
                T convert_to( T value, const Unit *target ) const {
                    return this->convert_to( value, *target );
                }

                template<typename T, ENABLE_IF( NCPA::types::is_iterable<T> )>
                T convert_to( const T& values, const Unit& target ) const {
                    T converted;
                    converted.reserve( values.size() );
                    for ( auto it = values.cbegin(); it != values.cend();
                          ++it ) {
                        converted.push_back( this->convert_to( *it, target ) );
                    }
                    return converted;
                }

                template<typename T, ENABLE_IF( NCPA::types::is_iterable<T> )>
                T convert_to( const T& values, const Unit *target ) const {
                    return this->convert_to( values, *target );
                }

                template<typename T = double,
                         ENABLE_IF( std::is_floating_point<T> )>
                void convert_to( size_t N, T *values, const Unit& target,
                                 T *& converted ) const {
                    for ( size_t i = 0; i < N; i++ ) {
                        converted[ i ]
                            = this->convert_to<T>( values[ i ], target );
                    }
                }

                template<typename T = double,
                         ENABLE_IF( std::is_floating_point<T> )>
                void convert_to( size_t N, T *values, const Unit *target,
                                 T *& converted ) const {
                    this->convert_to( N, values, *target, converted );
                }

                template<typename T = double,
                         ENABLE_IF( std::is_floating_point<T> )>
                T to_reference( T value ) const {
                    return ( value + (T)_prescale_offset ) * (T)_scale
                         + (T)_postscale_offset;
                }

                template<typename T = double,
                         ENABLE_IF( std::is_floating_point<T> )>
                T to_base_reference( T value ) const {
                    if ( is_base_reference() ) {
                        return value;
                    } else {
                        return _reference->to_base_reference<T>(
                            this->to_reference<T>( value ) );
                    }
                }

                template<typename T = double,
                         ENABLE_IF( std::is_floating_point<T> )>
                T from_reference( T value ) const {
                    return ( ( value - (T)_postscale_offset ) / (T)_scale )
                         - (T)_prescale_offset;
                }

                template<typename T = double,
                         ENABLE_IF( std::is_floating_point<T> )>
                T from_base_reference( T value ) const {
                    if ( is_base_reference() ) {
                        return value;
                    } else {
                        return from_reference(
                            _reference->from_base_reference( value ) );
                    }
                }

            private:
                std::string _name;
                std::vector<std::string> _aliases;
                double _scale;
                double _prescale_offset;
                double _postscale_offset;

                const Unit *_reference = nullptr;
        };

        // Standard units
        const Unit METERS( "m", { "meters" } );
        const Unit KELVIN( "K", { "Kelvin", "degK" } );
        const Unit METERS_PER_SECOND( "m/s", { "meters per second", "mps",
                                               "meters/second", "m/sec" } );
        const Unit PASCALS( "Pa", { "Pascals" } );
        const Unit KILOGRAMS_PER_CUBIC_METER( "kg/m3",
                                              { "kilograms per cubic meter",
                                                "kgpm3" } );
        const Unit DECIBELS_PER_METER( "dB/m",
                                       { "decibels per meter", "dB/meter" } );
        const Unit KILOGRAMS( "kg", { "kilograms" } );
        const Unit SECONDS( "s", { "seconds", "sec" } );

        // derived units
        const Unit KILOMETERS( "km", { "kilometers" }, &METERS, 1000.0 );
        const Unit MILLIMETERS( "mm", { "millimeters" }, &METERS, 0.001 );
        const Unit CELSIUS( "C", { "Celsius", "degrees C", "degC" }, &KELVIN,
                            1.0, 273.15 );
        const Unit FAHRENHEIT( "F", { "degrees F", "Fahrenheit", "degF" },
                               &CELSIUS, 5.0 / 9.0, -32.0 );
        const Unit KILOMETERS_PER_SECOND( "km/s",
                                          { "kilometers per second", "kmps",
                                            "km/sec", "kilometers/second" },
                                          &METERS_PER_SECOND, 1000.0 );
        const Unit HECTOPASCALS( "hPa", { "hectopascals" }, &PASCALS, 100.0 );
        const Unit MILLIBARS( "mbar", { "millibars" }, &HECTOPASCALS, 1.0 );
        const Unit ATMOSPHERES( "atm", { "atmospheres" }, &PASCALS, 101325.0 );
        const Unit GRAMS_PER_CUBIC_CENTIMETER(
            "g/cm3", { "grams per cubic centimeter", "gpcm3" },
            &KILOGRAMS_PER_CUBIC_METER, 1000.0 );
        const Unit NEPERS_PER_METER( "np/m",
                                     { "nepers/meter", "nepers per meter" },
                                     &DECIBELS_PER_METER, 8.685889638 );
        const Unit DECIBELS_PER_KILOMETER( "dB/km",
                                           { "dB/kilometer",
                                             "dB per kilometer" },
                                           &DECIBELS_PER_METER, 1000.0 );
        const Unit GRAMS( "g", { "grams" }, &KILOGRAMS, 0.001 );
        const Unit DAYS( "day", {}, &SECONDS, 86400.0 );
        const Unit MINUTES( "min", { "minutes" }, &SECONDS, 60.0 );
        const Unit HOURS( "hour", { "hr" }, &SECONDS, 3600.0 );

        namespace details {
            std::unordered_map<std::string, const Unit *> _units_map;
        }

        class Units {
            public:
                // iterable version, string lookup
                template<typename T, ENABLE_IF( NCPA::types::is_iterable<T> )>
                static T convert( T values, const std::string& from,
                                  const std::string& to ) {
                    return Units::from_string( from )->convert_to<T>(
                        values, *( Units::from_string( to ) ) );
                }

                // scalar version, string lookup
                template<typename T, ENABLE_IF( std::is_floating_point<T> )>
                static T convert( T value, const std::string& from,
                                  const std::string& to ) {
                    return Units::from_string( from )->convert_to<T>(
                        value, *( Units::from_string( to ) ) );
                }

                // array version, string lookup
                template<typename T, ENABLE_IF( std::is_floating_point<T> )>
                static void convert( size_t N, T *values,
                                     const std::string& from,
                                     const std::string& to, T *& converted ) {
                    Units::from_string( from )->convert_to<T>(
                        N, values, *( Units::from_string( to ) ), converted );
                }

                static void register_unit( const Unit& u ) {
                    _add_unit_to_map( u.name(), &u );
                    for ( auto it = u.aliases().cbegin();
                          it != u.aliases().cend(); ++it ) {
                        _add_unit_to_map( *it, &u );
                    }
                }

                static const Unit *from_string( const std::string& key ) {
                    _init_map();
                    return NCPA::units::details::_units_map.at(
                        NCPA::strings::to_upper( key ) );
                }

                static bool can_convert( const std::string& u1,
                                         const std::string& u2 ) {
                    return from_string( u1 )->is_convertible_to(
                        *from_string( u2 ) );
                }

                // static bool can_convert( const Unit& u1, const Unit& u2 ) {
                //     return u1.is_convertible_to( u2 );
                // }

                // static bool can_convert( const Unit *u1, const Unit *u2 ) {
                //     return u1->is_convertible_to( *u2 );
                // }


            private:
                static void _add_unit_to_map( const std::string& key,
                                              const Unit *u ) {
                    std::string capskey = NCPA::strings::to_upper( key );
                    if ( NCPA::units::details::_units_map.find( capskey )
                         == NCPA::units::details::_units_map.end() ) {
                        NCPA::units::details::_units_map[ capskey ] = u;
                    } else {
                        std::ostringstream oss;
                        oss << "Unit with name or alias " << key
                            << " is already registered";
                        throw std::invalid_argument( oss.str() );
                    }
                }

                static void _init_map() {
                    if ( NCPA::units::details::_units_map.size() == 0 ) {
                        register_unit( METERS );
                        register_unit( KELVIN );
                        register_unit( METERS_PER_SECOND );
                        register_unit( PASCALS );
                        register_unit( KILOGRAMS_PER_CUBIC_METER );
                        register_unit( DECIBELS_PER_METER );
                        register_unit( KILOGRAMS );
                        register_unit( SECONDS );

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
                    }
                }
        };

        template<typename T = double,
                 NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
        class ScalarWithUnits {
            public:
                ScalarWithUnits() : _value { 0.0 }, _units { nullptr } {}

                ScalarWithUnits( T value ) :
                    _value { value }, _units { nullptr } {}

                ScalarWithUnits( T value, const Unit& property_units ) :
                    _value { value }, _units { &property_units } {}

                // ScalarWithUnits( T value, const Unit& property_units ) :
                //     _value { value }, _units { &property_units } {}

                ScalarWithUnits( T value, const std::string& units ) :
                    _value { value },
                    _units { NCPA::units::Units::from_string( units ) } {}

                // copy constructor
                ScalarWithUnits( const ScalarWithUnits<T>& source ) :
                    _value { source._value }, _units { source._units } {}

                // move constructor
                ScalarWithUnits( ScalarWithUnits<T>&& that ) noexcept {
                    swap( *this, that );
                }

                // destructor
                virtual ~ScalarWithUnits() {}

                // swap
                friend void swap( ScalarWithUnits<T>& first,
                                  ScalarWithUnits<T>& second ) noexcept {
                    using std::swap;
                    swap( first._value, second._value );
                    swap( first._units, second._units );
                }

                virtual T get() const { return _value; }

                virtual const Unit *get_units() const { return _units; }

                virtual T get_as( const Unit& u ) const {
                    return _units->convert_to( _value, u );
                    // return NCPA::units::Units::convert( _value, _units, u );
                }

                virtual T get_as( Unit *u ) const {
                    return this->get_as( *u );
                    // return NCPA::units::Units::convert( _value, _units, u );
                }

                virtual T get_as( const std::string& units ) const {
                    return _units->convert_to(
                        _value, NCPA::units::Units::from_string( units ) );
                    // return NCPA::units::Units::convert( _value, _units, u );
                }

                virtual T get_as( const char *units ) const {
                    return this->get_as( std::string( units ) );
                    // return _units->convert_to(
                    //     _value, NCPA::units::Units::from_string( units ) );
                    // return NCPA::units::Units::convert( _value, _units, u );
                }

                virtual T as( const Unit& u ) const {
                    return this->get_as( u );
                }

                virtual T as( const Unit *u ) const { return this->as( *u ); }

                virtual T as( const std::string& units ) const {
                    return this->get_as( units );
                }

                virtual T as( const char *units ) const {
                    return this->get_as( units );
                }

                virtual void set_value( T newval ) { _value = newval; }

                virtual void set_units( const Unit& new_units ) {
                    _units = &new_units;
                }

                virtual void set_units( const Unit *new_units ) {
                    _units = new_units;
                }

                virtual void set_units( const std::string& units ) {
                    this->set_units(
                        NCPA::units::Units::from_string( units ) );
                }

                virtual void set_units( const char *units ) {
                    this->set_units( std::string( units ) );
                }

                virtual void set( T newval, const Unit& new_units ) {
                    set_value( newval );
                    set_units( new_units );
                }

                virtual void set( T newval, const Unit *new_units ) {
                    set_value( newval );
                    set_units( new_units );
                }

                virtual void set( T newval, const std::string& new_units ) {
                    this->set( newval,
                               NCPA::units::Units::from_string( new_units ) );
                }

                virtual void set( T newval, const char *new_units ) {
                    this->set( newval, std::string( new_units ) );
                }

                virtual void convert( const Unit& new_units ) {
                    this->convert_units( new_units );
                }

                virtual void convert( const Unit *new_units ) {
                    this->convert_units( new_units );
                }

                virtual void convert( const std::string& new_units ) {
                    this->convert_units( new_units );
                }

                virtual void convert( const char *new_units ) {
                    this->convert( std::string( new_units ) );
                }

                virtual void convert_units( const Unit& new_units ) {
                    // will throw invalid_conversion and leave original units
                    // unchanged if there's an error
                    if ( !new_units.equals( *_units ) ) {
                        _do_units_conversion( *_units, new_units );
                    }
                    _units = &new_units;
                }

                virtual void convert_units( const Unit *new_units ) {
                    // will throw invalid_conversion and leave original units
                    // unchanged if there's an error
                    convert_units( *new_units );
                }

                virtual void convert_units( const std::string& units ) {
                    this->convert_units(
                        NCPA::units::Units::from_string( units ) );
                }

                virtual void convert_units( const char *units ) {
                    this->convert_units( std::string( units ) );
                }

                ScalarWithUnits& operator=( ScalarWithUnits<T> other ) {
                    swap( *this, other );
                    return *this;
                }

                ScalarWithUnits operator=( T newval ) {
                    this->_value = newval;
                    return *this;
                }

                ScalarWithUnits operator+(
                    const ScalarWithUnits<T>& second ) const {
                    NCPA::units::ScalarWithUnits<T> temp1( *this );
                    temp1 += second;
                    return temp1;
                }

                ScalarWithUnits operator+( T second ) const {
                    NCPA::units::ScalarWithUnits<T> temp1( *this );
                    temp1._value += second;
                    return temp1;
                }

                ScalarWithUnits operator+=(
                    const ScalarWithUnits<T>& second ) {
                    if ( this->_units->is_convertible_to(
                             *( second._units ) ) ) {
                        this->_value += second.get_as( *( this->_units ) );
                        return *this;
                    } else {
                        throw invalid_conversion<>( this->_units,
                                                    second._units );
                    }
                }

                ScalarWithUnits operator+=( T val ) {
                    this->_value += val;
                    return *this;
                }

                ScalarWithUnits operator-(
                    const ScalarWithUnits<T>& second ) const {
                    NCPA::units::ScalarWithUnits<T> temp1( *this );
                    temp1 -= second;
                    return temp1;
                }

                ScalarWithUnits operator-( T second ) const {
                    NCPA::units::ScalarWithUnits<T> temp1( *this );
                    temp1._value -= second;
                    return temp1;
                }

                ScalarWithUnits operator-=(
                    const ScalarWithUnits<T>& second ) {
                    *this += -second;
                    return *this;
                }

                ScalarWithUnits operator-=( T second ) {
                    *this += -second;
                    return *this;
                }

                ScalarWithUnits operator*=( T factor ) {
                    this->_value *= factor;
                    return *this;
                }

                ScalarWithUnits operator/=( T factor ) {
                    this->_value /= factor;
                    return *this;
                }

                ScalarWithUnits operator*( T factor ) {
                    NCPA::units::ScalarWithUnits<T> temp( *this );
                    temp *= factor;
                    return temp;
                }

                ScalarWithUnits operator/( T factor ) {
                    NCPA::units::ScalarWithUnits<T> temp( *this );
                    temp /= factor;
                    return temp;
                }

                ScalarWithUnits operator-() const {
                    NCPA::units::ScalarWithUnits<T> temp1( *this );
                    temp1._value = -temp1._value;
                    return temp1;
                }

                T over( const ScalarWithUnits<T>& b ) const {
                    return this->_value / b.as( this->_units );
                }

                friend bool ::operator==( const ScalarWithUnits<T>& a,
                                          const ScalarWithUnits<T>& b );
                friend bool ::operator!=( const ScalarWithUnits<T>& a,
                                          const ScalarWithUnits<T>& b );
                friend bool ::operator>=( const ScalarWithUnits<T>& a,
                                          const ScalarWithUnits<T>& b );
                friend bool ::operator<=( const ScalarWithUnits<T>& a,
                                          const ScalarWithUnits<T>& b );
                friend bool ::operator>( const ScalarWithUnits<T>& a,
                                         const ScalarWithUnits<T>& b );
                friend bool ::operator<( const ScalarWithUnits<T>& a,
                                         const ScalarWithUnits<T>& b );
                friend std::ostream& ::operator<<(
                    std::ostream & output, const ScalarWithUnits<T>& D );

            protected:
                T _value;
                // std::stack< NCPA::Units& > _units;
                const Unit *_units = nullptr;

                void _do_units_conversion( const NCPA::units::Unit& fromUnits,
                                           const NCPA::units::Unit& toUnits ) {
                    // try to convert
                    double units_buffer = 0.0;

                    // throws invalid_conversion if conversion is undefined
                    units_buffer = fromUnits.convert_to( _value, toUnits );
                    // units_buffer = NCPA::units::Units::convert(
                    //     _value, fromUnits, toUnits );

                    // successful, so record the units change
                    _value = units_buffer;
                }
        };

        /**
         * VectorWithUnits class
         */
        template<typename T = double,
                 NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
        class VectorWithUnits : public std::vector<ScalarWithUnits<T>> {
            public:
                // constructors
                VectorWithUnits() : std::vector<ScalarWithUnits<T>>() {}

                VectorWithUnits( size_t n_points, const T *property_values,
                                 const Unit& property_units ) :
                    std::vector<ScalarWithUnits<T>>( n_points ) {
                    this->set( n_points, property_values, property_units );
                }

                VectorWithUnits( size_t n_points, const T *property_values,
                                 const std::string& property_units ) :
                    VectorWithUnits( n_points, property_values,
                                     Units::from_string( property_units ) ) {}

                VectorWithUnits( size_t n_points, const T *values,
                                 const char *units ) :
                    VectorWithUnits( n_points, values,
                                     *Units::from_string( units ) ) {}

                VectorWithUnits( size_t n_points,
                                 const ScalarWithUnits<T> *scalarvalues ) :
                    std::vector<ScalarWithUnits<T>>( n_points ) {
                    for ( size_t i = 0; i < n_points; i++ ) {
                        this->at( i ) = scalarvalues[ i ];
                    }
                    normalize_units();
                }

                VectorWithUnits( size_t n_points,
                                 const ScalarWithUnits<T>& singlevalue ) :
                    std::vector<ScalarWithUnits<T>>( n_points, singlevalue ) {}

                VectorWithUnits( size_t n_points, T singleValue,
                                 const Unit& units ) :
                    std::vector<ScalarWithUnits<T>>(
                        n_points, ScalarWithUnits<T>( singleValue, units ) ) {}

                VectorWithUnits( size_t n_points, T singleValue,
                                 const std::string& units ) :
                    std::vector<ScalarWithUnits<T>>(
                        n_points, ScalarWithUnits<T>( singleValue, units ) ) {}

                VectorWithUnits( size_t n_points, T singleValue,
                                 const char *units ) :
                    std::vector<ScalarWithUnits<T>>(
                        n_points, ScalarWithUnits<T>( singleValue, units ) ) {}

                VectorWithUnits( const std::vector<T>& values,
                                 const Unit& units ) :
                    std::vector<ScalarWithUnits<T>>( values.size() ) {
                    for ( size_t i = 0; i < values.size(); i++ ) {
                        this->at( i )
                            = ScalarWithUnits<T>( values[ i ], units );
                    }
                }

                VectorWithUnits( const std::vector<T>& values,
                                 const std::string& units ) :
                    std::vector<ScalarWithUnits<T>>( values.size() ) {
                    for ( size_t i = 0; i < values.size(); i++ ) {
                        this->at( i )
                            = ScalarWithUnits<T>( values[ i ], units );
                    }
                }

                VectorWithUnits( const std::vector<T>& values,
                                 const char *units ) :
                    std::vector<ScalarWithUnits<T>>( values.size() ) {
                    for ( size_t i = 0; i < values.size(); i++ ) {
                        this->at( i ) = ScalarWithUnits<T>(
                            values[ i ], *Units::from_string( units ) );
                    }
                }

                // copy constructor
                VectorWithUnits( const VectorWithUnits<T>& source ) :
                    std::vector<ScalarWithUnits<T>>( source ) {
                    normalize_units();
                }

                // move constructor
                VectorWithUnits( VectorWithUnits<T>&& source ) noexcept :
                    std::vector<ScalarWithUnits<T>>() {
                    ::swap( *this, source );
                }

                // destructor
                virtual ~VectorWithUnits() { this->clear(); }

                // assignment and swapping
                friend void ::swap( VectorWithUnits<T>& first,
                                    VectorWithUnits<T>& second ) noexcept;

                VectorWithUnits& operator=( VectorWithUnits<T> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                // methods
                virtual void as_array( ScalarWithUnits<T> *& buffer,
                                       bool normFirst = true ) {
                    if ( normFirst ) {
                        this->normalize_units();
                    }
                    if ( buffer == nullptr ) {
                        buffer = new ScalarWithUnits<T>[ this->size() ];
                    }
                    size_t i = 0;
                    for ( auto it = this->cbegin(); it != this->cend();
                          ++it ) {
                        buffer[ i++ ] = *it;
                    }
                }

                virtual void as_array( T *& buffer, bool normFirst = true ) {
                    if ( normFirst ) {
                        this->normalize_units();
                    } else if ( !this->is_normalized() ) {
                        throw std::logic_error( "Multiple units present in "
                                                "vector, normalize first!" );
                    }
                    if ( buffer == nullptr ) {
                        buffer = new T[ this->size() ];
                    }
                    this->get_values( buffer );
                }

                virtual std::vector<T> as_std_vector() const {
                    std::vector<T> v( this->size() );
                    auto it1 = this->cbegin();
                    auto it2 = v.begin();
                    for ( ; it1 != this->cend(); ++it1, ++it2 ) {
                        *it2 = it1->get();
                    }
                    return v;
                }

                virtual void convert_units( const Unit& new_units ) {
                    // will throw invalid_conversion and leave original units
                    // unchanged if there's an error if there's no change in
                    // units, don't bother with the calculation
                    const Unit *oldunits = this->get_units();
                    if ( !new_units.equals( *oldunits ) ) {
                        VectorWithUnits<T> buffer( *this );
                        for ( auto it = buffer.begin(); it != buffer.end();
                              ++it ) {
                            it->convert_units( new_units );
                        }
                        std::swap( *this, buffer );
                    }
                }

                virtual void convert_units( const std::string& new_units ) {
                    this->convert_units( *Units::from_string( new_units ) );
                }

                virtual void convert_units( const char *new_units ) {
                    this->convert_units( *Units::from_string( new_units ) );
                }

                // Fill the vector with identical values.  Does not resize.
                virtual void fill( T value, const Unit& units ) {
                    this->fill( ScalarWithUnits<T>( value, units ) );
                }

                virtual void fill( T value, const std::string& units ) {
                    this->fill( ScalarWithUnits<T>( value, units ) );
                }

                virtual void fill( T value, const char *units ) {
                    this->fill(
                        ScalarWithUnits<T>( value, std::string( units ) ) );
                }

                virtual void fill( const ScalarWithUnits<T>& value ) {
                    for ( auto it = this->begin(); it != this->end(); ++it ) {
                        *it = value;
                    }
                }

                virtual const Unit *get_units( bool normFirst = true ) {
                    if ( normFirst ) {
                        this->normalize_units();
                    } else if ( !this->is_normalized() ) {
                        throw std::logic_error( "Multiple units present in "
                                                "vector, normalize first!" );
                    }
                    if ( this->empty() ) {
                        return nullptr;
                    } else {
                        return this->begin()->get_units();
                    }
                }

                virtual void get_values( size_t& n, T *buffer,
                                         bool normFirst = true ) {
                    if ( normFirst ) {
                        this->normalize_units();
                    }
                    size_t i = 0;
                    for ( auto it = this->cbegin(); it != this->cend();
                          ++it ) {
                        buffer[ i++ ] = it->get();
                    }
                    n = this->size();
                }

                virtual void get_values( T *buffer, bool normFirst = true ) {
                    size_t n;
                    this->get_values( n, buffer, normFirst );
                }

                virtual bool is_normalized() const {
                    if ( !this->empty() ) {
                        const Unit *base = this->front().get_units();
                        for ( auto it = this->cbegin(); it != this->cend();
                              ++it ) {
                            // NCPA::units_t u = it->get_units();
                            if ( !( it->get_units()->equals( *base ) ) ) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual void normalize_units() {
                    if ( !this->empty() ) {
                        const Unit *base = this->front().get_units();
                        for ( auto it = this->begin() + 1; it != this->end();
                              ++it ) {
                            it->convert_units( *base );
                        }
                    }
                }

                virtual void set( size_t n_points, const T *values,
                                  const std::string& units ) {
                    this->set( n_points, values,
                               *Units::from_string( units ) );
                }

                virtual void set( size_t n_points, const T *values,
                                  const char *& units ) {
                    this->set( n_points, values,
                               *Units::from_string( units ) );
                }

                virtual void set( size_t n_points, const T *values,
                                  const Unit& units ) {
                    this->resize( n_points );
                    for ( size_t i = 0; i < n_points; i++ ) {
                        this->at( i )
                            = ScalarWithUnits<T>( values[ i ], units );
                    }
                }

                virtual void set( size_t n_points,
                                  const ScalarWithUnits<T> *values ) {
                    this->clear();
                    this->resize( n_points );
                    for ( size_t i = 0; i < n_points; i++ ) {
                        this->at( i ) = values[ i ];
                    }
                }

                virtual void set( const std::vector<T>& values,
                                  const Unit& units ) {
                    this->clear();
                    this->resize( values.size() );
                    for ( size_t i = 0; i < values.size(); i++ ) {
                        this->at( i )
                            = ScalarWithUnits<T>( values[ i ], units );
                    }
                }

                virtual void set( const std::vector<T>& values,
                                  const std::string& units ) {
                    this->set( values, *Units::from_string( units ) );
                }

                virtual void set( const std::vector<T>& values,
                                  const char *units ) {
                    this->set( values, *Units::from_string( units ) );
                }

                virtual void set_units( const Unit& new_units ) {
                    for ( auto it = this->begin(); it != this->end(); ++it ) {
                        it->set_units( new_units );
                    }
                }

                virtual void set_units( const std::string& new_units ) {
                    for ( auto it = this->begin(); it != this->end(); ++it ) {
                        it->set_units( new_units );
                    }
                }

                virtual void set_units( const char *new_units ) {
                    for ( auto it = this->begin(); it != this->end(); ++it ) {
                        it->set_units( new_units );
                    }
                }
        };

    }  // namespace units
}  // namespace NCPA

// template<typename T, NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
std::ostream& operator<<( std::ostream& output,
                          const NCPA::units::ScalarWithUnits<T>& D ) {
    NCPA::units::Unit *u = D.get_units();
    output << D.get() << " " << ( u == nullptr ? "<null>" : u->name() );
    return output;
}

// template<typename T, NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator==( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b ) {
    if ( a.get_units() == nullptr || b.get_units() == nullptr ) {
        throw NCPA::units::invalid_conversion<>( "One or both units is null" );
    } else if ( a.get_units()->is_convertible_to( *( b.get_units() ) ) ) {
        NCPA::units::ScalarWithUnits<T> bb = b;
        bb.convert_units( a._units );
        return ( a._value == b._value );
    } else {
        return false;
    }
}

// template<typename T, NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator!=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b ) {
    return !( a == b );
}

// template<typename T, NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator>( const NCPA::units::ScalarWithUnits<T>& a,
                const NCPA::units::ScalarWithUnits<T>& b ) {
    if ( a.get_units() == nullptr || b.get_units() == nullptr ) {
        throw NCPA::units::invalid_conversion<>( "One or both units is null" );
    } else if ( a.get_units()->is_convertible_to( *( b.get_units() ) ) ) {
        NCPA::units::ScalarWithUnits<T> bb = b;
        bb.convert_units( a._units );
        return ( a._value > b._value );
    } else {
        throw NCPA::units::invalid_conversion<>( a._units, b._units );
    }
}

// template<typename T, NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator<( const NCPA::units::ScalarWithUnits<T>& a,
                const NCPA::units::ScalarWithUnits<T>& b ) {
    if ( a.get_units() == nullptr || b.get_units() == nullptr ) {
        throw NCPA::units::invalid_conversion<>( "One or both units is null" );
    } else if ( a.get_units()->is_convertible_to( *( b.get_units() ) ) ) {
        NCPA::units::ScalarWithUnits<T> bb = b;
        bb.convert_units( a._units );
        return ( a._value < b._value );
    } else {
        throw NCPA::units::invalid_conversion<>( a._units, b._units );
    }
}

// template<typename T, NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator>=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b ) {
    return !( a < b );
}

// template<typename T, NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
template<typename T>
bool operator<=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b ) {
    return !( a > b );
}

/*
 * Friends
 */
template<typename T>
void swap( NCPA::units::VectorWithUnits<T>& a,
           NCPA::units::VectorWithUnits<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<std::vector<NCPA::units::ScalarWithUnits<T>>&>( a ),
          static_cast<std::vector<NCPA::units::ScalarWithUnits<T>>&>( b ) );
}
