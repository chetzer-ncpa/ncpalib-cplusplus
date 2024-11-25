#pragma once

#include "NCPA/strings.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/units/unit.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

// forward declarations
namespace NCPA {
    namespace units {
        template<typename T, ENABLE_FUNCTION_IF_REAL(T)>
        class ScalarWithUnits;
    }  // namespace units
}  // namespace NCPA

// template<typename T, ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
void swap( NCPA::units::ScalarWithUnits<T>&,
           NCPA::units::ScalarWithUnits<T>& ) noexcept;

// template<typename T, ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
bool operator==( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
bool operator!=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
bool operator>=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
bool operator<=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
bool operator>( const NCPA::units::ScalarWithUnits<T>& a,
                const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
bool operator<( const NCPA::units::ScalarWithUnits<T>& a,
                const NCPA::units::ScalarWithUnits<T>& b );

// template<typename T, ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
std::ostream& operator<<( std::ostream& output,
                          const NCPA::units::ScalarWithUnits<T>& D );

namespace NCPA {
    namespace units {
        template<typename T = double, typename std::enable_if<std::is_floating_point<T>::value, int>::type ENABLER>
                // ENABLE_IF_REAL(T)>
        class ScalarWithUnits {
            public:
                ScalarWithUnits() : _value { 0.0 }, _units { nullptr } {}

                ScalarWithUnits( T value ) :
                    _value { value }, _units { nullptr } {}

                ScalarWithUnits( T value, const Unit* property_units ) :
                    _value { value }, _units { property_units } {}

                ScalarWithUnits( T value, const Unit& property_units ) :
                    ScalarWithUnits( value, &property_units ) {}

                // ScalarWithUnits( T value, const Unit& property_units ) :
                //     _value { value }, _units { &property_units } {}

                ScalarWithUnits( T value, const std::string& units ) :
                    _value { value },
                    _units { NCPA::units::Units::from_string( units ) } {}

                ScalarWithUnits( T value, const char* units ) :
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

                virtual T get_as( const Unit *u ) const {
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

                virtual void set( T newval ) {
                    set_value( newval );
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

                friend bool ::operator== <>( const ScalarWithUnits<T>& a,
                                          const ScalarWithUnits<T>& b );
                friend bool ::operator!= <>( const ScalarWithUnits<T>& a,
                                          const ScalarWithUnits<T>& b );
                friend bool ::operator>= <>( const ScalarWithUnits<T>& a,
                                          const ScalarWithUnits<T>& b );
                friend bool ::operator<= <>( const ScalarWithUnits<T>& a,
                                          const ScalarWithUnits<T>& b );
                friend bool ::operator> <>( const ScalarWithUnits<T>& a,
                                         const ScalarWithUnits<T>& b );
                friend bool ::operator< <>( const ScalarWithUnits<T>& a,
                                         const ScalarWithUnits<T>& b );
                friend std::ostream& ::operator<< <>(
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

    }
}

// template<typename T, NO_DEFAULT_ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
std::ostream& operator<<( std::ostream& output,
                          const NCPA::units::ScalarWithUnits<T>& D ) {
    NCPA::units::Unit *u = D.get_units();
    output << D.get() << " " << ( u == nullptr ? "<null>" : u->name() );
    return output;
}

// template<typename T, NO_DEFAULT_ENABLE_FUNCTION_IF_REAL(T)>
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

// template<typename T, NO_DEFAULT_ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
bool operator!=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b ) {
    return !( a == b );
}

// template<typename T, NO_DEFAULT_ENABLE_FUNCTION_IF_REAL(T)>
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

// template<typename T, NO_DEFAULT_ENABLE_FUNCTION_IF_REAL(T)>
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

// template<typename T, NO_DEFAULT_ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
bool operator>=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b ) {
    return !( a < b );
}

// template<typename T, NO_DEFAULT_ENABLE_FUNCTION_IF_REAL(T)>
template<typename T>
bool operator<=( const NCPA::units::ScalarWithUnits<T>& a,
                 const NCPA::units::ScalarWithUnits<T>& b ) {
    return !( a > b );
}
