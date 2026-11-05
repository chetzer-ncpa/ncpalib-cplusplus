#pragma once

#include "NCPA/units/unit.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>

namespace NCPA {
    namespace units {
        namespace details {
            template<typename T,
                     typename std::enable_if<std::is_arithmetic<T>::value,
                                             bool>::type ENABLER = false>
            T _convert( const T& val, const Unit& Ufrom, const Unit& Uto ) {
                return Ufrom.convert_to( val, Uto );
            }
        }  // namespace details

        template<typename BaseType,typename DerivedType>
        class WithUnits {
            public:
                virtual ~WithUnits() = default;

                ScalarWithUnits& operator=( ScalarWithUnits<T> other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual DerivedType& derived( ) {
                    return static_cast<DerivedType&>(*this);
                }

                virtual const DerivedType& derived() const {
                    return static_cast<const DerivedType&>(*this);
                }

                virtual const BaseType& get() const { return _value; }

                virtual const Unit *get_units() const { return _units; }

                virtual BaseType get_as( const Unit& u ) const {
                    _check_units_pointer();
                    return details::_convert( _value, *_units, u );
                }

                virtual BaseType get_as( const std::string& units ) const {
                    _check_units_pointer();
                    return details::_convert(
                        _value, *_units,
                        *Units::from_string( units ) );
                }

                 virtual BaseType get_as( const char *units ) const {
                    return this->get_as( std::string( units ) );
                }

                virtual DerivedType as( const Unit& u ) const {
                    DerivedType converted( *this );
                    converted.convert_units( u );
                    return converted;
                }

                virtual DerivedType as( const Unit *u ) const {
                    return this->as( *u );
                }

                virtual DerivedType as(
                    const std::string& units ) const {
                    return this->as( Units::from_string( units ) );
                }

                virtual DerivedType as( const char *units ) const {
                    return this->as( Units::from_string( units ) );
                }

                virtual void convert( const Unit& new_units ) {
                    _check_units_pointer();
                    _value = details::_convert( _value, *_units, new_units );
                }

                virtual void convert( const Unit *new_units ) {
                    this->convert( *new_units );
                }

                virtual void convert( const std::string& new_units ) {
                    this->convert( *Units::from_string( new_units ) );
                }

                virtual void convert( const char *new_units ) {
                    this->convert( std::string( new_units ) );
                }

                virtual void set( const BaseType& newval ) { _value = newval; }

                virtual void set( const BaseType&  newval, const Unit& new_units ) {
                    set( newval);
                    set_units( new_units );
                }

                virtual void set( const BaseType&  newval, const Unit *new_units ) {
                    set( newval );
                    set_units( new_units );
                }

                virtual void set( const BaseType&  newval, const std::string& new_units ) {
                    this->set( newval,
                               NCPA::units::Units::from_string( new_units ) );
                }

                virtual void set( const BaseType&  newval, const char *new_units ) {
                    this->set( newval, std::string( new_units ) );
                }

                virtual void set_units( const Unit& new_units ) {
                    this->set_units(&new_units);
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





                

                ScalarWithUnits operator=( T newval ) {
                    this->set( newval );
                    return *this;
                }

                ScalarWithUnits operator+(
                    const DerivedType& second ) const {
                    DerivedType temp1( *this );
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
                    _check_units_pointer();
                    if (this->_units->is_convertible_to(
                            *( second._units ) )) {
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
                    return this->_value / b.get_as( this->_units );
                }

                friend bool operator==( const ScalarWithUnits<T>& a,
                                        const ScalarWithUnits<T>& b ) {
                    if (a.get_units() == nullptr || b.get_units() == nullptr) {
                        throw NCPA::units::invalid_conversion<>(
                            "One or both units is null" );
                    } else if (a.get_units()->is_convertible_to(
                                   *( b.get_units() ) )) {
                        NCPA::units::ScalarWithUnits<T> bb = b;
                        bb.convert_units( a._units );
                        return ( a._value == b._value );
                    } else {
                        return false;
                    }
                }

                friend bool operator!=( const ScalarWithUnits<T>& a,
                                        const ScalarWithUnits<T>& b ) {
                    return !( a == b );
                }

                friend bool operator>=( const ScalarWithUnits<T>& a,
                                        const ScalarWithUnits<T>& b ) {
                    return !( a < b );
                }

                friend bool operator<=( const ScalarWithUnits<T>& a,
                                        const ScalarWithUnits<T>& b ) {
                    return !( a > b );
                }

                friend bool operator>( const ScalarWithUnits<T>& a,
                                       const ScalarWithUnits<T>& b ) {
                    if (a.get_units() == nullptr || b.get_units() == nullptr) {
                        throw NCPA::units::invalid_conversion<>(
                            "One or both units is null" );
                    } else if (a.get_units()->is_convertible_to(
                                   *( b.get_units() ) )) {
                        NCPA::units::ScalarWithUnits<T> bb = b;
                        bb.convert_units( a._units );
                        return ( a._value > b._value );
                    } else {
                        throw NCPA::units::invalid_conversion<>( a._units,
                                                                 b._units );
                    }
                }

                friend bool operator<( const ScalarWithUnits<T>& a,
                                       const ScalarWithUnits<T>& b ) {
                    if (a.get_units() == nullptr || b.get_units() == nullptr) {
                        throw NCPA::units::invalid_conversion<>(
                            "One or both units is null" );
                    } else if (a.get_units()->is_convertible_to(
                                   *( b.get_units() ) )) {
                        NCPA::units::ScalarWithUnits<T> bb = b;
                        bb.convert_units( a._units );
                        return ( a._value < b._value );
                    } else {
                        throw NCPA::units::invalid_conversion<>( a._units,
                                                                 b._units );
                    }
                }

                friend std::ostream& operator<<(
                    std::ostream& output, const ScalarWithUnits<T>& D ) {
                    const Unit *u = D.get_units();
                    output << D.get() << " "
                           << ( u == nullptr ? "<null>" : u->name() );
                    return output;
                }











            protected:
                BaseType _value;
                const Unit *_units = nullptr;

                void _check_units_pointer() const {
                    if (_units == nullptr) {
                        throw std::logic_error( "No units have been set!" );
                    }
                }

                void _do_units_conversion( const NCPA::units::Unit& fromUnits,
                                           const NCPA::units::Unit& toUnits ) {
                    // try to convert
                    BaseType units_buffer
                        = fromUnits.convert_to( _value, toUnits );

                    // successful, so record the units change
                    _value = units_buffer;
                }
        };
    }  // namespace units
}  // namespace NCPA
