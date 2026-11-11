#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/units/declarations.hpp"
#include "NCPA/units/unit.hpp"

#include <type_traits>

namespace NCPA {
    namespace units {
        using namespace NCPA::arrays;

        template<typename BaseType>
        class WithUnits {
            public:
                WithUnits() : _value {}, _units { nullptr } {}

                WithUnits( BaseType value ) :
                    _value { std::move( value ) }, _units { nullptr } {}

                WithUnits( BaseType value, const Unit *u ) :
                    _value { std::move( value ) }, _units { u } {}

                WithUnits( BaseType value, const Unit& u ) :
                    _value { std::move( value ) }, _units { &u } {}

                WithUnits( BaseType value, std::string u ) :
                    _value { std::move( value ) }, _units { Units::from_string(u) } {}

                // template<typename... Args>
                // WithUnits( Args... args, const Unit *u ) :
                //     _value( args... ), _units { u } {}

                // template<typename... Args>
                // WithUnits( Args... args, const Unit& u ) :
                //     _value( args... ), _units { &u } {}

                template<typename T>
                WithUnits( const WithUnits<T>& other ) :
                    _value { other.get() }, _units { other.get_units() } {}

                WithUnits( const WithUnits& ) = default;

                // template<typename T,
                //           typename std::enable_if<
                //              std::is_same<T, BaseType>::value, int>::type =
                //              0>
                // WithUnits( const WithUnits<T>& other ) = default;

                WithUnits( WithUnits<BaseType>&& other ) noexcept = default;

                WithUnits<BaseType>& operator=( WithUnits<BaseType> other ) {
                    swap( *this, other );
                    return *this;
                }

                WithUnits<BaseType>& operator=( BaseType other ) {
                    _value = std::move( other );
                    return *this;
                }

                friend void swap( WithUnits<BaseType>& a,
                                  WithUnits<BaseType>& b ) noexcept {
                    using std::swap;
                    swap( a._value, b._value );
                    swap( a._units, b._units );
                }

                WithUnits<BaseType> as( const Unit& units ) const {
                    _check_units_pointer();
                    WithUnits<BaseType> converted( *this );
                    converted.convert_units( units );
                    return converted;
                }

                WithUnits<BaseType> as( const Unit* units ) const {
                    if (units == nullptr) {
                        throw UnitsNotSet("Can't convert to a null unit");
                    }
                    return this->as( *units );
                }

                WithUnits<BaseType> as( std::string units ) const {
                    return this->as( Units::from_string(units) );
                }

                void convert_units( const Unit& new_units ) {
                    // will throw invalid_conversion and leave original units
                    // unchanged if there's an error.  If there's no change in
                    // units, don't bother with the calculation
                    _check_units_pointer();
                    const Unit *oldunits = this->get_units();
                    if (!new_units.equals( *oldunits )) {
                        BaseType buffer
                            = _units->convert_to( _value, new_units );
                        _value = buffer;
                        _units = &new_units;
                    }
                }

                void convert_units( const Unit *new_units ) {
                    _check_units_pointer();
                    if (new_units == nullptr) {
                        throw UnitsNotSet(
                            "Passed null pointer for new units" );
                    }
                    this->convert_units( *new_units );
                }

                void convert_units( std::string new_units ) {
                    this->convert_units( *Units::from_string( new_units ) );
                }

                BaseType& get() { return _value; }

                const BaseType& get() const { return _value; }

                BaseType get_as( const Unit *units ) const {
                    _check_units_pointer();
                    if (units == nullptr) {
                        throw UnitsNotSet( "Target units are null" );
                    }
                    return this->get_units()->convert_to( _value, *units );
                }

                BaseType get_as( const Unit& units ) const {
                    return this->get_as( &units );
                }

                const Unit *get_units() const { return _units; }

                WithUnits<BaseType>& set( const BaseType& val ) {
                    return this->set_value( val );
                }

                WithUnits<BaseType>& set( const BaseType& val,
                                          const Unit *u ) {
                    this->set_value( val );
                    this->set_units( u );
                    return *this;
                }

                WithUnits<BaseType>& set( const BaseType& val,
                                          const Unit& u ) {
                    return this->set( val, &u );
                }

                WithUnits<BaseType>& set_units( const Unit *u ) {
                    _units = u;
                    return *this;
                }

                WithUnits<BaseType>& set_units( const Unit& u ) {
                    return this->set_units( &u );
                }

                WithUnits<BaseType>& set_units( std::string u ) {
                    return this->set_units( Units::from_string( u ) );
                }

                WithUnits<BaseType>& set_value( const BaseType& val ) {
                    _value = val;
                    return *this;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator+( const WithUnits<T>& second ) const {
                    WithUnits<T> temp1( *this );
                    temp1 += second;
                    return temp1;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator+( T second ) const {
                    WithUnits<T> temp1( *this );
                    temp1._value += second;
                    return temp1;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator+=( const WithUnits<T>& second ) {
                    _check_units_pointer();
                    if (this->_units->is_convertible_to(
                            *( second._units ) )) {
                        this->_value += second.get_as( *( this->_units ) );
                        return *this;
                    }
                    throw this->_units->invalid_conversion( *second._units );
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator+=( T val ) {
                    this->_value += val;
                    return *this;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator-( const WithUnits<T>& second ) const {
                    WithUnits<T> temp1( *this );
                    temp1 -= second;
                    return temp1;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator-( T second ) const {
                    WithUnits<T> temp1( *this );
                    temp1._value -= second;
                    return temp1;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator-=( const WithUnits<T>& second ) {
                    *this += -second;
                    return *this;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator-=( T second ) {
                    *this += -second;
                    return *this;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator*=( T factor ) {
                    this->_value *= factor;
                    return *this;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator/=( T factor ) {
                    this->_value /= factor;
                    return *this;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator*( T factor ) {
                    WithUnits<T> temp( *this );
                    temp *= factor;
                    return temp;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator/( T factor ) {
                    WithUnits<T> temp( *this );
                    temp /= factor;
                    return temp;
                }

                template<typename T                         = BaseType,
                         typename std::enable_if<std::is_arithmetic<T>::value,
                                                 int>::type = 0>
                WithUnits<T> operator-() const {
                    WithUnits<T> temp1( *this );
                    temp1._value = -temp1._value;
                    return temp1;
                }


            private:
                BaseType _value;
                const Unit *_units = nullptr;

                void _check_units_pointer() const {
                    if (_units == nullptr) {
                        throw UnitsNotSet();
                    }
                }
        };

        template<typename T>
        typename std::enable_if<std::is_arithmetic<T>::value, bool>::type
            operator==( const WithUnits<T>& a, const WithUnits<T>& b ) {
            if (a.get_units() == nullptr || b.get_units() == nullptr) {
                throw UnitsNotSet(
                    "One or both units in equality test is null" );
            } else if (a.get_units()->is_convertible_to(
                           *( b.get_units() ) )) {
                WithUnits<T> bb = b;
                bb.convert_units( a.get_units() );
                return ( a.get() == bb.get() );
            } else {
                return false;
            }
        }

        template<typename T>
        typename std::enable_if<std::is_arithmetic<T>::value, bool>::type
            operator!=( const WithUnits<T>& a, const WithUnits<T>& b ) {
            return !( a == b );
        }

        template<typename T>
        typename std::enable_if<std::is_arithmetic<T>::value, bool>::type
            operator>=( const WithUnits<T>& a, const WithUnits<T>& b ) {
            return !( a < b );
        }

        template<typename T>
        typename std::enable_if<std::is_arithmetic<T>::value, bool>::type
            operator<=( const WithUnits<T>& a, const WithUnits<T>& b ) {
            return !( a > b );
        }

        template<typename T>
        typename std::enable_if<std::is_arithmetic<T>::value, bool>::type
            operator>( const WithUnits<T>& a, const WithUnits<T>& b ) {
            if (a.get_units() == nullptr || b.get_units() == nullptr) {
                throw InvalidConversion( "One or both units is null in comparison" );
            } else if (a.get_units()->is_convertible_to(
                           *( b.get_units() ) )) {
                WithUnits<T> bb = b;
                bb.convert_units( a.get_units() );
                return ( a.get() > bb.get() );
            }
            throw a.get_units()->invalid_conversion( *b.get_units() );
        }

        template<typename T>
        typename std::enable_if<std::is_arithmetic<T>::value, bool>::type
            operator<( const WithUnits<T>& a, const WithUnits<T>& b ) {
            if (a.get_units() == nullptr || b.get_units() == nullptr) {
                throw InvalidConversion( "One or both units is null in comparison" );
            } else if (a.get_units()->is_convertible_to(
                           *( b.get_units() ) )) {
                WithUnits<T> bb = b;
                bb.convert_units( a.get_units() );
                return ( a.get() < bb.get() );
            }
            throw a.get_units()->invalid_conversion( *b.get_units() );
        }
    }  // namespace units
}  // namespace NCPA
