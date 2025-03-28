#pragma once
#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/strings.hpp"
#include "NCPA/units/scalar_with_units.hpp"
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
        template<typename T, ENABLE_FUNCTION_IF_REAL( T )>
        class VectorWithUnits;
    }  // namespace units
}  // namespace NCPA

template<typename T>
void swap( NCPA::units::VectorWithUnits<T>&,
           NCPA::units::VectorWithUnits<T>& ) noexcept;

namespace NCPA {
    namespace units {

        /**
         * VectorWithUnits class
         */
        template<typename T = double,
                 typename std::enable_if<std::is_floating_point<T>::value,
                                         int>::type ENABLER>
        class VectorWithUnits : public std::vector<T> {
            public:
                // constructors
                VectorWithUnits() : std::vector<T>(), _units { nullptr } {}

                VectorWithUnits( size_t n_points ) : VectorWithUnits<T>() {
                    this->resize( n_points );
                }

                VectorWithUnits( const std::vector<T>& values,
                                 const Unit *units ) :
                    VectorWithUnits<T>() {
                    set( values, units );
                }

                VectorWithUnits( size_t n_points, const Unit *units ) :
                    VectorWithUnits<T>( std::vector<T>( n_points ), units ) {}

                VectorWithUnits( size_t n_points, const char *units ) :
                    VectorWithUnits<T>( std::vector<T>( n_points ),
                                        Units::from_string( units ) ) {}

                VectorWithUnits( size_t n_points, const std::string& units ) :
                    VectorWithUnits<T>( std::vector<T>( n_points ),
                                        Units::from_string( units ) ) {}

                VectorWithUnits( const std::vector<T>& values,
                                 const Unit& units ) :
                    VectorWithUnits<T>( values, &units ) {}

                VectorWithUnits( size_t n_points, const T *property_values,
                                 const Unit *property_units ) :
                    VectorWithUnits<T>(
                        std::vector<T>( property_values,
                                        property_values + n_points ),
                        property_units ) {}

                VectorWithUnits( size_t n_points, const T *property_values,
                                 const Unit& property_units ) :
                    VectorWithUnits<T>( n_points, property_values,
                                        &property_units ) {}

                VectorWithUnits( size_t n_points, const T *property_values,
                                 const std::string& property_units ) :
                    VectorWithUnits<T>(
                        n_points, property_values,
                        Units::from_string( property_units ) ) {}

                VectorWithUnits( size_t n_points, const T *values,
                                 const char *property_units ) :
                    VectorWithUnits<T>(
                        n_points, values,
                        Units::from_string( property_units ) ) {}

                // VectorWithUnits( size_t n_points,
                //                  const ScalarWithUnits<T> *scalarvalues ) :
                //     std::vector<ScalarWithUnits<T>>( n_points ) {
                //     for ( size_t i = 0; i < n_points; i++ ) {
                //         this->at( i ) = scalarvalues[ i ];
                //     }
                //     normalize_units();
                // }

                VectorWithUnits( size_t n_points,
                                 const ScalarWithUnits<T>& singlevalue ) :
                    VectorWithUnits<T>(
                        std::vector<T>( n_points, singlevalue.get() ),
                        singlevalue.get_units() ) {}

                VectorWithUnits( size_t n_points, T singleValue,
                                 const Unit& units ) :
                    VectorWithUnits<T>(
                        std::vector<T>( n_points, singleValue ), &units ) {}

                VectorWithUnits( size_t n_points, T singleValue,
                                 const std::string& units ) :
                    VectorWithUnits<T>( n_points, singleValue,
                                        Units::from_string( units ) ) {}

                VectorWithUnits( size_t n_points, T singleValue,
                                 const char *units ) :
                    VectorWithUnits<T>( n_points, singleValue,
                                        Units::from_string( units ) ) {}

                VectorWithUnits( const std::vector<T>& values,
                                 const std::string& units ) :
                    VectorWithUnits<T>( values, Units::from_string( units ) ) {
                }

                VectorWithUnits( const std::vector<T>& values,
                                 const char *units ) :
                    VectorWithUnits<T>( values, Units::from_string( units ) ) {
                }

                // copy constructor
                VectorWithUnits( const VectorWithUnits<T>& source ) :
                    std::vector<T>( source ) {
                    _units = source._units;
                }

                // move constructor
                VectorWithUnits( VectorWithUnits<T>&& source ) noexcept :
                    std::vector<T>() {
                    ::swap( *this, source );
                }

                // destructor
                virtual ~VectorWithUnits() {}

                // assignment and swapping
                friend void ::swap<>( VectorWithUnits<T>& first,
                                      VectorWithUnits<T>& second ) noexcept;

                VectorWithUnits& operator=( VectorWithUnits<T> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                // methods
                virtual void as_array( ScalarWithUnits<T> *& buffer ) {
                    if ( buffer == nullptr ) {
                        buffer = new ScalarWithUnits<T>[ this->size() ];
                    }
                    size_t i = 0;
                    for ( auto it = this->cbegin(); it != this->cend();
                          ++it ) {
                        buffer[ i++ ] = ScalarWithUnits<T>( *it, *_units );
                    }
                }

                virtual void as_array( T *& buffer ) {
                    if ( buffer == nullptr ) {
                        buffer = new T[ this->size() ];
                    }
                    this->get_values( buffer );
                }

                // virtual std::vector<T> as_std_vector() const {
                //     std::vector<T> v( this->size() );
                //     auto it1 = this->cbegin();
                //     auto it2 = v.begin();
                //     for ( ; it1 != this->cend(); ++it1, ++it2 ) {
                //         *it2 = it1->get();
                //     }
                //     return v;
                // }

                virtual void convert_units( const Unit& new_units ) {
                    // will throw invalid_conversion and leave original units
                    // unchanged if there's an error.  If there's no change in
                    // units, don't bother with the calculation
                    const Unit *oldunits = this->get_units();
                    if ( !new_units.equals( *oldunits ) ) {
                        VectorWithUnits<T> buffer
                            = oldunits->convert_to( *this, new_units );
                        buffer.set_units( new_units );

                        std::swap( *this, buffer );
                        // *this = buffer;
                    }
                }

                virtual void convert_units( const std::string& new_units ) {
                    this->convert_units( *Units::from_string( new_units ) );
                }

                virtual void convert_units( const char *new_units ) {
                    this->convert_units( *Units::from_string( new_units ) );
                }

                // Fill the vector with identical values.  Does not resize.
                virtual void fill( T value, const Unit *units ) {
                    this->assign( this->size(), value );
                    _units = units;
                }

                virtual void fill( T value, const std::string& units ) {
                    this->fill( value, Units::from_string( units ) );
                }

                virtual void fill( T value, const char *units ) {
                    this->fill( value, Units::from_string( units ) );
                }

                virtual void fill( const ScalarWithUnits<T>& value ) {
                    this->fill( value.get(), value.get_units() );
                }

                virtual const Unit *get_units() const { return _units; }

                virtual void get_values( size_t& n, T *buffer ) {
                    if ( n == 0 ) {
                        if ( buffer != nullptr ) {
                            delete[] buffer;
                            buffer = nullptr;
                        }
                        n = this->size();
                    } else if ( n != this->size() ) {
                        std::ostringstream oss;
                        oss << "get_values: Size mismatch: vector has "
                            << this->size() << " elements, but " << n
                            << " elements were requested.";
                        throw std::invalid_argument( oss.str() );
                    }
                    if ( buffer == nullptr ) {
                        buffer = NCPA::arrays::zeros<T>( n );
                    }
                    std::copy( this->cbegin(), this->cend(), buffer );
                }

                virtual void get_values( T *buffer ) {
                    size_t n = this->size();
                    this->get_values( n, buffer );
                }

                virtual std::vector<T> get_values() const {
                    return std::vector<T>( *this );
                }

                virtual VectorWithUnits<T> as( const Unit &units ) const {
                    VectorWithUnits<T> converted( *this );
                    converted.convert_units( units );
                    return converted;
                }

                virtual VectorWithUnits<T> as( const Unit *units ) const {
                    return this->as( *units );
                }

                virtual VectorWithUnits<T> as( const char *units ) const {
                    return this->as( Units::from_string( units ) );
                }

                virtual VectorWithUnits<T> as( const std::string &units ) const {
                    return this->as( Units::from_string( units ) );
                }

                virtual T get( size_t ind ) const { return this->at( ind ); }

                virtual ScalarWithUnits<T> get_scalar( size_t ind ) const {
                    return ScalarWithUnits<T>( this->get( ind ),
                                               this->get_units() );
                }

                virtual T get_as( size_t ind, const Unit *units ) const {
                    return _units->convert_to( this->at( ind ), units );
                }

                virtual T get_as( size_t ind, const char *units ) const {
                    return _units->convert_to( this->at( ind ), units );
                }

                virtual T get_as( size_t ind,
                                  const std::string& units ) const {
                    return _units->convert_to( this->at( ind ), units );
                }

                virtual void set( size_t index, T value ) {
                    this->at( index ) = value;
                }

                virtual void set( size_t index,
                                  const ScalarWithUnits<T>& value ) {
                    this->at( index ) = value.get_as( this->get_units() );
                }

                virtual void set( const std::vector<T>& values,
                                  const Unit *units ) {
                    this->assign( values.cbegin(), values.cend() );
                    _units = units;
                }

                virtual void set( size_t n_points, const T *values,
                                  const std::string& units ) {
                    this->set( std::vector<T>( values, values + n_points ),
                               Units::from_string( units ) );
                }

                virtual void set( size_t n_points, const T *values,
                                  const char *& units ) {
                    this->set( std::vector<T>( values, values + n_points ),
                               Units::from_string( units ) );
                }

                virtual void set( size_t n_points, const T *values,
                                  const Unit *units ) {
                    this->set( std::vector<T>( values, values + n_points ),
                               units );
                }

                virtual void set( size_t n_points,
                                  const ScalarWithUnits<T> *values ) {
                    std::vector<T> v( n_points );
                    for ( size_t i = 0; i < n_points; i++ ) {
                        v[ i ] = values[ i ].get();
                    }
                    this->set( v, values[ 0 ].get_units() );
                }

                virtual void set( const std::vector<T>& values,
                                  const std::string& units ) {
                    this->set( values, Units::from_string( units ) );
                }

                virtual void set( const std::vector<T>& values,
                                  const char *units ) {
                    this->set( values, Units::from_string( units ) );
                }

                virtual void set_units( const Unit& new_units ) {
                    _units = &new_units;
                }

                virtual void set_units( const std::string& new_units ) {
                    _units = Units::from_string( new_units );
                }

                virtual void set_units( const char *new_units ) {
                    _units = Units::from_string( new_units );
                }

                explicit operator bool() const { return !this->empty(); }

                VectorWithUnits operator+(
                    const VectorWithUnits<T>& second ) const {
                    VectorWithUnits<T> temp1( *this );
                    temp1 += second;
                    return temp1;
                }

                VectorWithUnits operator+=(
                    const VectorWithUnits<T>& second ) {
                    if ( this->_units->is_convertible_to(
                             *( second._units ) ) ) {
                        for ( size_t i = 0; i < this->size(); i++ ) {
                            this->at( i )
                                += second.get_as( i, this->get_units() );
                        }
                        return *this;
                    } else {
                        throw invalid_conversion<>( this->_units,
                                                    second._units );
                    }
                }

            protected:
                const Unit *_units = nullptr;
        };
    }  // namespace units
}  // namespace NCPA

/*
 * Friends
 */
template<typename T>
void swap( NCPA::units::VectorWithUnits<T>& a,
           NCPA::units::VectorWithUnits<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<std::vector<T>&>( a ),
          static_cast<std::vector<T>&>( b ) );
    swap( a._units, b._units );
}
