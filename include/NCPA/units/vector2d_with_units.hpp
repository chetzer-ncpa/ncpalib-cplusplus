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

namespace NCPA {
    namespace units {

        /**
         * Vector2DWithUnits class
         */
        template<typename T = double>
        class Vector2DWithUnits : public NCPA::arrays::NDArray<2, T> {
                typedef NCPA::arrays::NDArray<2, T> parent_t;

            public:
                using parent_t::dim;
                using parent_t::fill;
                using parent_t::set;

                // constructors
                Vector2DWithUnits() : parent_t(), _units { nullptr } {}

                Vector2DWithUnits( size_t dim1, size_t dim2 ) :
                    parent_t( dim1, dim2 ), _units { nullptr } {}

                // Construct with no values
                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const Unit *units ) :
                    parent_t( dim1, dim2 ), _units { units } {}

                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const Unit& units ) :
                    parent_t( dim1, dim2 ), _units { &units } {}

                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const char *units ) :
                    Vector2DWithUnits<T>( dim1, dim2,
                                          Units::from_string( units ) ) {}

                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const std::string& units ) :
                    Vector2DWithUnits<T>( dim1, dim2,
                                          Units::from_string( units ) ) {}

                // construct with values from a vector2d_t
                Vector2DWithUnits( const parent_t& values,
                                   const Unit *units ) :
                    parent_t( values ), _units { units } {}

                Vector2DWithUnits( const parent_t& values,
                                   const Unit& units ) :
                    Vector2DWithUnits<T>( values, &units ) {}

                Vector2DWithUnits( const parent_t& values,
                                   const std::string& units ) :
                    Vector2DWithUnits<T>( values,
                                          Units::from_string( units ) ) {}

                Vector2DWithUnits( const parent_t& values,
                                   const char *units ) :
                    Vector2DWithUnits<T>( values,
                                          Units::from_string( units ) ) {}

                // construct with a constant value
                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const ScalarWithUnits<T>& singlevalue ) :
                    Vector2DWithUnits<T>( dim1, dim2, singlevalue.get(),
                                          singlevalue.get_units() ) {}

                Vector2DWithUnits( size_t dim1, size_t dim2, T singleValue,
                                   const Unit& units ) :
                    Vector2DWithUnits<T>( dim1, dim2, &units ) {
                    this->fill( singleValue );
                }

                Vector2DWithUnits( size_t dim1, size_t dim2, T singleValue,
                                   const Unit *units ) :
                    Vector2DWithUnits<T>( dim1, dim2, singleValue, units ) {}

                Vector2DWithUnits( size_t dim1, size_t dim2, T singleValue,
                                   const std::string& units ) :
                    Vector2DWithUnits<T>( dim1, dim2, singleValue,
                                          Units::from_string( units ) ) {}

                Vector2DWithUnits( size_t dim1, size_t dim2, T singleValue,
                                   const char *units ) :
                    Vector2DWithUnits<T>( dim1, dim2, singleValue,
                                          Units::from_string( units ) ) {}

                // copy constructor
                Vector2DWithUnits( const Vector2DWithUnits<T>& source )
                    = default;

                // move constructor
                Vector2DWithUnits( Vector2DWithUnits<T>&& source ) noexcept
                    = default;

                // destructor
                virtual ~Vector2DWithUnits() = default;

                // assignment and swapping
                friend void swap( Vector2DWithUnits<T>& a,
                                  Vector2DWithUnits<T>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<parent_t&>( a ),
                          static_cast<parent_t&>( b ) );
                    swap( a._units, b._units );
                }

                Vector2DWithUnits& operator=( Vector2DWithUnits<T> other ) {
                    swap( *this, other );
                    return *this;
                }

                // methods
                virtual void as_array( ScalarWithUnits<T> **& buffer ) {
                    size_t dim0 = this->dim( 0 ), dim1 = this->dim( 1 );
                    if (buffer == nullptr) {
                        buffer = NCPA::arrays::zeros<ScalarWithUnits<T>>(
                            dim0, dim1 );
                    }
                    for (size_t i = 0; i < dim0; ++i) {
                        for (size_t j = 0; j < dim1; ++j) {
                            buffer[ i ][ j ].set( this->at( i, j ),
                                                  this->get_units() );
                        }
                    }
                }

                virtual void as_array( T **& buffer ) {
                    if (buffer == nullptr) {
                        buffer = NCPA::arrays::zeros<T>( this->dim( 0 ),
                                                         this->dim( 1 ) );
                    }
                    this->get_values( buffer );
                }

                virtual void convert_units( const Unit& new_units ) {
                    // will throw invalid_conversion and leave original units
                    // unchanged if there's an error.  If there's no change in
                    // units, don't bother with the calculation
                    const Unit *oldunits = this->get_units();
                    if (!new_units.equals( *oldunits )) {
                        Vector2DWithUnits<T> buffer( *this );
                        for (auto& val : *this) {
                            val = oldunits->convert_to( val, new_units );
                        }
                        buffer.set_units( new_units );
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
                virtual void fill( T value, const Unit *units ) {
                    this->set( value );
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

                virtual void get_values( size_t& n1, size_t& n2, T **buffer ) {
                    if (n1 == 0 && n2 == 0) {
                        n1 = this->dim( 0 );
                        n2 = this->dim( 1 );
                    } else if (n1 != this->dim( 0 ) || n2 != this->dim( 1 )) {
                        std::ostringstream oss;
                        oss << "get_values: Size mismatch: vector has "
                            << this->dim( 0 ) << "x" << this->dim( 1 )
                            << " elements, but " << n1 << "x" << n2
                            << " elements were requested.";
                        throw std::invalid_argument( oss.str() );
                    }
                    if (buffer == nullptr) {
                        buffer = NCPA::arrays::zeros<T>( n1, n2 );
                    }
                    for (size_t i = 0; i < n1; ++i) {
                        for (size_t j = 0; j < n2; ++j) {
                            buffer[ i ][ j ] = this->at( i, j );
                        }
                    }
                }

                virtual void get_values( T **buffer ) {
                    size_t n1 = this->dim( 0 ), n2 = this->dim( 1 );
                    this->get_values( n1, n2, buffer );
                }

                virtual Vector2DWithUnits<T> as( const Unit& units ) const {
                    Vector2DWithUnits<T> converted( *this );
                    converted.convert_units( units );
                    return converted;
                }

                virtual Vector2DWithUnits<T> as( const Unit *units ) const {
                    return this->as( *units );
                }

                virtual Vector2DWithUnits<T> as( const char *units ) const {
                    return this->as( Units::from_string( units ) );
                }

                virtual Vector2DWithUnits<T> as(
                    const std::string& units ) const {
                    return this->as( Units::from_string( units ) );
                }

                template<typename U>
                Vector2DWithUnits<U> as() const {
                    Vector2DWithUnits<U> newv( this->dim( 0 ), this->dim( 1 ),
                                               this->get_units() );
                    for (size_t i = 0; i < this->dim( 0 ); ++i) {
                        for (size_t j = 0; j < this->dim( 1 ); ++j) {
                            newv.set( i, j, (U)this->get( i, j ) );
                        }
                    }
                    return newv;
                }

                virtual T get( size_t n1, size_t n2 ) const {
                    return this->at( n1, n2 );
                }

                virtual ScalarWithUnits<T> get_scalar( size_t n1,
                                                       size_t n2 ) const {
                    return ScalarWithUnits<T>( this->at( n1, n2 ),
                                               this->get_units() );
                }

                virtual T get_as( size_t n1, size_t n2,
                                  const Unit *units ) const {
                    return _units->convert_to( this->get( n1, n2 ), units );
                }

                virtual T get_as( size_t n1, size_t n2,
                                  const char *units ) const {
                    return _units->convert_to( this->get( n1, n2 ), units );
                }

                virtual T get_as( size_t n1, size_t n2,
                                  const std::string& units ) const {
                    return _units->convert_to( this->get( n1, n2 ), units );
                }

                virtual void set( size_t n1, size_t n2, T value ) {
                    this->at( n1, n2 ) = value;
                }

                virtual void set( size_t n1, size_t n2,
                                  const ScalarWithUnits<T>& value ) {
                    this->set( n1, n2, value.get_as( this->get_units() ) );
                }

                virtual void set( const parent_t& values, const Unit *units ) {
                    this->copy( values );
                    _units = units;
                }

                virtual void set( const parent_t& values, const Unit& units ) {
                    this->set( values, &units );
                }

                virtual void set( const parent_t& values,
                                  const std::string& units ) {
                    this->set( values, Units::from_string( units ) );
                }

                virtual void set( const parent_t& values, const char *units ) {
                    this->set( values, Units::from_string( units ) );
                }

                virtual void set( size_t n1, size_t n2, const T **values,
                                  const std::string& units ) {
                    this->resize( n1, n2 );
                    for (size_t i = 0; i < n1; ++i) {
                        for (size_t j = 0; j < n2; ++j) {
                            this->at( i, j ) = values[ i ][ j ];
                        }
                    }
                    this->set_units( Units::from_string( units ) );
                }

                virtual void set( size_t n1, size_t n2, const T **values,
                                  const char *& units ) {
                    this->set( n1, n2, values, Units::from_string( units ) );
                }

                virtual void set( size_t n1, size_t n2, const T **values,
                                  const Unit *units ) {
                    this->set( n1, n2, values, units );
                }

                virtual void set( size_t n1, size_t n2, const T **values,
                                  const Unit& units ) {
                    this->set( n1, n2, values, &units );
                }

                virtual void set( size_t n1, size_t n2,
                                  const ScalarWithUnits<T> **values ) {
                    this->resize( n1, n2 );
                    for (size_t i = 0; i < n1; i++) {
                        for (size_t j = 0; j < n2; ++j) {
                            this->at( i, j ) = values[ i ][ j ].get();
                        }
                    }
                    this->set_units( values[ 0 ][ 0 ].get_units() );
                }

                virtual void set_units( const Unit *new_units ) {
                    if (new_units != nullptr) {
                        this->set_units( *new_units );
                    }
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

                explicit operator bool() const {
                    // return !this->empty();
                    auto dims = this->shape();
                    return ( dims[ 0 ] > 0 && dims[ 1 ] > 0 );
                }

                Vector2DWithUnits operator+(
                    const Vector2DWithUnits<T>& second ) const {
                    Vector2DWithUnits<T> temp1( *this );
                    temp1 += second;
                    return temp1;
                }

                Vector2DWithUnits operator+=(
                    const Vector2DWithUnits<T>& second ) {
                    if (this->_units->is_convertible_to(
                            *( second._units ) )) {
                        for (size_t i = 0; i < this->size(); i++) {
                            this->at( i )
                                = this->at( i )
                                + second.get_as( i, this->get_units() );
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
