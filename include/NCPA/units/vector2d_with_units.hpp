#pragma once
#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/ndvector.hpp"
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
        class Vector2DWithUnits;
    }  // namespace units
}  // namespace NCPA

template<typename T>
void swap( NCPA::units::Vector2DWithUnits<T>&,
           NCPA::units::Vector2DWithUnits<T>& ) noexcept;

namespace NCPA {
    namespace units {

        /**
         * Vector2DWithUnits class
         */
        template<typename T = double,
                 typename std::enable_if<std::is_floating_point<T>::value,
                                         int>::type ENABLER>
        class Vector2DWithUnits : public NCPA::arrays::vector2d_t<T> {
            public:
                using NCPA::arrays::vector2d_t<T>::dim;

                // constructors
                Vector2DWithUnits() :
                    NCPA::arrays::vector2d_t<T>(), _units { nullptr } {}

                Vector2DWithUnits( size_t dim1, size_t dim2 ) :
                    NCPA::arrays::vector2d_t<T>( dim1, dim2 ),
                    _units { nullptr } {}

                // Construct with no values
                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const Unit *units ) :
                    NCPA::arrays::vector2d_t<T>( dim1, dim2 ),
                    _units { units } {}

                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const Unit& units ) :
                    NCPA::arrays::vector2d_t<T>( dim1, dim2 ),
                    _units { &units } {}

                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const char *units ) :
                    Vector2DWithUnits<T>(
                        NCPA::arrays::vector2d_t<T>( dim1, dim2 ),
                        Units::from_string( units ) ) {}

                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const std::string& units ) :
                    Vector2DWithUnits<T>(
                        NCPA::arrays::vector2d_t<T>( dim1, dim2 ),
                        Units::from_string( units ) ) {}

                // construct with values from a vector2d_t
                Vector2DWithUnits( const NCPA::arrays::vector2d_t<T>& values,
                                   const Unit *units ) :
                    NCPA::arrays::vector2d_t<T>( values ), _units { units } {}

                Vector2DWithUnits( const NCPA::arrays::vector2d_t<T>& values,
                                   const Unit& units ) :
                    Vector2DWithUnits<T>( values, &units ) {}

                Vector2DWithUnits( const NCPA::arrays::vector2d_t<T>& values,
                                   const std::string& units ) :
                    Vector2DWithUnits<T>( values,
                                          Units::from_string( units ) ) {}

                Vector2DWithUnits( const NCPA::arrays::vector2d_t<T>& values,
                                   const char *units ) :
                    Vector2DWithUnits<T>( values,
                                          Units::from_string( units ) ) {}

                // construct with a constant value
                Vector2DWithUnits( size_t dim1, size_t dim2,
                                   const ScalarWithUnits<T>& singlevalue ) :
                    Vector2DWithUnits<T>( NCPA::arrays::vector2d_t<T>(
                                              dim1, dim2, singlevalue.get() ),
                                          singlevalue.get_units() ) {}

                Vector2DWithUnits( size_t dim1, size_t dim2, T singleValue,
                                   const Unit& units ) :
                    Vector2DWithUnits<T>(
                        NCPA::arrays::vector2d_t<T>( dim1, dim2, singleValue ),
                        &units ) {}

                Vector2DWithUnits( size_t dim1, size_t dim2, T singleValue,
                                   const Unit *units ) :
                    Vector2DWithUnits<T>(
                        NCPA::arrays::vector2d_t<T>( dim1, dim2, singleValue ),
                        units ) {}

                Vector2DWithUnits( size_t dim1, size_t dim2, T singleValue,
                                   const std::string& units ) :
                    Vector2DWithUnits<T>(
                        NCPA::arrays::vector2d_t<T>( dim1, dim2, singleValue ),
                        Units::from_string( units ) ) {}

                Vector2DWithUnits( size_t dim1, size_t dim2, T singleValue,
                                   const char *units ) :
                    Vector2DWithUnits<T>(
                        NCPA::arrays::vector2d_t<T>( dim1, dim2, singleValue ),
                        Units::from_string( units ) ) {}

                // copy constructor
                Vector2DWithUnits( const Vector2DWithUnits<T>& source ) :
                    NCPA::arrays::vector2d_t<T>( source ) {
                    _units = source._units;
                }

                // move constructor
                Vector2DWithUnits( Vector2DWithUnits<T>&& source ) noexcept :
                    Vector2DWithUnits<T>() {
                    ::swap( *this, source );
                }

                // destructor
                virtual ~Vector2DWithUnits() {}

                // assignment and swapping
                friend void ::swap<>( Vector2DWithUnits<T>& first,
                                      Vector2DWithUnits<T>& second ) noexcept;

                Vector2DWithUnits& operator=( Vector2DWithUnits<T> other ) {
                    ::swap( *this, other );
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
                            buffer[ i ][ j ].set( this->at( i ).at( j ),
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
                        for (size_t i = 0; i < this->dim( 0 ); ++i) {
                            buffer.at( i ) = oldunits->convert_to(
                                this->at( i ), new_units );
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
                    for (auto it1 = this->begin(); it1 != this->end(); ++it1) {
                        for (auto it2 = it1->begin(); it2 != it1->end();
                             ++it2) {
                            *it2 = value;
                        }
                    }
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
                        this->size2d( n1, n2 );
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
                            buffer[ i ][ j ] = this->at( i ).at( j );
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
                    Vector2DWithUnits<U> newv( this->dim(0), this->dim(1), this->get_units() );
                    for (size_t i = 0; i < this->dim(0); ++i) {
                        for (size_t j = 0; j < this->dim(1); ++j) {
                            newv.set(i,j,(U)this->get(i,j));
                        }
                    }
                    return newv;
                }

                virtual T get( size_t n1, size_t n2 ) const {
                    return this->at( n1 ).at( n2 );
                }

                virtual ScalarWithUnits<T> get_scalar( size_t n1,
                                                       size_t n2 ) const {
                    return ScalarWithUnits<T>( this->get( n1, n2 ),
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
                    this->at( n1 ).at( n2 ) = value;
                }

                virtual void set( size_t n1, size_t n2,
                                  const ScalarWithUnits<T>& value ) {
                    this->set( n1, n2, value.get_as( this->get_units() ) );
                }

                virtual void set( const NCPA::arrays::vector2d_t<T>& values,
                                  const Unit *units ) {
                    this->assign( values.cbegin(), values.cend() );
                    _units = units;
                }

                virtual void set( const NCPA::arrays::vector2d_t<T>& values,
                                  const Unit& units ) {
                    this->assign( values.cbegin(), values.cend() );
                    _units = &units;
                }

                virtual void set( const NCPA::arrays::vector2d_t<T>& values,
                                  const std::string& units ) {
                    this->set( values, Units::from_string( units ) );
                }

                virtual void set( const NCPA::arrays::vector2d_t<T>& values,
                                  const char *units ) {
                    this->set( values, Units::from_string( units ) );
                }

                virtual void set( size_t n1, size_t n2, const T **values,
                                  const std::string& units ) {
                    this->set( NCPA::arrays::vector2d_t<T>( n1, n2, values ),
                               Units::from_string( units ) );
                }

                virtual void set( size_t n1, size_t n2, const T **values,
                                  const char *& units ) {
                    this->set( NCPA::arrays::vector2d_t<T>( n1, n2, values ),
                               Units::from_string( units ) );
                }

                virtual void set( size_t n1, size_t n2, const T **values,
                                  const Unit *units ) {
                    this->set( NCPA::arrays::vector2d_t<T>( n1, n2, values ),
                               units );
                }

                virtual void set( size_t n1, size_t n2, const T **values,
                                  const Unit& units ) {
                    this->set( NCPA::arrays::vector2d_t<T>( n1, n2, values ),
                               &units );
                }

                virtual void set( size_t n1, size_t n2,
                                  const ScalarWithUnits<T> **values ) {
                    NCPA::arrays::vector2d_t<T> v( n1, n2 );
                    for (size_t i = 0; i < n1; i++) {
                        for (size_t j = 0; j < n2; ++j) {
                            v[ i ][ j ] = values[ i ][ j ].get();
                        }
                    }
                    this->set( v, values[ 0 ][ 0 ].get_units() );
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
                    return (dims[0] > 0 && dims[1] > 0);
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

/*
 * Friends
 */
template<typename T>
void swap( NCPA::units::Vector2DWithUnits<T>& a,
           NCPA::units::Vector2DWithUnits<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::arrays::vector2d_t<T>&>( a ),
          static_cast<NCPA::arrays::vector2d_t<T>&>( b ) );
    swap( a._units, b._units );
}
