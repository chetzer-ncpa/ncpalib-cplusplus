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
        class Vector3DWithUnits;
    }  // namespace units
}  // namespace NCPA

template<typename T>
void swap( NCPA::units::Vector3DWithUnits<T>&,
           NCPA::units::Vector3DWithUnits<T>& ) noexcept;

namespace NCPA {
    namespace units {

        /**
         * Vector3DWithUnits class
         */
        template<typename T = double,
                 typename std::enable_if<std::is_floating_point<T>::value,
                                         int>::type ENABLER>
        class Vector3DWithUnits : public NCPA::arrays::vector3d_t<T> {
            public:
                using NCPA::arrays::vector3d_t<T>::dim;

                // constructors
                Vector3DWithUnits() :
                    NCPA::arrays::vector3d_t<T>(), _units { nullptr } {}

                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3 ) :
                    NCPA::arrays::vector3d_t<T>( dim1, dim2, dim3 ),
                    _units { nullptr } {}

                // Construct with no values
                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3,
                                   const Unit *units ) :
                    NCPA::arrays::vector3d_t<T>( dim1, dim2, dim3 ),
                    _units { units } {}

                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3,
                                   const Unit& units ) :
                    NCPA::arrays::vector3d_t<T>( dim1, dim2, dim3 ),
                    _units { &units } {}

                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3,
                                   const char *units ) :
                    Vector3DWithUnits<T>(
                        NCPA::arrays::vector3d_t<T>( dim1, dim2, dim3 ),
                        Units::from_string( units ) ) {}

                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3,
                                   const std::string& units ) :
                    Vector3DWithUnits<T>(
                        NCPA::arrays::vector3d_t<T>( dim1, dim2, dim3 ),
                        Units::from_string( units ) ) {}

                // construct with values from a vector3d_t
                Vector3DWithUnits( const NCPA::arrays::vector3d_t<T>& values,
                                   const Unit *units ) :
                    NCPA::arrays::vector3d_t<T>( values ), _units { units } {}

                Vector3DWithUnits( const NCPA::arrays::vector3d_t<T>& values,
                                   const Unit& units ) :
                    Vector3DWithUnits<T>( values, &units ) {}

                Vector3DWithUnits( const NCPA::arrays::vector3d_t<T>& values,
                                   const std::string& units ) :
                    Vector3DWithUnits<T>( values,
                                          Units::from_string( units ) ) {}

                Vector3DWithUnits( const NCPA::arrays::vector3d_t<T>& values,
                                   const char *units ) :
                    Vector3DWithUnits<T>( values,
                                          Units::from_string( units ) ) {}

                Vector3DWithUnits( const NCPA::arrays::ndvector<3, T>& values,
                                   const Unit *units ) :
                    NCPA::arrays::vector3d_t<T>( values ), _units { units } {}

                // construct with a constant value
                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3,
                                   const ScalarWithUnits<T>& singlevalue ) :
                    Vector3DWithUnits<T>(
                        NCPA::arrays::vector3d_t<T>( dim1, dim2, dim3,
                                                     singlevalue.get() ),
                        singlevalue.get_units() ) {}

                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3,
                                   T singleValue, const Unit& units ) :
                    Vector3DWithUnits<T>( NCPA::arrays::vector3d_t<T>(
                                              dim1, dim2, dim3, singleValue ),
                                          &units ) {}

                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3,
                                   T singleValue, const Unit *units ) :
                    Vector3DWithUnits<T>( NCPA::arrays::vector3d_t<T>(
                                              dim1, dim2, dim3, singleValue ),
                                          units ) {}

                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3,
                                   T singleValue, const std::string& units ) :
                    Vector3DWithUnits<T>( NCPA::arrays::vector3d_t<T>(
                                              dim1, dim2, dim3, singleValue ),
                                          Units::from_string( units ) ) {}

                Vector3DWithUnits( size_t dim1, size_t dim2, size_t dim3,
                                   T singleValue, const char *units ) :
                    Vector3DWithUnits<T>( NCPA::arrays::vector3d_t<T>(
                                              dim1, dim2, dim3, singleValue ),
                                          Units::from_string( units ) ) {}

                // copy constructor
                Vector3DWithUnits( const Vector3DWithUnits<T>& source ) :
                    NCPA::arrays::vector3d_t<T>( source ) {
                    _units = source._units;
                }

                // move constructor
                Vector3DWithUnits( Vector3DWithUnits<T>&& source ) noexcept :
                    Vector3DWithUnits<T>() {
                    ::swap( *this, source );
                }

                // destructor
                virtual ~Vector3DWithUnits() {}

                // assignment and swapping
                friend void ::swap<>( Vector3DWithUnits<T>& first,
                                      Vector3DWithUnits<T>& second ) noexcept;

                Vector3DWithUnits& operator=( Vector3DWithUnits<T> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                // methods
                virtual void as_array( ScalarWithUnits<T> ***& buffer ) {
                    size_t dim0 = this->dim( 0 ), dim1 = this->dim( 1 ),
                           dim2 = this->dim( 2 );
                    if (buffer == nullptr) {
                        buffer = NCPA::arrays::zeros<ScalarWithUnits<T>>(
                            dim0, dim1, dim2 );
                    }
                    for (size_t i = 0; i < dim0; ++i) {
                        for (size_t j = 0; j < dim1; ++j) {
                            for (size_t k = 0; k < dim2; ++k) {
                                buffer[ i ][ j ][ k ].set(
                                    this->at( i ).at( j ).at( k ),
                                    this->get_units() );
                            }
                        }
                    }
                }

                virtual void as_array( T ***& buffer ) {
                    if (buffer == nullptr) {
                        buffer = NCPA::arrays::zeros<T>(
                            this->dim( 0 ), this->dim( 1 ), this->dim( 2 ) );
                    }
                    this->get_values( buffer );
                }

                virtual void convert_units( const Unit& new_units ) {
                    // will throw invalid_conversion and leave original units
                    // unchanged if there's an error.  If there's no change in
                    // units, don't bother with the calculation
                    const Unit *oldunits = this->get_units();
                    if (!new_units.equals( *oldunits )) {
                        Vector3DWithUnits<T> buffer( *this );
                        for (size_t i = 0; i < this->dim( 0 ); ++i) {
                            for (size_t j = 0; j < this->dim( 1 ); ++j) {
                                buffer.at( i ).at( j ) = oldunits->convert_to(
                                    this->at( i ).at( j ), new_units );
                            }
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
                            for (auto it3 = it2->begin(); it3 != it2->end();
                                 ++it3) {
                                *it3 = value;
                            }
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

                virtual void get_values( size_t& n1, size_t& n2, size_t& n3,
                                         T ***buffer ) {
                    if (n1 == 0 && n2 == 0 && n3 == 0) {
                        this->size3d( n1, n2, n3 );
                    } else if (n1 != this->dim( 0 ) || n2 != this->dim( 1 )
                               || n3 != this->dim( 2 )) {
                        std::ostringstream oss;
                        oss << "get_values: Size mismatch: vector has "
                            << this->dim( 0 ) << "x" << this->dim( 1 ) << "x"
                            << this->dim( 2 ) << " elements, but " << n1 << "x"
                            << n2 << "x" << n3 << " elements were requested.";
                        throw std::invalid_argument( oss.str() );
                    }
                    if (buffer == nullptr) {
                        buffer = NCPA::arrays::zeros<T>( n1, n2, n3 );
                    }
                    for (size_t i = 0; i < n1; ++i) {
                        for (size_t j = 0; j < n2; ++j) {
                            for (size_t k = 0; k < n3; ++k) {
                                buffer[ i ][ j ][ k ]
                                    = this->at( i ).at( j ).at( k );
                            }
                        }
                    }
                }

                virtual void get_values( T ***buffer ) {
                    size_t n1 = this->dim( 0 ), n2 = this->dim( 1 ),
                           n3 = this->dim( 2 );
                    this->get_values( n1, n2, n3, buffer );
                }

                virtual Vector3DWithUnits<T> as( const Unit& units ) const {
                    Vector3DWithUnits<T> converted( *this );
                    converted.convert_units( units );
                    return converted;
                }

                virtual Vector3DWithUnits<T> as( const Unit *units ) const {
                    return this->as( *units );
                }

                virtual Vector3DWithUnits<T> as( const char *units ) const {
                    return this->as( Units::from_string( units ) );
                }

                virtual Vector3DWithUnits<T> as(
                    const std::string& units ) const {
                    return this->as( Units::from_string( units ) );
                }

                template<typename U>
                Vector3DWithUnits<U> as() const {
                    Vector3DWithUnits<U> newv( this->dim( 0 ), this->dim( 1 ),
                                               this->dim( 2 ),
                                               this->get_units() );
                    for (size_t i = 0; i < this->dim( 0 ); ++i) {
                        for (size_t j = 0; j < this->dim( 1 ); ++j) {
                            for (size_t k = 0; k < this->dim( 2 ); ++k) {
                                newv.set( i, j, k, (U)this->get( i, j, k ) );
                            }
                        }
                    }
                    return newv;
                }

                virtual T get( size_t n1, size_t n2, size_t n3 ) const {
                    return this->at( n1 ).at( n2 ).at( n3 );
                }

                virtual ScalarWithUnits<T> get_scalar( size_t n1, size_t n2,
                                                       size_t n3 ) const {
                    return ScalarWithUnits<T>( this->get( n1, n2, n3 ),
                                               this->get_units() );
                }

                virtual T get_as( size_t n1, size_t n2, size_t n3,
                                  const Unit *units ) const {
                    return _units->convert_to( this->get( n1, n2, n3 ),
                                               units );
                }

                virtual T get_as( size_t n1, size_t n2, size_t n3,
                                  const char *units ) const {
                    return _units->convert_to( this->get( n1, n2, n3 ),
                                               units );
                }

                virtual T get_as( size_t n1, size_t n2, size_t n3,
                                  const std::string& units ) const {
                    return _units->convert_to( this->get( n1, n2, n3 ),
                                               units );
                }

                virtual void set( size_t n1, size_t n2, size_t n3, T value ) {
                    this->at( n1 ).at( n2 ).at( n3 ) = value;
                }

                virtual void set( size_t n1, size_t n2, size_t n3,
                                  const ScalarWithUnits<T>& value ) {
                    this->set( n1, n2, n3, value.get_as( this->get_units() ) );
                }

                virtual void set( const NCPA::arrays::vector3d_t<T>& values,
                                  const Unit *units ) {
                    this->assign( values.cbegin(), values.cend() );
                    _units = units;
                }

                virtual void set( const NCPA::arrays::vector3d_t<T>& values,
                                  const Unit& units ) {
                    this->assign( values.cbegin(), values.cend() );
                    _units = &units;
                }

                virtual void set( const NCPA::arrays::vector3d_t<T>& values,
                                  const std::string& units ) {
                    this->set( values, Units::from_string( units ) );
                }

                virtual void set( const NCPA::arrays::vector3d_t<T>& values,
                                  const char *units ) {
                    this->set( values, Units::from_string( units ) );
                }

                virtual void set( size_t n1, size_t n2, size_t n3,
                                  const T ***values,
                                  const std::string& units ) {
                    this->set(
                        NCPA::arrays::vector3d_t<T>( n1, n2, n3, values ),
                        Units::from_string( units ) );
                }

                virtual void set( size_t n1, size_t n2, size_t n3,
                                  const T ***values, const char *& units ) {
                    this->set(
                        NCPA::arrays::vector3d_t<T>( n1, n2, n3, values ),
                        Units::from_string( units ) );
                }

                virtual void set( size_t n1, size_t n2, size_t n3,
                                  const T ***values, const Unit *units ) {
                    this->set(
                        NCPA::arrays::vector3d_t<T>( n1, n2, n3, values ),
                        units );
                }

                virtual void set( size_t n1, size_t n2, size_t n3,
                                  const T ***values, const Unit& units ) {
                    this->set(
                        NCPA::arrays::vector3d_t<T>( n1, n2, n3, values ),
                        &units );
                }

                virtual void set( size_t n1, size_t n2, size_t n3,
                                  const ScalarWithUnits<T> ***values ) {
                    NCPA::arrays::vector3d_t<T> v( n1, n2, n3 );
                    for (size_t i = 0; i < n1; i++) {
                        for (size_t j = 0; j < n2; ++j) {
                            for (size_t k = 0; k < n3; ++k) {
                                v[ i ][ j ][ k ] = values[ i ][ j ][ k ].get();
                            }
                        }
                    }
                    this->set( v, values[ 0 ][ 0 ][ 0 ].get_units() );
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

                Vector3DWithUnits operator+(
                    const Vector3DWithUnits<T>& second ) const {
                    Vector3DWithUnits<T> temp1( *this );
                    temp1 += second;
                    return temp1;
                }

                Vector3DWithUnits operator+=(
                    const Vector3DWithUnits<T>& second ) {
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
void swap( NCPA::units::Vector3DWithUnits<T>& a,
           NCPA::units::Vector3DWithUnits<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::arrays::vector3d_t<T>&>( a ),
          static_cast<NCPA::arrays::vector3d_t<T>&>( b ) );
    swap( a._units, b._units );
}
