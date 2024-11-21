#pragma once
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
        //  NO_DEFAULT_ENABLE_IF( std::is_floating_point<T> )>
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
                friend void ::swap<>( VectorWithUnits<T>& first,
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

                explicit operator bool() const { return !this->empty(); }
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
    swap( static_cast<std::vector<NCPA::units::ScalarWithUnits<T>>&>( a ),
          static_cast<std::vector<NCPA::units::ScalarWithUnits<T>>&>( b ) );
}
