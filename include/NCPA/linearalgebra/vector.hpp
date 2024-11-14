#pragma once


#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <initializer_list>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE(
            Vector, details::abstract_vector );
    }
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::Vector,
                                             ELEMENTTYPE );


template<typename ELEMENTTYPE>
NCPA::linear::Vector<ELEMENTTYPE> operator+(
    const NCPA::linear::Vector<ELEMENTTYPE>& c1,
    const NCPA::linear::Vector<ELEMENTTYPE>& c2 );
template<typename ELEMENTTYPE>
NCPA::linear::Vector<ELEMENTTYPE> operator-(
    const NCPA::linear::Vector<ELEMENTTYPE>& c1,
    const NCPA::linear::Vector<ELEMENTTYPE>& c2 );
template<typename ELEMENTTYPE>
NCPA::linear::Vector<ELEMENTTYPE> operator*(
    const NCPA::linear::Vector<ELEMENTTYPE>& c1,
    const NCPA::linear::Vector<ELEMENTTYPE>& c2 );

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class Vector<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC> {
            public:
                Vector() {}

                Vector( std::unique_ptr<details::abstract_vector<ELEMENTTYPE>>
                            ptr ) :
                    Vector<ELEMENTTYPE>() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                Vector( const Vector<ELEMENTTYPE>& other ) :
                    Vector<ELEMENTTYPE>() {
                    _ptr = std::move( other._ptr->clone() );
                }

                /**
                 * Move constructor.
                 * @param source The vector to assimilate.
                 */
                Vector( Vector<ELEMENTTYPE>&& source ) noexcept :
                    Vector<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                virtual ~Vector() {}

                friend void ::swap<ELEMENTTYPE>(
                    Vector<ELEMENTTYPE>& a, Vector<ELEMENTTYPE>& b ) noexcept;

                /**
                 * Assignment operator.
                 * @param other The vector to assign to this.
                 */
                Vector<ELEMENTTYPE>& operator=( Vector<ELEMENTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual size_t size() const {
                    return _ptr ? _ptr->size() : 0;
                };

                virtual ELEMENTTYPE& get( size_t n ) {
                    _check_ptr();
                    return _ptr->get( n );
                };

                virtual const ELEMENTTYPE& get( size_t n ) const {
                    _check_ptr();
                    return _ptr->get( n );
                }

                virtual std::vector<ELEMENTTYPE> as_std() const {
                    _check_ptr();
                    return _ptr->as_std();
                }

                virtual Vector<ELEMENTTYPE>& clear() {
                    if ( _ptr ) {
                        // _ptr->clear();
                        _ptr.reset();
                    }
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& resize( size_t n ) {
                    _check_ptr();
                    _ptr->resize( n );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& as_array( size_t n,
                                                       ELEMENTTYPE *& vals ) {
                    _check_ptr();
                    _ptr->as_array( n, vals );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& set( size_t n, ELEMENTTYPE val ) {
                    _check_ptr();
                    _ptr->set( n, val );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& set( size_t n,
                                                  ELEMENTTYPE *val ) {
                    _check_ptr();
                    _ptr->set( n, val );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& set(
                    const std::vector<ELEMENTTYPE>& v ) {
                    _check_ptr();
                    _ptr->set( v );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& scale( ELEMENTTYPE val ) {
                    _check_ptr();
                    _ptr->scale( val );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& scale(
                    const Vector<ELEMENTTYPE>& b ) {
                    _check_ptr();
                    _ptr->scale( *( b._ptr ) );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& add(
                    const Vector<ELEMENTTYPE>& b ) {
                    _check_ptr();
                    _ptr->add( *( b._ptr ) );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& add( ELEMENTTYPE b ) {
                    _check_ptr();
                    _ptr->add( b );
                    return *this;
                }

                virtual ELEMENTTYPE dot( const Vector<ELEMENTTYPE>& b ) const {
                    _check_ptr();
                    return _ptr->dot( *( b._ptr ) );
                }

                // implementations, not abstract
                virtual bool equals( const Vector<ELEMENTTYPE>& other ) const {
                    if ( !_ptr || !( other._ptr ) ) {
                        return false;
                    }
                    if ( size() != other.size() ) {
                        return false;
                    }
                    for ( auto i = 0; i < size(); i++ ) {
                        if ( get( i ) != other.get( i ) ) {
                            return false;
                        }
                    }
                    return true;
                }

                virtual Vector<ELEMENTTYPE> operator-() const {
                    _check_ptr();
                    return Vector<ELEMENTTYPE>( *this ).scale( -1.0 );
                }

                virtual Vector<ELEMENTTYPE>& operator+=(
                    const Vector<ELEMENTTYPE>& other ) {
                    _check_ptr();
                    this->add( other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator+=(
                    const ELEMENTTYPE& other ) {
                    _check_ptr();
                    this->add( other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator-=(
                    const Vector<ELEMENTTYPE>& other ) {
                    _check_ptr();
                    this->add( -other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator-=(
                    const ELEMENTTYPE& other ) {
                    _check_ptr();
                    this->add( -other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator*=(
                    const Vector<ELEMENTTYPE>& other ) {
                    _check_ptr();
                    this->scale( other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator*=(
                    const ELEMENTTYPE& other ) {
                    _check_ptr();
                    this->scale( other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator/=(
                    const ELEMENTTYPE& other ) {
                    _check_ptr();
                    this->scale( 1.0 / other );
                    return *this;
                }

                virtual ELEMENTTYPE& operator[]( size_t i ) {
                    _check_ptr();
                    return get( i );
                }

                virtual const ELEMENTTYPE& operator[]( size_t i ) const {
                    _check_ptr();
                    return get( i );
                }

                friend bool operator==( const Vector<ELEMENTTYPE>& a,
                                        const Vector<ELEMENTTYPE>& b ) {
                    return a.equals( b );
                }

                friend bool operator!=( const Vector<ELEMENTTYPE>& a,
                                        const Vector<ELEMENTTYPE>& b ) {
                    return !( a.equals( b ) );
                }

                // friend binary operators
                friend Vector<ELEMENTTYPE> operator+(
                    const Vector<ELEMENTTYPE>& c1,
                    const Vector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator-(
                    const Vector<ELEMENTTYPE>& c1,
                    const Vector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out -= c2;
                    return out;
                }

                friend NCPA::linear::Vector<ELEMENTTYPE> operator*(
                    const Vector<ELEMENTTYPE>& c1,
                    const Vector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out *= c2;
                    return out;
                }

            protected:
                void _check_ptr() const {
                    if ( !_ptr ) {
                        throw std::logic_error(
                            "Vector: Internal pointer has not been set!" );
                    }
                }

                void _check_same_size(
                    const Vector<ELEMENTTYPE>& other ) const {
                    if ( size() != other.size() ) {
                        throw std::invalid_argument(
                            "Vectors are not the same size!" );
                    }
                }

            private:
                std::unique_ptr<details::abstract_vector<ELEMENTTYPE>> _ptr;
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::Vector<T>& a,
                  NCPA::linear::Vector<T>& b ) noexcept {
    // using std::swap;
    a._ptr.swap( b._ptr );
}

// template<typename ELEMENTTYPE>
// NCPA::linear::Vector<ELEMENTTYPE> operator+(
//     const NCPA::linear::Vector<ELEMENTTYPE>& c1,
//     const NCPA::linear::Vector<ELEMENTTYPE>& c2 ) {
//     NCPA::linear::Vector<ELEMENTTYPE> out( c1 );
//     out += c2;
//     return out;
// }

// template<typename ELEMENTTYPE>
// NCPA::linear::Vector<ELEMENTTYPE> operator-(
//     const NCPA::linear::Vector<ELEMENTTYPE>& c1,
//     const NCPA::linear::Vector<ELEMENTTYPE>& c2 ) {
//     NCPA::linear::Vector<ELEMENTTYPE> out( c1 );
//     out -= c2;
//     return out;
// }

// template<typename ELEMENTTYPE>
// NCPA::linear::Vector<ELEMENTTYPE> operator*(
//     const NCPA::linear::Vector<ELEMENTTYPE>& c1,
//     const NCPA::linear::Vector<ELEMENTTYPE>& c2 ) {
//     NCPA::linear::Vector<ELEMENTTYPE> out( c1 );
//     out *= c2;
//     return out;
// }
