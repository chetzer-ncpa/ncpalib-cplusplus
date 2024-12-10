#pragma once


#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

// namespace NCPA {
//     namespace linear {
//         // NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE(
//         //     Vector, details::abstract_vector );
//         NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Vector );
//         NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( WrapperVector, Vector
//         );
//     }  // namespace linear
// }  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::Vector,
                                             ELEMENTTYPE );
template<typename ELEMENTTYPE>
std::ostream& operator<<( std::ostream& os,
                          const NCPA::linear::Vector<ELEMENTTYPE>& vec );

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::WrapperVector,
                                             ELEMENTTYPE );
template<typename ELEMENTTYPE>
std::ostream& operator<<(
    std::ostream& os, const NCPA::linear::WrapperVector<ELEMENTTYPE>& vec );

NCPA_LINEARALGEBRA_DECLARE_FRIEND_BINARY_OPERATORS( NCPA::linear::Vector,
                                                    ELEMENTTYPE )
NCPA_LINEARALGEBRA_DECLARE_FRIEND_BINARY_OPERATORS(
    NCPA::linear::WrapperVector, ELEMENTTYPE )

// template<typename ELEMENTTYPE>
// NCPA::linear::Vector<ELEMENTTYPE> operator+(
//     const NCPA::linear::Vector<ELEMENTTYPE>& c1,
//     const NCPA::linear::Vector<ELEMENTTYPE>& c2 );
// template<typename ELEMENTTYPE>
// NCPA::linear::Vector<ELEMENTTYPE> operator-(
//     const NCPA::linear::Vector<ELEMENTTYPE>& c1,
//     const NCPA::linear::Vector<ELEMENTTYPE>& c2 );
// template<typename ELEMENTTYPE>
// NCPA::linear::Vector<ELEMENTTYPE> operator*(
//     const NCPA::linear::Vector<ELEMENTTYPE>& c1,
//     const NCPA::linear::Vector<ELEMENTTYPE>& c2 );

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class Vector<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC> {
            public:
                friend class Matrix<ELEMENTTYPE>;

                Vector() {}

                Vector( std::unique_ptr<details::abstract_vector<ELEMENTTYPE>>
                            ptr ) :
                    Vector<ELEMENTTYPE>() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                Vector( const Vector<ELEMENTTYPE>& other ) :
                    Vector<ELEMENTTYPE>() {
                    _ptr = std::move( other.internal()->clone() );
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

                friend Vector<ELEMENTTYPE> operator*(
                    const Vector<ELEMENTTYPE>& a,
                    const Matrix<ELEMENTTYPE>& b );
                friend Vector<ELEMENTTYPE> operator*(
                    const Matrix<ELEMENTTYPE>& a,
                    const Vector<ELEMENTTYPE>& b );

                /**
                 * Assignment operator.
                 * @param other The vector to assign to this.
                 */
                Vector<ELEMENTTYPE>& operator=( Vector<ELEMENTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual size_t size() const {
                    return internal() ? internal()->size() : 0;
                };

                virtual bool is_zero() const {
                    return !(_ptr && !(_ptr->is_zero()));
                }

                virtual Vector<ELEMENTTYPE>& zero() {
                    check_pointer();
                    internal()->zero();
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& zero( size_t n ) {
                    check_pointer();
                    internal()->zero( n );
                    return *this;
                }

                virtual ELEMENTTYPE& get( size_t n ) {
                    check_pointer();
                    return internal()->get( n );
                };

                virtual const ELEMENTTYPE& get( size_t n ) const {
                    check_pointer();
                    return internal()->get( n );
                }

                virtual std::vector<ELEMENTTYPE> as_std() const {
                    check_pointer();
                    return internal()->as_std();
                }

                virtual Vector<ELEMENTTYPE>& clear() {
                    if ( _ptr ) {
                        // internal()clear();
                        _ptr.reset();
                    }
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& resize( size_t n ) {
                    check_pointer();
                    internal()->resize( n );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& as_array( size_t n,
                                                       ELEMENTTYPE *& vals ) {
                    check_pointer();
                    internal()->as_array( n, vals );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& set( size_t n, ELEMENTTYPE val ) {
                    check_pointer();
                    internal()->set( n, val );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& set( size_t n,
                                                  ELEMENTTYPE *val ) {
                    check_pointer();
                    internal()->set( n, val );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& set(
                    const std::vector<ELEMENTTYPE>& v ) {
                    check_pointer();
                    internal()->set( v );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& set( const ELEMENTTYPE& val ) {
                    check_pointer();
                    internal()->set( val );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& scale( ELEMENTTYPE val ) {
                    check_pointer();
                    internal()->scale( val );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& scale(
                    const Vector<ELEMENTTYPE>& b ) {
                    check_pointer();
                    internal()->scale( *( b.internal() ) );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& add(
                    const Vector<ELEMENTTYPE>& b ) {
                    check_pointer();
                    internal()->add( *( b.internal() ) );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& add( ELEMENTTYPE b ) {
                    check_pointer();
                    internal()->add( b );
                    return *this;
                }

                virtual ELEMENTTYPE dot( const Vector<ELEMENTTYPE>& b ) const {
                    check_pointer();
                    return internal()->dot( *( b.internal() ) );
                }

                // implementations, not abstract
                virtual bool equals( const Vector<ELEMENTTYPE>& other ) const {
                    if ( !internal() || !( other.internal() ) ) {
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
                    check_pointer();
                    return Vector<ELEMENTTYPE>( *this ).scale( -1.0 );
                }

                virtual Vector<ELEMENTTYPE>& operator+=(
                    const Vector<ELEMENTTYPE>& other ) {
                    check_pointer();
                    this->add( other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator+=(
                    const ELEMENTTYPE& other ) {
                    check_pointer();
                    this->add( other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator-=(
                    const Vector<ELEMENTTYPE>& other ) {
                    check_pointer();
                    this->add( -other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator-=(
                    const ELEMENTTYPE& other ) {
                    check_pointer();
                    this->add( -other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator*=(
                    const Vector<ELEMENTTYPE>& other ) {
                    check_pointer();
                    this->scale( other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator*=(
                    const ELEMENTTYPE& other ) {
                    check_pointer();
                    this->scale( other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator/=(
                    const ELEMENTTYPE& other ) {
                    check_pointer();
                    this->scale( 1.0 / other );
                    return *this;
                }

                virtual ELEMENTTYPE& operator[]( size_t i ) {
                    check_pointer();
                    return get( i );
                }

                virtual const ELEMENTTYPE& operator[]( size_t i ) const {
                    check_pointer();
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
                friend std::ostream& operator<<(
                    std::ostream& os, const Vector<ELEMENTTYPE>& vec ) {
                    os << "[ ";
                    for ( size_t r = 0; r < vec.size(); r++ ) {
                        if ( r > 0 ) {
                            os << ", ";
                        }
                        os << vec[ r ];
                    }
                    os << " ]";
                    return os;
                }

                friend Vector<ELEMENTTYPE> operator+(
                    const Vector<ELEMENTTYPE>& c1,
                    const Vector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator+(
                    const Vector<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator+(
                    ELEMENTTYPE c1, const Vector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c2 );
                    out += c1;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator-(
                    const Vector<ELEMENTTYPE>& c1,
                    const Vector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out -= c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator-(
                    const Vector<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
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

                friend Vector<ELEMENTTYPE> operator*(
                    const Vector<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out *= c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator*(
                    ELEMENTTYPE c1, const Vector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c2 );
                    out *= c1;
                    return out;
                }

                virtual details::abstract_vector<ELEMENTTYPE> *internal()
                    const {
                    return _ptr.get();
                }

                virtual void check_pointer() const {
                    if ( internal() == nullptr ) {
                        throw std::logic_error(
                            "Vector: Internal pointer has not been set!" );
                    }
                }

                explicit operator bool() const {
                    return ( _ptr ? true : false );
                }

            protected:
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

        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class WrapperVector<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public Vector<ELEMENTTYPE> {
            public:
                WrapperVector() : Vector<ELEMENTTYPE>(), _ptr { nullptr } {}

                WrapperVector( details::abstract_vector<ELEMENTTYPE>& vec ) {
                    _ptr = &vec;
                }

                WrapperVector( const WrapperVector<ELEMENTTYPE>& other ) :
                    WrapperVector<ELEMENTTYPE>() {
                    _ptr = other._ptr;
                }

                ~WrapperVector() {}

                friend void ::swap<ELEMENTTYPE>(
                    WrapperVector<ELEMENTTYPE>& a,
                    WrapperVector<ELEMENTTYPE>& b ) noexcept;

                /**
                 * Assignment operator.
                 * @param other The vector to assign to this.
                 */
                WrapperVector<ELEMENTTYPE>& operator=(
                    WrapperVector<ELEMENTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend bool operator==( const WrapperVector<ELEMENTTYPE>& a,
                                        const WrapperVector<ELEMENTTYPE>& b ) {
                    return a.equals( b );
                }

                friend bool operator!=( const WrapperVector<ELEMENTTYPE>& a,
                                        const WrapperVector<ELEMENTTYPE>& b ) {
                    return !( a.equals( b ) );
                }

                virtual details::abstract_vector<ELEMENTTYPE> *internal()
                    const override {
                    return _ptr;
                }

                virtual Vector<ELEMENTTYPE>& clear() override {
                    _ptr = nullptr;
                    return *static_cast<Vector<ELEMENTTYPE> *>( this );
                }

                virtual Vector<ELEMENTTYPE>& operator-=(
                    const Vector<ELEMENTTYPE>& other ) override {
                    this->check_pointer();
                    this->_check_same_size( other );
                    for ( size_t i = 0; i < this->size(); i++ ) {
                        internal()->set( i, this->get( i ) - other.get( i ) );
                    }
                    // this->add( -other );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE>& operator-=(
                    const ELEMENTTYPE& other ) override {
                    this->check_pointer();
                    this->add( -other );
                    return *this;
                }

                // friend binary operators
                friend std::ostream& operator<<(
                    std::ostream& os, const WrapperVector<ELEMENTTYPE>& vec ) {
                    os << "[ ";
                    for ( size_t r = 0; r < vec.size(); r++ ) {
                        if ( r > 0 ) {
                            os << ", ";
                        }
                        os << vec[ r ];
                    }
                    os << " ]";
                    return os;
                }

                friend Vector<ELEMENTTYPE> operator+(
                    const WrapperVector<ELEMENTTYPE>& c1,
                    const WrapperVector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator+(
                    const WrapperVector<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator+(
                    ELEMENTTYPE c1, const WrapperVector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c2 );
                    out += c1;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator-(
                    const WrapperVector<ELEMENTTYPE>& c1,
                    const WrapperVector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out -= c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator-(
                    const WrapperVector<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out -= c2;
                    return out;
                }

                friend NCPA::linear::Vector<ELEMENTTYPE> operator*(
                    const WrapperVector<ELEMENTTYPE>& c1,
                    const WrapperVector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out *= c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator*(
                    const WrapperVector<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    Vector<ELEMENTTYPE> out( c1 );
                    out *= c2;
                    return out;
                }

                friend Vector<ELEMENTTYPE> operator*(
                    ELEMENTTYPE c1, const WrapperVector<ELEMENTTYPE>& c2 ) {
                    Vector<ELEMENTTYPE> out( c2 );
                    out *= c1;
                    return out;
                }

                explicit operator bool() const { return ( _ptr != nullptr ); }

            private:
                details::abstract_vector<ELEMENTTYPE> *_ptr;
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::Vector<T>& a,
                  NCPA::linear::Vector<T>& b ) noexcept {
    // using std::swap;
    a._ptr.swap( b._ptr );
}

template<typename T>
static void swap( NCPA::linear::WrapperVector<T>& a,
                  NCPA::linear::WrapperVector<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::Vector<T>&>( a ),
            static_cast<NCPA::linear::Vector<T>&>( b ) );
    swap( a._ptr, b._ptr );
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
