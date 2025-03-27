#pragma once

// #include "NCPA/arrays.hpp"

// #include <array>
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace NCPA {
    namespace arrays {
        template<typename T>
        class _as;

        template<size_t N>
        class _dimarray;

        template<size_t N, typename T>
        class ndvector;
    }  // namespace arrays
}  // namespace NCPA

template<typename T>
static void swap( NCPA::arrays::_as<T>& a, NCPA::arrays::_as<T>& b ) noexcept;

template<size_t N>
static void swap( NCPA::arrays::_dimarray<N>& a,
                  NCPA::arrays::_dimarray<N>& b ) noexcept;

template<size_t N, typename T>
static void swap( NCPA::arrays::ndvector<N, T>& a,
                  NCPA::arrays::ndvector<N, T>& b ) noexcept;

template<typename T>
static void swap( NCPA::arrays::ndvector<1, T>& a,
                  NCPA::arrays::ndvector<1, T>& b ) noexcept;

template<typename T>
static void swap( NCPA::arrays::ndvector<0, T>& a,
                  NCPA::arrays::ndvector<0, T>& b ) noexcept;

namespace NCPA {
    namespace arrays {
        template<typename T>
        class _as {
            public:
                _as() : _val { 0 } {};

                _as( const T& rhs ) : _val( rhs ) {}

                virtual ~_as() {}

                operator T() const { return _val; }

                friend void ::swap( NCPA::arrays::_as<T>& a,
                                    NCPA::arrays::_as<T>& b ) noexcept;


            protected:
                T _val;
        };

        template<size_t N>
        class _dimarray : public std::vector<size_t> {
            public:
                _dimarray() : std::vector<size_t>( N ) {}

                _dimarray( const _dimarray<N>& v ) :
                    std::vector<size_t>( v ) {}

                _dimarray( const std::vector<size_t>& v ) :
                    std::vector<size_t>( v ) {
                    if (v.size() != N) {
                        std::ostringstream oss;
                        oss << "_dimarray(): For dimensionality " << N
                            << ", expected dimension vector of size " << N
                            << ", got size " << v.size();
                        throw std::invalid_argument( oss.str() );
                    }
                }

                _dimarray( _dimarray<N>&& source ) noexcept : _dimarray<N>() {
                    ::swap( *this, source );
                }

                _dimarray<N>& operator=( _dimarray<N> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                _dimarray<N>& operator=( std::vector<size_t> other ) {
                    if (other.size() != this->size()) {
                        throw std::invalid_argument(
                            "_dimarray.operator=(): Dimension vectors must be "
                            "the same size!" );
                    }
                    std::swap( static_cast<std::vector<size_t>&>( *this ),
                               other );
                    return *this;
                }

                virtual ~_dimarray() {}

                friend void ::swap<N>(
                    NCPA::arrays::_dimarray<N>& a,
                    NCPA::arrays::_dimarray<N>& b ) noexcept;

                _dimarray<N - 1> subdims() const {
                    _dimarray<N - 1> sub;
                    for (size_t i = 1; i < this->size(); ++i) {
                        sub[ i - 1 ] = this->at( i );
                    }
                    return sub;
                }
        };

        template<size_t N, typename T>
        class ndvector : public std::vector<ndvector<N - 1, T>> {
                typedef ndvector<N - 1, T> _subvector;

            public:
                using std::vector<ndvector<N - 1, T>>::operator[];

                ndvector() {}

                ndvector( const _dimarray<N>& dims ) {
                    _subvector sv( dims.subdims() );
                    this->assign( dims[ 0 ], sv );
                    // _v = std::vector<_subvector>( dims[ 0 ], sv );
                }

                ndvector( const std::vector<_subvector>& v ) {
                    // _v = v;
                    this->assign( v.cbegin(), v.cend() );
                }

                ndvector( const std::initializer_list<_subvector>& ls ) :
                    ndvector( std::vector<_subvector>( ls ) ) {}

                // ndvector( const ndvector<N, T>& other ) : ndvector<N, T>() {
                //     _v = other._v;
                // }

                // ndvector( ndvector<N, T>&& source ) noexcept :
                //     ndvector<N, T>() {
                //     ::swap( *this, source );
                // }

                // ndvector<N, T>& operator=( ndvector<N, T> other ) {
                //     ::swap( *this, other );
                //     return *this;
                // }

                ndvector<N, T>& operator=(
                    const std::initializer_list<_subvector>& other ) {
                    // _v = other;
                    // std::cout << "Making othervec" << std::endl;
                    ndvector<N, T> othervec( other );
                    // std::cout << "Swapping" << std::endl;
                    ::swap( *this, othervec );
                    return *this;
                }

                virtual ~ndvector() {}

                friend void ::swap<N, T>(
                    NCPA::arrays::ndvector<N, T>& a,
                    NCPA::arrays::ndvector<N, T>& b ) noexcept;

                // _subvector& operator[]( size_t n ) { return _v[ n ]; }

                // const _subvector& operator[]( size_t n ) const {
                //     return _v[ n ];
                // }

                // _subvector& at( size_t n ) { return _v.at( n ); }

                // const _subvector& at( size_t n ) const { return _v.at( n );
                // }

                T& operator[]( const std::initializer_list<size_t> inds ) {
                    return ( *this )[ _dimarray<N>( inds ) ];
                }

                const T& operator[](
                    const std::initializer_list<size_t> inds ) const {
                    return ( *this )[ _dimarray<N>( inds ) ];
                }

                T& operator[]( const _dimarray<N>& inds ) {
                    // if ( inds[ 0 ] > _v.size() ) {
                    //     _dimarray<N> dims = this->shape();
                    //     dims[ 0 ]         = inds[ 0 ];
                    //     this->reshape( dims );
                    // }
                    // return _v[ inds[ 0 ] ][ inds.subdims() ];
                    if (inds[ 0 ] > this->size()) {
                        _dimarray<N> dims = this->shape();
                        dims[ 0 ]         = inds[ 0 ];
                        this->reshape( dims );
                    }
                    return this->at( inds[ 0 ] )[ inds.subdims() ];
                    // return _v[ inds[ 0 ] ][ inds.subdims() ];
                }

                const T& operator[]( const _dimarray<N>& inds ) const {
                    return this->at( inds[ 0 ] )[ inds.subdims() ];
                    // return _v[ inds[ 0 ] ][ inds.subdims() ];
                }

                void set( const T& val ) {
                    for (auto it = this->begin(); it != this->end(); ++it) {
                        it->set( val );
                    }
                    // this->assign( this->size(), val );
                }

                void reshape( const _dimarray<N>& newshape ) {
                    this->resize( newshape[ 0 ] );
                    // _v.resize( newshape[ 0 ] );
                    for (size_t i = 0; i < newshape[ 0 ]; ++i) {
                        this->at( i ).reshape( newshape.subdims() );
                        // _v[ i ].reshape( newshape.subdims() );
                    }
                }

                void reshape( const _dimarray<N>& newshape, const T& val ) {
                    this->resize( newshape[ 0 ] );
                    // _v.resize( newshape[ 0 ] );
                    for (size_t i = 0; i < newshape[ 0 ]; ++i) {
                        // _v[ i ].reshape( newshape.subdims(), val );
                        this->at( i ).reshape( newshape.subdims(), val );
                    }
                }

                void reshape( const std::initializer_list<size_t>& newshape ) {
                    this->reshape( _dimarray<N>( newshape ) );
                }

                void reshape( const std::initializer_list<size_t>& newshape,
                              const T& val ) {
                    this->reshape( _dimarray<N>( newshape ), val );
                }

                _dimarray<N> shape() const {
                    _dimarray<N> dims;
                    if (this->size() > 0) {
                        dims[ 0 ]                 = this->size();
                        // dims[ 0 ]                 = _v.size();
                        _dimarray<N - 1> subshape = this->front().shape();
                        // _dimarray<N - 1> subshape = _v[ 0 ].shape();
                        for (size_t i = 1; i < N; ++i) {
                            dims[ i ] = subshape[ i - 1 ];
                        }
                    }
                    return dims;
                }

                // size_t size() const { return _v.size(); }

                // protected:
                // std::vector<_subvector> _v;
        };

        template<typename T>
        class ndvector<1, T> : public std::vector<T> {
            public:
                using std::vector<T>::operator[];
                using std::vector<T>::size;

                ndvector() : std::vector<T>() {}

                ndvector( size_t n ) : std::vector<T>( n ) {}

                ndvector( const ndvector<1, T>& other ) :
                    std::vector<T>( other ) {
        
                    }

                ndvector( ndvector<1, T>&& source ) noexcept :
                    ndvector<1, T>() {
                    ::swap( *this, source );
                }

                ndvector( const _dimarray<1>& dims ) :
                    ndvector<1, T>( dims[ 0 ] ) {
                    }

                ndvector<1, T>& operator=( ndvector<1, T> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                ndvector<1, T>& operator=( const std::vector<T>& other ) {
                    this->assign( other.cbegin(), other.cend() );
                    return *this;
                }

                virtual ~ndvector() {}

                friend void ::swap<T>(
                    NCPA::arrays::ndvector<1, T>& a,
                    NCPA::arrays::ndvector<1, T>& b ) noexcept;

                void set( const T& val ) { this->assign( this->size(), val ); }

                _dimarray<1> shape() const {
                    _dimarray<1> dims;
                    dims[ 0 ] = this->size();
                    return dims;
                }

                void reshape( const std::vector<size_t> newsize ) {
                    if (newsize.size() != 1) {
                        throw std::invalid_argument(
                            "ndvector<1>.reshape(): dimension vector must be "
                            "size 1" );
                    }
                    this->resize( newsize[ 0 ] );
                }

                void reshape( const std::vector<size_t> newsize,
                              const T& val ) {
                    if (newsize.size() != 1) {
                        throw std::invalid_argument(
                            "ndvector<1>.reshape(): dimension vector must be "
                            "size 1" );
                    }
                    this->assign( newsize[ 0 ], val );
                }

                T& operator[]( size_t n ) {
                    if (n > this->size()) {
                        this->resize( n );
                    }
                    return ( static_cast<std::vector<T>&>( *this )[ n ] );
                }

                const T& operator[]( size_t n ) const {
                    return (
                        static_cast<const std::vector<T>&>( *this )[ n ] );
                }

                T& operator[]( const std::vector<size_t>& inds ) {
                    return ( *this )[ inds[ 0 ] ];
                    // return ( static_cast<std::vector<T>&>( *this )[ inds[0]
                    // ] );
                }

                const T& operator[]( const std::vector<size_t>& inds ) const {
                    return ( *this )[ inds[ 0 ] ];
                    // return ( static_cast<const std::vector<T>&>( *this )[
                    // inds[0] ] );
                }
        };

        template<typename T>
        class ndvector<0, T> : public _as<T> {
            public:
                ndvector() : _as<T>() {}

                ndvector( T v ) : _as<T>( v ) {}

                ndvector( const ndvector<0, T>& other ) :
                    _as<T>( other._val ) {}

                ndvector( ndvector<1, T>&& source ) noexcept :
                    ndvector<0, T>() {
                    ::swap( *this, source );
                }

                ndvector<0, T>& operator=( ndvector<0, T> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                ndvector<0, T>& operator=( T other ) {
                    this->_val = other;
                    return *this;
                }

                virtual ~ndvector() {}

                friend void ::swap<T>(
                    NCPA::arrays::ndvector<0, T>& a,
                    NCPA::arrays::ndvector<0, T>& b ) noexcept;
        };

    }  // namespace arrays
}  // namespace NCPA

template<typename T>
static void swap( NCPA::arrays::_as<T>& a, NCPA::arrays::_as<T>& b ) noexcept {
    using std::swap;
    swap( a._val, b._val );
}

template<size_t N>
static void swap( NCPA::arrays::_dimarray<N>& a,
                  NCPA::arrays::_dimarray<N>& b ) noexcept {
    using std::swap;
    swap( dynamic_cast<std::vector<size_t>&>( a ),
          dynamic_cast<std::vector<size_t>&>( b ) );
}

template<size_t N, typename T>
static void swap( NCPA::arrays::ndvector<N, T>& a,
                  NCPA::arrays::ndvector<N, T>& b ) noexcept {
    using std::swap;
    swap( dynamic_cast<std::vector<NCPA::arrays::ndvector<N - 1, T>>&>( a ),
          dynamic_cast<std::vector<NCPA::arrays::ndvector<N - 1, T>>&>( b ) );

    // swap( a._v, b._v );
}

template<typename T>
static void swap( NCPA::arrays::ndvector<1, T>& a,
                  NCPA::arrays::ndvector<1, T>& b ) noexcept {
    using std::swap;
    swap( dynamic_cast<std::vector<T>&>( a ),
          dynamic_cast<std::vector<T>&>( b ) );
}

template<typename T>
static void swap( NCPA::arrays::ndvector<0, T>& a,
                  NCPA::arrays::ndvector<0, T>& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::arrays::_as<T>&>( a ),
            dynamic_cast<NCPA::arrays::_as<T>&>( b ) );
}
