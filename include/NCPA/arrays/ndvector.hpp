#pragma once

#include "NCPA/arrays/declarations.hpp"
#include "NCPA/arrays/exceptions.hpp"

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
        namespace details {
            template<typename T>
            class _as {
                public:
                    _as() : _val { 0 } {};

                    _as( const T& rhs ) : _val( rhs ) {}

                    virtual ~_as() {}

                    operator T() const { return _val; }

                    friend void swap( _as<T>& a, _as<T>& b ) noexcept {
                        using std::swap;
                        swap( a._val, b._val );
                    }

                protected:
                    T _val;
            };

            // template<size_t N>
            // using _dimarray = std::array<size_t, N>;
            typedef std::vector<size_t> _dimarray;

            // template<size_t N>
            static inline _dimarray subdims( _dimarray dims ) {
                _dimarray sub;
                if (dims.size() > 1) {
                    sub.assign( dims.begin() + 1, dims.end() );
                }
                return sub;
            }

            // template<size_t N>
            // class _dimarray : public std::vector<size_t> {
            //     public:
            //         _dimarray() : std::vector<size_t>( N ) {}

            //         _dimarray( const _dimarray& v ) :
            //             std::vector<size_t>( v ) {}

            //         _dimarray( const std::vector<size_t>& v ) :
            //             std::vector<size_t>( v ) {
            //             if (v.size() != N) {
            //                 std::ostringstream oss;
            //                 oss << "_dimarray(): For dimensionality " << N
            //                     << ", expected dimension vector of size " <<
            //                     N
            //                     << ", got size " << v.size();
            //                 throw std::invalid_argument( oss.str() );
            //             }
            //         }

            //         _dimarray( _dimarray&& source ) noexcept :
            //             _dimarray() {
            //             swap( *this, source );
            //         }

            //         _dimarray& operator=( _dimarray other ) {
            //             swap( *this, other );
            //             return *this;
            //         }

            //         _dimarray& operator=( std::vector<size_t> other ) {
            //             if (other.size() != this->size()) {
            //                 throw std::invalid_argument(
            //                     "_dimarray.operator=(): Dimension vectors "
            //                     "must be "
            //                     "the same size!" );
            //             }
            //             std::swap( static_cast<std::vector<size_t>&>( *this
            //             ),
            //                        other );
            //             return *this;
            //         }

            //         virtual ~_dimarray() {}

            //         friend void swap( _dimarray& a,
            //                           _dimarray& b ) noexcept {
            //             using std::swap;
            //             swap( dynamic_cast<std::vector<size_t>&>( a ),
            //                   dynamic_cast<std::vector<size_t>&>( b ) );
            //         }

            //         _dimarray subdims() const {
            //             _dimarray sub;
            //             for (size_t i = 1; i < this->size(); ++i) {
            //                 sub[ i - 1 ] = this->at( i );
            //             }
            //             return sub;
            //         }
            // };
        }  // namespace details

        template<size_t N, typename T>
        class ndvector : public std::vector<ndvector<N - 1, T>> {
                typedef ndvector<N - 1, T> _subvector;

            public:
                using std::vector<ndvector<N - 1, T>>::operator[];

                ndvector() {}

                ndvector( const details::_dimarray& dims ) {
                    _subvector sv( details::subdims( dims ) );
                    // _subvector sv( dims.subdims() );
                    this->assign( dims[ 0 ], sv );
                }

                ndvector( const std::initializer_list<size_t>& dims ) :
                    ndvector<N, T>( details::_dimarray { dims } ) {}

                ndvector( const std::vector<_subvector>& v ) {
                    // _v = v;
                    this->assign( v.cbegin(), v.cend() );
                }

                ndvector( const std::initializer_list<_subvector>& ls ) :
                    ndvector( std::vector<_subvector>( ls ) ) {}

                ndvector( const ndvector<N, T>& other ) :
                    std::vector<_subvector>( other ) {}

                template<typename U>
                ndvector( const ndvector<N, U>& other ) :
                    std::vector<_subvector>( other.size() ) {
                    for (size_t i = 0; i < this->size(); ++i) {
                        this->at( i ) = _subvector( other.at( i ) );
                    }
                }

                ndvector<N, T>& operator=(
                    const std::initializer_list<_subvector>& other ) {
                    ndvector<N, T> othervec( other );
                    swap( *this, othervec );
                    return *this;
                }

                virtual ~ndvector() {}

                friend void swap( ndvector<N, T>& a,
                                  ndvector<N, T>& b ) noexcept {
                    using std::swap;
                    swap(
                        dynamic_cast<std::vector<ndvector<N - 1, T>>&>( a ),
                        dynamic_cast<std::vector<ndvector<N - 1, T>>&>( b ) );
                    // swap( a._v, b._v );
                }

                T& operator[]( const std::initializer_list<size_t> inds ) {
                    return ( *this )[ details::_dimarray( inds ) ];
                }

                const T& operator[](
                    const std::initializer_list<size_t> inds ) const {
                    return ( *this )[ details::_dimarray( inds ) ];
                }

                T& operator[]( const details::_dimarray& inds ) {
                    if (inds[ 0 ] > this->size()) {
                        details::_dimarray dims = this->shape();
                        dims[ 0 ]               = inds[ 0 ];
                        this->reshape( dims );
                    }
                    return this->at( inds[ 0 ] )[ details::subdims( inds ) ];
                    // return this->at( inds[ 0 ] )[ inds.subdims() ];
                }

                const T& operator[]( const details::_dimarray& inds ) const {
                    return this->at( inds[ 0 ] )[ details::subdims( inds ) ];
                }

                _subvector& operator[]( size_t n ) { return this->at( n ); }

                const _subvector& operator[]( size_t n ) const {
                    return this->at( n );
                }

                void set( const T& val ) {
                    for (auto it = this->begin(); it != this->end(); ++it) {
                        it->set( val );
                    }
                    // this->assign( this->size(), val );
                }

                void reshape( const details::_dimarray& newshape ) {
                    this->resize( newshape[ 0 ] );
                    // _v.resize( newshape[ 0 ] );
                    for (size_t i = 0; i < newshape[ 0 ]; ++i) {
                        this->at( i ).reshape( details::subdims( newshape ) );
                        // _v[ i ].reshape( newshape.subdims() );
                    }
                }

                void reshape( const details::_dimarray& newshape,
                              const T& val ) {
                    this->resize( newshape[ 0 ] );
                    // _v.resize( newshape[ 0 ] );
                    for (size_t i = 0; i < newshape[ 0 ]; ++i) {
                        // _v[ i ].reshape( newshape.subdims(), val );
                        this->at( i ).reshape( details::subdims( newshape ),
                                               val );
                    }
                }

                void reshape( const std::initializer_list<size_t>& newshape ) {
                    this->reshape( details::_dimarray { newshape } );
                }

                void reshape( const std::initializer_list<size_t>& newshape,
                              const T& val ) {
                    this->reshape( details::_dimarray { newshape }, val );
                }

                details::_dimarray shape() const {
                    details::_dimarray dims( N, 0 );
                    if (this->size() > 0) {
                        dims[ 0 ]                   = this->size();
                        // dims[ 0 ]                 = _v.size();
                        details::_dimarray subshape = this->front().shape();
                        // details::_dimarray subshape = _v[ 0
                        // ].shape();
                        for (size_t i = 1; i < N; ++i) {
                            dims[ i ] = subshape[ i - 1 ];
                        }
                    } 
                    return dims;
                }

                void copy( const ndvector<N, T>& other ) {
                    *this = ndvector<N, T>( other );
                }

                template<typename U>
                void copy( const ndvector<N, U>& other ) {
                    this->resize( other.size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        this->at( i ).copy( other.at( i ) );
                    }
                }
        };

        template<typename T>
        class ndvector<1, T> : public std::vector<T> {
            public:
                using std::vector<T>::operator[];
                using std::vector<T>::size;

                ndvector() : std::vector<T>() {}

                ndvector( size_t n ) : std::vector<T>( n ) {}

                ndvector( const ndvector<1, T>& other ) :
                    std::vector<T>( other ) {}

                template<typename U>
                ndvector( const ndvector<1, U>& other ) :
                    std::vector<T>( other.size() ) {
                    for (size_t i = 0; i < other.size(); ++i) {
                        this->at( i ) = static_cast<T>( other.at( i ) );
                    }
                }

                ndvector( ndvector<1, T>&& source ) noexcept :
                    ndvector<1, T>() {
                    swap( *this, source );
                }

                ndvector( const details::_dimarray& dims ) :
                    ndvector<1, T>( dims[ 0 ] ) {}

                ndvector( const std::initializer_list<size_t>& dims ) :
                    ndvector<1, T>( details::_dimarray { *dims.begin() } ) {}

                ndvector<1, T>& operator=( ndvector<1, T> other ) {
                    swap( *this, other );
                    return *this;
                }

                ndvector<1, T>& operator=( const std::vector<T>& other ) {
                    this->assign( other.cbegin(), other.cend() );
                    return *this;
                }

                virtual ~ndvector() {}

                friend void swap( ndvector<1, T>& a,
                                  ndvector<1, T>& b ) noexcept {
                    using std::swap;
                    swap( dynamic_cast<std::vector<T>&>( a ),
                          dynamic_cast<std::vector<T>&>( b ) );
                }

                void set( const T& val ) { this->assign( this->size(), val ); }

                details::_dimarray shape() const {
                    return details::_dimarray { this->size() };
                    // dims[ 0 ] = this->size();
                    // return dims;
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
                    if (n >= this->size()) {
                        this->resize( n + 1 );
                    }
                    return ( static_cast<std::vector<T>&>( *this ).at( n ) );
                }

                const T& operator[]( size_t n ) const {
                    return (
                        static_cast<const std::vector<T>&>( *this ).at( n ) );
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

                void copy( const ndvector<1, T>& other ) {
                    *this = ndvector<1, T>( other );
                }

                template<typename U>
                void copy( const ndvector<1, U>& other ) {
                    this->resize( other.size() );
                    for (size_t i = 0; i < this->size(); ++i) {
                        this->at( i ) = static_cast<T>( other.at( i ) );
                    }
                }
        };

        template<typename T>
        class ndvector<0, T> : public details::_as<T> {
            public:
                ndvector() : details::_as<T>() {}

                ndvector( T v ) : details::_as<T>( v ) {}

                ndvector( const ndvector<0, T>& other ) :
                    details::_as<T>( other._val ) {}

                template<typename U>
                ndvector( const ndvector<0, U>& other ) :
                    details::_as<T>( static_cast<T>( other._val ) ) {}

                ndvector( ndvector<1, T>&& source ) noexcept :
                    ndvector<0, T>() {
                    swap( *this, source );
                }

                ndvector<0, T>& operator=( ndvector<0, T> other ) {
                    swap( *this, other );
                    return *this;
                }

                ndvector<0, T>& operator=( T other ) {
                    this->_val = other;
                    return *this;
                }

                virtual ~ndvector() {}

                friend void swap( NCPA::arrays::ndvector<0, T>& a,
                                  NCPA::arrays::ndvector<0, T>& b ) noexcept {
                    using std::swap;
                    swap( dynamic_cast<NCPA::arrays::details::_as<T>&>( a ),
                          dynamic_cast<NCPA::arrays::details::_as<T>&>( b ) );
                }

                void copy( const ndvector<0, T>& other ) {
                    this->_val = other._val;
                }

                template<typename U>
                void copy( const ndvector<0, U>& other ) {
                    this->_val = static_cast<T>( other._val );
                }
        };

        /**
         * @class vector2d_t
         * @brief A resizable 2D vector class.
         * @tparam T The type of data stored in the vector.
         */
        template<typename T>
        class vector2d_t : public ndvector<2, T> {
            public:
                using ndvector<2, T>::set;

                /**
                 * @brief Default constructor.
                 */
                vector2d_t() : ndvector<2, T>() {}

                /**
                 * @brief Constructor that initializes dimensions and optional
                 * default value.
                 * @param nx1 Size of the first dimension.
                 * @param nx2 Size of the second dimension.
                 * @param val Initial value for all elements.
                 */
                vector2d_t( size_t nx1, size_t nx2, const T& val = (T)0 ) :
                    ndvector<2, T>() {
                    this->resize2d( nx1, nx2, val );
                }

                /**
                 * @brief Constructor that initializes from a raw 2D array.
                 * @param nx1 Size of the first dimension.
                 * @param nx2 Size of the second dimension.
                 * @param vals Pointer to pointer of values to copy.
                 */
                vector2d_t( size_t nx1, size_t nx2, const T **vals ) :
                    vector2d_t( nx1, nx2 ) {
                    for (size_t i = 0; i < nx1; ++i) {
                        this->at( i ).assign( vals[ i ], vals[ i ] + nx2 );
                    }
                }

                /**
                 * @brief Resizes the 2D vector.
                 * @param nx1 New size for dimension 1.
                 * @param nx2 New size for dimension 2.
                 * @param val Value to fill in new elements.
                 */
                virtual void resize2d( size_t nx1, size_t nx2,
                                       const T& val = (T)0 ) {
                    this->reshape( { nx1, nx2 }, val );
                }

                /**
                 * @brief Retrieves the current sizes of the two dimensions.
                 * @param nx1 Reference to hold dimension 1 size.
                 * @param nx2 Reference to hold dimension 2 size.
                 */
                virtual void size2d( size_t& nx1, size_t& nx2 ) const {
                    auto dims = this->shape();
                    nx1       = dims[ 0 ];
                    nx2       = dims[ 1 ];
                }

                /**
                 * @brief Returns the size of a specific dimension.
                 * @param dimnum Index of the dimension (0 or 1).
                 * @return The size of the requested dimension.
                 */
                virtual size_t dim( size_t dimnum ) const {
                    return this->shape().at( dimnum );
                }

                /**
                 * @brief Fills the entire vector with a constant value.
                 * @param val The value to fill.
                 */
                virtual void fill( T val ) { this->set( val ); }
        };

        /**
         * @class vector3d_t
         * @brief A resizable 3D vector class.
         * @tparam T The type of data stored in the vector.
         */
        template<typename T>
        class vector3d_t : public ndvector<3, T> {
            public:
                /**
                 * @brief Default constructor.
                 */
                vector3d_t() : ndvector<3, T>() {}

                /**
                 * @brief Constructor that initializes dimensions and optional
                 * default value.
                 * @param nx1 Size of the first dimension.
                 * @param nx2 Size of the second dimension.
                 * @param nx3 Size of the third dimension.
                 * @param val Initial value for all elements.
                 */
                vector3d_t( size_t nx1, size_t nx2, size_t nx3,
                            const T& val = (T)0.0 ) :
                    ndvector<3, T>() {
                    this->resize3d( nx1, nx2, nx3, val );
                }

                /**
                 * @brief Constructor that initializes from a raw 3D array.
                 * @param nx1 Size of the first dimension.
                 * @param nx2 Size of the second dimension.
                 * @param nx3 Size of the third dimension.
                 * @param vals Pointer to 3D array of values.
                 */
                vector3d_t( size_t nx1, size_t nx2, size_t nx3,
                            const T ***vals ) :
                    vector3d_t( nx1, nx2, nx3 ) {
                    for (size_t i = 0; i < nx1; ++i) {
                        for (size_t j = 0; j < nx2; ++j) {
                            this->at( i ).at( j ).assign(
                                vals[ i ][ j ], vals[ i ][ j ] + nx3 );
                        }
                    }
                }

                /**
                 * @brief Copy constructor.
                 * @param other The vector to copy from.
                 */
                vector3d_t( const vector3d_t<T>& other ) :
                    ndvector<3, T>( other ) {}

                /**
                 * @brief Conversion constructor from an ndvector.
                 * @param other The base ndvector to copy.
                 */
                vector3d_t( const ndvector<3, T>& other ) :
                    ndvector<3, T>( other ) {}

                /**
                 * @brief Resizes the 3D vector.
                 * @param nx1 New size for dimension 1.
                 * @param nx2 New size for dimension 2.
                 * @param nx3 New size for dimension 3.
                 * @param val Value to fill in new elements.
                 */
                virtual void resize3d( size_t nx1, size_t nx2, size_t nx3,
                                       const T& val = (T)0.0 ) {
                    this->reshape( { nx1, nx2, nx3 }, val );
                }

                /**
                 * @brief Retrieves the current sizes of the three dimensions.
                 * @param nx1 Reference to hold dimension 1 size.
                 * @param nx2 Reference to hold dimension 2 size.
                 * @param nx3 Reference to hold dimension 3 size.
                 */
                virtual void size3d( size_t& nx1, size_t& nx2,
                                     size_t& nx3 ) const {
                    auto dims = this->shape();
                    nx1       = dims[ 0 ];
                    nx2       = dims[ 1 ];
                    nx3       = dims[ 2 ];
                }

                /**
                 * @brief Returns the size of a specific dimension.
                 * @param dimnum Index of the dimension (0, 1, or 2).
                 * @return The size of the requested dimension.
                 */
                virtual size_t dim( size_t dimnum ) const {
                    return this->shape()[ dimnum ];
                }
        };

    }  // namespace arrays
}  // namespace NCPA
