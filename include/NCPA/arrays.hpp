/**
 * @file arrays.hpp
 * @author Claus Hetzer
 * @date April 16, 2026
 * @version 1.0.0
 * @brief Header file providing N-dimensional array structures, views, and
 * utility functions.
 *
 * This file contains the implementation of multidimensional array
 * containers, including resizable vectors (2D/3D), array views for slicing,
 * and a variety of linear algebra and permutation utilities.
 */
#pragma once

#include "NCPA/constants.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/types.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

/**
 * @namespace NCPA
 * @brief The primary namespace for NCPA software components.
 */
namespace NCPA {
    namespace arrays {
        template<typename>
        class NDVector;
        template<typename>
        class NDVectorView;
        template<size_t, typename>
        class NDArray;
        template<size_t, typename>
        class NDArrayView;

        class InvalidRank : public std::invalid_argument {
            public:
                InvalidRank( const std::string& msg = "" ) :
                    std::invalid_argument( std::string( "InvalidRank" ) + ": "
                                           + msg ) {}
        };

        template<typename T>
        bool operator==( const std::vector<T>& a, const std::vector<T>& b ) {
            if (a.size() != b.size()) return false;
            auto ait = a.begin();
            auto bit = b.begin();
            while (ait != a.end() && bit != b.end()) {
                if (*ait != *bit) return false;
                ++ait;
                ++bit;
            }
            return true;
        }

        template<typename T>
        bool operator!=( const std::vector<T>& a, const std::vector<T>& b ) {
            return !( a == b );
        }

        namespace details {
            template<typename T>
            struct _ndstorage {
                    std::vector<T> data;
            };

            template<typename Shape>
            size_t total_size( const Shape& shape ) {
                size_t n = 1;
                for (size_t d : shape) {
                    n *= d;
                }
                return n;
            }

            template<typename Shape>
            void compute_strides( const Shape& shape, Shape& strides ) {
                strides                     = shape;
                strides[ shape.size() - 1 ] = 1;

                for (int i = static_cast<int>( shape.size() ) - 2; i >= 0;
                     --i) {
                    strides[ i ] = strides[ i + 1 ] * shape[ i + 1 ];
                }
            }

            template<typename ShapeA, typename ShapeB>
            size_t linear_index( const ShapeA& indices, const ShapeB& shape,
                                 const ShapeB& strides, size_t offset = 0 ) {
                size_t result = offset;

                if (indices.size() != shape.size()) {
                    throw std::invalid_argument( "wrong number of indices" );
                }

                for (size_t i = 0; i < shape.size(); ++i) {
                    if (indices[ i ] >= shape[ i ]) {
                        throw std::out_of_range( "index" );
                    }

                    result += indices[ i ] * strides[ i ];
                }

                return result;
            }

            template<typename Shape>
            size_t storage_index( size_t linear_pos, size_t offset,
                                  const Shape& shape, const Shape& strides ) {
                size_t idx = offset;
                size_t rem = linear_pos;

                for (int d = static_cast<int>( shape.size() ) - 1; d >= 0;
                     --d) {
                    size_t coord  = rem % shape[ d ];
                    rem          /= shape[ d ];

                    idx += coord * strides[ d ];
                }

                return idx;
            }

            template<typename ShapeIn, typename ShapeOut>
            void remove_dimension( const ShapeIn& shape,
                                   const ShapeIn& strides, size_t dim,
                                   ShapeOut& new_shape,
                                   ShapeOut& new_strides ) {
                size_t out = 0;

                for (size_t i = 0; i < shape.size(); ++i) {
                    if (i == dim) continue;

                    new_shape[ out ]   = shape[ i ];
                    new_strides[ out ] = strides[ i ];
                    ++out;
                }
            }

            template<typename ContainerA, typename ContainerB>
            void deal_into( const ContainerA& a, ContainerB& b ) {
                if (a.size() == b.size()) {
                    auto ait = a.begin();
                    auto bit = b.begin();
                    while (ait != a.end() && bit != b.end()) {
                        *bit = *ait;
                        ++ait;
                        ++bit;
                    }
                } else {
                    throw std::invalid_argument(
                        "Cannot deal N-D containers of different sizes" );
                }
            }

            template<typename ContainerA, typename ContainerB>
            void copy_to( const ContainerA& a, ContainerB& b ) {
                if (a.rank() == b.rank() || b.rank_mutable()) {
                    b.resize( a.shape() );
                    deal_into( a, b );
                } else {
                    throw std::invalid_argument(
                        "Cannot copy N-D containers of different ranks" );
                }
            }

            template<typename ContainerA, typename ContainerB>
            bool same_shape( const ContainerA& a, const ContainerB& b ) {
                if (a.rank() != b.rank()) return false;
                for (auto i = 0; i < a.rank(); ++i) {
                    if (a.dim(i) != b.dim(i)) return false;
                }
                return true;
            }

            template<typename Container>
            std::vector<typename Container::shape_type> indices(const Container& c) {
                auto shape = c.shape();
                auto digits = shape;
                digits.clear();
                std::vector<typename Container::shape_type> out;
                out.reserve( c.size() );
                do {
                    out.push_back( digits );
                }
                while (increment_digits( digits, shape ));
                return out;
            }
        }  // namespace details

        //***************************************************************************/
        // Definitions
        //***************************************************************************/
        /**
         * @brief Returns a pointer to the underlying data of a std::vector.
         * @tparam T Element type.
         * @param in The input vector.
         * @return Pointer to the first element.
         */
        template<typename T>
        T *as_array( std::vector<T>& in ) {
            return &in[ 0 ];
        }

        /**
         * @brief Returns a const pointer to the underlying data of a
         * std::vector.
         * @tparam T Element type.
         * @param in The input vector.
         * @return Const pointer to the first element.
         */
        template<typename T>
        const T *as_array( const std::vector<T>& in ) {
            return &in[ 0 ];
        }

        template<typename T, typename ShapeContainer>
        class NDViewIterator {
            public:
                typedef std::random_access_iterator_tag iterator_category;
                typedef T value_type;
                typedef std::ptrdiff_t difference_type;
                typedef T *pointer_t;
                typedef T& reference_t;

                template<typename, typename>
                friend class NDViewIterator;

                NDViewIterator( std::shared_ptr<details::_ndstorage<
                                    typename std::remove_const<T>::type>>
                                    storage,
                                size_t offset, const ShapeContainer& shape,
                                const ShapeContainer& strides,
                                size_t linear_pos ) :
                    _storage( storage ),
                    _offset( offset ),
                    _shape( &shape ),
                    _strides( &strides ),
                    _linear_pos( linear_pos ) {}

                template<
                    typename U,
                    typename std::enable_if<
                        std::is_const<T>::value
                            && std::is_same<
                                typename std::remove_const<T>::type, U>::value,
                        int>::type = 0>
                NDViewIterator(
                    const NDViewIterator<U, ShapeContainer>& other ) :
                    _storage( other._storage ),
                    _offset( other._offset ),
                    _shape( other._shape ),
                    _strides( other._strides ),
                    _linear_pos( other._linear_pos ) {}

                reference_t operator*() const {
                    return _storage->data[ storage_index() ];
                }

                pointer_t operator->() const {
                    return &_storage->data[ storage_index() ];
                }

                NDViewIterator& operator++() {
                    ++_linear_pos;
                    return *this;
                }

                NDViewIterator operator++( int ) {
                    NDViewIterator tmp( *this );
                    ++( *this );
                    return tmp;
                }

                NDViewIterator& operator--() {
                    --_linear_pos;
                    return *this;
                }

                NDViewIterator operator--( int ) {
                    NDViewIterator tmp( *this );
                    --( *this );
                    return tmp;
                }

                bool operator==( const NDViewIterator& rhs ) const {
                    return _storage.get() == rhs._storage.get()
                        && _offset == rhs._offset
                        && _linear_pos == rhs._linear_pos;
                }

                bool operator!=( const NDViewIterator& rhs ) const {
                    return !( *this == rhs );
                }

                template<typename U>
                bool operator==(
                    const NDViewIterator<U, ShapeContainer>& rhs ) const {
                    return _storage.get() == rhs._storage.get()
                        && _offset == rhs._offset
                        && _linear_pos == rhs._linear_pos;
                }

                template<typename U>
                bool operator!=(
                    const NDViewIterator<U, ShapeContainer>& rhs ) const {
                    return !( *this == rhs );
                }

                template<typename U>
                bool operator<(
                    const NDViewIterator<U, ShapeContainer>& rhs ) const {
                    assert( _storage.get() == rhs._storage.get() );
                    return _linear_pos < rhs._linear_pos;
                }

                template<typename U>
                bool operator>(
                    const NDViewIterator<U, ShapeContainer>& rhs ) const {
                    assert( _storage.get() == rhs._storage.get() );
                    return rhs < *this;
                }

                template<typename U>
                bool operator<=(
                    const NDViewIterator<U, ShapeContainer>& rhs ) const {
                    assert( _storage.get() == rhs._storage.get() );
                    return !( rhs < *this );
                }

                template<typename U>
                bool operator>=(
                    const NDViewIterator<U, ShapeContainer>& rhs ) const {
                    assert( _storage.get() == rhs._storage.get() );
                    return !( *this < rhs );
                }

                template<typename U>
                difference_type operator-(
                    const NDViewIterator<U, ShapeContainer>& rhs ) const {
                    return static_cast<difference_type>( _linear_pos )
                         - static_cast<difference_type>( rhs._linear_pos );
                }

                reference_t operator[]( difference_type n ) const {
                    return *( *this + n );
                }

                NDViewIterator& operator+=( difference_type n ) {
                    _linear_pos += n;
                    return *this;
                }

                NDViewIterator& operator-=( difference_type n ) {
                    _linear_pos -= n;
                    return *this;
                }

                NDViewIterator operator+( difference_type n ) const {
                    NDViewIterator tmp( *this );
                    tmp += n;
                    return tmp;
                }

                NDViewIterator operator-( difference_type n ) const {
                    NDViewIterator tmp( *this );
                    tmp -= n;
                    return tmp;
                }

            private:
                std::shared_ptr<
                    details::_ndstorage<typename std::remove_const<T>::type>>
                    _storage;

                difference_type _offset;
                const ShapeContainer *_shape;
                const ShapeContainer *_strides;
                size_t _linear_pos;

                size_t storage_index() const {
                    return details::storage_index( _linear_pos, _offset,
                                                   *_shape, *_strides );
                }
        };

        // template<typename T, typename ShapeType, typename ViewType,
        //          typename IteratorType, typename ConstIteratorType,
        //          typename ReverseIteratorType,
        //          typename ConstReverseIteratorType, typename DerivedType>
        // class NDBase {
        //         typedef NDBase<T, ShapeType, ViewType, IteratorType,
        //                        ConstIteratorType, ReverseIteratorType,
        //                        ConstReverseIteratorType, DerivedType>
        //             base_t;

        //     public:
        //         NDBase()                      = default;
        //         NDBase( const base_t& other ) = default;
        //         virtual ~NDBase()             = default;

        //         virtual size_t rank() const                              =
        //         0; virtual bool rank_mutable() const = 0; virtual size_t
        //         size() const                              = 0; virtual
        //         ViewType slice( size_t dim, size_t fixed_index ) = 0;
        //         virtual IteratorType begin()                             =
        //         0; virtual ConstIteratorType begin() const = 0; virtual
        //         ConstIteratorType cbegin() const                 = 0;
        //         virtual IteratorType end()                               =
        //         0; virtual ConstIteratorType end() const = 0; virtual
        //         ConstIteratorType cend() const                   = 0;
        //         virtual ReverseIteratorType rbegin()                     =
        //         0; virtual ConstReverseIteratorType rbegin() const = 0;
        //         virtual ReverseIteratorType rend()                       =
        //         0; virtual ConstReverseIteratorType rend() const = 0;

        //         DerivedType& clear() { return this->set( 0 ); }

        //         size_t dim( size_t n ) const {
        //             assert( n < rank() );
        //             return this->shape().at( n );
        //         }

        //         DerivedType& fill( T val ) { return this->set( val ); }

        //         std::vector<T> flatten() const {
        //             std::vector<T> target;
        //             std::copy( this->begin(), this->end(),
        //                        std::back_inserter( target ) );
        //             return target;
        //         }

        //         ViewType operator[]( size_t ind ) {
        //             return this->slice( 0, ind );
        //         }

        //         DerivedType& set( T val ) {
        //             for (auto it = begin(); it != end(); ++it) {
        //                 *it = val;
        //             }
        //             return static_cast<Derived&>( *this );
        //         }

        //         const ShapeType& shape() const { return _shape; }


        //     protected:
        //         std::shared_ptr<details::_ndstorage<T>> _storage;
        //         ShapeType _shape;
        //         ShapeType _strides;
        // };

        template<size_t N, typename T>
        class NDArray {
            public:
                typedef typename std::vector<T>::iterator iterator;
                typedef typename std::vector<T>::const_iterator const_iterator;
                typedef
                    typename std::vector<T>::reverse_iterator reverse_iterator;
                typedef typename std::vector<T>::const_reverse_iterator
                    const_reverse_iterator;
                typedef std::array<size_t, N> shape_type;
                typedef T element_t;

                template<typename... Dims>
                explicit NDArray( Dims... dims ) {
                    this->resize( dims... );
                }

                NDArray( std::vector<size_t> dims ) { this->resize( dims ); }

                NDArray( std::array<size_t, N> dims ) { this->resize( dims ); }

                NDArray( std::initializer_list<size_t> dims ) {
                    this->resize( dims );
                }

                NDArray( const NDArray<N, T>& other ) = default;

                NDArray( const NDArrayView<N, T>& other ) {
                    this->copy( other );
                }

                NDArray( const NDVector<T>& other ) { this->copy( other ); }

                NDArray( const NDVectorView<T>& other ) {
                    this->copy( other );
                }

                template<typename FromType>
                void copy( const FromType& other ) {
                    details::copy_to( other, *this );
                }

                template<typename... Indices>
                typename std::enable_if<( sizeof...( Indices ) > 1 ), T&>::type
                    at( Indices... indices ) {
                    static_assert( sizeof...( Indices ) == N,
                                   "Wrong number of indices" );

                    std::array<size_t, N> idx
                        = { { static_cast<size_t>( indices )... } };

                    return _storage->data.at(
                        details::linear_index( idx, _shape, _strides ) );
                }

                template<typename... Indices>
                typename std::enable_if<( sizeof...( Indices ) > 1 ),
                                        const T&>::type
                    at( Indices... indices ) const {
                    static_assert( sizeof...( Indices ) == N,
                                   "Wrong number of indices" );

                    std::array<size_t, N> idx
                        = { { static_cast<size_t>( indices )... } };

                    return _storage->data.at(
                        details::linear_index( idx, _shape, _strides ) );
                }

                template<size_t M = N>
                typename std::enable_if<M == 1, T&>::type at( size_t idx ) {
                    return _storage->data.at( details::linear_index(
                        shape_type { idx }, _shape, _strides ) );
                }

                template<size_t M = N>
                typename std::enable_if<M == 1, const T&>::type at(
                    size_t idx ) const {
                    return _storage->data.at( details::linear_index(
                        shape_type { idx }, _shape, _strides ) );
                }

                template<typename IntType>
                T& at( std::vector<IntType> idx ) {
                    // shape_type idx;
                    // for (auto i = 0; i < rank(); ++i) {
                    //     idx[ i ] = inds[ i ];
                    // }
                    if (idx.size() != rank()) {
                        throw InvalidRank(
                            "Array rank does not match number of indices" );
                    }
                    return _storage->data.at(
                        details::linear_index( idx, _shape, _strides ) );
                }

                template<typename IntType>
                const T& at( std::vector<IntType> idx ) const {
                    // shape_type idx;
                    // for (auto i = 0; i < rank(); ++i) {
                    //     idx[ i ] = inds[ i ];
                    // }
                    if (idx.size() != rank()) {
                        throw InvalidRank(
                            "Array rank does not match number of indices" );
                    }
                    return _storage->data.at(
                        details::linear_index( idx, _shape, _strides ) );
                }

                template<typename IntType>
                T& at( std::array<IntType, N> idx ) {
                    return _storage->data.at(
                        details::linear_index( idx, _shape, _strides ) );
                }

                template<typename IntType>
                const T& at( std::array<IntType, N> idx ) const {
                    return _storage->data.at(
                        details::linear_index( idx, _shape, _strides ) );
                }

                NDArray<N, T>& clear() { return this->set( 0 ); }

                size_t dim( size_t n ) const {
                    assert( n < rank() );
                    return _shape[ n ];
                }

                NDArray<N, T>& fill( T val ) { return this->set( val ); }

                std::vector<T> flatten() const {
                    std::vector<T> target;
                    std::copy( this->begin(), this->end(),
                               std::back_inserter( target ) );
                    return target;
                }

                template<typename... Indices>
                T get( Indices... indices ) const {
                    return this->at( indices... );
                }

                std::vector<shape_type> indices() const {
                    return details::indices(*this);
                }

                NDArrayView<N - 1, T> operator[]( size_t ind ) {
                    return this->slice( 0, ind );
                }

                template<typename... Indices>
                T& operator()( Indices... indices ) {
                    return this->at( indices... );
                }

                template<typename... Indices>
                const T& operator()( Indices... indices ) const {
                    return this->at( indices... );
                }

                size_t rank() const { return N; }

                bool rank_mutable() const { return false; }

                template<typename... Dims>
                NDArray<N, T>& reshape( Dims... dims ) {
                    return this->resize( dims... );
                }

                template<typename... Dims>
                NDArray<N, T>& resize( Dims... dims ) {
                    static_assert( sizeof...( Dims ) == N
                                       || sizeof...( Dims ) == 0,
                                   "Wrong number of dimensions" );
                    return this->resize( std::vector<size_t> {
                        static_cast<size_t>( dims )... } );
                }

                NDArray<N, T>& resize( std::initializer_list<size_t> dims ) {
                    return this->resize( std::vector<size_t> { dims } );
                }

                NDArray<N, T>& resize( std::array<size_t, N> dims ) {
                    // size_t total = 1;

                    for (size_t i = 0; i < N; ++i) {
                        _shape[ i ] = dims[ i ];
                        // total *= dims[i];
                    }
                    size_t total = details::total_size( _shape );
                    details::compute_strides( _shape, _strides );
                    // compute_strides();

                    _storage.reset( new details::_ndstorage<T> );
                    _storage->data.resize( total );
                    return *this;
                }

                NDArray<N, T>& resize( std::vector<size_t> dims ) {
                    if (dims.empty()) {
                        dims.resize( N );
                    }
                    if (dims.size() != N) {
                        throw std::out_of_range(
                            "Wrong number of dimensions" );
                    }
                    // size_t total = 1;

                    for (size_t i = 0; i < N; ++i) {
                        _shape[ i ] = dims[ i ];
                        // total *= dims[i];
                    }
                    size_t total = details::total_size( _shape );
                    details::compute_strides( _shape, _strides );
                    // compute_strides();

                    _storage.reset( new details::_ndstorage<T> );
                    _storage->data.resize( total );
                    return *this;
                }

                template<typename ContainerB>
                bool same_shape_as( const ContainerB& b ) const {
                    return details::same_shape( *this, b );
                }

                NDArray<N, T>& set( T val ) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it = val;
                    }
                    return *this;
                }

                NDArray<N, T>& set( const std::vector<T>& vals ) {
                    details::deal_into( vals, *this );
                    return *this;
                }

                const std::array<size_t, N>& shape() const { return _shape; }

                size_t size() const { return _storage->data.size(); }

                NDArrayView<N - 1, T> slice( size_t dim, size_t fixed_index ) {
                    if (dim >= N) throw std::out_of_range( "bad dimension" );

                    if (fixed_index >= _shape[ dim ])
                        throw std::out_of_range( "bad index" );

                    std::array<size_t, N - 1> shape;
                    std::array<size_t, N - 1> strides;

                    details::remove_dimension( _shape, _strides, dim, shape,
                                               strides );


                    size_t offset = fixed_index * _strides[ dim ];

                    return NDArrayView<N - 1, T>( _storage, offset, shape,
                                                  strides );
                }

                // iterators here
                iterator begin() { return _storage->data.begin(); }

                const_iterator begin() const { return _storage->data.begin(); }

                const_iterator cbegin() const {
                    return _storage->data.cbegin();
                }

                iterator end() { return _storage->data.end(); }

                const_iterator end() const { return _storage->data.end(); }

                const_iterator cend() const { return _storage->data.cend(); }

                reverse_iterator rbegin() { return _storage->data.rbegin(); }

                const_reverse_iterator rbegin() const {
                    return _storage->data.rbegin();
                }

                reverse_iterator rend() { return _storage->data.rend(); }

                const_reverse_iterator rend() const {
                    return _storage->data.rend();
                }

            private:
                friend class NDArrayView<N, T>;

                std::shared_ptr<details::_ndstorage<T>> _storage;

                shape_type _shape;
                shape_type _strides;
        };

        template<size_t N, typename T>
        class NDArrayView {
            public:
                using iterator = NDViewIterator<T, std::array<size_t, N>>;
                using const_iterator
                    = NDViewIterator<const T, std::array<size_t, N>>;
                typedef std::reverse_iterator<iterator> reverse_iterator;
                typedef std::reverse_iterator<const_iterator>
                    const_reverse_iterator;
                typedef std::array<size_t, N> shape_type;
                typedef T element_t;

                NDArrayView( std::shared_ptr<details::_ndstorage<T>> storage,
                             size_t offset, const std::array<size_t, N>& shape,
                             const std::array<size_t, N>& strides ) :
                    _storage( storage ),
                    _offset( offset ),
                    _shape( shape ),
                    _strides( strides ) {}

                template<typename... Indices>
                T& at( Indices... indices ) {
                    static_assert( sizeof...( Indices ) == N,
                                   "Wrong number of indices" );

                    std::array<size_t, N> idx
                        = { { static_cast<size_t>( indices )... } };

                    return _storage->data.at( details::linear_index(
                        idx, _shape, _strides, _offset ) );
                }

                template<size_t M = N>
                typename std::enable_if<M == 1, T&>::type at( size_t idx ) {
                    return _storage->data.at( details::linear_index(
                        shape_type { idx }, _shape, _strides, _offset ) );
                }

                template<typename... Indices>
                const T& at( Indices... indices ) const {
                    static_assert( sizeof...( Indices ) == N,
                                   "Wrong number of indices" );

                    std::array<size_t, N> idx
                        = { { static_cast<size_t>( indices )... } };

                    return _storage->data.at( details::linear_index(
                        idx, _shape, _strides, _offset ) );
                }

                template<size_t M = N>
                typename std::enable_if<M == 1, const T&>::type at(
                    size_t idx ) const {
                    return _storage->data.at( details::linear_index(
                        shape_type { idx }, _shape, _strides, _offset ) );
                }

                NDArrayView<N, T>& clear() { return this->set( 0 ); }

                size_t dim( size_t n ) const {
                    assert( n < rank() );
                    return _shape[ n ];
                }

                NDArrayView<N, T>& fill( T val ) { return this->set( val ); }

                template<typename... Indices>
                T get( Indices... indices ) const {
                    return this->at( indices... );
                }

                std::vector<shape_type> indices() const {
                    return details::indices(*this);
                }

                template<size_t M = N>
                typename std::enable_if<( M > 1 ), NDArrayView<N - 1, T>>::type
                    operator[]( size_t ind ) {
                    return this->slice( 0, ind );
                }

                template<size_t M = N>
                typename std::enable_if<M == 1, T&>::type operator[](
                    size_t i ) {
                    return this->at( i );
                }

                template<size_t M = N>
                typename std::enable_if<M == 1, const T&>::type operator[](
                    size_t i ) const {
                    return this->at( i );
                }

                template<typename... Indices>
                T& operator()( Indices... indices ) {
                    return this->at( indices... );
                }

                template<typename... Indices>
                const T& operator()( Indices... indices ) const {
                    return this->at( indices... );
                }

                size_t rank() const { return N; }

                bool rank_mutable() const { return false; }

                template<typename ContainerB>
                bool same_shape_as( const ContainerB& b ) const {
                    return details::same_shape( *this, b );
                }

                NDArrayView<N, T>& set( T val ) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it = val;
                    }
                    return *this;
                }

                NDArrayView<N, T>& set( const std::vector<T>& vals ) {
                    details::deal_into( vals, *this );
                    return *this;
                }

                const std::array<size_t, N>& shape() const { return _shape; }

                size_t size() const { return details::total_size( _shape ); }

                template<size_t M = N>
                typename std::enable_if<( M > 1 ), NDArrayView<N - 1, T>>::type
                    slice( size_t dim, size_t fixed_index ) {
                    if (dim >= N) {
                        throw std::out_of_range( "bad dimension" );
                    }

                    if (fixed_index >= _shape[ dim ]) {
                        throw std::out_of_range( "bad index" );
                    }

                    std::array<size_t, N - 1> shape;
                    std::array<size_t, N - 1> strides;
                    details::remove_dimension( _shape, _strides, dim, shape,
                                               strides );

                    return NDArrayView<N - 1, T>(
                        _storage, _offset + fixed_index * _strides[ dim ],
                        shape, strides );
                }

                iterator begin() {
                    return iterator( _storage, _offset, _shape, _strides, 0 );
                }

                iterator end() {
                    return iterator( _storage, _offset, _shape, _strides,
                                     size() );
                }

                const_iterator begin() const {
                    return const_iterator( _storage, _offset, _shape, _strides,
                                           0 );
                }

                const_iterator end() const {
                    return const_iterator( _storage, _offset, _shape, _strides,
                                           size() );
                }

                const_iterator cbegin() const { return begin(); }

                const_iterator cend() const { return end(); }

                reverse_iterator rbegin() { return reverse_iterator( end() ); }

                reverse_iterator rend() { return reverse_iterator( begin() ); }

                const_reverse_iterator rbegin() const {
                    return const_reverse_iterator( end() );
                }

                const_reverse_iterator rend() const {
                    return const_reverse_iterator( begin() );
                }

                const_reverse_iterator crbegin() const {
                    return const_reverse_iterator( end() );
                }

                const_reverse_iterator crend() const {
                    return const_reverse_iterator( begin() );
                }

            private:
                std::shared_ptr<details::_ndstorage<T>> _storage;

                size_t _offset;

                std::array<size_t, N> _shape;
                std::array<size_t, N> _strides;
        };

        template<typename T>
        class NDVectorView {
            public:
                using iterator = NDViewIterator<T, std::vector<size_t>>;
                using const_iterator
                    = NDViewIterator<const T, std::vector<size_t>>;
                typedef std::reverse_iterator<iterator> reverse_iterator;
                typedef std::reverse_iterator<const_iterator>
                    const_reverse_iterator;
                typedef T element_t;

                NDVectorView( std::shared_ptr<details::_ndstorage<T>> storage,
                              size_t offset, std::vector<size_t> shape,
                              std::vector<size_t> strides ) :
                    _storage( std::move( storage ) ),
                    _offset( offset ),
                    _shape( std::move( shape ) ),
                    _strides( std::move( strides ) ) {}

                size_t rank() const { return _shape.size(); }

                bool rank_mutable() const { return false; }

                T& at( std::initializer_list<size_t> idx ) {
                    return _storage->data.at( details::linear_index(
                        idx, _shape, _strides, _offset ) );
                }

                const T& at( std::initializer_list<size_t> idx ) const {
                    return _storage->data.at( details::linear_index(
                        idx, _shape, _strides, _offset ) );
                }

                T& at( const std::vector<size_t>& idx ) {
                    return _storage->data.at( details::linear_index(
                        idx, _shape, _strides, _offset ) );
                }

                const T& at( const std::vector<size_t>& idx ) const {
                    return _storage->data.at( details::linear_index(
                        idx, _shape, _strides, _offset ) );
                }

                template<typename... Indices>
                T& at( Indices... indices ) {
                    std::vector<size_t> idx { static_cast<size_t>(
                        indices )... };
                    return _storage->data.at( details::linear_index(
                        idx, _shape, _strides, _offset ) );
                }

                template<typename... Indices>
                const T& at( Indices... indices ) const {
                    std::vector<size_t> idx { static_cast<size_t>(
                        indices )... };
                    return _storage->data.at( details::linear_index(
                        idx, _shape, _strides, _offset ) );
                }

                NDVectorView<T>& clear() { return this->set( 0 ); }

                NDVectorView<T>& fill( T val ) { return this->set( val ); }

                template<typename... Indices>
                T get( Indices... indices ) const {
                    return this->at( indices... );
                }

                std::vector<shape_type> indices() const {
                    return details::indices(*this);
                }

                NDVectorView<T>& set( T val ) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it = val;
                    }
                    return *this;
                }

                NDVectorView<T>& set( const std::vector<T>& vals ) {
                    details::deal_into( vals, *this );
                    return *this;
                }

                NDVectorView slice( size_t dim, size_t fixed_index ) const {
                    if (dim >= rank()) {
                        throw std::out_of_range( "bad dimension" );
                    }

                    if (fixed_index >= _shape[ dim ]) {
                        throw std::out_of_range( "bad index" );
                    }

                    std::vector<size_t> new_shape;
                    std::vector<size_t> new_strides;

                    new_shape.reserve( rank() - 1 );
                    new_strides.reserve( rank() - 1 );

                    for (size_t i = 0; i < rank(); ++i) {
                        if (i != dim) {
                            new_shape.push_back( _shape[ i ] );
                            new_strides.push_back( _strides[ i ] );
                        }
                    }

                    return NDVectorView(
                        _storage, _offset + fixed_index * _strides[ dim ],
                        std::move( new_shape ), std::move( new_strides ) );
                }

                const std::vector<size_t>& shape() const { return _shape; }

                NDVectorView operator[]( size_t ind ) const {
                    return this->slice( 0, ind );
                }

                template<typename... Indices>
                T& operator()( Indices... indices ) {
                    return this->at( indices... );
                }

                template<typename... Indices>
                const T& operator()( Indices... indices ) const {
                    return this->at( indices... );
                }

                template<typename ContainerB>
                bool same_shape_as( const ContainerB& b ) const {
                    return details::same_shape( *this, b );
                }

                size_t size() const { return details::total_size( _shape ); }

                size_t dim( size_t n ) const {
                    assert( n < rank() );
                    return _shape[ n ];
                }

                iterator begin() {
                    return iterator( _storage, _offset, _shape, _strides, 0 );
                }

                iterator end() {
                    return iterator( _storage, _offset, _shape, _strides,
                                     size() );
                }

                const_iterator begin() const {
                    return const_iterator( _storage, _offset, _shape, _strides,
                                           0 );
                }

                const_iterator end() const {
                    return const_iterator( _storage, _offset, _shape, _strides,
                                           size() );
                }

                const_iterator cbegin() const { return begin(); }

                const_iterator cend() const { return end(); }

                reverse_iterator rbegin() { return reverse_iterator( end() ); }

                reverse_iterator rend() { return reverse_iterator( begin() ); }

                const_reverse_iterator rbegin() const {
                    return const_reverse_iterator( end() );
                }

                const_reverse_iterator rend() const {
                    return const_reverse_iterator( begin() );
                }

                const_reverse_iterator crbegin() const {
                    return const_reverse_iterator( end() );
                }

                const_reverse_iterator crend() const {
                    return const_reverse_iterator( begin() );
                }

            private:
                std::shared_ptr<details::_ndstorage<T>> _storage;

                size_t _offset;

                std::vector<size_t> _shape;
                std::vector<size_t> _strides;
        };

        template<typename T>
        class NDVector {
            public:
                typedef typename std::vector<T>::iterator iterator;
                typedef typename std::vector<T>::const_iterator const_iterator;
                typedef
                    typename std::vector<T>::reverse_iterator reverse_iterator;
                typedef typename std::vector<T>::const_reverse_iterator
                    const_reverse_iterator;
                typedef T element_t;

                template<typename... Dims>
                explicit NDVector( Dims... dims ) {
                    this->resize( dims... );
                }

                NDVector( std::vector<size_t> shape ) {
                    this->resize( shape );
                }

                NDVector( std::initializer_list<size_t> shape ) :
                    NDVector<T>( std::vector<size_t>( shape ) ) {}

                NDVector( const NDVector<T>& other ) = default;

                template<size_t N>
                NDVector( const NDArray<N, T>& other ) {
                    this->copy( other );
                }

                template<size_t N>
                NDVector( const NDArrayView<N, T>& other ) {
                    this->copy( other );
                }

                NDVector( const NDVectorView<T>& other ) {
                    this->copy( other );
                }

                T& at( const std::vector<size_t>& idx ) {
                    return _storage->data.at(
                        details::linear_index( idx, _shape, _strides ) );
                }

                const T& at( const std::vector<size_t>& idx ) const {
                    return _storage->data.at(
                        details::linear_index( idx, _shape, _strides ) );
                }

                T& at( std::initializer_list<size_t> idx ) {
                    return this->at( std::vector<size_t> { idx } );
                    // return _storage->data.at(
                    //     details::linear_index( idx, _shape, _strides ) );
                }

                const T& at( std::initializer_list<size_t> idx ) const {
                    return this->at( std::vector<size_t> { idx } );
                    // return _storage->data.at(
                    //     details::linear_index( idx, _shape, _strides ) );
                }

                // T& at( size_t idx ) { return this->at( vectorize( idx ) ); }

                // const T& at( size_t idx ) const {
                //     return this->at( vectorize( idx ) );
                // }

                template<typename... Indices>
                T& at( Indices... indices ) {
                    std::vector<size_t> idx { static_cast<size_t>(
                        indices )... };
                    return this->at( idx );
                    // return _storage->data.at(
                    //     details::linear_index( idx, _shape, _strides ) );
                }

                template<typename... Indices>
                const T& at( Indices... indices ) const {
                    std::vector<size_t> idx { static_cast<size_t>(
                        indices )... };
                    return this->at( idx );
                    // return _storage->data.at(
                    //     details::linear_index( idx, _shape, _strides ) );
                }

                NDVector<T>& clear() { return this->set( 0 ); }

                template<typename FromType>
                void copy( const FromType& other ) {
                    details::copy_to( other, *this );
                }

                T *data() { return _storage->data.data(); }

                const T *data() const { return _storage->data.data(); }

                size_t dim( size_t n ) const {
                    assert( n < rank() );
                    return _shape[ n ];
                }

                NDVector<T>& fill( T val ) { return this->set( val ); }

                template<typename... Indices>
                T get( Indices... indices ) const {
                    return this->at( indices... );
                }

                std::vector<shape_type> indices() const {
                    return details::indices(*this);
                }

                NDVectorView<T> operator[]( size_t ind ) const {
                    return this->slice( 0, ind );
                }

                template<typename... Indices>
                T& operator()( Indices... indices ) {
                    return this->at( indices... );
                }

                template<typename... Indices>
                const T& operator()( Indices... indices ) const {
                    return this->at( indices... );
                }

                size_t rank() const { return _shape.size(); }

                bool rank_mutable() const { return true; }

                template<typename... Dims>
                NDVector<T>& resize( Dims... dims ) {
                    return this->resize( std::vector<size_t> {
                        static_cast<size_t>( dims )... } );
                }

                NDVector<T>& resize( std::vector<size_t> shape ) {
                    if (shape.empty()) {
                        _shape.clear();
                        _strides.clear();
                        _storage = std::make_shared<details::_ndstorage<T>>();
                    } else if (_shape != shape) {
                        _shape = std::move( shape );
                        details::compute_strides( _shape, _strides );
                        size_t total = details::total_size( _shape );
                        _storage = std::make_shared<details::_ndstorage<T>>();
                        _storage->data.resize( total );
                    }

                    return *this;
                }

                NDVector<T>& resize( std::initializer_list<size_t> shape ) {
                    return this->resize( std::vector<size_t>( shape ) );
                }

                template<size_t N>
                NDVector<T>& resize( const std::array<size_t, N>& shape ) {
                    return this->resize(
                        std::vector<size_t>( shape.begin(), shape.end() ) );
                }

                template<typename ContainerB>
                bool same_shape_as( const ContainerB& b ) const {
                    return details::same_shape( *this, b );
                }

                NDVector<T>& set( T val ) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it = val;
                    }
                    return *this;
                }

                NDVectorView<T>& set( const std::vector<T>& vals ) {
                    details::deal_into( vals, *this );
                    return *this;
                }

                const std::vector<size_t>& shape() const { return _shape; }

                size_t size() const { return _storage->data.size(); }

                NDVectorView<T> slice( size_t dim, size_t fixed_index ) const {
                    if (dim >= rank()) {
                        throw std::out_of_range( "bad dimension" );
                    }

                    if (fixed_index >= _shape[ dim ]) {
                        throw std::out_of_range( "bad index" );
                    }

                    std::vector<size_t> shape;
                    std::vector<size_t> strides;

                    shape.reserve( rank() - 1 );
                    strides.reserve( rank() - 1 );

                    for (size_t i = 0; i < rank(); ++i) {
                        if (i != dim) {
                            shape.push_back( _shape[ i ] );
                            strides.push_back( _strides[ i ] );
                        }
                    }

                    return NDVectorView<T>(
                        _storage, fixed_index * _strides[ dim ],
                        std::move( shape ), std::move( strides ) );
                }

                // iterators here
                iterator begin() { return _storage->data.begin(); }

                const_iterator begin() const { return _storage->data.begin(); }

                const_iterator cbegin() const {
                    return _storage->data.cbegin();
                }

                iterator end() { return _storage->data.end(); }

                const_iterator end() const { return _storage->data.end(); }

                const_iterator cend() const { return _storage->data.cend(); }

                reverse_iterator rbegin() { return _storage->data.rbegin(); }

                const_reverse_iterator rbegin() const {
                    return _storage->data.rbegin();
                }

                reverse_iterator rend() { return _storage->data.rend(); }

                const_reverse_iterator rend() const {
                    return _storage->data.rend();
                }

            private:
                std::shared_ptr<details::_ndstorage<T>> _storage;

                std::vector<size_t> _shape;
                std::vector<size_t> _strides;
        };

        template<typename T>
        using vector2d_t = NDArray<2, T>;

        template<typename T>
        using vector3d_t = NDArray<3, T>;

        template<typename T, typename std::enable_if<
                                 std::is_arithmetic<T>::value, int>::type = 0>
        size_t vector_rank( const std::vector<T>& ) {
            return 1;
        }

        template<typename T, typename std::enable_if<
                                 !std::is_arithmetic<T>::value, int>::type = 0>
        size_t vector_rank( const std::vector<T>& vec ) {
            if (vec.empty())
                throw std::logic_error(
                    "Can't determine rank of empty nested vector!" );

            return vector_rank( vec.front() ) + 1;
        }

        template<typename T>
        std::vector<T> vectorify( const T& val ) {
            return std::vector<T>( 1, val );
        }

        /**
         * Dynamically allocates a new array and sets all elements to zero
         * before returning it.
         * @brief Returns a new array of all zeros.
         * @param n The size of the array.
         * @return A pointer to the newly-allocated array.
         */
        template<typename T>
        T *zeros( size_t n ) {
            T *out = new T[ n ]();
            return out;
        }

        /**
         * @brief Dynamically allocates and returns a zero-initialized
         * two-dimensional array.
         * @param nr The first dimension of the array.
         * @param nc The second dimension of the array.
         * @return A pointer to the newly-allocated array.
         */
        template<typename T>
        T **zeros( size_t nr, size_t nc ) {
            T **v;
            v = new T *[ nr ];
            for (size_t i = 0; i < nr; i++) {
                v[ i ] = new T[ nc ]();
            }
            return v;
        }

        /**
         * @brief Dynamically allocates and returns a zero-initialized
         * three-dimensional array.
         * @param nr The first dimension of the array.
         * @param nc The second dimension of the array.
         * @param nd The third dimension of the array.
         * @return A pointer to the newly-allocated array.
         */
        template<typename T>
        T ***zeros( size_t nr, size_t nc, size_t nd ) {
            T ***v;
            v = new T **[ nr ];
            for (size_t i = 0; i < nr; i++) {
                v[ i ] = new T *[ nc ];
                for (size_t j = 0; j < nc; j++) {
                    v[ i ][ j ] = new T[ nd ]();
                }
            }
            return v;
        }

        /**
         * @brief Frees a dynamically allocated 1-D array.
         * @param v The array to free.
         * @param nr The dimension of the array.
         * @param ENABLER SFINAE enabler for non-deletable types.
         */
        template<typename T>
        void free_array( T *& v, size_t nr,
                         ENABLE_FUNCTION_IF_NOT_DELETEABLE( T ) ) {
            delete[] v;
            v = nullptr;
        }

        /**
         * @brief Frees a dynamically allocated 1-D array of pointers.
         * @param v The array to free.
         * @param nr The dimension of the array.
         * @param ENABLER SFINAE enabler for deletable types.
         */
        template<typename T>
        void free_array( T *& v, size_t nr,
                         ENABLE_FUNCTION_IF_DELETEABLE( T ) ) {
            for (size_t i = 0; i < nr; i++) {
                delete v[ i ];
            }
            delete[] v;
            v = nullptr;
        }

        /**
         * @brief Frees a dynamically allocated 2-D array.
         * @param v The array to free.
         * @param nr The first dimension of the array.
         * @param nc The second dimension of the array.
         * @param ENABLER SFINAE enabler for non-deletable types.
         */
        template<typename T>
        void free_array( T **& v, size_t nr, size_t nc,
                         ENABLE_FUNCTION_IF_NOT_DELETEABLE( T ) ) {
            for (size_t i = 0; i < nr; i++) {
                delete[] v[ i ];
            }
            delete[] v;
            v = nullptr;
        }

        /**
         * @brief Frees a dynamically allocated 2-D array of pointers.
         * @param v The array to free.
         * @param nr The first dimension of the array.
         * @param nc The second dimension of the array.
         * @param ENABLER SFINAE enabler for deletable types.
         */
        template<typename T>
        void free_array( T **& v, size_t nr, size_t nc,
                         ENABLE_FUNCTION_IF_DELETEABLE( T ) ) {
            for (size_t i = 0; i < nr; i++) {
                for (size_t j = 0; j < nc; j++) {
                    delete v[ i ][ j ];
                }
                delete[] v[ i ];
            }
            delete[] v;
            v = nullptr;
        }

        /**
         * Frees a dynamically-allocated three-dimensional array.
         * @brief Frees a dynamically-allocated three-dimensional array and its
         * contents.
         * @param data The array to free.
         * @param nd1 The first dimension of the array.
         * @param nd2 The second dimension of the array.
         * @param nd3 The third dimension of the array.
         * @param ENABLER SFINAE enabler for non-deletable types.
         */
        template<typename T>
        void free_array( T ***& data, size_t nd1, size_t nd2, size_t nd3,
                         ENABLE_FUNCTION_IF_NOT_DELETEABLE( T ) ) {
            size_t i, j;
            for (i = 0; i < nd1; ++i) {
                if (data[ i ] != NULL) {
                    for (j = 0; j < nd2; ++j) {
                        delete[] data[ i ][ j ];
                    }
                    delete[] data[ i ];
                }
            }
            delete[] data;
            data = nullptr;
        }

        /**
         * @brief Frees a dynamically-allocated three-dimensional array of
         * pointers.
         * @param data The array to free.
         * @param nd1 The first dimension of the array.
         * @param nd2 The second dimension of the array.
         * @param nd3 The third dimension of the array.
         * @param ENABLER SFINAE enabler for deletable types.
         */
        template<typename T>
        void free_array( T ***& data, size_t nd1, size_t nd2, size_t nd3,
                         ENABLE_FUNCTION_IF_DELETEABLE( T ) ) {
            size_t i, j, k;
            for (i = 0; i < nd1; ++i) {
                if (data[ i ] != NULL) {
                    for (j = 0; j < nd2; ++j) {
                        for (k = 0; k < nd3; k++) {
                            delete data[ i ][ j ][ k ];
                        }
                        delete[] data[ i ][ j ];
                    }
                    delete[] data[ i ];
                }
            }
            delete[] data;
            data = nullptr;
        }

        /**
         * Circularly shifts the elements in an array X by K positions. If K is
         * positive, then the values of X are circularly shifted from the
         * beginning to the end. If K is negative, they are shifted from the
         * end to the beginning.
         * @brief Circularly shifts array elements.
         * @param X The array whose elements are to be shifted.
         * @param N The size of the array.
         * @param K The number of positions to shift.
         * @param out The shifted array output. Can be the same as input to
         * perform in-place.
         */
        template<typename T>
        void circshift( T *X, size_t N, int K, T *& out ) {
            while (std::abs( K ) > (int)N) {
                K -= ( (int)N ) * ( K < 0 ? -1 : 1 );
            }
            size_t i;
            T *tempvec = zeros<T>( N );
            if (K < 0) {
                size_t negshift = (size_t)( -K );
                for (i = 0; i < negshift; i++) {
                    tempvec[ i ] = X[ N - negshift + i ];
                }
                for (i = negshift; i < N; i++) {
                    tempvec[ i ] = X[ i - negshift ];
                }
            } else if (K > 0) {
                for (i = 0; i < (size_t)K; i++) {
                    tempvec[ N - K + i ] = X[ i ];
                }
                for (i = (size_t)K; i < N; i++) {
                    tempvec[ i - K ] = X[ i ];
                }
            } else {
                std::memcpy( tempvec, X, N * sizeof( T ) );
            }

            if (out == nullptr) {
                out = zeros<T>( N );
            }
            std::memcpy( out, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        /**
         * @brief Sets each element of an array to a constant value.
         * @param v The array to fill.
         * @param nr The dimension of the array.
         * @param val The value to set each element to.
         */
        template<typename T>
        void fill( T *v, size_t nr, const T& val ) {
            std::fill( v, v + nr, val );
        }

        /**
         * @brief Sets each element of a 2-D array to a constant value.
         * @param v The array to fill.
         * @param nr The first dimension of the array.
         * @param nc The second dimension of the array.
         * @param val The value to set each element to.
         */
        template<typename T>
        void fill( T **v, size_t nr, size_t nc, const T& val ) {
            for (size_t i = 0; i < nr; i++) {
                fill( v[ i ], nc, val );
            }
        }

        /**
         * @brief Sets each element of a three-dimensional array to a constant
         * value.
         * @param v The array to fill.
         * @param nr The first dimension of the array.
         * @param nc The second dimension of the array.
         * @param nz The third dimension of the array.
         * @param val The value to set each element to.
         */
        template<typename T>
        void fill( T ***v, size_t nr, size_t nc, size_t nz, const T& val ) {
            for (size_t i = 0; i < nr; i++) {
                fill( v[ i ], nc, nz, val );
            }
        }

        /**
         * @brief Copies one array to another.
         * @param from The array to copy.
         * @param n The number of values to copy.
         * @param to The array to copy into.
         */
        template<typename T>
        void copy( const T *from, size_t n, T *to ) {
            if (from != nullptr) {
                std::copy( from, from + n, to );
            } else if (to != nullptr) {
                free_array( to, n );
            }
        }

        /**
         * @brief Copies one 2-D array to another.
         * @param from The source array.
         * @param n1 The first dimension.
         * @param n2 The second dimension.
         * @param to The destination array.
         */
        template<typename T>
        void copy( const T **from, size_t n1, size_t n2, T **to ) {
            if (from != nullptr) {
                for (size_t i = 0; i < n1; i++) {
                    copy( from[ i ], n2, to[ i ] );
                }
            } else if (to != nullptr) {
                free_array( to, n1, n2 );
            }
        }

        /**
         * @brief Copies one 3-D array to another.
         * @param from The source array.
         * @param n1 The first dimension.
         * @param n2 The second dimension.
         * @param n3 The third dimension.
         * @param to The destination array.
         */
        template<typename T>
        void copy( const T ***from, size_t n1, size_t n2, size_t n3,
                   T ***to ) {
            if (from != nullptr) {
                for (size_t i = 0; i < n1; i++) {
                    copy( from[ i ], n2, n3, to[ i ] );
                }
            } else if (to != nullptr) {
                free_array( to, n1, n2, n3 );
            }
        }

        /**
         * Dynamically allocates a new array and sets the elements to
         * (a, a+1, a+2, ..., a+n-1).
         * @brief Returns a new array of index values.
         * @param n The size of the array.
         * @param a The offset from zero to use for each value. Defaults to 0.
         * @return A vector containing the index values.
         */
        template<typename T>
        std::vector<T> index_vector( size_t n, T a = 0 ) {
            std::vector<T> inds( n );
            for (size_t i = 0; i < n; i++) {
                inds[ i ] = (T)i + a;
            }
            return inds;
        }

        /**
         * @brief Reverses the order of the first N elements of an array.
         * @param in The input array.
         * @param N The number of elements of the array to reverse.
         * @param out Holds the output array. Can be the same as the input.
         */
        template<typename T>
        void reverse( T *in, size_t N, T *& out ) {
            T *tempvec = zeros<T>( N );
            size_t j;
            for (j = 0; j < N; j++) {
                tempvec[ N - 1 - j ] = in[ j ];
            }
            std::memcpy( out, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        template<typename ContainerType, typename ContainerType2>
        bool increment_digits( ContainerType& values,
                               const ContainerType2& sizes ) {
            if (values.empty()) return false;
            int ind = values.size() - 1;
            while (ind >= 0) {
                if (sizes.at( ind ) == 0) return false;
                ++values.at( ind );
                if (values.at( ind ) == sizes.at( ind )) {
                    values.at( ind ) = 0;
                    --ind;
                } else {
                    return true;
                }
            }
            return false;
        }

        /**
         * @brief Generates a sort permutation for a vector.
         * @tparam T Element type.
         * @tparam Compare Comparison functor type.
         * @param vec The vector to analyze.
         * @param compare Comparison functor.
         * @return Vector of indices representing the sorted order.
         */
        template<typename T, typename Compare>
        std::vector<std::size_t> sort_permutation( const std::vector<T>& vec,
                                                   Compare compare ) {
            std::vector<std::size_t> p( vec.size() );
            std::iota( p.begin(), p.end(), 0 );
            std::sort( p.begin(), p.end(),
                       [ & ]( std::size_t i, std::size_t j ) {
                           return compare( vec[ i ], vec[ j ] );
                       } );
            return p;
        }

        /**
         * @brief Generates a sort permutation based on two vectors (primary
         * and secondary keys).
         * @tparam T Element type.
         * @tparam Compare Comparison functor type.
         * @param vec1 Primary key vector.
         * @param vec2 Secondary key vector.
         * @param compare Comparison functor.
         * @return Vector of indices representing the sorted order.
         */
        template<typename T, typename Compare>
        std::vector<std::size_t> sort_permutation( const std::vector<T>& vec1,
                                                   const std::vector<T>& vec2,
                                                   Compare compare ) {
            if (vec1.size() != vec2.size()) {
                throw std::invalid_argument(
                    "Vectors must be the same size!" );
            }
            std::vector<std::size_t> p( vec1.size() );
            std::iota( p.begin(), p.end(), 0 );
            std::sort( p.begin(), p.end(),
                       [ & ]( std::size_t i, std::size_t j ) {
                           return compare( vec1[ i ], vec1[ j ] )
                               || ( vec1[ i ] == vec1[ j ]
                                    && compare( vec2[ i ], vec2[ j ] ) );
                       } );
            return p;
        }

        /**
         * @brief Returns permutation indices for increasing sort.
         * @param vec Input vector.
         * @return Permutation vector.
         */
        template<typename T>
        std::vector<std::size_t> sort_permutation_increasing(
            const std::vector<T>& vec ) {
            return sort_permutation<T>(
                vec, []( T const& a, T const& b ) { return a < b; } );
        }

        /**
         * @brief Returns permutation indices for increasing sort on two keys.
         * @param vec1 Primary vector.
         * @param vec2 Secondary vector.
         * @return Permutation vector.
         */
        template<typename T>
        std::vector<std::size_t> sort_permutation_increasing(
            const std::vector<T>& vec1, const std::vector<T>& vec2 ) {
            return sort_permutation<T>(
                vec1, vec2, []( T const& a, T const& b ) { return a < b; } );
        }

        /**
         * @brief Returns permutation indices for decreasing sort.
         * @param vec Input vector.
         * @return Permutation vector.
         */
        template<typename T>
        std::vector<std::size_t> sort_permutation_decreasing(
            const std::vector<T>& vec ) {
            return sort_permutation<T>(
                vec, []( T const& a, T const& b ) { return a > b; } );
        }

        /**
         * @brief Returns permutation indices for decreasing sort on two keys.
         * @param vec1 Primary vector.
         * @param vec2 Secondary vector.
         * @return Permutation vector.
         */
        template<typename T>
        std::vector<std::size_t> sort_permutation_decreasing(
            const std::vector<T>& vec1, const std::vector<T>& vec2 ) {
            return sort_permutation<T>(
                vec1, vec2, []( T const& a, T const& b ) { return a > b; } );
        }

        /**
         * @brief Applies a permutation to a vector, returning a new vector.
         * @param vec Source vector.
         * @param p Permutation vector.
         * @return The permuted vector.
         */
        template<typename T>
        std::vector<T> apply_permutation( const std::vector<T>& vec,
                                          const std::vector<std::size_t>& p ) {
            std::vector<T> sorted_vec( vec.size() );
            std::transform( p.begin(), p.end(), sorted_vec.begin(),
                            [ & ]( std::size_t i ) { return vec[ i ]; } );
            return sorted_vec;
        }

        /**
         * @brief Applies a permutation to a vector in-place.
         * @param vec Vector to permute.
         * @param p Permutation vector.
         */
        template<typename T>
        void apply_permutation_in_place( std::vector<T>& vec,
                                         const std::vector<std::size_t>& p ) {
            std::vector<bool> done( vec.size() );
            for (std::size_t i = 0; i < vec.size(); ++i) {
                if (done[ i ]) {
                    continue;
                }
                done[ i ]          = true;
                std::size_t prev_j = i;
                std::size_t j      = p[ i ];
                while (i != j) {
                    std::swap( vec[ prev_j ], vec[ j ] );
                    done[ j ] = true;
                    prev_j    = j;
                    j         = p[ j ];
                }
            }
        }

        /**
         * @brief Writes vector elements to an output stream.
         * @param os The output stream.
         * @param basevec The vector to write.
         * @param separator String to place between elements.
         */
        template<typename T>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::string& separator = ", " ) {
            for (size_t i = 0; i < basevec.size(); ++i) {
                if (i > 0) {
                    os << separator;
                }
                os << basevec[ i ];
            }
        }

        /**
         * @brief Writes a base vector and associated item vectors in columns.
         * @param os Output stream.
         * @param basevec Primary vector.
         * @param items Collection of dependent vectors.
         * @param separator Delimiter.
         */
        template<typename T, typename U>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::vector<std::vector<U>>& items,
                    const std::string& separator = " " ) {
            for (size_t i = 0; i < basevec.size(); i++) {
                os << basevec[ i ];
                for (size_t j = 0; j < items.size(); j++) {
                    os << separator << items[ j ][ i ];
                }
                os << std::endl;
            }
        }

        /**
         * @brief Writes two vectors as coordinate pairs for arithmetic types.
         * @param os Output stream.
         * @param basevec X-axis vector.
         * @param depvec Y-axis vector.
         * @param separator Delimiter.
         * @param ENABLER SFINAE for arithmetic types.
         */
        template<typename T, typename U>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::vector<U> depvec,
                    const std::string& separator = " ",
                    ENABLE_FUNCTION_IF_ARITHMETIC( U ) ) {
            for (size_t i = 0; i < basevec.size(); i++) {
                os << basevec[ i ] << separator << depvec[ i ] << std::endl;
            }
        }

        /**
         * @brief Writes two vectors as coordinate pairs for complex types
         * (splits real/imag).
         * @param os Output stream.
         * @param basevec X-axis vector.
         * @param depvec Complex Y-axis vector.
         * @param separator Delimiter.
         * @param ENABLER SFINAE for complex types.
         */
        template<typename T, typename U>
        void write( std::ostream& os, const std::vector<T> basevec,
                    const std::vector<U> depvec,
                    const std::string& separator = " ",
                    ENABLE_FUNCTION_IF_COMPLEX( U ) ) {
            for (size_t i = 0; i < basevec.size(); i++) {
                os << basevec[ i ] << separator << depvec[ i ].real()
                   << separator << depvec[ i ].imag() << std::endl;
            }
        }

        /**
         * @brief Performs element-wise vector addition.
         * @param v1 The first vector to add.
         * @param v2 The second vector to add.
         * @return The vector holding the sum.
         */
        template<typename T>
        T add_vectors( const T& v1, const T& v2,
                       ENABLE_FUNCTION_IF_ITERABLE( T ) ) {
            size_t N = std::min<size_t>( v1.size(), v2.size() );
            T v3     = v1.size() >= v2.size() ? v1 : v2;
            std::transform( v1.cbegin(), v1.cbegin() + N, v2.cbegin(),
                            v3.begin(), std::plus<typename T::value_type> {} );
            return v3;
        }

        /**
         * @brief Performs element-wise array addition.
         * @param N The number of points in the array.
         * @param v1 The first array to add.
         * @param v2 The second array to add.
         * @param v12 The array to hold the sum. Can be the same as either
         * input.
         */
        template<typename T>
        void add_arrays( size_t N, const T *v1, const T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for (size_t i = 0; i < N; i++) {
                tempvec[ i ] = v1[ i ] + v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        /**
         * Performs element-wise array division. If arrays are different
         * lengths, all values beyond the length of the shorter vector will be
         * zero.
         * @brief Performs element-wise array division.
         * @param v1 The array to divide.
         * @param v2 The array to divide by.
         * @return The array holding the quotient.
         */
        template<typename T>
        T divide_vectors( const T& v1, const T& v2,
                          ENABLE_FUNCTION_IF_ITERABLE( T ) ) {
            T v3 = T( std::max<size_t>( v1.size(), v2.size() ), 0.0 );
            std::transform( v1.cbegin(),
                            v1.cbegin()
                                + std::min<size_t>( v1.size(), v2.size() ),
                            v2.cbegin(), v3.begin(),
                            std::divides<typename T::value_type> {} );
            return v3;
        }

        /**
         * Divides one array by another element-wise, returning the
         * quotient in a supplied array.
         * @brief Performs element-wise array division.
         * @param N The number of points in the array.
         * @param v1 The array to divide.
         * @param v2 The array to divide by.
         * @param v12 The array to hold the quotient. Can be the same as input.
         */
        template<typename T>
        void divide_arrays( size_t N, const T *v1, const T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for (size_t i = 0; i < N; i++) {
                tempvec[ i ] = v1[ i ] / v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        /**
         * Performs element-wise array multiplication.
         * @brief Performs element-wise array multiplication.
         * @param v1 The first array to multiply.
         * @param v2 The second array to multiply.
         * @return The array holding the product.
         */
        template<typename T>
        T multiply_vectors( const T& v1, const T& v2,
                            ENABLE_FUNCTION_IF_ITERABLE( T ) ) {
            T v3 = T( std::max<size_t>( v1.size(), v2.size() ), 0.0 );
            std::transform( v1.cbegin(),
                            v1.cbegin()
                                + std::min<size_t>( v1.size(), v2.size() ),
                            v2.cbegin(), v3.begin(),
                            std::multiplies<typename T::value_type> {} );
            return v3;
        }

        /**
         * Multiplies two arrays together element-wise.
         * @brief Performs element-wise array multiplication.
         * @param N The number of points in the array.
         * @param v1 The first array to multiply.
         * @param v2 The second array to multiply.
         * @param v12 The array to hold the product.
         */
        template<typename T>
        void multiply_arrays( size_t N, const T *v1, const T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for (size_t i = 0; i < N; i++) {
                tempvec[ i ] = v1[ i ] * v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        /**
         * Scales an array of values by a constant value.
         * @brief Scales an array by a constant value.
         * @param N The number of points in the array.
         * @param in The array to scale.
         * @param factor The factor to scale by.
         * @param out The output scaled array.
         */
        template<typename T, typename U>
        void scale_array( size_t N, const U *in, T factor, U *& out ) {
            U *tempvec = NCPA::arrays::zeros<U>( N );
            for (size_t i = 0; i < N; i++) {
                tempvec[ i ] = in[ i ] * (U)factor;
            }
            std::memcpy( out, tempvec, N * sizeof( U ) );
            delete[] tempvec;
        }

        /**
         * Scales an iterable container.
         * @brief Performs container scaling.
         * @param v1 The array to multiply.
         * @param scalar The scalar to multiply by.
         * @param ENABLER SFINAE for iterable types.
         * @return The scaled container.
         */
        template<typename T, typename U>
        T scale_vector(
            const T& v1, U scalar,
            typename std::enable_if<NCPA::types::is_iterable_of<T, U>::value,
                                    int>::type ENABLER = 0 ) {
            if (scalar == NCPA::constants::one<U>()) {
                return v1;
            }
            T v3 = v1;
            std::transform( v3.begin(), v3.end(), v3.begin(),
                            [ scalar ]( U num ) { return num * scalar; } );
            return v3;
        }

        /**
         * @brief Scales an array of values in-place by a constant value.
         * @param N The number of points in the array.
         * @param in The array to scale.
         * @param factor The factor to scale by.
         */
        template<typename T, typename U>
        void scale_array( size_t N, U *in, T factor ) {
            if (factor != NCPA::constants::one<T>()) {
                for (size_t i = 0; i < N; i++) {
                    in[ i ] *= factor;
                }
            }
        }

        /**
         * Offsets an iterable container by a constant.
         * @brief Performs array offset.
         * @param v1 The container.
         * @param scalar The scalar to add.
         * @param ENABLER SFINAE for iterable types.
         * @return The offset container.
         */
        template<typename T, typename U>
        T offset_vector(
            const T& v1, U scalar,
            typename std::enable_if<NCPA::types::is_iterable_of<T, U>::value,
                                    int>::type ENABLER = 0 ) {
            if (NCPA::constants::is_zero<U>( scalar )) {
                return v1;
            }
            T v3 = v1;
            for (size_t i = 0; i < v1.size(); i++) {
                v3[ i ] += scalar;
            }
            return v3;
        }

        /**
         * @brief Offsets an array of values in-place by a constant value.
         * @param N The number of points in the array.
         * @param in The array to offset.
         * @param factor The amount to offset by.
         */
        template<typename T, typename U>
        void offset_array( size_t N, U *in, T factor ) {
            if (!NCPA::constants::is_zero<T>( factor )) {
                for (size_t i = 0; i < N; i++) {
                    in[ i ] += factor;
                }
            }
        }

        /**
         * @brief Performs element-wise container subtraction.
         * @param v1 The first container.
         * @param v2 The second container.
         * @param ENABLER SFINAE for iterable types.
         * @return Resulting container.
         */
        template<typename T>
        T subtract_vectors(
            const T& v1, const T& v2,
            typename std::enable_if<NCPA::types::is_iterable<T>::value,
                                    int>::type ENABLER = 0 ) {
            T v3 = scale_vector( v2, -1.0 );
            return add_vectors( v1, v3 );
        }

        /**
         * Subtracts one array from another element-wise.
         * @brief Performs element-wise array subtraction.
         * @param N The number of points in the array.
         * @param v1 The array to subtract from.
         * @param v2 The array to subtract.
         * @param v12 The array to hold the difference.
         */
        template<typename T>
        void subtract_arrays( size_t N, const T *v1, const T *v2, T *& v12 ) {
            T *tempvec = NCPA::arrays::zeros<T>( N );
            for (size_t i = 0; i < N; i++) {
                tempvec[ i ] = v1[ i ] - v2[ i ];
            }
            std::memcpy( v12, tempvec, N * sizeof( T ) );
            delete[] tempvec;
        }

        /**
         * @brief Casts elements of a vector from one type to another.
         * @tparam TOTYPE Target type.
         * @tparam FROMTYPE Source type.
         * @param vec Input vector.
         * @return New vector with casted elements.
         */
        template<typename TOTYPE, typename FROMTYPE>
        std::vector<TOTYPE> cast_vector( const std::vector<FROMTYPE>& vec ) {
            std::vector<TOTYPE> newvec( vec.size() );
            for (size_t i = 0; i < vec.size(); ++i) {
                newvec[ i ] = static_cast<TOTYPE>( vec.at( i ) );
            }
            return newvec;
        }

        template<typename ELEMENTTYPE, typename INDEXTYPE = size_t>
        std::vector<ELEMENTTYPE> extract_indices(
            const std::vector<ELEMENTTYPE>& source,
            const std::vector<INDEXTYPE>& indices ) {
            std::vector<ELEMENTTYPE> result( indices.size() );
            for (size_t i = 0; i < indices.size(); ++i) {
                result[ i ] = source.at( indices.at( i ) );
            }
            return result;
        }

        template<typename ELEMENTTYPE>
        std::vector<ELEMENTTYPE> vectorize( const ELEMENTTYPE& val,
                                            size_t n = 1 ) {
            return std::vector<ELEMENTTYPE>( n, val );
        }

        /**
         * @brief Operator overload for vector-vector addition.
         */
        template<typename T>
        std::vector<T> operator+( const std::vector<T>& a,
                                  const std::vector<T>& b ) {
            return NCPA::arrays::add_vectors<std::vector<T>>( a, b );
        }

        /**
         * @brief Operator overload for vector-scalar addition.
         */
        template<typename T>
        std::vector<T> operator+( const std::vector<T>& a, const T& b ) {
            return NCPA::arrays::offset_vector<std::vector<T>>( a, b );
        }

        /**
         * @brief Operator overload for vector-vector subtraction.
         */
        template<typename T>
        std::vector<T> operator-( const std::vector<T>& a,
                                  const std::vector<T>& b ) {
            return NCPA::arrays::subtract_vectors<std::vector<T>>( a, b );
        }

        /**
         * @brief Operator overload for vector-scalar subtraction.
         */
        template<typename T>
        std::vector<T> operator-( const std::vector<T>& a, const T& b ) {
            return NCPA::arrays::offset_vector<std::vector<T>>( a, -b );
        }

        /**
         * @brief Operator overload for element-wise vector-vector
         * multiplication.
         */
        template<typename T>
        std::vector<T> operator*( const std::vector<T>& a,
                                  const std::vector<T>& b ) {
            return NCPA::arrays::multiply_vectors<std::vector<T>>( a, b );
        }

        /**
         * @brief Operator overload for vector scaling.
         */
        template<typename T>
        std::vector<T> operator*( const std::vector<T>& a, const T& b ) {
            return NCPA::arrays::scale_vector<std::vector<T>>( a, b );
        }

        /**
         * @brief Operator overload for vector-scalar division.
         */
        template<typename T>
        std::vector<T> operator/( const std::vector<T>& a, const T& b ) {
            return NCPA::arrays::scale_vector<std::vector<T>>(
                a, std::pow( b, -1.0 ) );
        }

        /**
         * @brief Operator overload for unary vector negation.
         */
        template<typename T>
        std::vector<T> operator-( const std::vector<T>& a ) {
            return NCPA::arrays::scale_vector<std::vector<T>>(
                a, -NCPA::constants::one<T>() );
        }

        /**
         * @brief Stream operator for std::vector.
         */
        template<typename T>
        std::ostream& operator<<( std::ostream& os, const std::vector<T>& v ) {
            os << "{ ";
            for (auto it = v.cbegin(); it != v.cend(); ++it) {
                if (it != v.cbegin()) {
                    os << ", ";
                }
                os << *it;
            }
            os << " }";
            return os;
        }

        template<typename U, typename ShapeContainer>
        NDViewIterator<U, ShapeContainer> operator+(
            typename NDViewIterator<U, ShapeContainer>::difference_type n,
            NDViewIterator<U, ShapeContainer> it ) {
            it += n;
            return it;
        }


    }  // namespace arrays
}  // namespace NCPA
