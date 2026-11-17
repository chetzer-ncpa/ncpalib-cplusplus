#pragma once

#include "NCPA/arrays/declarations.hpp"
#include "NCPA/arrays/details.hpp"
#include "NCPA/arrays/exceptions.hpp"
#include "NCPA/arrays/SmartViewIterator.hpp"

#include <array>
#include <memory>
#include <vector>

namespace NCPA {
    namespace arrays {
        template<size_t N, typename T>
        class SmartArray {
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
                explicit SmartArray( Dims... dims ) {
                    this->resize( dims... );
                }

                SmartArray( std::vector<size_t> dims ) {
                    this->resize( dims );
                }

                SmartArray( std::array<size_t, N> dims ) {
                    this->resize( dims );
                }

                SmartArray( std::initializer_list<size_t> dims ) {
                    this->resize( dims );
                }

                SmartArray( const SmartArray<N, T>& other ) {
                    this->copy( other );
                }

                SmartArray( const SmartArrayView<N, T>& other ) {
                    this->copy( other );
                }

                SmartArray( const SmartVector<T>& other ) {
                    this->copy( other );
                }

                SmartArray( const SmartVectorView<T>& other ) {
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

                SmartArray<N, T>& clear() { return this->set( 0 ); }

                size_t dim( size_t n ) const {
                    assert( n < rank() );
                    return _shape[ n ];
                }

                SmartArray<N, T>& fill( T val ) { return this->set( val ); }

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
                    return details::indices( *this );
                }

                SmartArrayView<N - 1, T> operator[]( size_t ind ) {
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
                SmartArray<N, T>& reshape( Dims... dims ) {
                    return this->resize( dims... );
                }

                template<typename... Dims>
                SmartArray<N, T>& resize( Dims... dims ) {
                    static_assert( sizeof...( Dims ) == N
                                       || sizeof...( Dims ) == 0,
                                   "Wrong number of dimensions" );
                    return this->resize( std::vector<size_t> {
                        static_cast<size_t>( dims )... } );
                }

                SmartArray<N, T>& resize(
                    std::initializer_list<size_t> dims ) {
                    return this->resize( std::vector<size_t> { dims } );
                }

                SmartArray<N, T>& resize( std::array<size_t, N> dims ) {
                    // size_t total = 1;
                    size_t oldsize = this->size();
                    for (size_t i = 0; i < N; ++i) {
                        _shape[ i ] = dims[ i ];
                        // total *= dims[i];
                    }
                    size_t total = details::total_size( _shape );
                    details::compute_strides( _shape, _strides );
                    // compute_strides();
                    if (total != oldsize) {
                        _storage.reset( new details::_ndstorage<T> );
                        _storage->data.resize( total );
                    }
                    return *this;
                }

                SmartArray<N, T>& resize( std::vector<size_t> dims ) {
                    if (dims.empty()) {
                        dims.resize( N );
                    }
                    if (dims.size() != N) {
                        throw InvalidRank(
                            "SmartArray cannot be redimensioned" );
                    }
                    size_t oldsize = this->size();
                    // size_t total = 1;

                    for (size_t i = 0; i < N; ++i) {
                        _shape[ i ] = dims[ i ];
                        // total *= dims[i];
                    }
                    size_t total = details::total_size( _shape );
                    details::compute_strides( _shape, _strides );
                    // compute_strides();
                    if (total != oldsize) {
                        _storage.reset( new details::_ndstorage<T> );
                        _storage->data.resize( total );
                    }
                    return *this;
                }

                template<typename ContainerB>
                bool same_shape_as( const ContainerB& b ) const {
                    return details::same_shape( *this, b );
                }

                SmartArray<N, T>& set( T val ) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it = val;
                    }
                    return *this;
                }

                SmartArray<N, T>& set( const std::vector<T>& vals ) {
                    details::deal_into( vals, *this );
                    return *this;
                }

                const std::array<size_t, N>& shape() const { return _shape; }

                size_t size() const { return (_storage ? _storage->data.size() : 0); }

                SmartArrayView<N - 1, T> slice( size_t dim,
                                                size_t fixed_index ) {
                    if (dim >= N)
                        throw DimensionOutOfRange(
                            "Dimension too high to slice" );

                    if (fixed_index >= _shape[ dim ]) throw IndexOutOfRange();

                    std::array<size_t, N - 1> shape;
                    std::array<size_t, N - 1> strides;

                    details::remove_dimension( _shape, _strides, dim, shape,
                                               strides );

                    size_t offset = fixed_index * _strides[ dim ];

                    return SmartArrayView<N - 1, T>( _storage, offset, shape,
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
                friend class SmartArrayView<N, T>;

                std::shared_ptr<details::_ndstorage<T>> _storage;

                shape_type _shape;
                shape_type _strides;
        };

        template<size_t N, typename T>
        class SmartArrayView {
            public:
                using iterator = SmartViewIterator<T, std::array<size_t, N>>;
                using const_iterator
                    = SmartViewIterator<const T, std::array<size_t, N>>;
                typedef std::reverse_iterator<iterator> reverse_iterator;
                typedef std::reverse_iterator<const_iterator>
                    const_reverse_iterator;
                typedef std::array<size_t, N> shape_type;
                typedef T element_t;

                SmartArrayView(
                    std::shared_ptr<details::_ndstorage<T>> storage,
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

                SmartArrayView<N, T>& clear() { return this->set( 0 ); }

                size_t dim( size_t n ) const {
                    assert( n < rank() );
                    return _shape[ n ];
                }

                SmartArrayView<N, T>& fill( T val ) {
                    return this->set( val );
                }

                template<typename... Indices>
                T get( Indices... indices ) const {
                    return this->at( indices... );
                }

                std::vector<shape_type> indices() const {
                    return details::indices( *this );
                }

                template<size_t M = N>
                typename std::enable_if<( M > 1 ),
                                        SmartArrayView<N - 1, T>>::type
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

                SmartArrayView<N, T>& set( T val ) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it = val;
                    }
                    return *this;
                }

                SmartArrayView<N, T>& set( const std::vector<T>& vals ) {
                    details::deal_into( vals, *this );
                    return *this;
                }

                const std::array<size_t, N>& shape() const { return _shape; }

                size_t size() const { return details::total_size( _shape ); }

                template<size_t M = N>
                typename std::enable_if<( M > 1 ),
                                        SmartArrayView<N - 1, T>>::type
                    slice( size_t dim, size_t fixed_index ) {
                    if (dim >= N) {
                        throw DimensionOutOfRange(
                            "Dimension too high to slice" );
                    }

                    if (fixed_index >= _shape[ dim ]) {
                        throw IndexOutOfRange();
                    }

                    std::array<size_t, N - 1> shape;
                    std::array<size_t, N - 1> strides;
                    details::remove_dimension( _shape, _strides, dim, shape,
                                               strides );

                    return SmartArrayView<N - 1, T>(
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
    }  // namespace arrays
}  // namespace NCPA
