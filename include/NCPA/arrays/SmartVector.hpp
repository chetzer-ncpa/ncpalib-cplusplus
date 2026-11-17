#pragma once

#include "NCPA/arrays/declarations.hpp"
#include "NCPA/arrays/details.hpp"
#include "NCPA/arrays/exceptions.hpp"
#include "NCPA/arrays/SmartViewIterator.hpp"

#include <vector>
#include <memory>

namespace NCPA {
    namespace arrays {
        template<typename T>
        class SmartVectorView {
            public:
                using iterator = SmartViewIterator<T, std::vector<size_t>>;
                using const_iterator
                    = SmartViewIterator<const T, std::vector<size_t>>;
                typedef std::reverse_iterator<iterator> reverse_iterator;
                typedef std::reverse_iterator<const_iterator>
                    const_reverse_iterator;
                typedef T element_t;
                typedef std::vector<size_t> shape_type;

                SmartVectorView( std::shared_ptr<details::_ndstorage<T>> storage,
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

                SmartVectorView<T>& clear() { return this->set( 0 ); }

                SmartVectorView<T>& fill( T val ) { return this->set( val ); }

                template<typename... Indices>
                T get( Indices... indices ) const {
                    return this->at( indices... );
                }

                std::vector<shape_type> indices() const {
                    return details::indices(*this);
                }

                SmartVectorView<T>& set( T val ) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it = val;
                    }
                    return *this;
                }

                SmartVectorView<T>& set( const std::vector<T>& vals ) {
                    details::deal_into( vals, *this );
                    return *this;
                }

                SmartVectorView slice( size_t dim, size_t fixed_index ) const {
                    if (dim >= rank()) {
                        throw DimensionOutOfRange( "Dimension too high to slice" );
                    }

                    if (fixed_index >= _shape[ dim ]) {
                        throw IndexOutOfRange();
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

                    return SmartVectorView(
                        _storage, _offset + fixed_index * _strides[ dim ],
                        std::move( new_shape ), std::move( new_strides ) );
                }

                const std::vector<size_t>& shape() const { return _shape; }

                SmartVectorView operator[]( size_t ind ) const {
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
        class SmartVector {
            public:
                typedef typename std::vector<T>::iterator iterator;
                typedef typename std::vector<T>::const_iterator const_iterator;
                typedef
                    typename std::vector<T>::reverse_iterator reverse_iterator;
                typedef typename std::vector<T>::const_reverse_iterator
                    const_reverse_iterator;
                typedef T element_t;
                typedef std::vector<size_t> shape_type;

                template<typename... Dims>
                explicit SmartVector( Dims... dims ) {
                    this->resize( dims... );
                }

                SmartVector( std::vector<size_t> shape ) {
                    this->resize( shape );
                }

                SmartVector( std::initializer_list<size_t> shape ) :
                    SmartVector<T>( std::vector<size_t>( shape ) ) {}

                SmartVector( const SmartVector<T>& other ) {
                    this->copy( other );
                }

                template<size_t N>
                SmartVector( const SmartArray<N, T>& other ) {
                    this->copy( other );
                }

                template<size_t N>
                SmartVector( const SmartArrayView<N, T>& other ) {
                    this->copy( other );
                }

                SmartVector( const SmartVectorView<T>& other ) {
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

                SmartVector<T>& clear() { return this->set( 0 ); }

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

                SmartVector<T>& fill( T val ) { return this->set( val ); }

                template<typename... Indices>
                T get( Indices... indices ) const {
                    return this->at( indices... );
                }

                std::vector<shape_type> indices() const {
                    return details::indices(*this);
                }

                SmartVectorView<T> operator[]( size_t ind ) const {
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
                SmartVector<T>& resize( Dims... dims ) {
                    return this->resize( std::vector<size_t> {
                        static_cast<size_t>( dims )... } );
                }

                SmartVector<T>& resize( std::vector<size_t> shape ) {
                    if (shape.empty()) {
                        _shape.clear();
                        _strides.clear();
                        _storage = std::make_shared<details::_ndstorage<T>>();
                    } else if (_shape != shape) {
                        size_t oldsize = this->size();
                        _shape = std::move( shape );
                        details::compute_strides( _shape, _strides );
                        size_t total = details::total_size( _shape );
                        if (total != oldsize)
                       { _storage = std::make_shared<details::_ndstorage<T>>();
                        _storage->data.resize( total );}
                    }

                    return *this;
                }

                SmartVector<T>& resize( std::initializer_list<size_t> shape ) {
                    return this->resize( std::vector<size_t>( shape ) );
                }

                template<size_t N>
                SmartVector<T>& resize( const std::array<size_t, N>& shape ) {
                    return this->resize(
                        std::vector<size_t>( shape.begin(), shape.end() ) );
                }

                template<typename ContainerB>
                bool same_shape_as( const ContainerB& b ) const {
                    return details::same_shape( *this, b );
                }

                SmartVector<T>& set( T val ) {
                    for (auto it = begin(); it != end(); ++it) {
                        *it = val;
                    }
                    return *this;
                }

                SmartVectorView<T>& set( const std::vector<T>& vals ) {
                    details::deal_into( vals, *this );
                    return *this;
                }

                const std::vector<size_t>& shape() const { return _shape; }

                size_t size() const { return _storage ? _storage->data.size() : 0; }

                SmartVectorView<T> slice( size_t dim, size_t fixed_index ) const {
                    if (dim >= rank()) {
                        throw DimensionOutOfRange( "Dimension too high to slice" );
                    }

                    if (fixed_index >= _shape[ dim ]) {
                        throw IndexOutOfRange();
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

                    return SmartVectorView<T>(
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
    }
}