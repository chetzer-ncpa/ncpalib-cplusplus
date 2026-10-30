#pragma once

#include <array>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace NCPA {
    namespace arrays {
        template<typename T>
        struct _ndstorage {
                std::vector<T> data;
        };

        template<size_t N, typename T>
        class ndarrayview;

        template<size_t N, typename T>
        class ndarray {
            public:
                template<typename... Dims>
                explicit ndarray( Dims... dims ) {
                    static_assert( sizeof...( Dims ) == N,
                                   "Wrong number of dimensions" );

                    size_t tmp[ N ] = { static_cast<size_t>( dims )... };

                    size_t total = 1;

                    for (size_t i = 0; i < N; ++i) {
                        _shape[ i ]  = tmp[ i ];
                        total       *= tmp[ i ];
                    }

                    compute_strides();

                    _storage.reset( new _ndstorage<T> );
                    _storage->data.resize( total );
                }

                template<typename... Indices>
                T& at( Indices... indices ) {
                    static_assert( sizeof...( Indices ) == N,
                                   "Wrong number of indices" );

                    std::array<size_t, N> idx
                        = { { static_cast<size_t>( indices )... } };

                    return _storage->data.at( linear_index( idx ) );
                }

                template<typename... Indices>
                const T& at( Indices... indices ) const {
                    static_assert( sizeof...( Indices ) == N,
                                   "Wrong number of indices" );

                    std::array<size_t, N> idx
                        = { { static_cast<size_t>( indices )... } };

                    return _storage->data.at( linear_index( idx ) );
                }

                ndarrayview<N - 1, T> slice( size_t dim, size_t fixed_index ) {
                    if (dim >= N) throw std::out_of_range( "bad dimension" );

                    if (fixed_index >= _shape[ dim ])
                        throw std::out_of_range( "bad index" );

                    std::array<size_t, N - 1> shape;
                    std::array<size_t, N - 1> strides;

                    size_t slice_ind = 0;

                    for (size_t i = 0; i < N; ++i) {
                        if (i != dim) {
                            shape[ slice_ind ]   = _shape[ i ];
                            strides[ slice_ind ] = _strides[ i ];
                            ++slice_ind;
                        }
                    }

                    size_t offset = fixed_index * _strides[ dim ];

                    return ndarrayview<N - 1, T>( _storage, offset, shape,
                                                  strides );
                }

                // std::vector<T>& vector() { return _storage->data; }

                // const std::vector<T>& vector() const { return _storage->data; }

                // T *data() { return _storage->data.data(); }

                // const T *data() const { return _storage->data.data(); }


            private:
                friend class ndarrayview<N, T>;

                std::shared_ptr<_ndstorage<T>> _storage;

                std::array<size_t, N> _shape;
                std::array<size_t, N> _strides;

                void compute_strides() {
                    _strides[ N - 1 ] = 1;

                    for (int i = (int)N - 2; i >= 0; --i) {
                        _strides[ i ] = _strides[ i + 1 ] * _shape[ i + 1 ];
                    }
                }

                size_t linear_index( const std::array<size_t, N>& idx ) const {
                    size_t result = 0;

                    for (size_t i = 0; i < N; ++i) {
                        if (idx[ i ] >= _shape[ i ]) {
                            throw std::out_of_range( "index" );
                        }
                        result += idx[ i ] * _strides[ i ];
                    }

                    return result;
                }
        };

        template<size_t N, typename T>
        class ndarrayview {
            public:
                ndarrayview( std::shared_ptr<_ndstorage<T>> storage,
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

                    return _storage->data.at( linear_index( idx ) );
                }

                template<typename... Indices>
                const T& at( Indices... indices ) const {
                    static_assert( sizeof...( Indices ) == N,
                                   "Wrong number of indices" );

                    std::array<size_t, N> idx
                        = { { static_cast<size_t>( indices )... } };

                    return _storage->data.at( linear_index( idx ) );
                }

                ndarrayview<N - 1, T> slice( size_t dim, size_t fixed_index ) {
                    if (dim >= N) {
                        throw std::out_of_range( "bad dimension" );
                    }

                    if (fixed_index >= _shape[ dim ]) {
                        throw std::out_of_range( "bad index" );
                    }

                    std::array<size_t, N - 1> shape;
                    std::array<size_t, N - 1> strides;

                    size_t slice_ind = 0;

                    for (size_t i = 0; i < N; ++i) {
                        if (i != dim) {
                            shape[ slice_ind ]   = _shape[ i ];
                            strides[ slice_ind ] = _strides[ i ];
                            ++slice_ind;
                        }
                    }

                    return ndarrayview<N - 1, T>(
                        _storage, _offset + fixed_index * _strides[ dim ],
                        shape, strides );
                }

                std::vector<T>& vector() { return _storage->data; }

                const std::vector<T>& vector() const { return _storage->data; }

                T *data() { return _storage->data.data(); }

                const T *data() const { return _storage->data.data(); }


            private:
                std::shared_ptr<_ndstorage<T>> _storage;

                size_t _offset;

                std::array<size_t, N> _shape;
                std::array<size_t, N> _strides;

                size_t linear_index( const std::array<size_t, N>& idx ) const {
                    size_t result = _offset;

                    for (size_t i = 0; i < N; ++i) {
                        if (idx[ i ] >= _shape[ i ]) {
                            throw std::out_of_range( "index" );
                        }

                        result += idx[ i ] * _strides[ i ];
                    }

                    return result;
                }
        };

        template<typename T>
        class ndarrayview<1, T> {
            public:
                ndarrayview( std::shared_ptr<_ndstorage<T>> storage,
                             size_t offset, const std::array<size_t, 1>& shape,
                             const std::array<size_t, 1>& strides ) :
                    _storage( storage ),
                    _offset( offset ),
                    _shape( shape ),
                    _strides( strides ) {}

                T& at( size_t i ) {
                    if (i >= _shape[ 0 ]) {
                        throw std::out_of_range( "index" );
                    }
                    return _storage->data.at( _offset + i * _strides[ 0 ] );
                }

                const T& at( size_t i ) const {
                    if (i >= _shape[ 0 ]) {
                        throw std::out_of_range( "index" );
                    }

                    return _storage->data.at( _offset + i * _strides[ 0 ] );
                }

                // std::vector<T>& vector() { return _storage->data; }

                // const std::vector<T>& vector() const { return _storage->data; }

                // T *data() { return _storage->data.data(); }

                // const T *data() const { return _storage->data.data(); }


            private:
                std::shared_ptr<_ndstorage<T>> _storage;

                size_t _offset;

                std::array<size_t, 1> _shape;
                std::array<size_t, 1> _strides;
        };

        template<typename T>
        class ndvectorview {
            public:
                ndvectorview( std::shared_ptr<_ndstorage<T>> storage,
                              size_t offset, std::vector<size_t> shape,
                              std::vector<size_t> strides ) :
                    _storage( std::move( storage ) ),
                    _offset( offset ),
                    _shape( std::move( shape ) ),
                    _strides( std::move( strides ) ) {}

                size_t rank() const { return _shape.size(); }

                T& at( std::initializer_list<size_t> indices ) {
                    return _storage->data.at( linear_index( indices ) );
                }

                const T& at( std::initializer_list<size_t> indices ) const {
                    return _storage->data.at( linear_index( indices ) );
                }

                T& at( const std::vector<size_t>& indices ) {
                    return _storage->data.at( linear_index( indices ) );
                }

                const T& at( const std::vector<size_t>& indices ) const {
                    return _storage->data.at( linear_index( indices ) );
                }

                ndvectorview slice( size_t dim, size_t fixed_index ) const {
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

                    return ndvectorview(
                        _storage, _offset + fixed_index * _strides[ dim ],
                        std::move( new_shape ), std::move( new_strides ) );
                }

                const std::vector<size_t>& shape() const { return _shape; }

                // std::vector<T>& vector() { return _storage->data; }

                // const std::vector<T>& vector() const { return _storage->data; }

                // T *data() { return _storage->data.data(); }

                // const T *data() const { return _storage->data.data(); }


            private:
                std::shared_ptr<_ndstorage<T>> _storage;

                size_t _offset;

                std::vector<size_t> _shape;
                std::vector<size_t> _strides;

                size_t linear_index(
                    std::initializer_list<size_t> idx ) const {
                    if (idx.size() != rank()) {
                        throw std::invalid_argument(
                            "wrong number of indices" );
                    }

                    size_t result = _offset;

                    size_t i = 0;
                    for (size_t v : idx) {
                        if (v >= _shape[ i ]) {
                            throw std::out_of_range( "index" );
                        }

                        result += v * _strides[ i ];
                        ++i;
                    }

                    return result;
                }

                size_t linear_index( const std::vector<size_t>& idx ) const {
                    if (idx.size() != rank()) {
                        throw std::invalid_argument(
                            "wrong number of indices" );
                    }

                    size_t result = _offset;

                    for (size_t i = 0; i < rank(); ++i) {
                        if (idx[ i ] >= _shape[ i ]) {
                            throw std::out_of_range( "index" );
                        }

                        result += idx[ i ] * _strides[ i ];
                    }

                    return result;
                }
        };

        template<typename T>
        class ndvector {
            public:
                explicit ndvector( std::vector<size_t> shape ) :
                    _shape( std::move( shape ) ) {
                    if (_shape.empty()) {
                        throw std::invalid_argument( "rank must be > 0" );
                    }

                    compute_strides();

                    size_t total = 1;

                    for (size_t d : _shape) {
                        total *= d;
                    }

                    _storage = std::make_shared<_ndstorage<T>>();
                    _storage->data.resize( total );
                }

                ndvector( std::initializer_list<size_t> shape ) :
                    ndvector<T>( std::vector<size_t>( shape ) ) {}

                size_t rank() const { return _shape.size(); }

                T& at( std::initializer_list<size_t> indices ) {
                    return _storage->data.at( linear_index( indices ) );
                }

                const T& at( std::initializer_list<size_t> indices ) const {
                    return _storage->data.at( linear_index( indices ) );
                }

                T& at( const std::vector<size_t>& indices ) {
                    return _storage->data.at( linear_index( indices ) );
                }

                const T& at( const std::vector<size_t>& indices ) const {
                    return _storage->data.at( linear_index( indices ) );
                }

                template<typename... Indices>
                T& at( Indices... indices ) {
                    std::vector<size_t> idx { static_cast<size_t>(
                        indices )... };
                    return _storage->data.at( linear_index( idx ) );
                }

                template<typename... Indices>
                const T& at( Indices... indices ) const {
                    std::vector<size_t> idx { static_cast<size_t>(
                        indices )... };
                    return _storage->data.at( linear_index( idx ) );
                }

                ndvectorview<T> slice( size_t dim, size_t fixed_index ) const {
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

                    return ndvectorview<T>(
                        _storage, fixed_index * _strides[ dim ],
                        std::move( shape ), std::move( strides ) );
                }

                const std::vector<size_t>& shape() const { return _shape; }

                std::vector<T>& vector() { return _storage->data; }

                const std::vector<T>& vector() const { return _storage->data; }

                T *data() { return _storage->data.data(); }

                const T *data() const { return _storage->data.data(); }

            private:
                std::shared_ptr<_ndstorage<T>> _storage;

                std::vector<size_t> _shape;
                std::vector<size_t> _strides;

                void compute_strides() {
                    _strides.resize( rank() );

                    _strides.back() = 1;

                    for (int i = static_cast<int>( rank() ) - 2; i >= 0; --i) {
                        _strides[ i ] = _strides[ i + 1 ] * _shape[ i + 1 ];
                    }
                }

                size_t linear_index(
                    std::initializer_list<size_t> idx ) const {
                    if (idx.size() != rank()) {
                        throw std::invalid_argument(
                            "wrong number of indices" );
                    }

                    size_t result = 0;

                    size_t i = 0;
                    for (size_t v : idx) {
                        if (v >= _shape[ i ]) {
                            throw std::out_of_range( "index" );
                        }

                        result += v * _strides[ i ];
                        ++i;
                    }

                    return result;
                }

                size_t linear_index( const std::vector<size_t>& idx ) const {
                    if (idx.size() != rank()) {
                        throw std::invalid_argument(
                            "wrong number of indices" );
                    }

                    size_t result = 0;

                    for (size_t i = 0; i < rank(); ++i) {
                        if (idx[ i ] >= _shape[ i ]) {
                            throw std::out_of_range( "index" );
                        }

                        result += idx[ i ] * _strides[ i ];
                    }

                    return result;
                }
        };

        template<typename T>
        using vector3d_t = ndarray<3, T>;

        template<typename T>
        using vector2d_t = ndarray<2, T>;
    }  // namespace arrays
}  // namespace NCPA
