#pragma once

#include "NCPA/arrays/declarations.hpp"
#include "NCPA/arrays/details.hpp"
#include <memory>
#include <cstddef>
#include <type_traits>

namespace NCPA {
    namespace arrays {
        template<typename T, typename ShapeContainer>
        class SmartViewIterator {
            public:
                typedef std::random_access_iterator_tag iterator_category;
                typedef T value_type;
                typedef std::ptrdiff_t difference_type;
                typedef T *pointer_t;
                typedef T& reference_t;

                template<typename, typename>
                friend class SmartViewIterator;

                SmartViewIterator( std::shared_ptr<details::_ndstorage<
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
                SmartViewIterator(
                    const SmartViewIterator<U, ShapeContainer>& other ) :
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

                SmartViewIterator& operator++() {
                    ++_linear_pos;
                    return *this;
                }

                SmartViewIterator operator++( int ) {
                    SmartViewIterator tmp( *this );
                    ++( *this );
                    return tmp;
                }

                SmartViewIterator& operator--() {
                    --_linear_pos;
                    return *this;
                }

                SmartViewIterator operator--( int ) {
                    SmartViewIterator tmp( *this );
                    --( *this );
                    return tmp;
                }

                bool operator==( const SmartViewIterator& rhs ) const {
                    return _storage.get() == rhs._storage.get()
                        && _offset == rhs._offset
                        && _linear_pos == rhs._linear_pos;
                }

                bool operator!=( const SmartViewIterator& rhs ) const {
                    return !( *this == rhs );
                }

                template<typename U>
                bool operator==(
                    const SmartViewIterator<U, ShapeContainer>& rhs ) const {
                    return _storage.get() == rhs._storage.get()
                        && _offset == rhs._offset
                        && _linear_pos == rhs._linear_pos;
                }

                template<typename U>
                bool operator!=(
                    const SmartViewIterator<U, ShapeContainer>& rhs ) const {
                    return !( *this == rhs );
                }

                template<typename U>
                bool operator<(
                    const SmartViewIterator<U, ShapeContainer>& rhs ) const {
                    assert( _storage.get() == rhs._storage.get() );
                    return _linear_pos < rhs._linear_pos;
                }

                template<typename U>
                bool operator>(
                    const SmartViewIterator<U, ShapeContainer>& rhs ) const {
                    assert( _storage.get() == rhs._storage.get() );
                    return rhs < *this;
                }

                template<typename U>
                bool operator<=(
                    const SmartViewIterator<U, ShapeContainer>& rhs ) const {
                    assert( _storage.get() == rhs._storage.get() );
                    return !( rhs < *this );
                }

                template<typename U>
                bool operator>=(
                    const SmartViewIterator<U, ShapeContainer>& rhs ) const {
                    assert( _storage.get() == rhs._storage.get() );
                    return !( *this < rhs );
                }

                template<typename U>
                difference_type operator-(
                    const SmartViewIterator<U, ShapeContainer>& rhs ) const {
                    return static_cast<difference_type>( _linear_pos )
                         - static_cast<difference_type>( rhs._linear_pos );
                }

                reference_t operator[]( difference_type n ) const {
                    return *( *this + n );
                }

                SmartViewIterator& operator+=( difference_type n ) {
                    _linear_pos += n;
                    return *this;
                }

                SmartViewIterator& operator-=( difference_type n ) {
                    _linear_pos -= n;
                    return *this;
                }

                SmartViewIterator operator+( difference_type n ) const {
                    SmartViewIterator tmp( *this );
                    tmp += n;
                    return tmp;
                }

                SmartViewIterator operator-( difference_type n ) const {
                    SmartViewIterator tmp( *this );
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
    }
}