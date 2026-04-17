#pragma once

#include "NCPA/processing/AbstractDataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"

template<typename T>
void swap( NCPA::processing::DataWrapper<T>& a,
           NCPA::processing::DataWrapper<T>& b ) noexcept;

namespace NCPA {
    namespace processing {
        template<typename T>
        class DataWrapper
            : public AbstractDataWrapper,
              public std::enable_shared_from_this<DataWrapper<T>> {
            public:
                DataWrapper() { 
                    // _internal = std::make_shared<T>();
                 }

                DataWrapper( const T in ) {
                    _internal = std::make_shared<T>( in );
                }

                // DataWrapper( const DataWrapper<T>& input ) :
                // DataWrapper<T>() {
                //     _internal = std::make_unique<T>( input.get() );
                // }

                DataWrapper( std::shared_ptr<T>& input ) { _internal = input; }

                DataWrapper( DataWrapper<T>&& input ) noexcept :
                    DataWrapper<T>() {
                    ::swap( *this, input );
                }

                virtual ~DataWrapper() {}

                friend void swap<>( DataWrapper<T>& a,
                                    DataWrapper<T>& b ) noexcept;

                DataWrapper<T>& operator=( DataWrapper<T> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                DataWrapper<T>& set( const T& input ) {
                    _internal = std::make_shared<T>( input );
                    return *this;
                }

                DataWrapper<T>& set( std::shared_ptr<T> input ) {
                    _internal = input;
                    return *this;
                }

                T& get() {
                    if (_internal) {
                        return *_internal;
                    } else {
                        throw std::logic_error(
                            "DataWrapper: Nothing has been "
                            "assigned to internal pointer" );
                    }
                }

                const T& get() const {
                    if (_internal) {
                        return *_internal;
                    } else {
                        throw std::logic_error(
                            "DataWrapper: Nothing has been "
                            "assigned to internal pointer" );
                    }
                }

                std::shared_ptr<T> ptr() {
                    return _internal;
                    // return _internal.get();
                    // return std::shared_ptr<T>( this->shared_from_this(),
                    //                            _internal.get() );
                }

                const std::shared_ptr<T> ptr() const {
                    return _internal;
                    // return _internal.get();
                    // return std::shared_ptr<T>( this->shared_from_this(),
                    //                            _internal.get() );
                }

                explicit operator bool() const { return (bool)_internal; }

            private:
                // std::unique_ptr<T> _internal;
                std::shared_ptr<T> _internal;
        };
    }  // namespace processing
}  // namespace NCPA

template<typename T>
void swap( NCPA::processing::DataWrapper<T>& a,
           NCPA::processing::DataWrapper<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::processing::AbstractDataWrapper&>( a ),
          static_cast<NCPA::processing::AbstractDataWrapper&>( b ) );
    swap( a._internal, b._internal );
}
