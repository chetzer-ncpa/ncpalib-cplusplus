#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"

#include <memory>

namespace NCPA {
    namespace processing {

        template<typename T>
        void swap( NCPA::processing::GenericPacket<T>& a,
                   NCPA::processing::GenericPacket<T>& b ) noexcept;

        template<typename T>
        class GenericPacket : public InputPacket {
            public:
                GenericPacket() {}

                GenericPacket( const T in ) {}

                GenericPacket( const GenericPacket<T>& input ) :
                    GenericPacket<T>() {
                    _internal = std::unique_ptr<T>( input.get() );
                    // _internal = std::make_unique<T>( input.get() );
                }

                GenericPacket( GenericPacket<T>&& input ) noexcept :
                    GenericPacket<T>() {
                    ::swap( *this, input );
                }

                virtual ~GenericPacket() {}

                friend void swap<>( GenericPacket<T>& a,
                                    GenericPacket<T>& b ) noexcept;

                GenericPacket<T>& operator=( GenericPacket<T> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                GenericPacket<T>& set( const T& input ) {
                    // _internal = std::make_unique<T>( input );
                    _internal = std::unique_ptr<T>( input );
                    return *this;
                }

                T& get() {
                    if (_internal) {
                        return *_internal;
                    } else {
                        throw std::logic_error(
                            "GenericPacket: Nothing has been "
                            "assigned to internal pointer" );
                    }
                }

                const T& get() const {
                    if (_internal) {
                        return *_internal;
                    } else {
                        throw std::logic_error(
                            "GenericPacket: Nothing has been "
                            "assigned to internal pointer" );
                    }
                }

                T *ptr() { return _internal.get(); }

                static input_ptr_t build() {
                    return input_ptr_t( new GenericPacket<T>() );
                }

                static input_ptr_t build( const T in ) {
                    return input_ptr_t( new GenericPacket<T>( in ) );
                }

                static input_ptr_t build( const GenericPacket<T>& input ) {
                    return input_ptr_t( new GenericPacket<T>( input ) );
                }

            private:
                std::unique_ptr<T> _internal;
        };

        template<typename T>
        void swap( NCPA::processing::GenericPacket<T>& a,
                   NCPA::processing::GenericPacket<T>& b ) noexcept {
            using std::swap;
            swap( dynamic_cast<NCPA::processing::InputPacket&>( a ),
                  dynamic_cast<NCPA::processing::InputPacket&>( b ) );
            swap( a._internal, b._internal );
        }
    }  // namespace processing
}  // namespace NCPA
