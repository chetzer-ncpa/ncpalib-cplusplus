#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"

#include <memory>

template<typename T>
void swap( NCPA::processing::GenericPacket<T>& a,
           NCPA::processing::GenericPacket<T>& b ) noexcept;

namespace NCPA::processing {

    template<typename T>
    class GenericPacket : public InputPacket {
        public:
            GenericPacket() { _internal = std::make_unique<T>(); }

            GenericPacket( const T in ) { _internal = std::make_unique<T>( in ); }

            GenericPacket( const GenericPacket<T>& input ) : GenericPacket<T>() {
                _internal = std::make_unique<T>( input.get() );
            }

            GenericPacket( GenericPacket<T>&& input ) noexcept : GenericPacket<T>() {
                ::swap( *this, input );
            }

            virtual ~GenericPacket() {}

            friend void swap<>( GenericPacket<T>& a, GenericPacket<T>& b ) noexcept;

            GenericPacket<T>& operator=( GenericPacket<T> other ) {
                ::swap( *this, other );
                return *this;
            }

            GenericPacket<T>& set( const T& input ) {
                _internal = std::make_unique<T>( input );
                return *this;
            }

            T& get() {
                if (_internal) {
                    return *_internal;
                } else {
                    throw std::logic_error( "GenericPacket: Nothing has been "
                                            "assigned to internal pointer" );
                }
            }

            const T& get() const {
                if (_internal) {
                    return *_internal;
                } else {
                    throw std::logic_error( "GenericPacket: Nothing has been "
                                            "assigned to internal pointer" );
                }
            }

            T *ptr() { return _internal.get(); }

            static std::unique_ptr<InputPacket> build() {
                return std::unique_ptr<InputPacket>(
                    new GenericPacket<T>() );
            }

            static std::unique_ptr<InputPacket> build(const T in) {
                return std::unique_ptr<InputPacket>(
                    new GenericPacket<T>(in) );
            }

            static std::unique_ptr<InputPacket> build( const GenericPacket<T>& input ) {
                return std::unique_ptr<InputPacket>(
                    new GenericPacket<T>(input) );
            }

        private:
            std::unique_ptr<T> _internal;
    };
}  // namespace NCPA::processing

template<typename T>
void ::swap( NCPA::processing::GenericPacket<T>& a,
             NCPA::processing::GenericPacket<T>& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::processing::InputPacket&>( a ),
          dynamic_cast<NCPA::processing::InputPacket&>( b ) );
    swap( a._internal, b._internal );
}
