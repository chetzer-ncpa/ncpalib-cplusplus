#pragma once

#include "NCPA/processing/packets/ResponsePacket.hpp"

#include <memory>

namespace NCPA {
    namespace processing {

        template<typename T>
        void swap( NCPA::processing::ProductPacket<T>& a,
                   NCPA::processing::ProductPacket<T>& b ) noexcept;

        template<typename T>
        class ProductPacket
            : public ResponsePacket,
              public std::enable_shared_from_this<ProductPacket<T>> {
            public:
                ProductPacket() : ResponsePacket( response_id_t::PRODUCT ) {
                    // _internal = std::make_shared<T>();
                }

                ProductPacket( const T in ) :
                    ResponsePacket( response_id_t::PRODUCT ) {
                    // _internal = std::make_shared<T>( in );
                    _internal.set( in );
                }

                ProductPacket( const T in, const std::string& tag ) :
                    ResponsePacket( response_id_t::PRODUCT, tag ) {
                    // _internal = std::make_shared<T>( in );
                    _internal.set( in );
                }

                ProductPacket( std::shared_ptr<T> in ) :
                    ResponsePacket( response_id_t::PRODUCT ) {
                    // _internal = in;
                    _internal.set( in );
                }

                ProductPacket( std::shared_ptr<T> in,
                               const std::string& tag ) :
                    ResponsePacket( response_id_t::PRODUCT, tag ) {
                    // _internal = in;
                    _internal.set( in );
                }

                // ProductPacket( const ProductPacket<T>& input ) :
                //     ProductPacket<T>() {
                //     _internal = std::make_unique<T>( input.get() );
                // }

                ProductPacket( ProductPacket<T>&& input ) noexcept :
                    ProductPacket<T>() {
                    ::swap( *this, input );
                }

                virtual ~ProductPacket() {}

                friend void swap<>( ProductPacket<T>& a,
                                    ProductPacket<T>& b ) noexcept;

                // ProductPacket<T>& operator=( ProductPacket<T> other ) {
                //     ::swap( *this, other );
                //     return *this;
                // }

                ProductPacket<T>& set( const T& input ) {
                    _internal.set( input );
                    return *this;
                }

                ProductPacket<T>& set( std::shared_ptr<T> input ) {
                    _internal = input;
                    return *this;
                }

                T& get() {
                    if (_internal) {
                        return _internal.get();
                    } else {
                        throw std::logic_error(
                            "ProductPacket: Nothing has been "
                            "assigned to internal pointer" );
                    }
                }

                const T& get() const {
                    if (_internal) {
                        return _internal.get();
                    } else {
                        throw std::logic_error(
                            "ProductPacket: Nothing has been "
                            "assigned to internal pointer" );
                    }
                }

                std::shared_ptr<T> *ptr() { return _internal.ptr(); }

            private:
                DataWrapper<T> _internal;
        };

        template<typename T>
        const T& get_product( const ResponsePacket& packet ) {
            auto packet_ptr
                = dynamic_cast<const ProductPacket<T> *>( &packet );
            if (packet_ptr) {
                return packet_ptr->get();
            } else {
                throw std::invalid_argument(
                    "Packet is not a ProductPacket!" );
            }
        }

        template<typename T>
        void swap( NCPA::processing::ProductPacket<T>& a,
                   NCPA::processing::ProductPacket<T>& b ) noexcept {
            using std::swap;
            swap( dynamic_cast<NCPA::processing::ResponsePacket&>( a ),
                  dynamic_cast<NCPA::processing::ResponsePacket&>( b ) );
            swap( a._internal, b._internal );
        }
    }  // namespace processing
}  // namespace NCPA
