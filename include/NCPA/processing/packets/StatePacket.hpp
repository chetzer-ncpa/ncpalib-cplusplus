#pragma once

#include "NCPA/processing/declarations.hpp"
// #include "NCPA/processing/packets.hpp"
#include "NCPA/processing/packets/ResponsePacket.hpp"
#include "NCPA/stateful.hpp"

#include <memory>

namespace NCPA {
    namespace processing {
        template<typename T>
        class StatePacket
            : public ProductPacket<T>,
              public std::enable_shared_from_this<StatePacket<T>> {
            public:
                StatePacket() : ProductPacket<T>() {
                    this->id() = response_id_t::STATE;
                }

                StatePacket( const T& in ) : ProductPacket<T>( in ) {}

                StatePacket( const Stateful<T>& s ) : StatePacket<T>( s.state() ) {}

                // StatePacket( const AbstractProcessingStep& in ) :
                //     ResponsePacket( response_id_t::STATE, in.tag() ),
                //     _ptr { &in } {}

                StatePacket( const StatePacket<T>& other ) :
                    ProductPacket<T>( other ) {}

                StatePacket( StatePacket<T>&& input ) noexcept : StatePacket<T>() {
                    swap( *this, input );
                }

                virtual ~StatePacket() {}

                friend void swap( StatePacket<T>& a, StatePacket<T>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<ProductPacket<T>&>( a ),
                          static_cast<ProductPacket<T>&>( b ) );
                }

                StatePacket& operator=( StatePacket<T> other ) {
                    swap( *this, other );
                    return *this;
                }
        };

        // const AbstractProcessingStep& get_state_pointer(
        //     const ResponsePacket& packet ) {
        //     auto packet_ptr = dynamic_cast<const StatePacket *>( &packet );
        //     if (packet_ptr) {
        //         return packet_ptr->get();
        //     } else {
        //         throw std::invalid_argument( "Packet is not a StatePacket!" );
        //     }
        // }

        // template<class T>
        // const T& get_state_pointer_as( const ResponsePacket& packet ) {
        //     T *ptr = dynamic_cast<const T *>( &get_state_pointer( packet ) );
        //     if (ptr == nullptr) {
        //         throw std::invalid_argument(
        //             "Can't cast step pointer to requested type!" );
        //     }
        //     return *ptr;
        // }

        // void swap( StatePacket& a, StatePacket& b ) noexcept {
        //     using std::swap;
        //     swap( dynamic_cast<ResponsePacket&>( a ),
        //             dynamic_cast<ResponsePacket&>( b ) );
        //     swap( a._ptr, b._ptr );
        // }
    }  // namespace processing
}  // namespace NCPA
