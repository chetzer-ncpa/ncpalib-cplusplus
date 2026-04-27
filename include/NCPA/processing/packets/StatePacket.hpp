#pragma once

#include "NCPA/processing/AbstractProcessingStep.hpp"
#include "NCPA/processing/packets.hpp"
#include "NCPA/processing/packets/ResponsePacket.hpp"

#include <memory>

namespace NCPA {
    namespace processing {

        class StatePacket : public ResponsePacket,
                            public std::enable_shared_from_this<StatePacket> {
            public:
                StatePacket() :
                    ResponsePacket( response_id_t::STATE ), _ptr { nullptr } {}

                StatePacket( const AbstractProcessingStep& in ) :
                    ResponsePacket( response_id_t::STATE, in.tag() ),
                    _ptr { &in } {}

                StatePacket( const AbstractProcessingStep *in ) :
                    ResponsePacket( response_id_t::STATE,
                                    ( in ? in->tag() : "" ) ),
                    _ptr { in } {}

                StatePacket( const StatePacket& other ) :
                    ResponsePacket( other ), _ptr { other._ptr } {}

                StatePacket( StatePacket&& input ) noexcept : StatePacket() {
                    swap( *this, input );
                }

                virtual ~StatePacket() {}

                friend void swap( StatePacket& a, StatePacket& b ) noexcept;

                StatePacket& operator=( StatePacket other ) {
                    swap( *this, other );
                    return *this;
                }

                StatePacket& set( const AbstractProcessingStep& step ) {
                    _ptr = &step;
                    return *this;
                }

                StatePacket& set(
                    const std::unique_ptr<AbstractProcessingStep>& input ) {
                    _ptr = input.get();
                    return *this;
                }

                const AbstractProcessingStep& get() const {
                    if (_ptr) {
                        return *_ptr;
                    } else {
                        throw std::logic_error(
                            "StatePacket: Nothing has been "
                            "assigned to internal pointer" );
                    }
                }

            private:
                const AbstractProcessingStep *_ptr;
        };

        const AbstractProcessingStep& get_state_pointer(
            const ResponsePacket& packet ) {
            auto packet_ptr = dynamic_cast<const StatePacket *>( &packet );
            if (packet_ptr) {
                return packet_ptr->get();
            } else {
                throw std::invalid_argument( "Packet is not a StatePacket!" );
            }
        }

        void swap( StatePacket& a, StatePacket& b ) noexcept {
            using std::swap;
            swap( dynamic_cast<ResponsePacket&>( a ),
                  dynamic_cast<ResponsePacket&>( b ) );
            swap( a._ptr, b._ptr );
        }
    }  // namespace processing
}  // namespace NCPA
