#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/ProcessingStep.hpp"
#include "NCPA/stateful.hpp"

#include <type_traits>

namespace NCPA {
    namespace processing {

        template<typename inout_t, typename state_t>
        class StatefulStep : public ProcessingStep<inout_t, inout_t>,
                             public Stateful<state_t> {
            public:
                using Stateful<state_t>::state;

                StatefulStep() = default;

                StatefulStep( const std::string& key ) :
                    ProcessingStep<inout_t, inout_t>( key ),
                    Stateful<state_t>() {}

                virtual ~StatefulStep() = default;
                StatefulStep( const StatefulStep<inout_t, state_t>& other )
                    = default;
                StatefulStep( StatefulStep<inout_t, state_t>&& other ) noexcept
                    = default;

                void swap( StatefulStep<inout_t, state_t>& a,
                           StatefulStep<inout_t, state_t>& b ) noexcept {
                    using std::swap;
                    swap(
                        static_cast<ProcessingStep<inout_t, inout_t>&>( a ),
                        static_cast<ProcessingStep<inout_t, inout_t>&>( b ) );
                    swap( static_cast<Stateful<state_t>&>( a ),
                          static_cast<Stateful<state_t>&>( b ) );
                }

            protected:
                virtual response_ptr_t _build_state_packet() const override {
                    return response_ptr_t( new StatePacket( this ) );
                }

                virtual packet_processing_result_t
                    _process_state_request_packet(
                        const StateRequestPacket *packet_ptr,
                        std::vector<std::string>& message ) override {
                    CHECK_PACKET_POINTER_NOT_NULL( packet_ptr )
                    if (packet_ptr->tag() == this->tag()) {
                        return packet_processing_result_t::SUCCESS_PRODUCT;
                    } else {
                        return packet_processing_result_t::
                            PACKET_NOT_APPLICABLE;
                    }
                }
        };
    }  // namespace processing
}  // namespace NCPA
