#pragma once

#include "NCPA/processing/AbstractProcessingStep.hpp"
#include "NCPA/processing/declarations.hpp"
#include "NCPA/stateful.hpp"

#include <type_traits>


namespace NCPA {
    namespace processing {

        template<typename this_t, typename state_t>
        class StatefulProcessingStep : public virtual AbstractProcessingStep,
                                       public Stateful<state_t> {
                // static_assert( std::__is_base_of<int, int>::value, "check" );

                // static_assert(
                //     std::is_base_of<AbstractProcessingStep, this_t>::value,
                //     "this_t must derive from AbstractProcessingStep" );

            public:
                StatefulProcessingStep() {}

                virtual ~StatefulProcessingStep() {}

            protected:
                virtual response_ptr_t _build_state_packet() const override {
                    return response_ptr_t( new StatePacket( this ) );
                }

                virtual packet_processing_result_t
                    _process_state_request_packet(
                        const StateRequestPacket *packet_ptr,
                        std::string& message ) override {
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
