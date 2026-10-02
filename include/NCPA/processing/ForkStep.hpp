#pragma once

#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/packets.hpp"
#include "NCPA/processing/ProcessingStep.hpp"

namespace NCPA {
    namespace processing {

        template<typename intype, typename outtype>
        class ForkStep : public ProcessingStep<intype, outtype> {
            public:
                ForkStep() = default;
                ForkStep( AbstractProcessingStep *if_true,
                          AbstractProcessingStep *if_false );
                virtual ~ForkStep                                  = default;
                ForkStep( const ForkStep<intype, outtype>& other ) = default;
                ForkStep( ForkStep<intype, outtype>&& other ) noexcept
                    = default;

                friend void swap( ForkStep<intype, outtype>& a,
                                  ForkStep<intype, outtype>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<ProcessingStep<intype, outtype>&>( a ),
                          static_cast<ProcessingStep<intype, outtype>&>( b ) );
                }

                virtual bool test() const = 0;

                response_ptr_t process_data_packet(
                    InputPacket& packet ) override {
                    NCPA_DEBUG << this->tag() << ": processing data packet"
                               << std::endl;
                    std::string msg;
                    switch (this->_process_data_packet( packet, msg )) {
                        case packet_processing_result_t::PACKET_NOT_APPLICABLE:
                            return this->pass_to_next(
                                packet, this->response(
                                            response_id_t::ERROR,
                                            "Data packet reached last step "
                                            "without being handled!" ) );
                        case packet_processing_result_t::ERROR:
                            return this->response( response_id_t::ERROR, msg );
                        case packet_processing_result_t::FAILURE_NO_PRODUCT:
                        case packet_processing_result_t::FAILURE_PRODUCT:
                            return this->if_false( packet );
                        case packet_processing_result_t::SUCCESS_NO_PRODUCT:
                        case packet_processing_result_t::SUCCESS_PRODUCT:
                            return this->if_true( packet );
                        default:
                            return return_code_unsupported(
                                "process_data_packet" );
                    }
                }

                virtual response_ptr_t if_true( InputPacket& packet ) {
                    if (_branch_if_true != nullptr) {
                        return _branch_if_true->process( packet );
                    } else {
                        return this->_build_product_packet_from_input();
                    }
                }

                virtual response_ptr_t if_false( InputPacket& packet ) {
                    if (_branch_if_false != nullptr) {
                        return _branch_if_false->process( packet );
                    } else {
                        return this->_build_product_packet_from_input();
                    }
                }

            protected:
                AbstractProcessingStep *_branch_if_true  = nullptr;
                AbstractProcessingStep *_branch_if_false = nullptr;

                packet_processing_result_t _process_input() override {
                    return (
                        this->test()
                            ? packet_processing_result_t::SUCCESS_NO_PRODUCT
                            : packet_processing_result_t::FAILURE_NO_PRODUCT );
                }
        };
    }  // namespace processing
}  // namespace NCPA
