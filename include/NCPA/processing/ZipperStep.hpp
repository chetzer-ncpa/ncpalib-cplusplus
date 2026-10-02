#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/ProcessingStep.hpp"

#include <tuple>

namespace NCPA {
    namespace processing {

        template<typename intype, typename outtype1, typename outtype2>
        class ZipperStep
            : public ProcessingStep<intype, std::tuple<outtype1, outtype2>> {
                using parent_t
                    = ProcessingStep<intype, std::tuple<outtype1, outtype2>>;

            public:
                explicit ZipperStep(
                    ProcessingStep<intype, outtype1>& step1,
                    ProcessingStep<intype, outtype2>& step2 ) :
                    parent_t() {
                    _step1 = &step1;
                    _step2 = &step2;
                }

                ZipperStep()          = default;
                virtual ~ZipperStep() = default;
                ZipperStep(
                    ZipperStep<intype, outtype1, outtype2>&& other ) noexcept
                    = default;


            protected:
                ProcessingStep<intype, outtype1> *_step1 = nullptr;
                ProcessingStep<intype, outtype2> *_step2 = nullptr;

                virtual void _define_parameters() override {}

                virtual packet_processing_result_t _process_input() override {
                    if (_step1 == nullptr || _step2 == nullptr) {
                        throw std::logic_error(
                            "One or more unset branches in ZipperStep" );
                    }
                    this->_flags.clear();
                    _step1._next = nullptr;
                    _step2._next = nullptr;
                    response_ptr_t resp1
                        = _step1->process( this->_ditto_input_packet() );
                    response_ptr_t resp2
                        = _step2->process( this->_ditto_input_packet() );
                    _product.set( std::make_tuple(
                        dynamic_cast<ProductPacket<outtype1> *>( resp1.get() )
                            ->get(),
                        dynamic_cast<ProductPacket<outtype2> *>( resp2.get() )
                            ->get() ) );
                    return packet_processing_result_t::SUCCESS_PRODUCT;
                }
        };
    }  // namespace processing
}  // namespace NCPA
