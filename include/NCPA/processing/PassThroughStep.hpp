#pragma once

#include "NCPA/DAMES/declarations.hpp"

#include <deque>

namespace NCPA::DAMES {
    using namespace NCPA::processing;

    template<typename contents_t>
    class PassThroughStep : public ProcessingStep<contents_t, contents_t> {
        public:
            PassThroughStep() : PassThroughStep( "" ) {}

            PassThroughStep( const std::string& tag ) :
                StatefulStep<contents_t, std::deque<contents_t>>( tag ) {
                _define_parameters();
            }

            virtual ~PassThroughStep() = default;
            PassThroughStep( const PassThroughStep<contents_t>& other )
                = default;
            PassThroughStep( PassThroughStep<contents_t>&& other ) noexcept
                = default;

        protected:
            virtual void _define_parameters() override {}

            virtual packet_processing_result_t _process_input() override {
                this->_flags.clear();
                _product.set( this->_input.get() );
                return packet_processing_result_t::SUCCESS_PRODUCT;
            }
    };


}  // namespace NCPA::DAMES
