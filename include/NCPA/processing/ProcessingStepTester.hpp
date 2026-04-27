#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets.hpp"
#include "NCPA/processing/ProcessingStep.hpp"

#include <vector>

namespace NCPA {
    namespace processing {
        class ProcessingStepTester {
            public:
                ProcessingStepTester() {}

                template<typename INTYPE, typename OUTTYPE>
                response_ptr_t test( ProcessingStep<INTYPE, OUTTYPE> *step_ptr,
                              const std::vector<INTYPE>& input,
                              const OUTTYPE& expected ) {
                    response_ptr_t resp;
                    for (auto in : input) {
                        DataPacket<INTYPE> packet( in );
                        resp = step_ptr->process( packet );
                    }
                    return resp;
                }
        };
    }  // namespace processing
}  // namespace NCPA
