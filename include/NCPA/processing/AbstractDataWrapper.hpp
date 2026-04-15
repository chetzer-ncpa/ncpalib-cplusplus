#pragma once

#include "NCPA/processing/declarations.hpp"


void swap( NCPA::processing::AbstractDataWrapper& a,
           NCPA::processing::AbstractDataWrapper& b ) noexcept;

namespace NCPA {
    namespace processing {
        class AbstractDataWrapper {
            public:
                AbstractDataWrapper() {}

                virtual ~AbstractDataWrapper() {}

        };
    }  // namespace processing
}  // namespace NCPA

void swap( NCPA::processing::AbstractDataWrapper& a,
           NCPA::processing::AbstractDataWrapper& b ) noexcept {}
