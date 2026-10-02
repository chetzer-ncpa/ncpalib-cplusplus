#pragma once

#include "NCPA/processing/declarations.hpp"

namespace NCPA {
    namespace processing {
        class AbstractDataWrapper {
            public:
                AbstractDataWrapper()                             = default;
                virtual ~AbstractDataWrapper()                    = default;
                AbstractDataWrapper( const AbstractDataWrapper& ) = default;
                AbstractDataWrapper( AbstractDataWrapper&& ) noexcept
                    = default;

                friend void swap( AbstractDataWrapper& a,
                                  AbstractDataWrapper& b ) noexcept {}
        };
    }  // namespace processing
}  // namespace NCPA
