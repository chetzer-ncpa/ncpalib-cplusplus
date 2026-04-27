#pragma once

#include "NCPA/processing/declarations.hpp"

namespace NCPA {
    namespace processing {
        class _step_state {
            public:
                _step_state() {}

                _step_state( const _step_state& other ) {}

                _step_state( _step_state&& other ) noexcept : _step_state() {
                    swap( *this, other );
                }

                virtual ~_step_state() {}

                friend void swap( _step_state& a, _step_state& b ) noexcept;
        };

        void swap( _step_state& a, _step_state& b ) noexcept {}
    }  // namespace processing
}  // namespace NCPA
