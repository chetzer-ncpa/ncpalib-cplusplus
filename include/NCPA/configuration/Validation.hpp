#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <functional>
#include <string>
#include <utility>

namespace NCPA {
    namespace config {
        class Validation : public Cloneable<Validation> {
            public:
                Validation() {}

                virtual ~Validation() {}

                Validation( const Validation& other ) {}

                Validation( Validation&& other ) noexcept {
                    swap( *this, other );
                }

                friend void swap( Validation& a, Validation& b ) noexcept {}
        };
    }  // namespace config
}  // namespace NCPA
