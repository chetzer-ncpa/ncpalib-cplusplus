/**
 * NCPA pointers library
 * @version 1.0.0
 * @author Claus Hetzer
 * @date 2026-04-14
 *
 * Implements a generic base for cloneable members.
 */
#pragma once

#define NCPA_CLONE_METHOD( THISTYPE, PARENTTYPE )                     \
    std::unique_ptr<PARENTTYPE> clone() const override {              \
        return std::unique_ptr<PARENTTYPE>(                           \
            new THISTYPE( dynamic_cast<const THISTYPE&>( *this ) ) ); \
    }

#include <memory>
#include <utility>

namespace NCPA {
    template<typename BASE>
    class Cloneable {
        public:
            virtual ~Cloneable() {}

            virtual std::unique_ptr<BASE> clone() const = 0;
    };
}  // namespace NCPA
