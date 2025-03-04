#pragma once

#include "NCPA/atmosphere/declarations.hpp"

#include <memory>

namespace NCPA {
    namespace atmos {
        class abstract_atmospheric_property;
    }
}  // namespace NCPA

static void swap( NCPA::atmos::abstract_atmospheric_property&,
                  NCPA::atmos::abstract_atmospheric_property& ) noexcept;

namespace NCPA {
    namespace atmos {
        class abstract_atmospheric_property {
            public:
                virtual ~abstract_atmospheric_property() {}
                virtual size_t dimensions() const = 0;
                virtual std::unique_ptr<abstract_atmospheric_property> clone()
                    const
                    = 0;
                virtual abstract_atmospheric_property& copy(
                    const abstract_atmospheric_property& source )
                    = 0;
        };

        typedef std::unique_ptr<abstract_atmospheric_property> _atm_prop_ptr_t;
    }  // namespace atmos
}  // namespace NCPA
