#pragma once

#include "NCPA/atmosphere/declarations.hpp"

#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace atmos {
        class abstract_atmospheric_property;
    }
}  // namespace NCPA

#define RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY \
    return static_cast<abstract_atmospheric_property&>( *this );

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
                virtual bool strict() const { return _strict; }
                virtual abstract_atmospheric_property& set_strict( bool tf ) {
                    _strict = tf;
                    return *this;
                }
            protected:
                bool _strict = false;
                virtual void _validate_axis( size_t n ) const {
                    if (n >= this->dimensions()) {
                        std::ostringstream oss;
                        oss << "Invalid axis " << n;
                        throw std::range_error( oss.str() );
                    }
                }
        };

        typedef std::unique_ptr<abstract_atmospheric_property> _atm_prop_ptr_t;
    }  // namespace atmos
}  // namespace NCPA
