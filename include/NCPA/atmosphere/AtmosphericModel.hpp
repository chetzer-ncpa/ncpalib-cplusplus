#pragma once

static void swap( NCPA::atmos::AtmosphericModel&,
                  NCPA::atmos::AtmosphericModel& ) noexcept;

namespace NCPA {
    namespace atmos {
        class AtmosphericModel {
            public:
                virtual ~AtmosphericModel() = default;
                virtual bool is_stratified() const = 0;
        };

    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphericModel& a,
                  NCPA::atmos::AtmosphericModel& b ) noexcept {}
