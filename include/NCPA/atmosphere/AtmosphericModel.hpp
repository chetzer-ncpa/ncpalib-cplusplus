#pragma once

namespace NCPA {
    namespace atmos {
        class AtmosphericModel {
            public:
                virtual ~AtmosphericModel() = default;

                friend void swap( AtmosphericModel& a,
                                  AtmosphericModel& b ) noexcept {}

                virtual bool is_stratified() const = 0;
        };

    }  // namespace atmos
}  // namespace NCPA
