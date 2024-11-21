#pragma once

#include "NCPA/units.hpp"

namespace NCPA {
    namespace atmos {
        typedef NCPA::units::ScalarWithUnits<double> scalar_t;
        typedef NCPA::units::VectorWithUnits<double> vector_t;
        typedef const NCPA::units::Unit* units_ptr_t;
    }
}