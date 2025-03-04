#pragma once

#include "NCPA/interpolation.hpp"
#include "NCPA/units.hpp"

#include <unordered_map>
#include <vector>

namespace NCPA {
    namespace atmos {

        class constants {
            public:
                constexpr static double GAMMA() { return 1.4; }

                constexpr static double R() { return 287.0; }
        };

        enum class dimension_t { X, Y, Z, VERTICAL, RADIAL, AZIMUTHAL, POLAR };
        enum class atmospheric_property_1d_t { TUPLE };
        enum class atmospheric_property_2d_t {
            STRATIFIED,
            PIECEWISE_STRATIFIED,
            GRID
        };
        enum class atmosphere_1d_t { TUPLE };
        enum class atmosphere_2d_t { GRID, STRATIFIED, PIECEWISE_STRATIFIED };
        // enum class extrapolation_t { ZERO, CONSTANT, LINEAR };

        typedef NCPA::units::ScalarWithUnits<double> scalar_u_t;
        typedef NCPA::units::VectorWithUnits<double> vector_u_t;
        typedef NCPA::units::Vector2DWithUnits<double> vector2d_u_t;
        typedef const NCPA::units::Unit *units_ptr_t;

        class abstract_atmosphere_1d;
        class tuple_atmosphere_1d;
        class abstract_atmosphere_2d;
        class stratified_atmosphere_2d;
        class grid_atmosphere_2d;
        class abstract_atmosphere_3d;

        // atmospheric properties
        class abstract_atmospheric_property;
        class AtmosphericProperty1D;
        // class ValuePairAtmosphericProperty1D;
        // class VerticalProperty1D;
        class AtmosphericProperty2D;
        class AtmosphericProperty3D;

        // public API classes
        class AtmosphericModel;
        class Atmosphere1D;
        class Atmosphere2D;
        class Atmosphere3D;

        class AtmosphereFactory;
    }  // namespace atmos
}  // namespace NCPA
