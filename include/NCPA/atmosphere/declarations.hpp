#pragma once

#include "NCPA/interpolation.hpp"
#include "NCPA/units.hpp"

#include <unordered_map>
#include <vector>

#ifdef NCPA_INTERPOLATION_GSL_STEFFEN_SPLINE_AVAILABLE
#  if NCPA_INTERPOLATION_GSL_STEFFEN_SPLINE_AVAILABLE
#    define NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR \
        NCPA::interpolation::interpolator_1d_type_t::GSL_STEFFEN
#  else
#    define NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR \
        NCPA::interpolation::interpolator_1d_type_t::LANL_CUBIC
#  endif
#else
#  define NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR \
      NCPA::interpolation::interpolator_1d_type_t::LANL_CUBIC
#endif

#ifndef NCPA_ATMOSPHERE_DEFAULT_2D_INTERPOLATOR
#  define NCPA_ATMOSPHERE_DEFAULT_2D_INTERPOLATOR \
      NCPA::interpolation::interpolator_2d_type_t::LANL_BICUBIC
#endif

#ifndef NCPA_ATMOSPHERE_DEFAULT_3D_INTERPOLATOR
#  define NCPA_ATMOSPHERE_DEFAULT_3D_INTERPOLATOR \
      NCPA::interpolation::interpolator_3d_type_t::LANL_HYBRID
#endif

#ifndef NCPA_ATMOSPHERE_DEFAULT_DISTANCE_UNITS
#  define NCPA_ATMOSPHERE_DEFAULT_DISTANCE_UNITS NCPA::units::KILOMETERS
#endif

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
        enum class atmospheric_property_3d_t { STRATIFIED, GRID };
        enum class atmosphere_1d_t { TUPLE };
        enum class atmosphere_2d_t { GRID, STRATIFIED, PIECEWISE_STRATIFIED };
        enum class atmosphere_3d_t { GRID, STRATIFIED };

        enum class reader_1d_t { NCPAPROP };
        enum class reader_2d_t {
            NCPAPROP,
            NCPAPROP_STRATIFIED,
            NCPAPROP_PIECEWISE_STRATIFIED
        };
        enum class reader_3d_t { NCPAPROP, NCPAPROP_STRATIFIED };
        // enum class extrapolation_t { ZERO, CONSTANT, LINEAR };

        // Public API
        class AtmosphericProperty1D;
        class AtmosphericProperty2D;
        class AtmosphericProperty3D;
        class AtmosphericModel;
        class Atmosphere1D;
        class Atmosphere2D;
        class Atmosphere3D;
        class AtmosphereReader1D;
        class AtmosphereReader2D;
        class AtmosphereReader3D;
        class AtmosphereFactory;


        // internal declarations and convenience typedefs
        typedef NCPA::units::ScalarWithUnits<double> scalar_u_t;
        typedef NCPA::units::VectorWithUnits<double> vector_u_t;
        typedef NCPA::units::Vector2DWithUnits<double> vector2d_u_t;
        typedef NCPA::units::Vector3DWithUnits<double> vector3d_u_t;
        typedef NCPA::units::units_ptr_t units_ptr_t;

        // atmospheres
        class abstract_atmosphere_1d;
        class tuple_atmosphere_1d;
        class abstract_atmosphere_2d;
        class stratified_atmosphere_2d;
        class piecewise_stratified_atmosphere_2d;
        class grid_atmosphere_2d;
        class abstract_atmosphere_3d;
        class grid_atmosphere_3d;
        class stratified_atmosphere_3d;

        // readers
        class _abstract_atmosphere_reader_1d;
        class ncpaprop_atmosphere_reader_1d;
        class _abstract_atmosphere_reader_2d;
        class ncpaprop_atmosphere_reader_2d;
        class ncpaprop_atmosphere_reader_stratified_2d;
        class _abstract_atmosphere_reader_3d;
        class ncpaprop_atmosphere_reader_3d;
        class ncpaprop_atmosphere_reader_stratified_3d;

        // atmospheric properties
        class abstract_atmospheric_property;
        class abstract_atmospheric_property_1d;
        class tuple_atmospheric_property_1d;
        class abstract_atmospheric_property_2d;
        class grid_atmospheric_property_2d;
        class stratified_atmospheric_property_2d;
        class piecewise_stratified_atmospheric_property_2d;
        class abstract_atmospheric_property_3d;
        class grid_atmospheric_property_3d;
        class stratified_atmospheric_property_3d;

        typedef std::unique_ptr<abstract_atmospheric_property> _atm_prop_ptr_t;
        typedef std::unique_ptr<abstract_atmospheric_property_1d>
            _atm_prop_1d_ptr_t;
        typedef std::unique_ptr<abstract_atmospheric_property_2d>
            _atm_prop_2d_ptr_t;
        typedef std::unique_ptr<abstract_atmospheric_property_3d>
            _atm_prop_3d_ptr_t;
    }  // namespace atmos
}  // namespace NCPA
