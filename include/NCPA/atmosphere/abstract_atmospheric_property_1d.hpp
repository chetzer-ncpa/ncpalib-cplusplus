#pragma once

/**
 * @file
 * @brief
 */

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <memory>
#include <stdexcept>

/**
 * @brief
 * @param a
 * @param b
 */
static void swap( NCPA::atmos::abstract_atmospheric_property_1d& a,
                  NCPA::atmos::abstract_atmospheric_property_1d& b ) noexcept;

namespace NCPA {
    namespace atmos {

        /**
         * @brief
         */
        class abstract_atmospheric_property_1d
            : public abstract_atmospheric_property {
            public:
                /**
                 * @brief
                 */
                abstract_atmospheric_property_1d() {}

                /**
                 * @brief
                 */
                virtual ~abstract_atmospheric_property_1d() {}

                /**
                 * @brief
                 * @param a
                 * @param b
                 */
                friend void ::swap(
                    abstract_atmospheric_property_1d& a,
                    abstract_atmospheric_property_1d& b ) noexcept;

                /**
                 * @brief
                 * @return size_t
                 */
                virtual size_t size() const = 0;

                /**
                 * @brief
                 * @return std::unique_ptr<abstract_atmospheric_property_1d>
                 */
                virtual std::unique_ptr<abstract_atmospheric_property_1d>
                    clone1d() const = 0;

                /**
                 * @brief
                 * @return abstract_atmospheric_property_1d&
                 */
                virtual abstract_atmospheric_property_1d& clear() = 0;

                /**
                 * @brief
                 * @param interp_type
                 * @return abstract_atmospheric_property_1d&
                 */
                virtual abstract_atmospheric_property_1d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;

                /**
                 * @brief
                 * @param extrap_type
                 * @return abstract_atmospheric_property_1d&
                 */
                virtual abstract_atmospheric_property_1d& set_extrapolator(
                    NCPA::interpolation::extrapolator_1d_type_t extrap_type )
                    = 0;

                /**
                 * @brief
                 * @return NCPA::interpolation::Interpolator1D<double, double>&
                 */
                virtual NCPA::interpolation::Interpolator1D<double, double>&
                    interpolator() = 0;

                /**
                 * @brief
                 * @param z
                 * @param units
                 * @return abstract_atmospheric_property_1d&
                 */
                virtual abstract_atmospheric_property_1d& set(
                    const vector_u_t& z, const units_ptr_t units ) = 0;

                /**
                 * @brief
                 * @param z
                 * @param p
                 * @return abstract_atmospheric_property_1d&
                 */
                virtual abstract_atmospheric_property_1d& set(
                    const vector_u_t& z, const vector_u_t& p ) = 0;

                /**
                 * @brief
                 * @return vector_u_t&
                 */
                virtual vector_u_t& values() = 0;

                /**
                 * @brief
                 * @return const vector_u_t&
                 */
                virtual const vector_u_t& values() const = 0;

                /**
                 * @brief
                 * @return vector_u_t&
                 */
                virtual vector_u_t& axis() = 0;

                /**
                 * @brief
                 * @return const vector_u_t&
                 */
                virtual const vector_u_t& axis() const = 0;

                /**
                 * @brief
                 * @param altitude
                 * @return double
                 */
                virtual double get( double altitude ) = 0;

                /**
                 * @brief
                 * @param altitude
                 * @return double
                 */
                virtual double get_first_derivative( double altitude ) = 0;

                /**
                 * @brief
                 * @param altitude
                 * @return double
                 */
                virtual double get_second_derivative( double altitude ) = 0;

                /**
                 * @brief
                 * @return const units_ptr_t
                 */
                virtual const units_ptr_t get_units() const = 0;

                /**
                 * @brief
                 * @return const units_ptr_t
                 */
                virtual const units_ptr_t get_axis_units() const = 0;

                /**
                 * @brief
                 * @param u
                 * @return abstract_atmospheric_property_1d&
                 */
                virtual abstract_atmospheric_property_1d& convert_units(
                    const NCPA::units::Unit& u ) = 0;

                /**
                 * @brief
                 * @param u
                 * @return abstract_atmospheric_property_1d&
                 */
                virtual abstract_atmospheric_property_1d& convert_axis_units(
                    const NCPA::units::Unit& u ) = 0;

                /**
                 * @brief
                 * @param new_z
                 * @return abstract_atmospheric_property_1d&
                 */
                virtual abstract_atmospheric_property_1d& resample(
                    vector_u_t new_z ) = 0;

                /**
                 * @brief
                 * @return size_t
                 */
                virtual size_t dimensions() const override { return 1; }
        };
    }  // namespace atmos
}  // namespace NCPA

/**
 * @brief
 * @param a
 * @param b
 */
static void swap( NCPA::atmos::abstract_atmospheric_property_1d& a,
                  NCPA::atmos::abstract_atmospheric_property_1d& b ) noexcept {
}
