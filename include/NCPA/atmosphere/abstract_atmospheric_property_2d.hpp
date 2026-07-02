#pragma once

/**
 * @file
 * @brief
 */

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/exceptions.hpp"
#include "NCPA/interpolation.hpp"

#include <array>
#include <stdexcept>
#include <unordered_map>

/**
 * @brief
 * @param a
 * @param b
 */
static void swap( NCPA::atmos::abstract_atmospheric_property_2d& a,
                  NCPA::atmos::abstract_atmospheric_property_2d& b ) noexcept;

namespace NCPA {
    namespace atmos {

        /**
         * @brief
         */
        class abstract_atmospheric_property_2d
            : public abstract_atmospheric_property {
            public:
                /**
                 * @brief
                 */
                virtual ~abstract_atmospheric_property_2d() {}

                /**
                 * @brief
                 * @param a
                 * @param b
                 */
                friend void ::swap(
                    abstract_atmospheric_property_2d& a,
                    abstract_atmospheric_property_2d& b ) noexcept;

                /**
                 * @brief
                 * @return size_t
                 */
                virtual size_t dimensions() const override { return 2; }

                /**
                 * @brief
                 * @return std::unique_ptr<abstract_atmospheric_property_2d>
                 */
                virtual std::unique_ptr<abstract_atmospheric_property_2d>
                    clone2d() const = 0;

                /**
                 * @brief
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& clear() = 0;

                /**
                 * @brief
                 * @param ranges
                 * @param atmos1ds
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ranges,
                    const std::vector<
                        const abstract_atmospheric_property_1d *>& atmos1ds )
                    = 0;

                /**
                 * @brief
                 * @param range
                 * @param atmos1d
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d ) = 0;

                /**
                 * @brief
                 * @param atmos1d
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set(
                    const abstract_atmospheric_property_1d& atmos1d ) = 0;

                /**
                 * @brief
                 * @param ax1
                 * @param ax2
                 * @param units
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const units_ptr_t units ) = 0;

                /**
                 * @brief
                 * @param range
                 * @param atmos1d
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& append(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d ) = 0;

                /**
                 * @brief
                 * @param range
                 * @return _atm_prop_1d_ptr_t
                 */
                virtual _atm_prop_1d_ptr_t extract( double range ) = 0;

                /**
                 * @brief
                 * @param d
                 * @return size_t
                 */
                virtual size_t size( size_t d ) const = 0;

                /**
                 * @brief
                 * @param interp_type
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    = 0;

                /**
                 * @brief
                 * @param interp_type
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;

                /**
                 * @brief
                 * @param dim
                 * @param ax
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set(
                    size_t dim, const vector_u_t& ax ) = 0;

                /**
                 * @brief
                 * @param ax
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set(
                    const vector2d_u_t& ax ) = 0;

                /**
                 * @brief
                 * @param ax1
                 * @param ax2
                 * @param vals
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector2d_u_t& vals ) = 0;

                /**
                 * @brief
                 * @param n
                 * @return vector_u_t&
                 */
                virtual vector_u_t& axis( size_t n ) = 0;

                /**
                 * @brief
                 * @param n
                 * @return const vector_u_t&
                 */
                virtual const vector_u_t& axis( size_t n ) const = 0;

                /**
                 * @brief
                 * @param dim
                 * @param minax
                 * @param maxax
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& set_limits(
                    size_t dim, double minax, double maxax ) = 0;

                /**
                 * @brief
                 * @param dim
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& reset_limits(
                    size_t dim ) = 0;

                /**
                 * @brief
                 * @param dim
                 * @return const std::pair<double, double>
                 */
                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const = 0;


                /**
                 * @brief
                 * @return vector2d_u_t&
                 */
                virtual vector2d_u_t& values() = 0;

                /**
                 * @brief
                 * @return const vector2d_u_t&
                 */
                virtual const vector2d_u_t& values() const = 0;

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @return double
                 */
                virtual double get( double val1, double val2 ) = 0;

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @return double
                 */
                virtual double f( double val1, double val2 ) {
                    return this->get( val1, val2 );
                }

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param rel
                 * @return double
                 */
                virtual double get_first_derivative( double val1, double val2,
                                                     size_t rel ) = 0;

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param rel
                 * @return double
                 */
                virtual double df( double val1, double val2, size_t rel ) {
                    return this->get_first_derivative( val1, val2, rel );
                }

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param rel1
                 * @param rel2
                 * @return double
                 */
                virtual double get_second_derivative( double val1, double val2,
                                                      size_t rel1,
                                                      size_t rel2 ) = 0;

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param rel1
                 * @param rel2
                 * @return double
                 */
                virtual double ddf( double val1, double val2, size_t rel1,
                                    size_t rel2 ) {
                    return this->get_second_derivative( val1, val2, rel1,
                                                        rel2 );
                }

                /**
                 * @brief
                 * @return const units_ptr_t
                 */
                virtual const units_ptr_t get_units() const = 0;

                /**
                 * @brief
                 * @param n
                 * @return const units_ptr_t
                 */
                virtual const units_ptr_t get_axis_units( size_t n ) const = 0;

                /**
                 * @brief
                 * @param u
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& convert_units(
                    const NCPA::units::Unit& u ) = 0;

                /**
                 * @brief
                 * @param n
                 * @param u
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) = 0;

                /**
                 * @brief
                 * @param new_r
                 * @param new_z
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& resample(
                    const vector_u_t& new_r, const vector_u_t& new_z ) = 0;

                /**
                 * @brief
                 * @param dim
                 * @param new_ax
                 * @return abstract_atmospheric_property_2d&
                 */
                virtual abstract_atmospheric_property_2d& resample(
                    size_t dim, const vector_u_t& new_ax ) {
                    if (dim == 0) {
                        return this->resample( new_ax, this->axis( 1 ) );
                    } else {
                        return this->resample( this->axis( 0 ), new_ax );
                    }
                }
        };


    }  // namespace atmos
}  // namespace NCPA

/**
 * @brief
 * @param a
 * @param b
 */
static void swap( NCPA::atmos::abstract_atmospheric_property_2d& a,
                  NCPA::atmos::abstract_atmospheric_property_2d& b ) noexcept {
}
