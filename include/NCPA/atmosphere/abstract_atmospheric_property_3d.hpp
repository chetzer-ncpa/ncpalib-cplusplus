#pragma once

/**
 * @file
 * @brief
 */

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
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
static void swap( NCPA::atmos::abstract_atmospheric_property_3d& a,
                  NCPA::atmos::abstract_atmospheric_property_3d& b ) noexcept;

namespace NCPA {
    namespace atmos {

        /**
         * @brief
         */
        class abstract_atmospheric_property_3d
            : public abstract_atmospheric_property {
            public:
                /**
                 * @brief
                 */
                virtual ~abstract_atmospheric_property_3d() {}

                /**
                 * @brief
                 * @param a
                 * @param b
                 */
                friend void ::swap(
                    abstract_atmospheric_property_3d& a,
                    abstract_atmospheric_property_3d& b ) noexcept;

                /**
                 * @brief
                 * @param dim
                 * @param dimval
                 * @param newslice
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& append(
                    size_t dim, double dimval, const vector2d_u_t& newslice )
                    = 0;

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
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& clear() = 0;

                /**
                 * @brief
                 * @return std::unique_ptr<abstract_atmospheric_property_3d>
                 */
                virtual std::unique_ptr<abstract_atmospheric_property_3d>
                    clone3d() const = 0;

                /**
                 * @brief
                 * @param n
                 * @param u
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) = 0;

                /**
                 * @brief
                 * @param u
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& convert_units(
                    const NCPA::units::Unit& u ) = 0;

                /**
                 * @brief
                 * @param n
                 * @return size_t
                 */
                virtual size_t dim( size_t n ) const = 0;

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param val3
                 * @return double
                 */
                virtual double get( double val1, double val2, double val3 )
                    = 0;

                /**
                 * @brief
                 * @param n
                 * @return const units_ptr_t
                 */
                virtual const units_ptr_t get_axis_units( size_t n ) const = 0;

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param val3
                 * @param rel
                 * @return double
                 */
                virtual double get_first_derivative( double val1, double val2,
                                                     double val3, size_t rel )
                    = 0;

                /**
                 * @brief
                 * @param dim
                 * @return const std::pair<double, double>
                 */
                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const = 0;

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param val3
                 * @param rel1
                 * @param rel2
                 * @return double
                 */
                virtual double get_second_derivative( double val1, double val2,
                                                      double val3, size_t rel1,
                                                      size_t rel2 ) = 0;

                /**
                 * @brief
                 * @return const units_ptr_t
                 */
                virtual const units_ptr_t get_units() const = 0;

                /**
                 * @brief
                 * @param ax1
                 * @param ax2
                 * @param ax3
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& resample(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector_u_t& ax3 ) = 0;

                /**
                 * @brief
                 * @param dim
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& reset_limits(
                    size_t dim ) = 0;

                /**
                 * @brief
                 * @param ax
                 * @param axvals
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& set(
                    size_t ax, const vector_u_t& axvals ) = 0;

                /**
                 * @brief
                 * @param ax
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& set(
                    const vector3d_u_t& ax ) = 0;

                /**
                 * @brief
                 * @param ax1
                 * @param ax2
                 * @param ax3
                 * @param vals
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector_u_t& ax3, const vector3d_u_t& vals ) = 0;

                /**
                 * @brief
                 * @param nx1
                 * @param nx2
                 * @param nx3
                 * @param units
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& set(
                    size_t nx1, size_t nx2, size_t nx3,
                    const units_ptr_t units ) = 0;

                /**
                 * @brief
                 * @param interp_type
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type )
                    = 0;

                /**
                 * @brief
                 * @param interp_type
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    = 0;

                /**
                 * @brief
                 * @param interp_type
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;

                /**
                 * @brief
                 * @param dim
                 * @param minax
                 * @param maxax
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& set_limits(
                    size_t dim, double minax, double maxax ) = 0;

                /**
                 * @brief
                 * @return std::vector<size_t>
                 */
                virtual std::vector<size_t> shape() const = 0;

                /**
                 * @brief
                 * @return vector3d_u_t&
                 */
                virtual vector3d_u_t& values() = 0;

                /**
                 * @brief
                 * @return const vector3d_u_t&
                 */
                virtual const vector3d_u_t& values() const = 0;

                /**
                 * @brief
                 * @param v1
                 * @param v2
                 * @param v3
                 * @return vector3d_u_t
                 */
                virtual vector3d_u_t get( const std::vector<double>& v1,
                                          const std::vector<double>& v2,
                                          const std::vector<double>& v3 ) {
                    vector3d_u_t v_out( v1.size(), v2.size(), v3.size(),
                                        this->get_units() );
                    for (auto i = 0; i < v1.size(); ++i) {
                        for (auto j = 0; j < v2.size(); ++j) {
                            for (auto k = 0; k < v3.size(); ++k) {
                                v_out[ i ][ j ][ k ]
                                    = this->get( v1[ i ], v2[ j ], v3[ k ] );
                            }
                        }
                    }
                    return v_out;
                }

                /**
                 * @brief
                 * @param v1
                 * @param v2
                 * @param v3
                 * @return vector3d_u_t
                 */
                virtual vector3d_u_t get( const vector_u_t& v1,
                                          const vector_u_t& v2,
                                          const vector_u_t& v3 ) {
                    return this->get( v1.as( this->get_axis_units( 0 ) ),
                                      v2.as( this->get_axis_units( 1 ) ),
                                      v3.as( this->get_axis_units( 2 ) ) );
                }

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param val3
                 * @param rel1
                 * @param rel2
                 * @return double
                 */
                virtual double ddf( double val1, double val2, double val3,
                                    size_t rel1, size_t rel2 ) {
                    return this->get_second_derivative( val1, val2, val3, rel1,
                                                        rel2 );
                }

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param val3
                 * @param rel
                 * @return double
                 */
                virtual double df( double val1, double val2, double val3,
                                   size_t rel ) {
                    return this->get_first_derivative( val1, val2, val3, rel );
                }

                /**
                 * @brief
                 * @return size_t
                 */
                virtual size_t dimensions() const override { return 3; }

                /**
                 * @brief
                 * @param val1
                 * @param val2
                 * @param val3
                 * @return double
                 */
                virtual double f( double val1, double val2, double val3 ) {
                    return this->get( val1, val2, val3 );
                }

                /**
                 * @brief
                 * @param dim
                 * @param new_ax
                 * @return abstract_atmospheric_property_3d&
                 */
                virtual abstract_atmospheric_property_3d& resample(
                    size_t dim, const vector_u_t& new_ax ) {
                    switch (dim) {
                        case 0:
                            return this->resample( new_ax, this->axis( 1 ),
                                                   this->axis( 2 ) );
                            break;
                        case 1:
                            return this->resample( this->axis( 0 ), new_ax,
                                                   this->axis( 2 ) );
                            break;
                        case 2:
                            return this->resample( this->axis( 0 ),
                                                   this->axis( 1 ), new_ax );
                            break;
                        default:
                            std::ostringstream oss;
                            oss << "abstract_atmopheric_property_3d.resample()"
                                   ": Invalid axis "
                                << dim;
                            throw std::range_error( oss.str() );
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
static void swap( NCPA::atmos::abstract_atmospheric_property_3d& a,
                  NCPA::atmos::abstract_atmospheric_property_3d& b ) noexcept {
}
