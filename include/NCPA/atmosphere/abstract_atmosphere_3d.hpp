#pragma once

/**
 * @file
 * @brief
 */

#include "NCPA/atmosphere/AtmosphericProperty3D.hpp"
#include "NCPA/atmosphere/calculations.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace NCPA {
    namespace atmos {
        class abstract_atmosphere_3d;
    }  // namespace atmos
}  // namespace NCPA

/**
 * @brief
 * @param a
 * @param b
 */
static void swap( NCPA::atmos::abstract_atmosphere_3d&,
                  NCPA::atmos::abstract_atmosphere_3d& ) noexcept;

namespace NCPA {
    namespace atmos {

        /**
         * @brief
         */
        class abstract_atmosphere_3d {
            public:
                /**
                 * @brief
                 */
                virtual ~abstract_atmosphere_3d() {}

                /**
                 * @brief
                 * @param a
                 * @param b
                 */
                friend void ::swap( abstract_atmosphere_3d& a,
                                    abstract_atmosphere_3d& b ) noexcept;

                /**
                 * @brief
                 * @param atmos1d
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& set(
                    abstract_atmosphere_1d& atmos1d ) = 0;

                /**
                 * @brief
                 * @param components
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& set(
                    NCPA::arrays::ndvector<2, abstract_atmosphere_1d *>&
                        components ) = 0;

                /**
                 * @brief
                 * @param ax1
                 * @param ax2
                 * @param components
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    NCPA::arrays::ndvector<2, abstract_atmosphere_1d *>&
                        components ) = 0;

                /**
                 * @brief
                 * @param dim
                 * @return size_t
                 */
                virtual size_t size( size_t dim ) const = 0;

                /**
                 * @brief
                 * @param interp_type
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type )
                    = 0;

                /**
                 * @brief
                 * @param interp_type
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    = 0;

                /**
                 * @brief
                 * @param interp_type
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;

                /**
                 * @brief
                 * @param axis
                 * @param vals
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& set_axis( size_t axis,
                                                          vector_u_t vals )
                    = 0;

                /**
                 * @brief
                 * @param key
                 * @param property
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& add_property(
                    const std::string& key,
                    const AtmosphericProperty3D& property ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param property
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& add_property(
                    const std::string& key, const vector3d_u_t& property ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param property
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& add_property(
                    const std::string& key, const vector2d_u_t& property ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param property
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& add_property(
                    const std::string& key, const scalar_u_t& property ) = 0;

                /**
                 * @brief
                 * @param key
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& remove_property(
                    const std::string& key ) = 0;

                /**
                 * @brief
                 * @param old_key
                 * @param new_key
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& copy_property(
                    const std::string& old_key, const std::string& new_key )
                    = 0;

                /**
                 * @brief
                 * @return std::unique_ptr<abstract_atmosphere_3d>
                 */
                virtual std::unique_ptr<abstract_atmosphere_3d> clone() const
                    = 0;

                /**
                 * @brief
                 * @param n
                 * @return vector_u_t
                 */
                virtual vector_u_t get_axis_vector( size_t n ) = 0;

                /**
                 * @brief
                 * @param n
                 * @return vector_u_t
                 */
                virtual vector_u_t get_axis_vector( size_t n ) const = 0;

                /**
                 * @brief
                 * @param key
                 * @return vector3d_u_t
                 */
                virtual vector3d_u_t get_values( const std::string& key ) = 0;

                /**
                 * @brief
                 * @param key
                 * @return vector3d_u_t
                 */
                virtual vector3d_u_t get_values( const std::string& key ) const
                    = 0;

                /**
                 * @brief
                 * @param key
                 * @param x1
                 * @param x2
                 * @return double
                 */
                virtual double get( const std::string& key, double x1,
                                    double x2 ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param x1
                 * @param x2
                 * @return vector2d_u_t
                 */
                virtual vector2d_u_t get( const std::string& key,
                                          const std::vector<double>& x1,
                                          const std::vector<double>& x2 ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param x1
                 * @param x2
                 * @param wrt
                 * @return double
                 */
                virtual double get_first_derivative( const std::string& key,
                                                     double x1, double x2,
                                                     size_t wrt ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param x1
                 * @param x2
                 * @param wrt1
                 * @param wrt2
                 * @return double
                 */
                virtual double get_second_derivative( const std::string& key,
                                                      double x1, double x2,
                                                      size_t wrt1,
                                                      size_t wrt2 ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param x1
                 * @param x2
                 * @param x3
                 * @return double
                 */
                virtual double get( const std::string& key, double x1,
                                    double x2, double x3 ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param v1
                 * @param v2
                 * @param v3
                 * @return vector3d_u_t
                 */
                virtual vector3d_u_t get( const std::string& key,
                                          const std::vector<double>& v1,
                                          const std::vector<double>& v2,
                                          const std::vector<double>& v3 ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param x1
                 * @param x2
                 * @param x3
                 * @param wrt
                 * @return double
                 */
                virtual double get_first_derivative( const std::string& key,
                                                     double x1, double x2,
                                                     double x3, size_t wrt )
                    = 0;

                /**
                 * @brief
                 * @param key
                 * @param x1
                 * @param x2
                 * @param x3
                 * @param wrt1
                 * @param wrt2
                 * @return double
                 */
                virtual double get_second_derivative( const std::string& key,
                                                      double x1, double x2,
                                                      double x3, size_t wrt1,
                                                      size_t wrt2 ) = 0;

                /**
                 * @brief
                 * @param n
                 * @return units_ptr_t
                 */
                virtual units_ptr_t get_axis_units( size_t n ) const = 0;

                /**
                 * @brief
                 * @param key
                 * @return units_ptr_t
                 */
                virtual units_ptr_t get_units( const std::string& key ) const
                    = 0;

                /**
                 * @brief
                 * @param n
                 * @return double
                 */
                virtual double get_minimum_axis( size_t n ) const = 0;

                /**
                 * @brief
                 * @param n
                 * @return double
                 */
                virtual double get_maximum_axis( size_t n ) const = 0;

                /**
                 * @brief
                 * @param n
                 * @param new_units
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& convert_axis_units(
                    size_t n, units_ptr_t new_units ) = 0;

                /**
                 * @brief
                 * @param key
                 * @param new_units
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& convert_units(
                    const std::string& key, units_ptr_t new_units ) = 0;

                /**
                 * @brief
                 * @param axis
                 * @param new_d
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& resample( size_t axis,
                                                          double new_d ) = 0;

                /**
                 * @brief
                 * @param axis
                 * @param new_z
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& resample( size_t axis,
                                                          vector_u_t new_z )
                    = 0;

                /**
                 * @brief
                 * @param new_ax1
                 * @param new_ax2
                 * @param new_ax3
                 * @return abstract_atmosphere_3d&
                 */
                virtual abstract_atmosphere_3d& resample( vector_u_t new_ax1,
                                                          vector_u_t new_ax2,
                                                          vector_u_t new_ax3 )
                    = 0;

                /**
                 * @brief
                 * @return std::vector<std::string>
                 */
                virtual std::vector<std::string> get_keys() const = 0;

                /**
                 * @brief
                 * @return std::vector<std::string>
                 */
                virtual std::vector<std::string> get_vector_keys() const = 0;

                /**
                 * @brief
                 * @return std::vector<std::string>
                 */
                virtual std::vector<std::string> get_scalar_keys() const = 0;

                /**
                 * @brief
                 * @param key
                 * @return true
                 * @return false
                 */
                virtual bool contains_scalar( const std::string& key ) const
                    = 0;

                /**
                 * @brief
                 * @param key
                 * @return true
                 * @return false
                 */
                virtual bool contains_vector( const std::string& key ) const
                    = 0;

                /**
                 * @brief
                 * @param key
                 * @return true
                 * @return false
                 */
                virtual bool contains_key( const std::string& key ) const = 0;

                /**
                 * @brief
                 * @param x1_1
                 * @param x2_1
                 * @param x1_2
                 * @param x2_2
                 * @return true
                 * @return false
                 */
                virtual bool same( scalar_u_t x1_1, scalar_u_t x2_1,
                                   scalar_u_t x1_2, scalar_u_t x2_2 ) const
                    = 0;

                /**
                 * @brief
                 * @param x1_1
                 * @param x2_1
                 * @param x1_2
                 * @param x2_2
                 * @return true
                 * @return false
                 */
                virtual bool same( double x1_1, double x2_1, double x1_2,
                                   double x2_2 ) const = 0;

                /**
                 * @brief
                 * @return true
                 * @return false
                 */
                virtual bool is_stratified() const = 0;


                /**
                 * @brief
                 * @param os
                 */
                virtual void print( std::ostream& os ) = 0;

                /**
                 * @brief
                 * @param n
                 * @throws std::range_error
                 */
                virtual void validate_axis( size_t n ) const {
                    if (n >= 3) {
                        throw std::range_error(
                            "Only axes 0, 1 and 2 are valid "
                            "for 3-D atmosphere!" );
                    }
                }
        };

        /**
         * @brief
         */
        typedef std::unique_ptr<abstract_atmosphere_3d> _atm_3d_ptr_t;
    }  // namespace atmos
}  // namespace NCPA

/**
 * @brief
 * @param a
 * @param b
 */
static void swap( NCPA::atmos::abstract_atmosphere_3d& a,
                  NCPA::atmos::abstract_atmosphere_3d& b ) noexcept {}
