#pragma once

/**
 * @file
 * @brief 
 */

#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/calculations.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace NCPA {
    namespace atmos {
        class abstract_atmosphere_2d;
    }  // namespace atmos
}  // namespace NCPA

/**
 * @brief 
 */
#define RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D \
    return static_cast<abstract_atmosphere_2d&>( *this );

/**
 * @brief 
 * @param a 
 * @param b 
 */
static void swap( NCPA::atmos::abstract_atmosphere_2d&,
                  NCPA::atmos::abstract_atmosphere_2d& ) noexcept;

namespace NCPA {
    namespace atmos {

        /**
         * @brief 
         */
        class abstract_atmosphere_2d {
            public:
                /**
                 * @brief 
                 */
                virtual ~abstract_atmosphere_2d() {}

                /**
                 * @brief 
                 * @param a 
                 * @param b 
                 */
                friend void ::swap( abstract_atmosphere_2d& a,
                                    abstract_atmosphere_2d& b ) noexcept;

                /**
                 * @brief 
                 * @param atmos1d 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& set(
                    abstract_atmosphere_1d& atmos1d )
                    = 0;

                /**
                 * @brief 
                 * @param ranges 
                 * @param components 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& set(
                    const vector_u_t ranges,
                    std::vector<abstract_atmosphere_1d *> components )
                    = 0;

                /**
                 * @brief 
                 * @param range 
                 * @param atmos1d 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& append(
                    scalar_u_t range, abstract_atmosphere_1d& atmos1d )
                    = 0;

                /**
                 * @brief 
                 * @param dim 
                 * @return size_t 
                 */
                virtual size_t size( size_t dim ) const = 0;

                /**
                 * @brief 
                 * @param interp_type 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    = 0;

                /**
                 * @brief 
                 * @param interp_type 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;

                /**
                 * @brief 
                 * @param axis 
                 * @param vals 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& set_axis( size_t axis,
                                                          vector_u_t vals )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param property 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const AtmosphericProperty2D& property )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param property 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key, const vector2d_u_t& property )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param property 
                 * @param index 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key, const vector_u_t& property,
                    const vector_u_t& index )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param property 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key, const vector_u_t& property )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param property 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key, const scalar_u_t& property )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& remove_property(
                    const std::string& key )
                    = 0;

                /**
                 * @brief 
                 * @param old_key 
                 * @param new_key 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& copy_property(
                    const std::string& old_key, const std::string& new_key )
                    = 0;

                /**
                 * @brief 
                 * @return std::unique_ptr<abstract_atmosphere_2d> 
                 */
                virtual std::unique_ptr<abstract_atmosphere_2d> clone() const
                    = 0;

                /**
                 * @brief 
                 * @param range 
                 * @return std::unique_ptr<abstract_atmosphere_1d> 
                 */
                virtual std::unique_ptr<abstract_atmosphere_1d> extract(
                    double range )
                    = 0;

                /**
                 * @brief 
                 * @param n 
                 * @return vector_u_t 
                 */
                virtual vector_u_t get_axis_vector( size_t n )       = 0;

                /**
                 * @brief 
                 * @param n 
                 * @return vector_u_t 
                 */
                virtual vector_u_t get_axis_vector( size_t n ) const = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param r 
                 * @return double 
                 */
                virtual double get( const std::string& key, double r )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param r 
                 * @param z 
                 * @return double 
                 */
                virtual double get( const std::string& key, double r,
                                    double z )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param r 
                 * @param z 
                 * @param wrt 
                 * @return double 
                 */
                virtual double get_first_derivative( const std::string& key,
                                                     double r, double z,
                                                     size_t wrt )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param r 
                 * @param z 
                 * @param wrt1 
                 * @param wrt2 
                 * @return double 
                 */
                virtual double get_second_derivative( const std::string& key,
                                                      double r, double z,
                                                      size_t wrt1,
                                                      size_t wrt2 )
                    = 0;

                /**
                 * @brief 
                 * @return true 
                 * @return false 
                 */
                virtual bool is_stratified() const = 0;

                /**
                 * @brief 
                 * @param r1 
                 * @param r2 
                 * @return true 
                 * @return false 
                 */
                virtual bool same( scalar_u_t r1, scalar_u_t r2 ) const = 0;

                /**
                 * @brief 
                 * @param r1 
                 * @param r2 
                 * @return true 
                 * @return false 
                 */
                virtual bool same( double r1, double r2 ) const         = 0;

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
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& convert_axis_units(
                    size_t n, units_ptr_t new_units )
                    = 0;

                /**
                 * @brief 
                 * @param key 
                 * @param new_units 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& convert_units(
                    const std::string& key, units_ptr_t new_units )
                    = 0;

                /**
                 * @brief 
                 * @param axis 
                 * @param new_d 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& resample( size_t axis,
                                                          double new_d )
                    = 0;

                /**
                 * @brief 
                 * @param axis 
                 * @param new_z 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& resample( size_t axis,
                                                          vector_u_t new_z )
                    = 0;

                /**
                 * @brief 
                 * @param new_ax1 
                 * @param new_ax2 
                 * @return abstract_atmosphere_2d& 
                 */
                virtual abstract_atmosphere_2d& resample( vector_u_t new_ax1,
                                                          vector_u_t new_ax2 )
                    = 0;

                /**
                 * @brief 
                 * @return std::vector<std::string> 
                 */
                virtual std::vector<std::string> get_keys() const        = 0;

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
                 * @param os 
                 */
                virtual void print( std::ostream& os ) = 0;

                /**
                 * @brief 
                 * @param n 
                 * @throws std::range_error 
                 */
                virtual void validate_axis( size_t n ) const {
                    if (n >= 2) {
                        throw std::range_error( "Only axis 0 and 1 are valid "
                                                "for 2-D atmosphere!" );
                    }
                }
        };

        /**
         * @brief 
         */
        typedef std::unique_ptr<abstract_atmosphere_2d> _atm_2d_ptr_t;
    }  // namespace atmos
}  // namespace NCPA

/**
 * @brief 
 * @param a 
 * @param b 
 */
static void swap( NCPA::atmos::abstract_atmosphere_2d& a,
                  NCPA::atmos::abstract_atmosphere_2d& b ) noexcept {}