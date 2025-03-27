#pragma once

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

#define RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D \
    return static_cast<abstract_atmosphere_2d&>( *this );

static void swap( NCPA::atmos::abstract_atmosphere_2d&,
                  NCPA::atmos::abstract_atmosphere_2d& ) noexcept;

namespace NCPA {
    namespace atmos {

        class abstract_atmosphere_2d {
            public:
                virtual ~abstract_atmosphere_2d() {}

                friend void ::swap( abstract_atmosphere_2d& a,
                                    abstract_atmosphere_2d& b ) noexcept;

                virtual abstract_atmosphere_2d& set(
                    abstract_atmosphere_1d& atmos1d )
                    = 0;

                virtual abstract_atmosphere_2d& set(
                    const vector_u_t ranges,
                    std::vector<abstract_atmosphere_1d *> components )
                    = 0;

                virtual abstract_atmosphere_2d& append(
                    scalar_u_t range, abstract_atmosphere_1d& atmos1d )
                    = 0;

                virtual size_t size( size_t dim ) const = 0;

                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    = 0;
                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;

                virtual abstract_atmosphere_2d& set_axis( size_t axis,
                                                          vector_u_t vals )
                    = 0;

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const AtmosphericProperty2D& property )
                    = 0;
                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key, const vector2d_u_t& property )
                    = 0;
                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key, const vector_u_t& property,
                    const vector_u_t& index )
                    = 0;
                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key, const vector_u_t& property )
                    = 0;
                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key, const scalar_u_t& property )
                    = 0;
                virtual abstract_atmosphere_2d& remove_property(
                    const std::string& key )
                    = 0;
                virtual abstract_atmosphere_2d& copy_property(
                    const std::string& old_key, const std::string& new_key )
                    = 0;

                virtual std::unique_ptr<abstract_atmosphere_2d> clone() const
                    = 0;

                // virtual AtmosphericProperty2D& get_property(
                //     const std::string& key )
                //     = 0;
                // virtual const AtmosphericProperty2D& get_property(
                //     const std::string& key ) const
                //     = 0;

                virtual std::unique_ptr<abstract_atmosphere_1d> extract(
                    double range )
                    = 0;

                virtual vector_u_t get_axis_vector( size_t n )       = 0;
                virtual vector_u_t get_axis_vector( size_t n ) const = 0;

                virtual double get( const std::string& key, double r )
                    = 0;  // scalars
                virtual double get( const std::string& key, double r,
                                    double z )
                    = 0;
                virtual double get_first_derivative( const std::string& key,
                                                     double r, double z,
                                                     size_t wrt )
                    = 0;
                virtual double get_second_derivative( const std::string& key,
                                                      double r, double z,
                                                      size_t wrt1,
                                                      size_t wrt2 )
                    = 0;
                virtual bool same( scalar_u_t r1, scalar_u_t r2 ) const = 0;
                virtual bool same( double r1, double r2 ) const         = 0;

                virtual units_ptr_t get_axis_units( size_t n ) const = 0;
                virtual units_ptr_t get_units( const std::string& key ) const
                    = 0;

                virtual double get_minimum_axis( size_t n ) const = 0;
                virtual double get_maximum_axis( size_t n ) const = 0;

                virtual abstract_atmosphere_2d& convert_axis_units(
                    size_t n, units_ptr_t new_units )
                    = 0;
                virtual abstract_atmosphere_2d& convert_units(
                    const std::string& key, units_ptr_t new_units )
                    = 0;
                virtual abstract_atmosphere_2d& resample( size_t axis,
                                                          double new_d )
                    = 0;
                virtual abstract_atmosphere_2d& resample( size_t axis,
                                                          vector_u_t new_z )
                    = 0;
                virtual abstract_atmosphere_2d& resample( vector_u_t new_ax1,
                                                          vector_u_t new_ax2 )
                    = 0;

                virtual std::vector<std::string> get_keys() const        = 0;
                virtual std::vector<std::string> get_vector_keys() const = 0;
                virtual std::vector<std::string> get_scalar_keys() const = 0;
                virtual bool contains_scalar( const std::string& key ) const
                    = 0;
                virtual bool contains_vector( const std::string& key ) const
                    = 0;
                virtual bool contains_key( const std::string& key ) const = 0;

                virtual void print( std::ostream& os ) = 0;

                virtual void validate_axis( size_t n ) const {
                    if (n >= 2) {
                        throw std::range_error( "Only axis 0 and 1 are valid "
                                                "for 2-D atmosphere!" );
                    }
                }
        };

        typedef std::unique_ptr<abstract_atmosphere_2d> _atm_2d_ptr_t;
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract_atmosphere_2d& a,
                  NCPA::atmos::abstract_atmosphere_2d& b ) noexcept {}
