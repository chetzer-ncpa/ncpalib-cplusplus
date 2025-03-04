#pragma once

#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/calculations.hpp"
#include "NCPA/atmosphere/declarations.hpp"
// #include "NCPA/atmosphere/types.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <string>
#include <vector>


static void swap( NCPA::atmos::abstract::atmosphere_1d&,
                  NCPA::atmos::abstract::atmosphere_1d& ) noexcept;

namespace NCPA {
    namespace atmos {
        namespace abstract {
            class atmosphere_1d {
                public:
                    virtual ~atmosphere_1d() {}

                    friend void ::swap( atmosphere_1d& a,
                                        atmosphere_1d& b ) noexcept;

                    virtual size_t nz() const = 0;

                    virtual atmosphere_1d& set_interpolator(
                        NCPA::interpolation::interpolator_1d_type_t interp_type )
                        = 0;

                    virtual atmosphere_1d& set_altitude_vector(
                        vector_u_t z )
                        = 0;
                    virtual atmosphere_1d& add_property(
                        const std::string& key,
                        const AtmosphericProperty1D& property )
                        = 0;
                    virtual atmosphere_1d& add_property(
                        const std::string& key, const vector_u_t& property )
                        = 0;
                    virtual atmosphere_1d& add_property(
                        const std::string& key, const scalar_u_t& property )
                        = 0;
                    virtual atmosphere_1d& remove_property(
                        const std::string& key )
                        = 0;
                    virtual atmosphere_1d& copy_property(
                        const std::string& old_key,
                        const std::string& new_key )
                        = 0;

                    virtual std::unique_ptr<atmosphere_1d> clone()
                        const
                        = 0;

                    virtual AtmosphericProperty1D& get_property(
                        const std::string& key )
                        = 0;
                    virtual const AtmosphericProperty1D& get_property(
                        const std::string& key ) const = 0;
                    
                    virtual vector_u_t& get_altitude_vector() = 0;
                    virtual const vector_u_t& get_altitude_vector() const = 0;

                    virtual double get( const std::string& key ) const
                        = 0;  // scalars
                    virtual double get( const std::string& key,
                                        double altitude )
                        = 0;
                    virtual double get_first_derivative(
                        const std::string& key, double altitude )
                        = 0;
                    virtual double get_second_derivative(
                        const std::string& key, double altitude )
                        = 0;

                    virtual units_ptr_t get_altitude_units() const = 0;
                    virtual units_ptr_t get_property_units(
                        const std::string& key ) const
                        = 0;

                    virtual double get_minimum_altitude() const = 0;
                    virtual double get_maximum_altitude() const = 0;

                    virtual atmosphere_1d& convert_altitude_units(
                        units_ptr_t new_units )
                        = 0;
                    virtual atmosphere_1d& convert_property_units(
                        const std::string& key, units_ptr_t new_units )
                        = 0;
                    virtual atmosphere_1d& resample( double new_dz )
                        = 0;
                    virtual atmosphere_1d& resample( vector_u_t new_z )
                        = 0;

                    virtual std::vector<std::string> get_keys() const = 0;
                    virtual std::vector<std::string> get_vector_keys() const
                        = 0;
                    virtual std::vector<std::string> get_scalar_keys() const
                        = 0;
                    virtual bool contains_scalar(
                        const std::string& key ) const
                        = 0;
                    virtual bool contains_vector(
                        const std::string& key ) const
                        = 0;
                    virtual bool contains_key( const std::string& key ) const
                        = 0;

                    virtual void print( std::ostream& os ) = 0;
            };
        }
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract::atmosphere_1d& a,
                  NCPA::atmos::abstract::atmosphere_1d& b ) noexcept {}
