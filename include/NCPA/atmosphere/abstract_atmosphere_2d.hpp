#pragma once

#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/calculations.hpp"
#include "NCPA/atmosphere/types.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <string>
#include <vector>

namespace NCPA {
    namespace atmos {
            class abstract_atmosphere_2d;
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract_atmosphere_2d&,
                  NCPA::atmos::abstract_atmosphere_2d& ) noexcept;

namespace NCPA {
    namespace atmos {
        
            class abstract_atmosphere_2d {
                public:
                    virtual ~abstract_atmosphere_2d() {}

                    friend void ::swap( abstract_atmosphere_2d& a,
                                        abstract_atmosphere_2d& b ) noexcept;

                    virtual size_t nz() const = 0;

                    virtual abstract_atmosphere_2d& set_interpolator(
                        NCPA::interpolation::interpolator_2d_type_t interp_type )
                        = 0;

                    virtual abstract_atmosphere_2d& set_altitude_vector(
                        vector_t z )
                        = 0;
                    virtual abstract_atmosphere_2d& add_property(
                        const std::string& key,
                        const AtmosphericProperty1D& property )
                        = 0;
                    virtual abstract_atmosphere_2d& add_property(
                        const std::string& key, const vector_t& property )
                        = 0;
                    virtual abstract_atmosphere_2d& add_property(
                        const std::string& key, const scalar_t& property )
                        = 0;
                    virtual abstract_atmosphere_2d& remove_property(
                        const std::string& key )
                        = 0;
                    virtual abstract_atmosphere_2d& copy_property(
                        const std::string& old_key,
                        const std::string& new_key )
                        = 0;

                    virtual std::unique_ptr<abstract_atmosphere_2d> clone()
                        const
                        = 0;

                    virtual AtmosphericProperty1D& get_property(
                        const std::string& key )
                        = 0;
                    virtual const AtmosphericProperty1D& get_property(
                        const std::string& key ) const = 0;
                    
                    virtual vector_t& get_altitude_vector() = 0;
                    virtual const vector_t& get_altitude_vector() const = 0;

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

                    virtual abstract_atmosphere_2d& convert_altitude_units(
                        units_ptr_t new_units )
                        = 0;
                    virtual abstract_atmosphere_2d& convert_property_units(
                        const std::string& key, units_ptr_t new_units )
                        = 0;
                    virtual abstract_atmosphere_2d& resample( double new_dz )
                        = 0;
                    virtual abstract_atmosphere_2d& resample( vector_t new_z )
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
        
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract_atmosphere_2d& a,
                  NCPA::atmos::abstract_atmosphere_2d& b ) noexcept {}
