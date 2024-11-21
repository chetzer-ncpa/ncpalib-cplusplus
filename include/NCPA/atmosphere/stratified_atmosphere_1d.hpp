#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_1d.hpp"
#include "NCPA/atmosphere/atmospheric_property_1d.hpp"
#include "NCPA/atmosphere/calculations.hpp"
#include "NCPA/atmosphere/types.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/units.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace NCPA {
    namespace atmos {
        class stratified_atmosphere_1d
            : public details::abstract_atmosphere_1d {
            public:
                virtual ~stratified_atmosphere_1d() {}

                virtual size_t nz() const override { return _z.size(); }

                virtual abstract_atmosphere_1d& set_altitude_vector(
                    vector_t z )
                    = 0;
                virtual abstract_atmosphere_1d& add_property(
                    const std::string& key, vector_t property )
                    = 0;
                virtual abstract_atmosphere_1d& add_property(
                    const std::string& key, scalar_t property )
                    = 0;
                virtual void remove_property( const std::string& key ) = 0;
                virtual abstract_atmosphere_1d& copy_property(
                    const std::string& old_key, const std::string& new_key )
                    = 0;

                virtual scalar_t get( const std::string& key ) const
                    = 0;  // scalars
                virtual scalar_t get( const std::string& key,
                                      scalar_t altitude ) const
                    = 0;
                virtual scalar_t get_first_derivative(
                    const std::string& key, scalar_t altitude ) const
                    = 0;
                virtual scalar_t get_second_derivative(
                    const std::string& key, scalar_t altitude ) const
                    = 0;
                virtual units_ptr_t get_altitude_units() const = 0;
                virtual units_ptr_t get_property_units(
                    const std::string& key ) const
                    = 0;

                virtual scalar_t get_minimum_altitude() const = 0;
                virtual scalar_t get_maximum_altitude() const = 0;

                virtual abstract_atmosphere_1d& convert_altitude_units(
                    units_ptr_t new_units )
                    = 0;
                virtual abstract_atmosphere_1d& convert_property_units(
                    const std::string& key, units_ptr_t new_units )
                    = 0;
                virtual abstract_atmosphere_1d& resample( double new_dz ) = 0;

                virtual std::vector<std::string> get_keys() const        = 0;
                virtual std::vector<std::string> get_vector_keys() const = 0;
                virtual std::vector<std::string> get_scalar_keys() const = 0;
                virtual bool contains_scalar( const std::string& key ) const
                    = 0;
                virtual bool contains_vector( const std::string& key ) const
                    = 0;
                virtual bool contains_key( const std::string& key ) const = 0;

            private:
                std::unordered_map<std::string, AtmosphericProperty1D> _properties;
                std::unordered_map<std::string, scalar_t>
                    _scalar_properties;

                // double *z_;
                // size_t nz_;
                // std::stack< units_t > z_units_;
                vector_t _z;
                // std::vector<std::string> _headerlines;

                // void cleanup();
                // void do_units_conversion_( size_t n_points, double *inplace,
                //                            NCPA::units_t fromUnits,
                //                            NCPA::units_t toUnits );
        };
    }  // namespace atmos
}  // namespace NCPA
