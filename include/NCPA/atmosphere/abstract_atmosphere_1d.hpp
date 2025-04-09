#pragma once

#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/calculations.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <string>
#include <vector>

// namespace NCPA {
//     namespace atmos {
//             class abstract_atmosphere_1d;
//     }  // namespace atmos
// }  // namespace NCPA

#define RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D \
    return static_cast<abstract_atmosphere_1d&>( *this );

static void swap( NCPA::atmos::abstract_atmosphere_1d&,
                  NCPA::atmos::abstract_atmosphere_1d& ) noexcept;

namespace NCPA {
    namespace atmos {

        /*
        TO DO:
        Remove set_altitude_vector.
        */
        class abstract_atmosphere_1d {
            public:
                virtual ~abstract_atmosphere_1d() {}

                friend void ::swap( abstract_atmosphere_1d& a,
                                    abstract_atmosphere_1d& b ) noexcept;

                virtual size_t size() const = 0;
                virtual abstract_atmosphere_1d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;
                virtual abstract_atmosphere_1d& set_axis( vector_u_t z ) = 0;
                virtual abstract_atmosphere_1d& add_property(
                    const std::string& key,
                    const AtmosphericProperty1D& property )
                    = 0;
                // virtual abstract_atmosphere_1d& add_property(
                //     const std::string& key, const vector_u_t& property )
                //     = 0;
                virtual abstract_atmosphere_1d& add_property(
                    const std::string& key, const scalar_u_t& property )
                    = 0;
                virtual abstract_atmosphere_1d& remove_property(
                    const std::string& key )
                    = 0;
                virtual abstract_atmosphere_1d& copy_property(
                    const std::string& old_key, const std::string& new_key )
                    = 0;
                virtual std::unique_ptr<abstract_atmosphere_1d> clone() const
                    = 0;
                virtual AtmosphericProperty1D& get_property(
                    const std::string& key )
                    = 0;
                virtual const AtmosphericProperty1D& get_property(
                    const std::string& key ) const
                    = 0;
                virtual vector_u_t get_axis_vector() const = 0;
                virtual vector_u_t get_vector(const std::string& key)            = 0;
                virtual double get( const std::string& key ) const
                    = 0;  // scalars
                virtual double get( const std::string& key, double altitude )
                    = 0;
                virtual double get_first_derivative( const std::string& key,
                                                     double altitude )
                    = 0;
                virtual double get_second_derivative( const std::string& key,
                                                      double altitude )
                    = 0;
                virtual units_ptr_t get_axis_units() const = 0;
                virtual units_ptr_t get_units( const std::string& key ) const
                    = 0;
                virtual double get_minimum_axis() const = 0;
                virtual double get_maximum_axis() const = 0;
                virtual abstract_atmosphere_1d& convert_axis_units(
                    units_ptr_t new_units )
                    = 0;
                virtual abstract_atmosphere_1d& convert_units(
                    const std::string& key, units_ptr_t new_units )
                    = 0;
                virtual abstract_atmosphere_1d& resample( double new_dz ) = 0;
                virtual abstract_atmosphere_1d& resample( vector_u_t new_z )
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

                virtual std::vector<AtmosphericProperty1D> properties() const {
                    std::vector<AtmosphericProperty1D> props;
                    auto keys = this->get_vector_keys();
                    for (auto keyit = keys.cbegin(); keyit != keys.cend(); ++keyit) {
                        props.emplace_back( this->get_property( *keyit ) );
                    }
                    return props;
                }

                virtual std::vector<scalar_u_t> scalars() const {
                    std::vector<scalar_u_t> sc;
                    auto keys = this->get_scalar_keys();
                    for (auto keyit = keys.cbegin(); keyit != keys.cend(); ++keyit) {
                        sc.emplace_back( this->get( *keyit ), this->get_units( *keyit ) );
                    }
                    return sc;
                }
        };

        typedef std::unique_ptr<abstract_atmosphere_1d> _atm_1d_ptr_t;
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract_atmosphere_1d& a,
                  NCPA::atmos::abstract_atmosphere_1d& b ) noexcept {}
