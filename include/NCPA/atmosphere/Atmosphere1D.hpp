#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_1d.hpp"
#include "NCPA/atmosphere/AtmosphericModel.hpp"
#include "NCPA/atmosphere/types.hpp"
#include "NCPA/units.hpp"
#include "NCPA/interpolation.hpp"

#include <memory>
#include <string>

namespace NCPA {
    namespace atmos {
        class Atmosphere1D;
    }
}  // namespace NCPA

static void swap( NCPA::atmos::Atmosphere1D&, NCPA::atmos::Atmosphere1D& ) noexcept;

namespace NCPA {
    namespace atmos {
        class Atmosphere1D : public AtmosphericModel {
            public:
                Atmosphere1D() : AtmosphericModel() {}

                Atmosphere1D(
                    std::unique_ptr<details::abstract_atmosphere_1d> ptr ) :
                    Atmosphere1D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                Atmosphere1D( const Atmosphere1D& other ) : Atmosphere1D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                Atmosphere1D( Atmosphere1D&& source ) noexcept :
                    Atmosphere1D() {
                    ::swap( *this, source );
                }

                virtual ~Atmosphere1D() {}

                friend void ::swap( Atmosphere1D& a,
                                    Atmosphere1D& b ) noexcept;

                Atmosphere1D& operator=( Atmosphere1D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual Atmosphere1D set_interpolator(
                    NCPA::interpolation::interpolator_type_t interp_type ) {
                        check_pointer();
                        if (NCPA::interpolation::InterpolatorFactory::can_build( interp_type )) {
                            _ptr->set_interpolator( interp_type );
                        } else {
                            throw std::logic_error( "Selected interpolator type not available" );
                        }
                        return *this;
                    }

                virtual Atmosphere1D& add_property(
                    const std::string& key,
                    const AtmosphericProperty1D& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere1D& add_property(
                    const std::string& key, const vector_t& property ) {
                    check_pointer();
                     _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere1D& add_property(
                    const std::string& key, const scalar_t& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere1D remove_property(
                    const std::string& key ) {
                    check_pointer();
                    _ptr->remove_property( key );
                    return *this;
                }

                virtual Atmosphere1D& copy_property(
                    const std::string& old_key, const std::string& new_key ) {
                    check_pointer();
                    _ptr->copy_property( old_key, new_key );
                    return *this;
                }

                virtual AtmosphericProperty1D& get_property(
                    const std::string& key ) const {
                    check_pointer();
                    return _ptr->get_property( key );
                }

                virtual double get( const std::string& key ) {
                    check_pointer();
                    return _ptr->get( key );
                }

                virtual double get( const std::string& key, double altitude ) {
                    check_pointer();
                    return _ptr->get( key, altitude );
                }

                virtual double get( const std::string& key, const NCPA::units::ScalarWithUnits<double>& altitude ) {
                    check_pointer();
                    return _ptr->get( key, altitude.get_as( this->get_altitude_units() ) );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double altitude ) {
                    check_pointer();
                    return _ptr->get_first_derivative( key, altitude );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     const NCPA::units::ScalarWithUnits<double>& altitude ) {
                    check_pointer();
                    return _ptr->get_first_derivative( key, altitude.get_as( this->get_altitude_units() ) );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double altitude ) {
                    check_pointer();
                    return _ptr->get_second_derivative( key, altitude );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      const NCPA::units::ScalarWithUnits<double>& altitude ) {
                    check_pointer();
                    return _ptr->get_second_derivative( key, altitude.get_as( this->get_altitude_units() ) );
                }

                virtual units_ptr_t get_altitude_units() {
                    check_pointer();
                    return _ptr->get_altitude_units();
                }

                virtual units_ptr_t get_property_units(
                    const std::string& key ) {
                    check_pointer();
                    return _ptr->get_property_units( key );
                }

                virtual double get_minimum_altitude() const {
                    check_pointer();
                    return _ptr->get_minimum_altitude();
                }

                virtual double get_maximum_altitude() const {
                    check_pointer();
                    return _ptr->get_maximum_altitude();
                }

                virtual Atmosphere1D& convert_altitude_units(
                    units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_altitude_units( new_units );
                    return *this;
                }

                virtual Atmosphere1D& convert_property_units(
                    const std::string& key, units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_property_units( key, new_units );
                    return *this;
                }

                virtual Atmosphere1D& resample( double new_dz ) {
                    check_pointer();
                    _ptr->resample( new_dz );
                    return *this;
                }

                virtual Atmosphere1D& resample( vector_t new_z ) {
                    check_pointer();
                    _ptr->resample( new_z );
                    return *this;
                }

                virtual std::vector<std::string> get_keys() const {
                    check_pointer();
                    return _ptr->get_keys();
                }

                virtual std::vector<std::string> get_vector_keys() const {
                    check_pointer();
                    return _ptr->get_vector_keys();
                }

                virtual std::vector<std::string> get_scalar_keys() const {
                    check_pointer();
                    return _ptr->get_scalar_keys();
                }

                virtual bool contains_scalar( const std::string& key ) const {
                    check_pointer();
                    return _ptr->contains_scalar( key );
                }

                virtual bool contains_vector( const std::string& key ) const {
                    check_pointer();
                    return _ptr->contains_vector( key );
                }

                virtual bool contains_key( const std::string& key ) const {
                    check_pointer();
                    return _ptr->contains_key( key );
                }

                virtual void check_pointer() const {
                    if ( !_ptr ) {
                        throw std::logic_error(
                            "Atmosphere1D: Internal pointer not set!" );
                    }
                }


                explicit operator bool() const { return (_ptr ? true : false); }

            private:
                std::unique_ptr<details::abstract_atmosphere_1d> _ptr;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::Atmosphere1D& a,
             NCPA::atmos::Atmosphere1D& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::atmos::AtmosphericModel&>( a ),
          static_cast<NCPA::atmos::AtmosphericModel&>( b ) );
    a._ptr.swap( b._ptr );
}
