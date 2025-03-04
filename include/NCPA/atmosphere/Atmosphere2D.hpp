#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_2d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/AtmosphericModel.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/units.hpp"

#include <memory>
#include <string>

static void swap( NCPA::atmos::Atmosphere2D&,
                  NCPA::atmos::Atmosphere2D& ) noexcept;

namespace NCPA {
    namespace atmos {
        class Atmosphere2D : public AtmosphericModel {
            public:
                Atmosphere2D() : AtmosphericModel() {}

                Atmosphere2D( std::unique_ptr<abstract_atmosphere_2d> ptr ) :
                    Atmosphere2D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                Atmosphere2D( const Atmosphere2D& other ) : Atmosphere2D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                Atmosphere2D( Atmosphere2D&& source ) noexcept :
                    Atmosphere2D() {
                    ::swap( *this, source );
                }

                virtual ~Atmosphere2D() {}

                friend void ::swap( Atmosphere2D& a,
                                    Atmosphere2D& b ) noexcept;

                Atmosphere2D& operator=( Atmosphere2D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                // Atmosphere2D& set( const vector_u_t ranges, const std::vector<const abstract_atmosphere_1d*> atmos1ds ) {
                //     check_pointer();
                //     _ptr->set( ranges, atmos1ds );
                //     return *this;
                // }

                Atmosphere2D& set(  Atmosphere1D& atmos1d ) {
                    check_pointer();
                    _ptr->set( *atmos1d.internal() );
                    // return this->set( scalar_u_t( 0.0, this->get_axis_units(0) ), atmos1d );
                    return *this;
                }

                Atmosphere2D& append( const scalar_u_t& range, const Atmosphere1D& atmos1d ) {
                    check_pointer();
                    _ptr->append( range, *atmos1d.internal()->clone().get() );
                    return *this;
                }

                virtual Atmosphere2D set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type ) {
                    check_pointer();
                    if ( NCPA::interpolation::InterpolatorFactory<
                             double, double>::can_build( interp_type ) ) {
                        _ptr->set_interpolator( interp_type );
                    } else {
                        throw std::logic_error(
                            "Selected interpolator type not available" );
                    }
                    return *this;
                }

                virtual Atmosphere2D set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type ) {
                    check_pointer();
                    if ( NCPA::interpolation::InterpolatorFactory<
                             double, double>::can_build( interp_type ) ) {
                        _ptr->set_interpolator( interp_type );
                    } else {
                        throw std::logic_error(
                            "Selected interpolator type not available" );
                    }
                    return *this;
                }

                virtual Atmosphere2D& add_property(
                    const std::string& key,
                    const AtmosphericProperty2D& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere2D& add_property(
                    const std::string& key, const vector2d_u_t& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere2D& add_property(
                    const std::string& key, const scalar_u_t& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere2D remove_property(
                    const std::string& key ) {
                    check_pointer();
                    _ptr->remove_property( key );
                    return *this;
                }

                virtual Atmosphere2D& copy_property(
                    const std::string& old_key, const std::string& new_key ) {
                    check_pointer();
                    _ptr->copy_property( old_key, new_key );
                    return *this;
                }

                // virtual AtmosphericProperty2D& get_property(
                //     const std::string& key ) const {
                //     check_pointer();
                //     return _ptr->get_property( key );
                // }

                virtual vector_u_t get_axis_vector( size_t dim ) {
                    check_pointer();
                    return _ptr->get_axis_vector( dim );
                }

                virtual double get( const std::string& key ) {
                    check_pointer();
                    return _ptr->get( key, 0.0 );
                }

                virtual double get( const std::string& key, double range,
                                    double altitude ) {
                    check_pointer();
                    return _ptr->get( key, range, altitude );
                }

                virtual double get(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& range,
                    const NCPA::units::ScalarWithUnits<double>& altitude ) {
                    check_pointer();
                    return _ptr->get(
                        key, range.get_as( this->get_axis_units( 0 ) ),
                        altitude.get_as( this->get_axis_units( 1 ) ) );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double range,
                                                     double altitude,
                                                     size_t wrt1 ) {
                    check_pointer();
                    return _ptr->get_first_derivative( key, range, altitude,
                                                       wrt1 );
                }

                virtual double get_first_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& range,
                    const NCPA::units::ScalarWithUnits<double>& altitude,
                    size_t wrt1 ) {
                    check_pointer();
                    return _ptr->get_first_derivative(
                        key, range.get_as( this->get_axis_units( 0 ) ),
                        altitude.get_as( this->get_axis_units( 1 ) ), wrt1 );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double range,
                                                      double altitude,
                                                      size_t wrt1,
                                                      size_t wrt2 ) {
                    check_pointer();
                    return _ptr->get_second_derivative( key, range, altitude,
                                                        wrt1, wrt2 );
                }

                virtual double get_second_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& range,
                    const NCPA::units::ScalarWithUnits<double>& altitude,
                    size_t wrt1, size_t wrt2 ) {
                    check_pointer();
                    return _ptr->get_second_derivative(
                        key, range.get_as( this->get_axis_units( 0 ) ),
                        altitude.get_as( this->get_axis_units( 1 ) ), wrt1,
                        wrt2 );
                }

                virtual units_ptr_t get_axis_units( size_t dim ) {
                    check_pointer();
                    return _ptr->get_axis_units( dim );
                }

                virtual units_ptr_t get_units( const std::string& key ) {
                    check_pointer();
                    return _ptr->get_units( key );
                }

                virtual double get_minimum_axis( size_t dim ) const {
                    check_pointer();
                    return _ptr->get_minimum_axis( dim );
                }

                virtual double get_maximum_axis( size_t dim ) const {
                    check_pointer();
                    return _ptr->get_maximum_axis( dim );
                }

                virtual Atmosphere2D& convert_axis_units(
                    size_t dim, units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_axis_units( dim, new_units );
                    return *this;
                }

                virtual Atmosphere2D& convert_units( const std::string& key,
                                                     units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_units( key, new_units );
                    return *this;
                }

                virtual Atmosphere2D& resample( size_t dim, double new_dz ) {
                    check_pointer();
                    _ptr->resample( dim, new_dz );
                    return *this;
                }

                virtual Atmosphere2D& resample( size_t dim, vector_u_t new_z ) {
                    check_pointer();
                    _ptr->resample( dim, new_z );
                    return *this;
                }

                virtual Atmosphere2D& resample( vector_u_t new_ax1, vector_u_t new_ax2 ) {
                    check_pointer();
                    _ptr->resample( new_ax1, new_ax2 );
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
                            "Atmosphere2D: Internal pointer not set!" );
                    }
                }

                // friend binary operators
                friend std::ostream& operator<<( std::ostream& os,
                                                 const Atmosphere2D& atm ) {
                    if ( atm ) {
                        atm._ptr->print( os );
                    }
                    return os;
                }

                explicit operator bool() const {
                    return ( _ptr ? true : false );
                }

            private:
                std::unique_ptr<abstract_atmosphere_2d> _ptr;
        };


    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::Atmosphere2D& a,
                  NCPA::atmos::Atmosphere2D& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::atmos::AtmosphericModel&>( a ),
          static_cast<NCPA::atmos::AtmosphericModel&>( b ) );
    a._ptr.swap( b._ptr );
}
