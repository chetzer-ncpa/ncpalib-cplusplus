#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_3d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/AtmosphericModel.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/grid_atmosphere_3d.hpp"
#include "NCPA/atmosphere/stratified_atmosphere_3d.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/units.hpp"

#include <cfloat>
#include <memory>
#include <string>

static void swap( NCPA::atmos::Atmosphere3D&,
                  NCPA::atmos::Atmosphere3D& ) noexcept;

namespace NCPA {
    namespace atmos {
        class Atmosphere3D : public AtmosphericModel {
            public:
                Atmosphere3D() : AtmosphericModel() {}

                Atmosphere3D( std::unique_ptr<abstract_atmosphere_3d> ptr ) :
                    Atmosphere3D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                Atmosphere3D( const Atmosphere3D& other ) : Atmosphere3D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                Atmosphere3D( Atmosphere3D&& source ) noexcept :
                    Atmosphere3D() {
                    ::swap( *this, source );
                }

                virtual ~Atmosphere3D() {}

                friend void ::swap( Atmosphere3D& a,
                                    Atmosphere3D& b ) noexcept;

                Atmosphere3D& operator=( Atmosphere3D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                Atmosphere3D& set( Atmosphere1D& atmos1d ) {
                    check_pointer();
                    _ptr->set( *atmos1d.internal() );
                    return *this;
                }

                Atmosphere3D& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    NCPA::arrays::ndvector<2, Atmosphere1D>& components ) {
                    check_pointer();
                    NCPA::arrays::ndvector<2, abstract_atmosphere_1d *>
                        vectors;
                    vectors.reshape( components.shape() );
                    for (size_t i = 0; i < ax1.size(); ++i) {
                        for (size_t j = 0; j < ax2.size(); ++j) {
                            vectors[ i ][ j ]
                                = components[ i ][ j ].internal();
                        }
                    }
                    _ptr->set( ax1, ax2, vectors );
                    return *this;
                }

                Atmosphere3D& set(
                    NCPA::arrays::ndvector<2, Atmosphere1D>& components ) {
                    check_pointer();
                    NCPA::arrays::ndvector<2, abstract_atmosphere_1d *>
                        vectors;
                    auto dims = components.shape();
                    vectors.reshape( dims );

                    for (size_t i = 0; i < dims[ 0 ]; ++i) {
                        for (size_t j = 0; j < dims[ 1 ]; ++j) {
                            vectors[ i ][ j ]
                                = components[ i ][ j ].internal();
                        }
                    }
                    _ptr->set( vectors );
                    return *this;
                }

                virtual Atmosphere3D set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type ) {
                    check_pointer();
                    if (NCPA::interpolation::InterpolatorFactory<
                            double, double>::can_build( interp_type )) {
                        _ptr->set_interpolator( interp_type );
                    } else {
                        throw std::logic_error(
                            "Selected interpolator type not available" );
                    }
                    return *this;
                }

                virtual Atmosphere3D set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type ) {
                    check_pointer();
                    if (NCPA::interpolation::InterpolatorFactory<
                            double, double>::can_build( interp_type )) {
                        _ptr->set_interpolator( interp_type );
                    } else {
                        throw std::logic_error(
                            "Selected interpolator type not available" );
                    }
                    return *this;
                }

                virtual Atmosphere3D set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type ) {
                    check_pointer();
                    if (NCPA::interpolation::InterpolatorFactory<
                            double, double>::can_build( interp_type )) {
                        _ptr->set_interpolator( interp_type );
                    } else {
                        throw std::logic_error(
                            "Selected interpolator type not available" );
                    }
                    return *this;
                }

                virtual Atmosphere3D& add_property(
                    const std::string& key,
                    const AtmosphericProperty3D& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere3D& add_property(
                    const std::string& key, const vector3d_u_t& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere3D& add_property(
                    const std::string& key, const vector2d_u_t& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere3D& add_property(
                    const std::string& key, const scalar_u_t& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere3D remove_property(
                    const std::string& key ) {
                    check_pointer();
                    _ptr->remove_property( key );
                    return *this;
                }

                virtual Atmosphere3D& copy_property(
                    const std::string& old_key, const std::string& new_key ) {
                    check_pointer();
                    _ptr->copy_property( old_key, new_key );
                    return *this;
                }

                virtual vector_u_t get_axis_vector( size_t dim ) {
                    check_pointer();
                    return _ptr->get_axis_vector( dim );
                }

                virtual vector3d_u_t get_values( const std::string& key ) {
                    check_pointer();
                    return _ptr->get_values(key);
                }

                virtual double get( const std::string& key ) {
                    check_pointer();
                    return _ptr->get( key, 0.0, 0.0 );
                }

                virtual double get( const std::string& key, double x1,
                                    double x2 ) {
                    check_pointer();
                    return _ptr->get( key, x1, x2 );
                }

                virtual double get( const std::string& key, double x1,
                                    double x2, double x3 ) {
                    check_pointer();
                    return _ptr->get( key, x1, x2, x3 );
                }

                virtual vector3d_u_t get( const std::string& key,
                                          const std::vector<double>& v1,
                                          const std::vector<double>& v2,
                                          const std::vector<double>& v3 ) {
                    check_pointer();
                    return _ptr->get( key, v1, v2, v3 );
                }

                virtual vector2d_u_t get( const std::string& key,
                                          const std::vector<double>& v1,
                                          const std::vector<double>& v2 ) {
                    check_pointer();
                    return _ptr->get( key, v1, v2 );
                }

                virtual double get(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& x1,
                    const NCPA::units::ScalarWithUnits<double>& x2 ) {
                    check_pointer();
                    return _ptr->get( key,
                                      x1.get_as( this->get_axis_units( 0 ) ),
                                      x2.get_as( this->get_axis_units( 1 ) ) );
                }

                virtual double get(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& x1,
                    const NCPA::units::ScalarWithUnits<double>& x2,
                    const NCPA::units::ScalarWithUnits<double>& x3 ) {
                    check_pointer();
                    return _ptr->get( key,
                                      x1.get_as( this->get_axis_units( 0 ) ),
                                      x2.get_as( this->get_axis_units( 1 ) ),
                                      x3.get_as( this->get_axis_units( 2 ) ) );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double x1, double x2,
                                                     size_t wrt1 ) {
                    check_pointer();
                    return _ptr->get_first_derivative( key, x1, x2, wrt1 );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double x1, double x2,
                                                     double x3, size_t wrt1 ) {
                    check_pointer();
                    return _ptr->get_first_derivative( key, x1, x2, x3, wrt1 );
                }

                virtual double get_first_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& x1,
                    const NCPA::units::ScalarWithUnits<double>& x2,
                    size_t wrt1 ) {
                    check_pointer();
                    return _ptr->get_first_derivative(
                        key, x1.get_as( this->get_axis_units( 0 ) ),
                        x2.get_as( this->get_axis_units( 1 ) ), wrt1 );
                }

                virtual double get_first_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& x1,
                    const NCPA::units::ScalarWithUnits<double>& x2,
                    const NCPA::units::ScalarWithUnits<double>& x3,
                    size_t wrt1 ) {
                    check_pointer();
                    return _ptr->get_first_derivative(
                        key, x1.get_as( this->get_axis_units( 0 ) ),
                        x2.get_as( this->get_axis_units( 1 ) ),
                        x3.get_as( this->get_axis_units( 2 ) ), wrt1 );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double x1, double x2,
                                                      size_t wrt1,
                                                      size_t wrt2 ) {
                    check_pointer();
                    return _ptr->get_second_derivative( key, x1, x2, wrt1,
                                                        wrt2 );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double x1, double x2,
                                                      double x3, size_t wrt1,
                                                      size_t wrt2 ) {
                    check_pointer();
                    return _ptr->get_second_derivative( key, x1, x2, x3, wrt1,
                                                        wrt2 );
                }

                virtual double get_second_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& x1,
                    const NCPA::units::ScalarWithUnits<double>& x2,
                    size_t wrt1, size_t wrt2 ) {
                    check_pointer();
                    return _ptr->get_second_derivative(
                        key, x1.get_as( this->get_axis_units( 0 ) ),
                        x2.get_as( this->get_axis_units( 1 ) ), wrt1, wrt2 );
                }

                virtual double get_second_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& x1,
                    const NCPA::units::ScalarWithUnits<double>& x2,
                    const NCPA::units::ScalarWithUnits<double>& x3,
                    size_t wrt1, size_t wrt2 ) {
                    check_pointer();
                    return _ptr->get_second_derivative(
                        key, x1.get_as( this->get_axis_units( 0 ) ),
                        x2.get_as( this->get_axis_units( 1 ) ),
                        x3.get_as( this->get_axis_units( 2 ) ), wrt1, wrt2 );
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

                virtual Atmosphere3D& convert_axis_units(
                    size_t dim, units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_axis_units( dim, new_units );
                    return *this;
                }

                virtual Atmosphere3D& convert_units( const std::string& key,
                                                     units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_units( key, new_units );
                    return *this;
                }

                virtual Atmosphere3D& resample( size_t dim, double new_dz ) {
                    check_pointer();
                    _ptr->resample( dim, new_dz );
                    return *this;
                }

                virtual Atmosphere3D& resample( size_t dim,
                                                vector_u_t new_z ) {
                    check_pointer();
                    _ptr->resample( dim, new_z );
                    return *this;
                }

                virtual Atmosphere3D& resample( vector_u_t new_ax1,
                                                vector_u_t new_ax2,
                                                vector_u_t new_ax3 ) {
                    check_pointer();
                    _ptr->resample( new_ax1, new_ax2, new_ax3 );
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

                virtual bool is_stratified() const override {
                    check_pointer();
                    return _ptr->is_stratified();
                }

                virtual bool same( scalar_u_t x1_1, scalar_u_t x2_1,
                                   scalar_u_t x1_2, scalar_u_t x2_2 ) const {
                    check_pointer();
                    return _ptr->same( x1_1, x2_1, x1_2, x2_2 );
                }

                virtual bool same( double x1_1, double x2_1, double x1_2,
                                   double x2_2 ) const {
                    check_pointer();
                    return _ptr->same( x1_1, x2_1, x1_2, x2_2 );
                }

                virtual void check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "Atmosphere3D: Internal pointer not set!" );
                    }
                }

                // friend binary operators
                friend std::ostream& operator<<( std::ostream& os,
                                                 const Atmosphere3D& atm ) {
                    if (atm) {
                        atm._ptr->print( os );
                    }
                    return os;
                }

                explicit operator bool() const {
                    return ( _ptr ? true : false );
                }

            private:
                std::unique_ptr<abstract_atmosphere_3d> _ptr;
        };


    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::Atmosphere3D& a,
                  NCPA::atmos::Atmosphere3D& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::atmos::AtmosphericModel&>( a ),
          static_cast<NCPA::atmos::AtmosphericModel&>( b ) );
    a._ptr.swap( b._ptr );
}
