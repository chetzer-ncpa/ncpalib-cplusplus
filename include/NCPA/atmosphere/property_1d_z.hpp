#pragma once

// #include "NCPA/atmosphere/types.hpp"
#include "NCPA/atmosphere/abstract/atmospheric_property_1d.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <memory>
#include <stdexcept>


static void swap( NCPA::atmos::VerticalProperty1D& a,
                  NCPA::atmos::VerticalProperty1D& b ) noexcept;

bool operator==( const NCPA::atmos::VerticalProperty1D& a,
                 const NCPA::atmos::VerticalProperty1D& b );

bool operator!=( const NCPA::atmos::VerticalProperty1D& a,
                 const NCPA::atmos::VerticalProperty1D& b );

namespace NCPA {
    namespace atmos {
        class VerticalProperty1D
            : public abstract::atmospheric_property_1d {
            public:
                VerticalProperty1D() {
                    this->set_interpolator(
                        NCPA::interpolation::interpolator_1d_type_t::
                            NEAREST_NEIGHBOR );
                }

                VerticalProperty1D( const vector_u_t& z,
                                       const vector_u_t& p ) :
                    VerticalProperty1D() {
                    set( z, p );
                }

                VerticalProperty1D( const VerticalProperty1D& source ) :
                    VerticalProperty1D() {
                    set( source._z, source._x );
                }

                VerticalProperty1D(
                    VerticalProperty1D&& source ) noexcept :
                    VerticalProperty1D() {
                    ::swap( *this, source );
                }

                virtual ~VerticalProperty1D() {}

                VerticalProperty1D& operator=(
                    VerticalProperty1D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( VerticalProperty1D& a,
                                    VerticalProperty1D& b ) noexcept;

                virtual size_t size() const override { return _z.size(); }

                // virtual VerticalProperty1D& set_interpolator(
                //     NCPA::interpolation::interpolator_1d_type_t interp_type
                //     ) override { _spline =
                //     NCPA::interpolation::InterpolatorFactory<
                //         double, double>::build( interp_type );
                //     _init_spline();
                //     return *this;
                // }

                virtual double get( double x1 ) override {
                    return this->spline().eval_f( x1 );
                }

                virtual double get_first_derivative( double x1 ) override {
                    return this->spline().eval_df( x1 );
                }

                virtual double get_second_derivative( double x1 ) override {
                    return this->spline().eval_ddf( x1 );
                }

                virtual atmospheric_property_1d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t
                        interp_type ) override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    this->setup_spline();
                    return *this;
                }

                virtual abstract::atmospheric_property_1d& set(
                    const vector_u_t& z, const vector_u_t& p ) override {
                    if ( z.size() != p.size() ) {
                        throw std::invalid_argument(
                            "Altitude and property vectors are not the same "
                            "size!" );
                    }
                    _z = z;
                    _x = p;
                    setup_spline();
                    RETURN_THIS_AS_NCPA_ATMOS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D;
                    // return *this;
                }

                virtual vector_u_t& vector() override { return _x; }

                virtual const vector_u_t& vector() const override {
                    return _x;
                }

                // virtual double get( double altitude ) override {
                //     return _spline.eval_f( altitude );
                // }

                // virtual double get_first_derivative(
                //     double altitude ) override {
                //     return _spline.eval_df( altitude );
                // }

                // virtual double get_second_derivative(
                //     double altitude ) override {
                //     return _spline.eval_ddf( altitude );
                // }

                virtual const units_ptr_t get_units() const override {
                    return _x.get_units();
                }

                virtual const units_ptr_t get_independent_units(
                    size_t dimension ) const override {
                    return _z.get_units();
                }

                virtual const units_ptr_t get_altitude_units() const {
                    return _z.get_units();
                }

                virtual abstract::atmospheric_property& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _x.convert_units( u );
                    // _init_spline();
                    return this->setup_spline();
                    // RETURN_THIS_AS_NCPA_ATMOS_ABSTRACT_ATMOSPHERIC_PROPERTY;
                    // return *this;
                }

                virtual abstract::atmospheric_property&
                    convert_independent_units(
                        size_t dimension,
                        const NCPA::units::Unit& u ) override {
                    _x.convert_units( u );
                    return this->setup_spline();
                    // _init_spline();
                    // RETURN_THIS_AS_NCPA_ATMOS_ABSTRACT_ATMOSPHERIC_PROPERTY;
                    // return *this;
                }

                virtual VerticalProperty1D& convert_altitude_units(
                    const NCPA::units::Unit& u ) {
                    // auto parent = this->convert_independent_units( 0, u );
                    // _init_spline();
                    return dynamic_cast<VerticalProperty1D&>(
                        this->convert_independent_units( 0, u )
                            .setup_spline() );
                    // return *this;
                }

                virtual abstract::atmospheric_property_1d& resample(
                    vector_u_t new_z ) override {
                    if ( !_spline ) {
                        throw std::logic_error(
                            "Cannot resample, interpolator has not been "
                            "set." );
                    }
                    _z.convert_units( *new_z.get_units() );
                    vector_u_t new_x( new_z.size() );
                    new_x.set_units( *_x.get_units() );
                    for ( size_t i = 0; i < new_z.size(); i++ ) {
                        new_x[ i ] = this->get( new_z[ i ] );
                    }
                    _x = new_x;
                    _z = new_z;
                    this->setup_spline();
                    RETURN_THIS_AS_NCPA_ATMOS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D;
                }

                virtual abstract::atmospheric_property& setup_spline()
                    override {
                    if ( _z && _x && this->spline() ) {
                        this->spline().init( size() );
                        this->spline().fill( _z, _x );
                    }
                    RETURN_THIS_AS_NCPA_ATMOS_ABSTRACT_ATMOSPHERIC_PROPERTY;
                }

            protected:
                vector_u_t _z, _x;
                NCPA::interpolation::Interpolator1D<double, double>
                        _spline;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::VerticalProperty1D& a,
                  NCPA::atmos::VerticalProperty1D& b ) noexcept {
    using std::swap;
    ::swap(
        dynamic_cast<NCPA::atmos::abstract::atmospheric_property_1d&>( a ),
        dynamic_cast<NCPA::atmos::abstract::atmospheric_property_1d&>( b ) );
    swap( a._z, b._z );
    swap( a._x, b._x );
    swap( a._spline, b._spline );
}
