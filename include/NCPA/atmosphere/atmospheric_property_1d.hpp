#pragma once

#include "NCPA/atmosphere/types.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace atmos {
        class AtmosphericProperty1D;
    }
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphericProperty1D& a,
                  NCPA::atmos::AtmosphericProperty1D& b ) noexcept;

bool operator==( const NCPA::atmos::AtmosphericProperty1D& a,
                 const NCPA::atmos::AtmosphericProperty1D& b );

bool operator!=( const NCPA::atmos::AtmosphericProperty1D& a,
                 const NCPA::atmos::AtmosphericProperty1D& b );

namespace NCPA {
    namespace atmos {
        class AtmosphericProperty1D {
            public:
                AtmosphericProperty1D() {
                    set_interpolator(
                        NCPA::interpolation::interpolator_type_t::
                            NEAREST_NEIGHBOR );
                }

                AtmosphericProperty1D( const vector_t& z, const vector_t& p ) :
                    AtmosphericProperty1D() {
                    set( z, p );
                }

                AtmosphericProperty1D( const AtmosphericProperty1D& source ) :
                    AtmosphericProperty1D() {
                    set( source._z, source._x );
                }

                AtmosphericProperty1D(
                    AtmosphericProperty1D&& source ) noexcept :
                    AtmosphericProperty1D() {
                    ::swap( *this, source );
                }

                virtual ~AtmosphericProperty1D() {}

                AtmosphericProperty1D& operator=(
                    AtmosphericProperty1D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( AtmosphericProperty1D& a,
                                    AtmosphericProperty1D& b ) noexcept;

                virtual size_t size() const { return _z.size(); }

                virtual AtmosphericProperty1D& set_interpolator(
                    NCPA::interpolation::interpolator_type_t interp_type ) {
                    _spline = NCPA::interpolation::InterpolatorFactory::build<
                        double, double>( interp_type );
                    _init_spline();
                    return *this;
                }

                virtual AtmosphericProperty1D& set( const vector_t& z,
                                                    const vector_t& p ) {
                    if ( z.size() != p.size() ) {
                        throw std::invalid_argument(
                            "Altitude and property vectors are not the same "
                            "size!" );
                    }
                    _z = z;
                    _x = p;
                    _init_spline();
                    return *this;
                }

                virtual vector_t& vector() {
                    return _x;
                }

                virtual double get( double altitude ) {
                    return _spline.eval_f( altitude );
                }

                virtual double get_first_derivative( double altitude ) {
                    return _spline.eval_df( altitude );
                }

                virtual double get_second_derivative( double altitude ) {
                    return _spline.eval_ddf( altitude );
                }

                virtual const units_ptr_t get_units() const {
                    return _x.get_units();
                }

                virtual const units_ptr_t get_altitude_units() const {
                    return _z.get_units();
                }

                virtual AtmosphericProperty1D& convert_units(
                    const NCPA::units::Unit& u ) {
                    _x.convert_units( u );
                    _init_spline();
                    return *this;
                }

                virtual AtmosphericProperty1D& convert_altitude_units(
                    const NCPA::units::Unit& u ) {
                    _z.convert_units( u );
                    _init_spline();
                    return *this;
                }

                virtual AtmosphericProperty1D& resample( vector_t new_z ) {
                    if ( !_spline ) {
                        throw std::logic_error(
                            "Cannot resample, interpolator has not been "
                            "set." );
                    }
                    _z.convert_units( *new_z.get_units() );
                    vector_t new_x( new_z.size() );
                    new_x.set_units( *_x.get_units() );
                    for ( size_t i = 0; i < new_z.size(); i++ ) {
                        new_x[ i ] = this->get( new_z[ i ] );
                    }
                    _x = new_x;
                    _z = new_z;
                    _init_spline();
                    return *this;
                }

            protected:
                void _init_spline() {
                    if ( _z && _x && _spline ) {
                        _spline.init( size() );
                        _spline.fill( _z, _x );
                    }
                }

            private:
                vector_t _z, _x;
                NCPA::interpolation::Interpolator1D<double, double> _spline;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphericProperty1D& a,
                  NCPA::atmos::AtmosphericProperty1D& b ) noexcept {
    using std::swap;
    swap( a._z, b._z );
    swap( a._x, b._x );
    swap( a._spline, b._spline );
}
