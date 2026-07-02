#pragma once

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/abstract_atmospheric_property_1d.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <memory>
#include <stdexcept>

static void swap( NCPA::atmos::tuple_atmospheric_property_1d& a,
                  NCPA::atmos::tuple_atmospheric_property_1d& b ) noexcept;

namespace NCPA {
    namespace atmos {
        class tuple_atmospheric_property_1d
            : public abstract_atmospheric_property_1d {
            public:
                tuple_atmospheric_property_1d() {
                    set_interpolator(
                        NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR );
                    set_extrapolator( NCPA::interpolation::
                                          extrapolator_1d_type_t::CONSTANT );
                }

                tuple_atmospheric_property_1d( const vector_u_t& z,
                                               const vector_u_t& p ) :
                    tuple_atmospheric_property_1d() {
                    set( z, p );
                }

                tuple_atmospheric_property_1d(
                    const tuple_atmospheric_property_1d& source ) :
                    tuple_atmospheric_property_1d() {
                    set( source._z, source._x );
                    set_interpolator( source._spline.interptype() );
                    _init_spline();
                }

                virtual tuple_atmospheric_property_1d& clear() override {
                    _z.clear();
                    _x.clear();
                    _spline.clear();
                    return *this;
                }

                tuple_atmospheric_property_1d(
                    tuple_atmospheric_property_1d&& source ) noexcept :
                    tuple_atmospheric_property_1d() {
                    ::swap( *this, source );
                }

                virtual ~tuple_atmospheric_property_1d() {}

                tuple_atmospheric_property_1d& operator=(
                    tuple_atmospheric_property_1d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap(
                    tuple_atmospheric_property_1d& a,
                    tuple_atmospheric_property_1d& b ) noexcept;

                virtual std::unique_ptr<abstract_atmospheric_property> clone()
                    const override {
                    return std::unique_ptr<abstract_atmospheric_property>(
                        new tuple_atmospheric_property_1d( *this ) );
                }

                virtual std::unique_ptr<abstract_atmospheric_property_1d>
                    clone1d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_1d>(
                        new tuple_atmospheric_property_1d( *this ) );
                }

                virtual size_t size() const override { return _z.size(); }

                virtual tuple_atmospheric_property_1d& set(
                    const vector_u_t& z, const units_ptr_t units ) override {
                    this->clear();
                    _z = z;
                    _x = vector_u_t( _z.size(), units );
                    this->_init_spline();
                    return *this;
                }

                virtual tuple_atmospheric_property_1d& copy(
                    const abstract_atmospheric_property& source ) override {
                    if (source.dimensions() == 1) {
                        const abstract_atmospheric_property_1d *srcptr
                            = dynamic_cast<
                                const abstract_atmospheric_property_1d *>(
                                &source );
                        set( srcptr->axis(), srcptr->values() );
                    } else {
                        throw std::invalid_argument(
                            "tuple_atmospheric_property_1d.copy(): Cannot "
                            "copy property with more than 1 dimension" );
                    }
                    return *this;
                }

                virtual tuple_atmospheric_property_1d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    _init_spline();
                    return *this;
                }

                virtual tuple_atmospheric_property_1d& set_extrapolator(
                    NCPA::interpolation::extrapolator_1d_type_t extrap_type )
                    override {
                    this->interpolator().set_extrapolation( extrap_type );
                    return *this;
                }

                virtual NCPA::interpolation::Interpolator1D<double, double>&
                    interpolator() override {
                    return _spline;
                }

                virtual tuple_atmospheric_property_1d& set(
                    const vector_u_t& z, const vector_u_t& p ) override {
                    if (z.size() != p.size()) {
                        throw std::invalid_argument(
                            "Altitude and property vectors are not the same "
                            "size!" );
                    }
                    _z = z;
                    _x = p;
                    _init_spline();
                    return *this;
                }

                virtual vector_u_t& values() override { return _x; }

                virtual const vector_u_t& values() const override {
                    return _x;
                }

                virtual vector_u_t& axis() override { return _z; }

                virtual const vector_u_t& axis() const override { return _z; }

                virtual double get( double altitude ) override {
                    return _spline.eval_f( altitude );
                }

                virtual double get_first_derivative(
                    double altitude ) override {
                    return _spline.eval_df( altitude );
                }

                virtual double get_second_derivative(
                    double altitude ) override {
                    return _spline.eval_ddf( altitude );
                }

                virtual const units_ptr_t get_units() const override {
                    return _x.get_units();
                }

                virtual const units_ptr_t get_axis_units() const override {
                    return _z.get_units();
                }

                virtual tuple_atmospheric_property_1d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _x.convert_units( u );
                    _init_spline();
                    return *this;
                }

                virtual tuple_atmospheric_property_1d& convert_axis_units(
                    const NCPA::units::Unit& u ) override {
                    _z.convert_units( u );
                    _init_spline();
                    return *this;
                }

                virtual tuple_atmospheric_property_1d& resample(
                    vector_u_t new_z ) override {
                    if (!_spline) {
                        throw std::logic_error(
                            "Cannot resample, interpolator has not been "
                            "set." );
                    }
                    _z.convert_units( *new_z.get_units() );
                    vector_u_t new_x( new_z.size() );
                    new_x.set_units( *_x.get_units() );
                    for (size_t i = 0; i < new_z.size(); i++) {
                        new_x[ i ] = this->get( new_z[ i ] );
                    }
                    _x = new_x;
                    _z = new_z;
                    _init_spline();
                    return *this;
                }

            protected:
                void _init_spline() {
                    if (_z && _x && _spline) {
                        _spline.init( size() ).fill( _z, _x ).ready();
                    }
                }

            protected:
                vector_u_t _z, _x;
                NCPA::interpolation::Interpolator1D<double, double> _spline;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::tuple_atmospheric_property_1d& a,
                  NCPA::atmos::tuple_atmospheric_property_1d& b ) noexcept {
    using std::swap;
    ::swap(
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_1d&>( a ),
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_1d&>( b ) );
    swap( a._z, b._z );
    swap( a._x, b._x );
    swap( a._spline, b._spline );
}
