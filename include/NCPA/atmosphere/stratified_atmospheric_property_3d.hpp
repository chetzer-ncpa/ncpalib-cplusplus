#pragma once

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/abstract_atmospheric_property_3d.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/exceptions.hpp"
#include "NCPA/interpolation.hpp"

#include <array>
#include <stdexcept>
#include <unordered_map>

static void swap(
    NCPA::atmos::stratified_atmospheric_property_3d& a,
    NCPA::atmos::stratified_atmospheric_property_3d& b ) noexcept;

namespace NCPA {
    namespace atmos {
        class stratified_atmospheric_property_3d
            : public abstract_atmospheric_property_3d {
            public:
                using abstract_atmospheric_property_3d::get;
                using abstract_atmospheric_property_3d::resample;

                stratified_atmospheric_property_3d() {
                    _dummy = vector_u_t( { 0 }, NCPA::units::KILOMETERS );
                    set_interpolator(
                        NCPA::interpolation::interpolator_1d_type_t::
                            NEAREST_NEIGHBOR );
                }

                stratified_atmospheric_property_3d( const vector_u_t& z,
                                                    const vector_u_t& p ) :
                    stratified_atmospheric_property_3d() {
                    if (p.size() != z.size()) {
                        throw std::range_error( "Values dimensions do not "
                                                "match axes dimensions!" );
                    }
                    _z    = z;
                    _vals = p;
                    this->_init_spline();
                    _copy_to_3d();
                }

                stratified_atmospheric_property_3d( const vector_u_t& ax1,
                                                    const vector_u_t& ax2,
                                                    const vector_u_t& ax3,
                                                    const vector3d_u_t& p ) :
                    stratified_atmospheric_property_3d() {
                    _z    = ax3;
                    _vals = vector_u_t( p[ 0 ][ 0 ], p.get_units() );
                    _copy_to_3d();
                }

                stratified_atmospheric_property_3d(
                    const stratified_atmospheric_property_3d& source ) :
                    stratified_atmospheric_property_3d() {
                    _z        = source._z;
                    _vals     = source._vals;
                    _z_limits = source._z_limits;
                    _spline   = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( source._spline.interptype() );
                    _as_3d = source._as_3d;
                    _dummy = source._dummy;
                }

                stratified_atmospheric_property_3d(
                    stratified_atmospheric_property_3d&& source ) noexcept :
                    stratified_atmospheric_property_3d() {
                    ::swap( *this, source );
                }

                virtual ~stratified_atmospheric_property_3d() {}

                stratified_atmospheric_property_3d& operator=(
                    stratified_atmospheric_property_3d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap(
                    stratified_atmospheric_property_3d& a,
                    stratified_atmospheric_property_3d& b ) noexcept;

                virtual std::unique_ptr<abstract_atmospheric_property> clone()
                    const override {
                    return std::unique_ptr<abstract_atmospheric_property>(
                        new stratified_atmospheric_property_3d( *this ) );
                }

                virtual stratified_atmospheric_property_3d& clear() override {
                    _z.clear();
                    _vals.clear();
                    _z_limits = std::pair<double, double>( { 0, 0 } );
                    _spline.clear();
                    _as_3d.clear();
                    return *this;
                }

                virtual stratified_atmospheric_property_3d& copy(
                    const abstract_atmospheric_property& source ) override {
                    vector_u_t dummy( 1, &NCPA::units::KILOMETERS );
                    vector3d_u_t strat;
                    const abstract_atmospheric_property_1d *source1d;
                    const abstract_atmospheric_property_3d *source3d;
                    switch (source.dimensions()) {
                        case 1:
                            source1d = static_cast<
                                const abstract_atmospheric_property_1d *>(
                                &source );
                            _z    = source1d->axis();
                            _vals = source1d->values();
                            _copy_to_3d();
                            _init_spline();
                            break;
                        case 3:
                            source3d = static_cast<
                                const abstract_atmospheric_property_3d *>(
                                &source );
                            this->set( dummy, dummy, source3d->axis( 2 ),
                                       source3d->values() );
                            break;
                        default:
                            throw std::range_error(
                                "stratified_atmospheric_property_3d.copy(): "
                                "Supplied property must be 1-D or 3-D" );
                    }
                    return *this;
                }

                virtual size_t dim( size_t n ) const override {
                    _validate_axis( n );
                    return ( n == 2 ? _z.size() : 1 );
                }

                virtual stratified_atmospheric_property_3d& set(
                    size_t ax, const vector_u_t& axvals ) override {
                    if (ax == 2) {
                        if (axvals.size() != _vals.size()) {
                            std::ostringstream oss;
                            oss << "stratified_atmospheric_property_3d.set(): "
                                   "size "
                                   "mismatch between provided vector"
                                   " of size "
                                << axvals.size() << " and existing axis " << ax
                                << " of size " << _vals.size();
                            throw std::range_error( oss.str() );
                        }
                        _z = axvals;
                        reset_limits( ax );
                        _init_spline();
                        _copy_to_3d();
                    }
                    return *this;
                }

                virtual stratified_atmospheric_property_3d& set(
                    const vector3d_u_t& ax ) override {
                    if (ax.dim( 2 ) != _vals.size()) {
                        std::ostringstream oss;
                        oss << "stratospheric_atmospheric_property_3d.set(): "
                               "Size "
                               "mismatch between existing property vertical "
                               "dimension "
                            << _vals.size()
                            << " and supplied property vertical dimension "
                            << ax.dim( 2 );
                        throw std::range_error( oss.str() );
                    }
                    _vals = vector_u_t( ax[ 0 ][ 0 ], ax.get_units() );
                    _init_spline();
                    _copy_to_3d();
                    return *this;
                }

                virtual stratified_atmospheric_property_3d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector_u_t& ax3,
                    const vector3d_u_t& vals ) override {
                    this->clear();
                    _z    = ax3;
                    _vals = vector_u_t( vals[ 0 ][ 0 ], vals.get_units() );
                    reset_limits( 2 );
                    _init_spline();
                    _copy_to_3d();
                    return *this;
                }

                virtual stratified_atmospheric_property_3d& set(
                    size_t nx1, size_t nx2, size_t nx3,
                    const units_ptr_t units ) override {
                    this->clear();
                    _z    = vector_u_t( nx3 );
                    _vals = vector_u_t( nx3, units );
                    _init_spline();
                    _copy_to_3d();
                    return *this;
                }

                virtual stratified_atmospheric_property_3d& append(
                    size_t dim, double dimval,
                    const vector2d_u_t& newslice ) override {
                    if (dim != 2) {
                        throw std::range_error(
                            "stratified_atmospheric_property_3d.append(): "
                            "Cannot append to stratified atmosphere except in "
                            "third dimension!" );
                    }
                    if (dimval >= _z.back()) {
                        std::ostringstream oss;
                        oss << "stratified_atmospheric_property_3d.append(): "
                               "New coordinate "
                            << dimval
                            << " is within existing limits, cannot append";
                        throw std::range_error( oss.str() );
                    }
                    auto converted = newslice;
                    converted.convert_units( *_vals.get_units() );
                    _z.push_back( dimval );
                    _vals.push_back( converted[ 0 ][ 0 ] );
                    _copy_to_3d();
                    _init_spline();
                    return *this;
                }

                virtual std::unique_ptr<abstract_atmospheric_property_3d>
                    clone3d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_3d>(
                        new stratified_atmospheric_property_3d( *this ) );
                }

                virtual std::vector<size_t> shape() const override {
                    return std::vector<size_t> { 1, 1, _z.size() };
                }

                virtual stratified_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Stratified property requires a 1-D interpolator" );
                }

                virtual stratified_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Stratified property requires a 1-D interpolator" );
                }

                virtual stratified_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    this->_init_spline();
                    return *this;
                }

                virtual vector_u_t& axis( size_t n ) override {
                    this->_validate_axis( n );
                    return ( n == 2 ? _z : _dummy );
                }

                virtual const vector_u_t& axis( size_t n ) const override {
                    this->_validate_axis( n );
                    return ( n == 2 ? _z : _dummy );
                }

                virtual stratified_atmospheric_property_3d& set_limits(
                    size_t dim, double minax, double maxax ) override {
                    this->_validate_axis( dim );
                    if (dim == 2) {
                        if (_z_limits.first < minax
                            || _z_limits.second > maxax) {
                            std::ostringstream oss;
                            oss << "Error setting axis " << dim
                                << " limits to [" << minax << "," << maxax
                                << "]: cannot contract existing axes.";
                            throw std::range_error( oss.str() );
                        }
                        _z_limits.first  = minax;
                        _z_limits.second = maxax;
                    }
                    return *this;
                }

                virtual stratified_atmospheric_property_3d& reset_limits(
                    size_t dim ) override {
                    this->_validate_axis( dim );
                    if (dim == 2) {
                        _z_limits.first  = _z.front();
                        _z_limits.second = _z.back();
                    }
                    return *this;
                }

                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const override {
                    this->_validate_axis( dim );
                    return ( dim == 2
                                 ? _z_limits
                                 : std::pair<double, double> { 0.0, 0.0 } );
                }

                virtual vector3d_u_t& values() override { return _as_3d; }

                virtual const vector3d_u_t& values() const override {
                    return _as_3d;
                }

                virtual double get( double val1, double val2,
                                    double val3 ) override {
                    this->_snap_to_limits( val3 );
                    return _spline.eval_f( val3 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     double val3,
                                                     size_t rel ) override {
                    _validate_axis( rel );
                    this->_snap_to_limits( val3 );
                    return ( rel == 2 ? _spline.eval_df( val3 ) : 0.0 );
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      double val3, size_t rel1,
                                                      size_t rel2 ) override {
                    _validate_axis( rel1 );
                    _validate_axis( rel2 );
                    this->_snap_to_limits( val3 );
                    return ( rel1 == 2 && rel2 == 2 ? _spline.eval_ddf( val3 )
                                                    : 0.0 );
                }

                virtual const units_ptr_t get_units() const override {
                    return _vals.get_units();
                }

                virtual const units_ptr_t get_axis_units(
                    size_t n ) const override {
                    _validate_axis( n );
                    return ( n == 2 ? _z.get_units() : _dummy.get_units() );
                }

                virtual stratified_atmospheric_property_3d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _vals.convert_units( u );
                    _init_spline();
                    return *this;
                }

                virtual stratified_atmospheric_property_3d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) override {
                    _validate_axis( n );
                    if (n == 2) {
                        _z.convert_units( u );
                        this->reset_limits( 2 );
                        _init_spline();
                    } else {
                        _dummy.convert_units( u );
                    }
                    return *this;
                }

                virtual stratified_atmospheric_property_3d& resample(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector_u_t& ax3 ) override {
                    if (!_spline) {
                        throw std::logic_error(
                            "Cannot resample, interpolator has not been "
                            "set." );
                    }
                    _dummy.convert_units( *ax1.get_units() );
                    _z.convert_units( *ax3.get_units() );

                    vector_u_t new_x( ax3.size() );
                    new_x.set_units( *_vals.get_units() );

                    for (size_t k = 0; k < ax3.size(); k++) {
                        new_x[ k ] = this->get( ax1[ 0 ], ax2[ 0 ], ax3[ k ] );
                    }
                    _z    = ax3;
                    _vals = new_x;
                    reset_limits( 2 );
                    _init_spline();
                    _copy_to_3d();
                    return *this;
                }


            protected:
                void _init_spline() {
                    if (_z && _vals && _spline) {
                        _spline.init( _z.size() ).fill( _z, _vals ).ready();
                    }
                }

                void _copy_to_3d() {
                    _as_3d
                        = vector3d_u_t( 1, 1, _z.size(), _vals.get_units() );
                    _as_3d[ 0 ][ 0 ].assign( _vals.cbegin(), _vals.cend() );
                }

                void _check_limits( double val1, double val2 ) const {
                    if (val2 < _z_limits.first || val2 > _z_limits.second) {
                        std::ostringstream oss;
                        oss << "Requested second dimension value " << val1
                            << " outside dimension 2 limits ["
                            << _z_limits.first << "," << _z_limits.second
                            << "]";
                        throw std::range_error( oss.str() );
                    }
                }

                void _snap_to_limits( double& val3 ) const {
                    if (this->strict()) {
                        if (val3 < _z_limits.first
                            || val3 > _z_limits.second) {
                            std::ostringstream oss;
                            oss << "Requested third dimension value " << val3
                                << " outside dimension 3 limits ["
                                << _z_limits.first << "," << _z_limits.second
                                << "]";
                            throw std::range_error( oss.str() );
                        }
                    }
                    val3 = std::min( std::max( val3, _z.front() ), _z.back() );
                }

                vector_u_t _z;
                vector_u_t _vals;
                std::pair<double, double> _z_limits;
                NCPA::interpolation::Interpolator1D<double, double> _spline;
                vector_u_t _dummy;
                vector3d_u_t _as_3d;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap(
    NCPA::atmos::stratified_atmospheric_property_3d& a,
    NCPA::atmos::stratified_atmospheric_property_3d& b ) noexcept {
    using std::swap;
    ::swap(
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_3d&>( a ),
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_3d&>( a ) );
    swap( a._z, b._z );
    swap( a._vals, b._vals );
    swap( a._spline, b._spline );
    swap( a._z_limits, b._z_limits );
    swap( a._dummy, b._dummy );
    swap( a._as_3d, b._as_3d );
}
