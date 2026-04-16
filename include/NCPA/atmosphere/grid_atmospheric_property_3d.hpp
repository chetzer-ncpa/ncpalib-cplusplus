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

static void swap( NCPA::atmos::grid_atmospheric_property_3d& a,
                  NCPA::atmos::grid_atmospheric_property_3d& b ) noexcept;

namespace NCPA {
    namespace atmos {
        class grid_atmospheric_property_3d
            : public abstract_atmospheric_property_3d {
            public:
                using abstract_atmospheric_property_3d::get;
                using abstract_atmospheric_property_3d::resample;

                grid_atmospheric_property_3d() {
                    set_interpolator(
                        NCPA::interpolation::interpolator_3d_type_t::
                            LANL_HYBRID );
                }

                grid_atmospheric_property_3d( const vector_u_t& ax1,
                                              const vector_u_t& ax2,
                                              const vector_u_t& ax3,
                                              const vector3d_u_t& p ) :
                    grid_atmospheric_property_3d() {
                    if (p.dim( 0 ) != ax1.size() || p.dim( 1 ) != ax2.size()
                        || p.dim( 2 ) != ax3.size()) {
                        throw std::range_error( "grid_atmospheric_property_3d("
                                                "): Values dimensions do not "
                                                "match axes dimensions!" );
                    }
                    this->set( ax1, ax2, ax3, p );
                }

                grid_atmospheric_property_3d(
                    const grid_atmospheric_property_3d& source ) :
                    grid_atmospheric_property_3d() {
                    _axes        = source._axes;
                    _vals        = source._vals;
                    _axis_limits = source._axis_limits;
                    set_interpolator( source._spline.interptype() );
                    this->_init_spline();
                }

                grid_atmospheric_property_3d(
                    grid_atmospheric_property_3d&& source ) noexcept :
                    grid_atmospheric_property_3d() {
                    ::swap( *this, source );
                }

                virtual ~grid_atmospheric_property_3d() {}

                grid_atmospheric_property_3d& operator=(
                    grid_atmospheric_property_3d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( grid_atmospheric_property_3d& a,
                                    grid_atmospheric_property_3d& b ) noexcept;

                virtual std::unique_ptr<abstract_atmospheric_property> clone()
                    const override {
                    return std::unique_ptr<abstract_atmospheric_property>(
                        new grid_atmospheric_property_3d( *this ) );
                }

                virtual grid_atmospheric_property_3d& copy(
                    const abstract_atmospheric_property& source ) override {
                    vector_u_t dummy;
                    vector3d_u_t strat;
                    const abstract_atmospheric_property_3d *source3d;
                    const abstract_atmospheric_property_1d *source1d;
                    switch (source.dimensions()) {
                        case 3:
                            source3d = static_cast<
                                const abstract_atmospheric_property_3d *>(
                                &source );
                            this->set(
                                source3d->axis( 0 ), source3d->axis( 1 ),
                                source3d->axis( 2 ), source3d->values() );
                            break;
                        case 1:
                            source1d = static_cast<
                                const abstract_atmospheric_property_1d *>(
                                &source );
                            dummy = vector_u_t( 1, &NCPA::units::KILOMETERS );
                            strat
                                = vector3d_u_t( 1, 1, source1d->axis().size(),
                                                source1d->get_units() );
                            strat[ 0 ][ 0 ] = source1d->values();
                            this->set( dummy, dummy, source1d->axis(), strat );
                            break;
                        default:
                            throw std::range_error(
                                "grid_atmospheric_property_3d.copy(): "
                                "Supplied property must be 1-D or 3-D" );
                    }
                    return *this;
                }

                virtual grid_atmospheric_property_3d& clear() override {
                    _axes[ 0 ].clear();
                    _axes[ 1 ].clear();
                    _axes[ 2 ].clear();
                    _vals.clear();
                    _axis_limits[ 0 ] = std::pair<double, double>( { 0, 0 } );
                    _axis_limits[ 1 ] = std::pair<double, double>( { 0, 0 } );
                    _axis_limits[ 2 ] = std::pair<double, double>( { 0, 0 } );
                    _spline.clear();
                    return *this;
                }

                virtual size_t dim( size_t n ) const override {
                    return _vals.dim( n );
                }

                virtual grid_atmospheric_property_3d& set(
                    size_t ax, const vector_u_t& axvals ) override {
                    if (axvals.size() != _vals.dim( ax )) {
                        std::ostringstream oss;
                        oss << "grid_atmospheric_property_3d.set(): size "
                               "mismatch between provided vector"
                               " of size "
                            << axvals.size() << " and existing axis " << ax
                            << " of size " << _vals.dim( ax );
                        throw std::range_error( oss.str() );
                    }
                    _axes[ ax ] = axvals;
                    reset_limits( ax );
                    _init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_3d& set(
                    const vector3d_u_t& ax ) override {
                    if (ax.shape() != _vals.shape()) {
                        std::ostringstream oss;
                        oss << "grid_atmospheric_property_3d.set(): Size "
                               "mismatch between existing property size "
                            << _vals.dim( 0 ) << "x" << _vals.dim( 1 ) << "x"
                            << _vals.dim( 2 ) << " and supplied property size "
                            << ax.dim( 0 ) << "x" << ax.dim( 1 ) << "x"
                            << ax.dim( 2 );
                        throw std::range_error( oss.str() );
                    }
                    _vals = ax;
                    _init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_3d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector_u_t& ax3,
                    const vector3d_u_t& vals ) override {
                    this->clear();
                    _axes[ 0 ] = ax1;
                    _axes[ 1 ] = ax2;
                    _axes[ 2 ] = ax3;
                    _vals      = vals;
                    reset_limits( 0 );
                    reset_limits( 1 );
                    reset_limits( 2 );
                    _init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_3d& set(
                    size_t nx1, size_t nx2, size_t nx3,
                    const units_ptr_t units ) override {
                    this->clear();
                    _axes[ 0 ] = vector_u_t( nx1 );
                    _axes[ 1 ] = vector_u_t( nx2 );
                    _axes[ 2 ] = vector_u_t( nx3 );
                    _vals      = vector3d_u_t( nx1, nx2, nx3, units );
                    reset_limits( 0 );
                    reset_limits( 1 );
                    reset_limits( 2 );
                    _init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_3d& append(
                    size_t dim, double dimval,
                    const vector2d_u_t& newslice ) override {
                    std::ostringstream oss;
                    size_t oldax1, oldax2;
                    switch (dim) {
                        case 0:
                            oldax1 = 1;
                            oldax2 = 2;
                            break;
                        case 1:
                            oldax1 = 0;
                            oldax2 = 2;
                            break;
                        case 3:
                            oldax1 = 0;
                            oldax2 = 1;
                            break;
                        default:
                            oss << "grid_atmospheric_property_3d.append(): "
                                   "Invalid axis "
                                << dim << ", must be 0, 1, or 2";
                            throw std::range_error( oss.str() );
                    }
                    if (newslice.dim( 0 ) != _axes[ oldax1 ].size()
                        || newslice.dim( 1 ) != _axes[ oldax2 ].size()) {
                        oss << "Size mismatch between dimensions of new slice "
                            << newslice.dim( 0 ) << "x" << newslice.dim( 1 )
                            << " and existing slices "
                            << _axes[ oldax1 ].size() << "x"
                            << _axes[ oldax2 ].size();
                        throw std::range_error( oss.str() );
                    }
                    if (_axes[ dim ].back() >= dimval) {
                        oss << "grid_atmospheric_property_3d.append(): "
                               "Cannot append, supplied coordinate "
                            << dimval << " for axis " << dim
                            << " is less than current maximum value "
                            << _axes[ dim ].back();
                        throw std::range_error( oss.str() );
                    }
                    vector2d_u_t converted = newslice;
                    converted.convert_units( *_vals.get_units() );

                    auto shape = _vals.shape();
                    shape[ dim ]++;
                    _vals.reshape( shape );
                    auto coords = shape;
                    for (size_t i = 0; i < shape[ oldax1 ]; ++i) {
                        for (size_t j = 0; j < shape[ oldax2 ]; ++j) {
                            coords[ oldax1 ] = i;
                            coords[ oldax2 ] = j;
                            _vals[ coords ]  = converted[ i ][ j ];
                        }
                    }
                    reset_limits( dim );
                    _init_spline();
                    return *this;
                }

                virtual std::unique_ptr<abstract_atmospheric_property_3d>
                    clone3d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_3d>(
                        new grid_atmospheric_property_3d( *this ) );
                }

                virtual std::vector<size_t> shape() const override {
                    return _vals.shape();
                }

                virtual grid_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    this->_init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Grid property requires a 3-D interpolator" );
                }

                virtual grid_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Grid property requires a 3-D interpolator" );
                }

                virtual vector_u_t& axis( size_t n ) override {
                    return _axes[ n ];
                }

                virtual const vector_u_t& axis( size_t n ) const override {
                    return _axes[ n ];
                }

                virtual grid_atmospheric_property_3d& set_limits(
                    size_t dim, double minax, double maxax ) override {
                    if (_axis_limits[ dim ].first < minax
                        || _axis_limits[ dim ].second > maxax) {
                        std::ostringstream oss;
                        oss << "Error setting axis " << dim << " limits to ["
                            << minax << "," << maxax
                            << "]: cannot contract existing axes.";
                        throw std::range_error( oss.str() );
                    }
                    _axis_limits[ dim ].first  = minax;
                    _axis_limits[ dim ].second = maxax;
                    return *this;
                }

                virtual grid_atmospheric_property_3d& reset_limits(
                    size_t dim ) override {
                    _axis_limits[ dim ].first  = _axes[ dim ].front();
                    _axis_limits[ dim ].second = _axes[ dim ].back();
                    return *this;
                }

                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const override {
                    return _axis_limits[ dim ];
                }

                virtual vector3d_u_t& values() override { return _vals; }

                virtual const vector3d_u_t& values() const override {
                    return _vals;
                }

                virtual double get( double val1, double val2,
                                    double val3 ) override {
                    this->_snap_to_limits( val1, val2, val3 );
                    return _spline.eval_f( val1, val2, val3 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     double val3,
                                                     size_t rel ) override {
                    this->_snap_to_limits( val1, val2, val3 );
                    return _spline.eval_df( val1, val2, val3, rel );
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      double val3, size_t rel1,
                                                      size_t rel2 ) override {
                    this->_snap_to_limits( val1, val2, val3 );
                    return _spline.eval_ddf( val1, val2, val3, rel1, rel2 );
                }

                virtual const units_ptr_t get_units() const override {
                    return _vals.get_units();
                }

                virtual const units_ptr_t get_axis_units(
                    size_t n ) const override {
                    return _axes[ n ].get_units();
                }

                virtual grid_atmospheric_property_3d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _vals.convert_units( u );
                    _init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_3d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) override {
                    _axes[ n ].convert_units( u );
                    this->reset_limits( n );
                    _init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_3d& resample(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector_u_t& ax3 ) override {
                    if (!_spline) {
                        throw std::logic_error(
                            "Cannot resample, interpolator has not been "
                            "set." );
                    }
                    _axes[ 0 ].convert_units( *ax1.get_units() );
                    _axes[ 1 ].convert_units( *ax2.get_units() );
                    _axes[ 2 ].convert_units( *ax3.get_units() );

                    vector3d_u_t new_x( ax1.size(), ax2.size(), ax3.size() );
                    new_x.set_units( *_vals.get_units() );
                    for (size_t i = 0; i < ax1.size(); i++) {
                        for (size_t j = 0; j < ax2.size(); j++) {
                            for (size_t k = 0; k < ax3.size(); k++) {
                                new_x[ i ][ j ][ k ] = this->get(
                                    ax1[ i ], ax2[ j ], ax3[ k ] );
                            }
                        }
                    }
                    this->set( ax1, ax2, ax3, new_x );
                    _init_spline();
                    return *this;
                }

            protected:
                void _snap_to_limits( double& val1, double& val2,
                                      double& val3 ) {
                    if (this->strict()) {
                        if (val1 < _axis_limits[ 0 ].first
                            || val1 > _axis_limits[ 0 ].second) {
                            std::ostringstream oss;
                            oss << "Requested first dimension value " << val1
                                << " outside dimension 1 limits ["
                                << _axis_limits[ 0 ].first << ","
                                << _axis_limits[ 0 ].second << "]";
                            throw std::range_error( oss.str() );
                        }
                        if (val2 < _axis_limits[ 1 ].first
                            || val2 > _axis_limits[ 1 ].second) {
                            std::ostringstream oss;
                            oss << "Requested second dimension value " << val2
                                << " outside dimension 2 limits ["
                                << _axis_limits[ 1 ].first << ","
                                << _axis_limits[ 1 ].second << "]";
                            throw std::range_error( oss.str() );
                        }
                        if (val3 < _axis_limits[ 2 ].first
                            || val3 > _axis_limits[ 2 ].second) {
                            std::ostringstream oss;
                            oss << "Requested third dimension value " << val3
                                << " outside dimension 3 limits ["
                                << _axis_limits[ 2 ].first << ","
                                << _axis_limits[ 2 ].second << "]";
                            throw std::range_error( oss.str() );
                        }
                    }
                    val1 = std::min( std::max( val1, _axes[ 0 ].front() ),
                                     _axes[ 0 ].back() );
                    val2 = std::min( std::max( val2, _axes[ 1 ].front() ),
                                     _axes[ 1 ].back() );
                    val3 = std::min( std::max( val3, _axes[ 2 ].front() ),
                                     _axes[ 2 ].back() );
                }

                void _init_spline() {
                    if (_axes[ 0 ] && _axes[ 1 ] && _axes[ 2 ] && _vals
                        && _spline) {
                        _spline
                            .init( this->dim( 0 ), this->dim( 1 ),
                                   this->dim( 2 ) )
                            .fill( _axes[ 0 ], _axes[ 1 ], _axes[ 2 ], _vals )
                            .ready();
                    }
                }

                std::array<vector_u_t, 3> _axes;
                std::array<std::pair<double, double>, 3> _axis_limits;
                vector3d_u_t _vals;
                NCPA::interpolation::Interpolator3D<double, double>
                    _spline;  // stored as [x][y][z]
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::grid_atmospheric_property_3d& a,
                  NCPA::atmos::grid_atmospheric_property_3d& b ) noexcept {
    using std::swap;
    ::swap(
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_3d&>( a ),
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_3d&>( a ) );
    swap( a._axes, b._axes );
    swap( a._vals, b._vals );
    swap( a._spline, b._spline );
    swap( a._axis_limits, b._axis_limits );
}
