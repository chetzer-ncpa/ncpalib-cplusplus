#pragma once

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/exceptions.hpp"
#include "NCPA/interpolation.hpp"

#include <array>
#include <stdexcept>
#include <unordered_map>

#define RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D \
    return static_cast<abstract_atmospheric_property_3d&>( *this );

namespace NCPA {
    namespace atmos {
        class abstract_atmospheric_property_3d;
        class grid_atmospheric_property_3d;
        class stratified_atmospheric_property_3d;
        // class cylindrical_atmospheric_property_3d;
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract_atmospheric_property_3d& a,
                  NCPA::atmos::abstract_atmospheric_property_3d& b ) noexcept;
static void swap( NCPA::atmos::grid_atmospheric_property_3d& a,
                  NCPA::atmos::grid_atmospheric_property_3d& b ) noexcept;
static void swap(
    NCPA::atmos::stratified_atmospheric_property_3d& a,
    NCPA::atmos::stratified_atmospheric_property_3d& b ) noexcept;
// static void swap(
//     NCPA::atmos::cylindrical_atmospheric_property_3d& a,
//     NCPA::atmos::cylindrical_atmospheric_property_3d& b ) noexcept;
static void swap( NCPA::atmos::AtmosphericProperty3D& a,
                  NCPA::atmos::AtmosphericProperty3D& b ) noexcept;

namespace NCPA {
    namespace atmos {

        class abstract_atmospheric_property_3d
            : public abstract_atmospheric_property {
            public:
                virtual ~abstract_atmospheric_property_3d() {}

                friend void ::swap(
                    abstract_atmospheric_property_3d& a,
                    abstract_atmospheric_property_3d& b ) noexcept;

                virtual abstract_atmospheric_property_3d& append(
                    size_t dim, double dimval, const vector2d_u_t& newslice )
                    = 0;
                virtual vector_u_t& axis( size_t n )              = 0;
                virtual const vector_u_t& axis( size_t n ) const  = 0;
                virtual abstract_atmospheric_property_3d& clear() = 0;
                virtual std::unique_ptr<abstract_atmospheric_property_3d>
                    clone3d() const = 0;
                virtual abstract_atmospheric_property_3d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u )
                    = 0;
                virtual abstract_atmospheric_property_3d& convert_units(
                    const NCPA::units::Unit& u )
                    = 0;
                virtual size_t dim( size_t n ) const = 0;
                virtual double get( double val1, double val2, double val3 )
                    = 0;
                virtual const units_ptr_t get_axis_units( size_t n ) const = 0;
                virtual double get_first_derivative( double val1, double val2,
                                                     double val3, size_t rel )
                    = 0;
                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const
                    = 0;
                virtual double get_second_derivative( double val1, double val2,
                                                      double val3, size_t rel1,
                                                      size_t rel2 )
                    = 0;
                virtual const units_ptr_t get_units() const = 0;
                virtual abstract_atmospheric_property_3d& resample(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector_u_t& ax3 )
                    = 0;
                virtual abstract_atmospheric_property_3d& reset_limits(
                    size_t dim )
                    = 0;
                virtual abstract_atmospheric_property_3d& set(
                    size_t ax, const vector_u_t& axvals )
                    = 0;
                virtual abstract_atmospheric_property_3d& set(
                    const vector3d_u_t& ax )
                    = 0;
                virtual abstract_atmospheric_property_3d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector_u_t& ax3, const vector3d_u_t& vals )
                    = 0;
                virtual abstract_atmospheric_property_3d& set(
                    size_t nx1, size_t nx2, size_t nx3,
                    const units_ptr_t units )
                    = 0;
                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type )
                    = 0;
                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    = 0;
                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;
                virtual abstract_atmospheric_property_3d& set_limits(
                    size_t dim, double minax, double maxax )
                    = 0;
                virtual std::vector<size_t> shape() const  = 0;
                virtual vector3d_u_t& values()             = 0;
                virtual const vector3d_u_t& values() const = 0;

                // convenience defines and fixed-value methods
                virtual double ddf( double val1, double val2, double val3,
                                    size_t rel1, size_t rel2 ) {
                    return this->get_second_derivative( val1, val2, val3, rel1,
                                                        rel2 );
                }

                virtual double df( double val1, double val2, double val3,
                                   size_t rel ) {
                    return this->get_first_derivative( val1, val2, val3, rel );
                }

                virtual size_t dimensions() const override { return 3; }

                virtual double f( double val1, double val2, double val3 ) {
                    return this->get( val1, val2, val3 );
                }

                virtual abstract_atmospheric_property_3d& resample(
                    size_t dim, const vector_u_t& new_ax ) {
                    switch (dim) {
                        case 0:
                            return this->resample( new_ax, this->axis( 1 ),
                                                   this->axis( 2 ) );
                            break;
                        case 1:
                            return this->resample( this->axis( 0 ), new_ax,
                                                   this->axis( 2 ) );
                            break;
                        case 2:
                            return this->resample( this->axis( 0 ),
                                                   this->axis( 1 ), new_ax );
                            break;
                        default:
                            std::ostringstream oss;
                            oss << "abstract_atmopheric_property_3d.resample()"
                                   ": Invalid axis "
                                << dim;
                            throw std::range_error( oss.str() );
                    }
                }
        };

        typedef std::unique_ptr<abstract_atmospheric_property_3d>
            _atm_prop_3d_ptr_t;

        class grid_atmospheric_property_3d
            : public abstract_atmospheric_property_3d {
            public:
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

                DECLARE_BOILERPLATE_METHODS( grid_atmospheric_property_3d,
                                             abstract_atmospheric_property )

                virtual abstract_atmospheric_property& copy(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY;
                }

                virtual abstract_atmospheric_property_3d& clear() override {
                    _axes[ 0 ].clear();
                    _axes[ 1 ].clear();
                    _axes[ 2 ].clear();
                    _vals.clear();
                    _axis_limits[ 0 ] = std::pair<double, double>( { 0, 0 } );
                    _axis_limits[ 1 ] = std::pair<double, double>( { 0, 0 } );
                    _axis_limits[ 2 ] = std::pair<double, double>( { 0, 0 } );
                    _spline.clear();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual size_t dim( size_t n ) const override {
                    return _vals.dim( n );
                }

                virtual abstract_atmospheric_property_3d& set(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& set(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& set(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& set(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& append(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual std::unique_ptr<abstract_atmospheric_property_3d>
                    clone3d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_3d>(
                        new grid_atmospheric_property_3d( *this ) );
                }

                // DECLARE_BOILERPLATE_METHODS( grid_atmospheric_property_3d,
                //                              abstract_atmospheric_property )

                virtual std::vector<size_t> shape() const override {
                    return _vals.shape();
                }

                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    this->_init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Grid property requires a 3-D interpolator" );
                }

                virtual abstract_atmospheric_property_3d& set_interpolator(
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

                virtual abstract_atmospheric_property_3d& set_limits(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& reset_limits(
                    size_t dim ) override {
                    _axis_limits[ dim ].first  = _axes[ dim ].front();
                    _axis_limits[ dim ].second = _axes[ dim ].back();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
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

                virtual abstract_atmospheric_property_3d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _vals.convert_units( u );
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D
                }

                virtual abstract_atmospheric_property_3d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) override {
                    _axes[ n ].convert_units( u );
                    this->reset_limits( n );
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D
                }

                virtual abstract_atmospheric_property_3d& resample(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D
                }

            protected:
                void _snap_to_limits( double& val1, double& val2,
                                      double& val3 ) {
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

        // class cylindrical_atmospheric_property_3d
        //     : public grid_atmospheric_property_3d {
        //     public:
        //         // using abstract_atmospheric_property_3d::resample;

        //         cylindrical_atmospheric_property_3d() :
        //             grid_atmospheric_property_3d() {}

        //         cylindrical_atmospheric_property_3d( const vector_u_t& ax1,
        //                                              const vector_u_t& ax2,
        //                                              const vector_u_t& ax3,
        //                                              const vector3d_u_t& p )
        //                                              :
        //             grid_atmospheric_property_3d( ax1, ax2, ax3, p ) {}

        //         cylindrical_atmospheric_property_3d(
        //             const cylindrical_atmospheric_property_3d& source ) :
        //             grid_atmospheric_property_3d( source ) {}

        //         DECLARE_BOILERPLATE_METHODS(
        //             cylindrical_atmospheric_property_3d,
        //             abstract_atmospheric_property )

        //         virtual abstract_atmospheric_property& copy(
        //             const abstract_atmospheric_property& source ) override {
        //             vector_u_t dummy;
        //             vector3d_u_t strat;
        //             switch ( source.dimensions() ) {
        //                 case 3:
        //                     const abstract_atmospheric_property_3d *source3d
        //                         = static_cast<
        //                             const abstract_atmospheric_property_3d
        //                             *>( &source );
        //                     this->set(
        //                         source3d->axis( 0 ), source3d->axis( 1 ),
        //                         source3d->axis( 2 ), source3d->values() );
        //                     break;
        //                 case 1:
        //                     const abstract_atmospheric_property_1d *source1d
        //                         = static_cast<
        //                             const abstract_atmospheric_property_1d
        //                             *>( &source );
        //                     dummy = vector_u_t( 1, &NCPA::units::KILOMETERS
        //                     ); strat
        //                         = vector3d_u_t( 1, 1,
        //                         source1d->axis().size(),
        //                                         source1d->get_units() );
        //                     strat[ 0 ][ 0 ] = source1d->values();
        //                     this->set( dummy, dummy, source1d->axis(), strat
        //                     ); break;
        //                 default:
        //                     throw std::range_error(
        //                         "cylindrical_atmospheric_property_3d.copy():
        //                         " "Supplied property must be 1-D or 3-D" );
        //             }
        //             RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY;
        //         }

        //         virtual double get_first_derivative( double val1, double
        //         val2,
        //                                              double val3,
        //                                              size_t rel ) override {
        //             this->_snap_to_limits( val1, val2, val3 );
        //             double factor = 1.0;
        //             if (rel == 1) {
        //                 factor = val1;
        //             }
        //             return factor * this->_spline.eval_df( val1, val2, val3,
        //             rel );
        //         }

        //         virtual double get_second_derivative( double val1, double
        //         val2,
        //                                               double val3, size_t
        //                                               rel1, size_t rel2 )
        //                                               override {
        //             this->_snap_to_limits( val1, val2, val3 );
        //             double factor = 1.0;
        //             if (rel1 == 1) {
        //                 factor *= val1;
        //             }
        //             if (rel2 == 1) {
        //                 factor *= val1;
        //             }
        //             return factor * this->_spline.eval_ddf( val1, val2,
        //             val3, rel1, rel2 );
        //         }
        // };

        class stratified_atmospheric_property_3d
            : public abstract_atmospheric_property_3d {
            public:
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

                DECLARE_BOILERPLATE_METHODS(
                    stratified_atmospheric_property_3d,
                    abstract_atmospheric_property )

                virtual abstract_atmospheric_property_3d& clear() override {
                    _z.clear();
                    _vals.clear();
                    _z_limits = std::pair<double, double>( { 0, 0 } );
                    _spline.clear();
                    _as_3d.clear();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property& copy(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY;
                }

                virtual size_t dim( size_t n ) const override {
                    _validate_axis( n );
                    return ( n == 2 ? _z.size() : 1 );
                }

                virtual abstract_atmospheric_property_3d& set(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& set(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector_u_t& ax3,
                    const vector3d_u_t& vals ) override {
                    this->clear();
                    _z    = ax3;
                    _vals = vector_u_t( vals[ 0 ][ 0 ], vals.get_units() );
                    reset_limits( 2 );
                    _init_spline();
                    _copy_to_3d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& set(
                    size_t nx1, size_t nx2, size_t nx3,
                    const units_ptr_t units ) override {
                    this->clear();
                    _z    = vector_u_t( nx3 );
                    _vals = vector_u_t( nx3, units );
                    _init_spline();
                    _copy_to_3d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& append(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual std::unique_ptr<abstract_atmospheric_property_3d>
                    clone3d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_3d>(
                        new stratified_atmospheric_property_3d( *this ) );
                }

                // DECLARE_BOILERPLATE_METHODS( grid_atmospheric_property_3d,
                //                              abstract_atmospheric_property )

                virtual std::vector<size_t> shape() const override {
                    return std::vector<size_t> { 1, 1, _z.size() };
                }

                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Stratified property requires a 1-D interpolator" );
                }

                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Stratified property requires a 1-D interpolator" );
                }

                virtual abstract_atmospheric_property_3d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    this->_init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual vector_u_t& axis( size_t n ) override {
                    this->_validate_axis( n );
                    return ( n == 2 ? _z : _dummy );
                }

                virtual const vector_u_t& axis( size_t n ) const override {
                    this->_validate_axis( n );
                    return ( n == 2 ? _z : _dummy );
                }

                virtual abstract_atmospheric_property_3d& set_limits(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
                }

                virtual abstract_atmospheric_property_3d& reset_limits(
                    size_t dim ) override {
                    this->_validate_axis( dim );
                    if (dim == 2) {
                        _z_limits.first  = _z.front();
                        _z_limits.second = _z.back();
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D;
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

                virtual abstract_atmospheric_property_3d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _vals.convert_units( u );
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D
                }

                virtual abstract_atmospheric_property_3d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) override {
                    _validate_axis( n );
                    if (n == 2) {
                        _z.convert_units( u );
                        this->reset_limits( 2 );
                        _init_spline();
                    } else {
                        _dummy.convert_units( u );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D
                }

                virtual abstract_atmospheric_property_3d& resample(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_3D
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
                    if (val3 < _z_limits.first || val3 > _z_limits.second) {
                        std::ostringstream oss;
                        oss << "Requested third dimension value " << val3
                            << " outside dimension 3 limits ["
                            << _z_limits.first << "," << _z_limits.second
                            << "]";
                        throw std::range_error( oss.str() );
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

static void swap( NCPA::atmos::abstract_atmospheric_property_3d& a,
                  NCPA::atmos::abstract_atmospheric_property_3d& b ) noexcept {
}

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

// static void swap(
//     NCPA::atmos::cylindrical_atmospheric_property_3d& a,
//     NCPA::atmos::cylindrical_atmospheric_property_3d& b ) noexcept {
//     using std::swap;
//     ::swap( dynamic_cast<NCPA::atmos::grid_atmospheric_property_3d&>( a ),
//             dynamic_cast<NCPA::atmos::grid_atmospheric_property_3d&>( a ) );
// }

namespace NCPA {
    namespace atmos {

        class AtmosphericProperty3D {
            public:
                AtmosphericProperty3D() {}

                AtmosphericProperty3D( _atm_prop_3d_ptr_t engine ) {
                    _ptr = std::move( engine->clone3d() );
                }

                AtmosphericProperty3D( const AtmosphericProperty3D& other ) :
                    AtmosphericProperty3D() {
                    _ptr = std::move( other._ptr->clone3d() );
                }

                AtmosphericProperty3D(
                    const abstract_atmospheric_property& prop ) {
                    if (prop.dimensions() == 3) {
                        _ptr = std::move(
                            dynamic_cast<
                                const abstract_atmospheric_property_3d&>(
                                prop )
                                .clone3d() );
                    } else {
                        throw std::invalid_argument(
                            "AtmosphericProperty3D: Can only ingest 3-D "
                            "atmospheric properties" );
                    }
                }

                AtmosphericProperty3D& set(
                    const vector_u_t& ax1, const vector_u_t ax2,
                    const vector_u_t ax3,
                    const NCPA::arrays::ndvector<
                        2, const AtmosphericProperty1D *>& props ) {
                    _check_pointer();
                    if (props.shape()[ 0 ] != ax1.size()
                        || props.shape()[ 1 ] != ax2.size()) {
                        throw std::range_error(
                            "AtmosphericProperty3D.set(): Size mismatch "
                            "between supplied axes and supplied 2-D vector of "
                            "properties" );
                    }
                    vector3d_u_t grid( ax1.size(), ax2.size(), ax3.size(),
                                       props[ 0 ][ 0 ]->get_units() );
                    for (size_t i = 0; i < ax1.size(); ++i) {
                        for (size_t j = 0; j < ax2.size(); ++j) {
                            AtmosphericProperty1D prop( *props[ i ][ j ] );
                            prop.convert_units( *grid.get_units() );
                            prop.convert_axis_units( *ax3.get_units() );
                            // prop.resample( ax3 );
                            for (size_t k = 0; k < ax3.size(); ++k) {
                                grid[ i ][ j ][ k ] = prop.get( ax3[ k ] );
                            }
                        }
                    }
                    _ptr->set( ax1, ax2, ax3, grid );
                    return *this;
                }

                AtmosphericProperty3D& append(
                    size_t axis, const scalar_u_t& range,
                    const AtmosphericProperty2D& prop ) {
                    _check_pointer();
                    _ptr->append( axis,
                                  range.get_as( _ptr->get_axis_units( axis ) ),
                                  prop.internal()->values() );
                    return *this;
                }

                DECLARE_WRAPPER_BOILERPLATE_METHODS( AtmosphericProperty3D )

                virtual size_t size( size_t d ) const {
                    _check_pointer();
                    return _ptr->dim( d );
                }

                virtual AtmosphericProperty3D& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type ) {
                    _check_pointer();
                    _ptr->set_interpolator( interp_type );
                    return *this;
                }

                virtual AtmosphericProperty3D& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type ) {
                    _check_pointer();
                    _ptr->set_interpolator( interp_type );
                    return *this;
                }

                virtual AtmosphericProperty3D& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type ) {
                    _check_pointer();
                    _ptr->set_interpolator( interp_type );
                    return *this;
                }

                virtual AtmosphericProperty3D& set( size_t dim,
                                                    const vector_u_t& ax ) {
                    _check_pointer();
                    _ptr->set( dim, ax );
                    return *this;
                }

                virtual AtmosphericProperty3D& set( const vector3d_u_t& ax ) {
                    _check_pointer();
                    _ptr->set( ax );
                    return *this;
                }

                virtual AtmosphericProperty3D& set( const vector_u_t& ax1,
                                                    const vector_u_t& ax2,
                                                    const vector_u_t& ax3,
                                                    const vector3d_u_t vals ) {
                    _check_pointer();
                    _ptr->set( ax1, ax2, ax3, vals );
                    return *this;
                }

                virtual vector_u_t& axis( size_t n ) {
                    _check_pointer();
                    return _ptr->axis( n );
                }

                virtual const vector_u_t& axis( size_t n ) const {
                    _check_pointer();
                    return _ptr->axis( n );
                }

                virtual vector3d_u_t& values() {
                    _check_pointer();
                    return _ptr->values();
                }

                virtual const vector3d_u_t& values() const {
                    _check_pointer();
                    return _ptr->values();
                }

                virtual double get( double val1, double val2, double val3 ) {
                    _check_pointer();
                    return _ptr->get( val1, val2, val3 );
                }

                virtual double f( double val1, double val2, double val3 ) {
                    _check_pointer();
                    return _ptr->f( val1, val2, val3 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     double val3,
                                                     size_t rel ) {
                    _check_pointer();
                    return _ptr->get_first_derivative( val1, val2, val3, rel );
                }

                virtual double df( double val1, double val2, double val3,
                                   size_t rel ) {
                    _check_pointer();
                    return _ptr->df( val1, val2, val3, rel );
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      double val3, size_t rel1,
                                                      size_t rel2 ) {
                    _check_pointer();
                    return _ptr->get_second_derivative( val1, val2, val3, rel1,
                                                        rel2 );
                }

                virtual double ddf( double val1, double val2, double val3,
                                    size_t rel1, size_t rel2 ) {
                    _check_pointer();
                    return _ptr->ddf( val1, val2, val3, rel1, rel2 );
                }

                virtual const units_ptr_t get_units() const {
                    _check_pointer();
                    return _ptr->get_units();
                }

                virtual const units_ptr_t get_axis_units( size_t n ) const {
                    _check_pointer();
                    return _ptr->get_axis_units( n );
                }

                virtual AtmosphericProperty3D& convert_units(
                    const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_units( u );
                    return *this;
                }

                virtual AtmosphericProperty3D& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_axis_units( n, u );
                    return *this;
                }

                virtual AtmosphericProperty3D& resample(
                    const vector_u_t& new_ax1, const vector_u_t& new_ax2,
                    const vector_u_t& new_ax3 ) {
                    _check_pointer();
                    _ptr->resample( new_ax1, new_ax2, new_ax3 );
                    return *this;
                }

                virtual AtmosphericProperty3D& resample(
                    size_t axis, const vector_u_t& new_ax ) {
                    _check_pointer();
                    _ptr->resample( axis, new_ax );
                    return *this;
                }

                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const {
                    _check_pointer();
                    return _ptr->get_limits( dim );
                }

                virtual AtmosphericProperty3D& set_limits( size_t dim,
                                                           double minax,
                                                           double maxax ) {
                    _check_pointer();
                    _ptr->set_limits( dim, minax, maxax );
                    return *this;
                }

                virtual _atm_prop_3d_ptr_t& internal() { return _ptr; }

                virtual const _atm_prop_3d_ptr_t& internal() const {
                    return _ptr;
                }


            protected:
                _atm_prop_3d_ptr_t _ptr;

                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphericProperty3D: no pointer has been set" );
                    }
                }
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphericProperty3D& a,
                  NCPA::atmos::AtmosphericProperty3D& b ) noexcept {
    using std::swap;
    swap( a._ptr, b._ptr );
}
