#pragma once

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/abstract_atmospheric_property_2d.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/tuple_atmospheric_property_1d.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/exceptions.hpp"
#include "NCPA/interpolation.hpp"

#include <array>
#include <stdexcept>
#include <unordered_map>

static void swap( NCPA::atmos::grid_atmospheric_property_2d& a,
                  NCPA::atmos::grid_atmospheric_property_2d& b ) noexcept;

namespace NCPA {
    namespace atmos {
        class grid_atmospheric_property_2d
            : public abstract_atmospheric_property_2d {
            public:
                using abstract_atmospheric_property_2d::resample;

                grid_atmospheric_property_2d() {
                    set_interpolator(
                        NCPA::interpolation::interpolator_2d_type_t::
                            NEAREST_NEIGHBOR );
                }

                grid_atmospheric_property_2d( const vector_u_t& r,
                                              const vector_u_t& z,
                                              const vector2d_u_t& p ) :
                    grid_atmospheric_property_2d() {
                    if (p.dim( 0 ) != r.size() || p.dim( 1 ) != z.size()) {
                        throw std::range_error( "Values dimensions do not "
                                                "match axes dimensions!" );
                    }
                    this->set( 0, r );
                    this->set( 1, z );
                    this->set( p );
                    this->_init_spline();
                }

                grid_atmospheric_property_2d(
                    const grid_atmospheric_property_2d& source ) :
                    grid_atmospheric_property_2d() {
                    _axes        = source._axes;
                    _vals        = source._vals;
                    _axis_limits = source._axis_limits;
                    set_interpolator( source._spline.interptype() );
                    _init_spline();
                }

                virtual grid_atmospheric_property_2d& clear() override {
                    _axes[ 0 ].clear();
                    _axes[ 1 ].clear();
                    _vals.clear();
                    _axis_limits[ 0 ] = std::pair<double, double>( { 0, 0 } );
                    _axis_limits[ 1 ] = std::pair<double, double>( { 0, 0 } );
                    _spline.clear();
                    return *this;
                    ;
                }

                virtual grid_atmospheric_property_2d& set(
                    const vector_u_t& ranges,
                    const std::vector<const abstract_atmospheric_property_1d
                                          *>& atmos1ds ) override {
                    if (ranges.size() != atmos1ds.size()) {
                        throw std::invalid_argument(
                            "Range and component vectors must be the same "
                            "size!" );
                    }
                    this->clear();
                    _axes[ 0 ] = ranges;
                    _axes[ 1 ] = atmos1ds[ 0 ]->axis();
                    _vals = vector2d_u_t( _axes[ 0 ].size(), _axes[ 1 ].size(),
                                          atmos1ds[ 0 ]->get_units() );
                    for (size_t i = 0; i < ranges.size(); ++i) {
                        auto clone = atmos1ds[ i ]->clone1d();
                        abstract_atmospheric_property_1d *temp = clone.get();
                        temp->resample( _axes[ 1 ] );
                        _vals[ i ] = temp->values();
                    }
                    this->reset_limits( 0 );
                    this->reset_limits( 1 );
                    this->_init_spline();
                    return *this;
                    ;
                }

                virtual grid_atmospheric_property_2d& set(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    return this->set(
                        vector_u_t( 1, range ),
                        std::vector<const abstract_atmospheric_property_1d *>(
                            1, &atmos1d ) );
                }

                virtual grid_atmospheric_property_2d& set(
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    return this->set(
                        scalar_u_t( 0.0, _axes[ 0 ].get_units() ), atmos1d );
                }

                virtual grid_atmospheric_property_2d& append(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    _axes[ 0 ].push_back(
                        range.get_as( _axes[ 0 ].get_units() ) );
                    auto clone                             = atmos1d.clone1d();
                    abstract_atmospheric_property_1d *temp = clone.get();
                    temp->resample( _axes[ 1 ] );
                    auto shape = _vals.shape();
                    shape[ 0 ]++;
                    _vals.reshape( shape );
                    _vals[ shape[ 0 ] - 1 ] = temp->values();
                    this->reset_limits( 0 );
                    this->_init_spline();
                    return *this;
                    ;
                }

                grid_atmospheric_property_2d(
                    grid_atmospheric_property_2d&& source ) noexcept :
                    grid_atmospheric_property_2d() {
                    ::swap( *this, source );
                }

                virtual ~grid_atmospheric_property_2d() {}

                grid_atmospheric_property_2d& operator=(
                    grid_atmospheric_property_2d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( grid_atmospheric_property_2d& a,
                                    grid_atmospheric_property_2d& b ) noexcept;

                virtual std::unique_ptr<abstract_atmospheric_property> clone()
                    const override {
                    return std::unique_ptr<abstract_atmospheric_property>(
                        new grid_atmospheric_property_2d( *this ) );
                }

                virtual grid_atmospheric_property_2d& copy(
                    const abstract_atmospheric_property& source ) override {
                    if (source.dimensions() == 2) {
                        const abstract_atmospheric_property_2d *srcptr
                            = dynamic_cast<
                                const abstract_atmospheric_property_2d *>(
                                &source );
                        this->set( srcptr->axis( 0 ), srcptr->axis( 1 ),
                                   srcptr->values() );
                    } else if (source.dimensions() == 1) {
                        const abstract_atmospheric_property_1d *srcptr
                            = dynamic_cast<
                                const abstract_atmospheric_property_1d *>(
                                &source );

                        const vector_u_t srcvec = srcptr->values();
                        vector2d_u_t values2d( 1, srcvec.size(),
                                               srcvec.get_units() );
                        for (size_t i = 0; i < srcptr->size(); ++i) {
                            values2d.set( 0, i, srcvec[ i ] );
                        }
                        this->set( vector_u_t( std::vector<double>( { 0.0 } ),
                                               &NCPA::units::METERS ),
                                   srcptr->axis(), values2d );
                    } else {
                        throw std::invalid_argument(
                            "grid_atmospheric_property_2d.copy(): Cannot "
                            "copy property with more than 2 dimension" );
                    }
                    return *this;
                    ;
                }

                virtual std::unique_ptr<abstract_atmospheric_property_2d>
                    clone2d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_2d>(
                        new grid_atmospheric_property_2d( *this ) );
                }

                virtual size_t size( size_t d ) const override {
                    return _axes[ d ].size();
                }

                virtual grid_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    this->_init_spline();
                    return *this;
                    ;
                }

                virtual grid_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Grid property requires a 2-D interpolator" );
                }

                virtual grid_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const units_ptr_t units ) override {
                    _axes[ 0 ] = ax1;
                    _axes[ 1 ] = ax2;
                    _vals      = vector2d_u_t( ax1.size(), ax2.size(), units );
                    reset_limits( 0 );
                    reset_limits( 1 );
                    this->_init_spline();
                    return *this;
                    ;
                }

                virtual grid_atmospheric_property_2d& set(
                    size_t dim, const vector_u_t& ax ) override {
                    _axes[ dim ] = ax;
                    if (_vals.dim( dim ) != ax.size()) {
                        auto shape   = _vals.shape();
                        shape[ dim ] = ax.size();
                        _vals.reshape( shape );
                    }
                    this->reset_limits( dim );
                    this->_init_spline();
                    return *this;
                    ;
                };

                virtual grid_atmospheric_property_2d& set(
                    const vector2d_u_t& ax ) override {
                    if (ax.dim( 0 ) != _axes[ 0 ].size()
                        || ax.dim( 1 ) != _axes[ 1 ].size()) {
                        std::ostringstream oss;
                        oss << "2-D vector size " << ax.dim( 0 ) << "x"
                            << ax.dim( 1 )
                            << " does not match existing axes sizes "
                            << _axes[ 0 ].size() << "x" << _axes[ 1 ].size();
                        throw std::logic_error( oss.str() );
                    }
                    _vals = ax;
                    this->_init_spline();
                    return *this;
                    ;
                };

                virtual grid_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector2d_u_t& vals ) override {
                    this->clear();
                    this->set( 0, ax1 );
                    this->set( 1, ax2 );
                    this->set( vals );
                    return *this;
                    ;
                }

                virtual vector_u_t& axis( size_t n ) override {
                    return _axes[ n ];
                }

                virtual const vector_u_t& axis( size_t n ) const override {
                    return _axes[ n ];
                }

                virtual vector2d_u_t& values() override { return _vals; }

                virtual const vector2d_u_t& values() const override {
                    return _vals;
                }

                virtual double get( double val1, double val2 ) override {
                    this->_snap_to_limits( val1, val2 );
                    return _spline.eval_f( val1, val2 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     size_t rel ) override {
                    this->_snap_to_limits( val1, val2 );
                    return _spline.eval_df( val1, val2, rel );
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      size_t rel1,
                                                      size_t rel2 ) override {
                    this->_snap_to_limits( val1, val2 );
                    return _spline.eval_ddf( val1, val2, rel1, rel2 );
                }

                virtual _atm_prop_1d_ptr_t extract( double range ) override {
                    vector_u_t z = _axes[ 1 ];
                    vector_u_t p( z.size(), this->get_units() );
                    for (size_t i = 0; i < z.size(); ++i) {
                        p[ i ] = this->get( range, z[ i ] );
                    }
                    return _atm_prop_1d_ptr_t(
                        new tuple_atmospheric_property_1d( z, p ) );
                }

                virtual const units_ptr_t get_units() const override {
                    return _vals.get_units();
                }

                virtual const units_ptr_t get_axis_units(
                    size_t n ) const override {
                    return _axes[ n ].get_units();
                }

                virtual grid_atmospheric_property_2d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _vals.convert_units( u );
                    _init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_2d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) override {
                    _axes[ n ].convert_units( u );
                    this->reset_limits( n );
                    _init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_2d& resample(
                    const vector_u_t& new_r,
                    const vector_u_t& new_z ) override {
                    if (!_spline) {
                        throw std::logic_error(
                            "Cannot resample, interpolator has not been "
                            "set." );
                    }
                    _axes[ 1 ].convert_units( *new_z.get_units() );
                    _axes[ 0 ].convert_units( *new_r.get_units() );
                    vector2d_u_t new_x( new_r.size(), new_z.size() );
                    new_x.set_units( *_vals.get_units() );
                    for (size_t i = 0; i < new_r.size(); i++) {
                        for (size_t j = 0; j < new_z.size(); j++) {
                            new_x[ i ][ j ]
                                = this->get( new_r[ i ], new_z[ j ] );
                        }
                    }
                    this->set( 0, new_r );
                    this->set( 1, new_z );
                    this->set( new_x );
                    _init_spline();
                    return *this;
                }

                virtual grid_atmospheric_property_2d& set_limits(
                    size_t dim, double minax, double maxax ) override {
                    std::pair<double, double> *dimlimits
                        = &_axis_limits[ dim ];
                    if (dimlimits->first < minax
                        || dimlimits->second > maxax) {
                        std::ostringstream oss;
                        oss << "Error setting axis " << dim << " limits to ["
                            << minax << "," << maxax
                            << "]: cannot contract existing axes.";
                        throw std::range_error( oss.str() );
                    }
                    dimlimits->first  = minax;
                    dimlimits->second = maxax;
                    return *this;
                }

                virtual grid_atmospheric_property_2d& reset_limits(
                    size_t dim ) override {
                    _axis_limits[ dim ].first  = _axes[ dim ].front();
                    _axis_limits[ dim ].second = _axes[ dim ].back();
                    return *this;
                    ;
                }

                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const override {
                    return _axis_limits[ dim ];
                }

            protected:
                void _snap_to_limits( double& val1, double& val2 ) {
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
                        oss << "Requested second dimension value " << val1
                            << " outside dimension 2 limits ["
                            << _axis_limits[ 1 ].first << ","
                            << _axis_limits[ 1 ].second << "]";
                        throw std::range_error( oss.str() );
                    }
                    val1 = std::min( std::max( val1, _axes[ 0 ].front() ),
                                     _axes[ 0 ].back() );
                    val2 = std::min( std::max( val2, _axes[ 1 ].front() ),
                                     _axes[ 1 ].back() );
                }

                void _init_spline() {
                    if (_axes.size() == 2 && _axes[ 0 ] && _axes[ 1 ] && _vals
                        && _spline) {
                        _spline.init( this->size( 0 ), this->size( 1 ) )
                            .fill( _axes[ 0 ], _axes[ 1 ], _vals )
                            .ready();
                    }
                }

                std::array<vector_u_t, 2> _axes;
                std::array<std::pair<double, double>, 2> _axis_limits;
                vector2d_u_t _vals;
                NCPA::interpolation::Interpolator2D<double, double>
                    _spline;  // stored as [r][z]
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::grid_atmospheric_property_2d& a,
                  NCPA::atmos::grid_atmospheric_property_2d& b ) noexcept {
    using std::swap;
    ::swap(
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_2d&>( a ),
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_2d&>( a ) );
    swap( a._axes, b._axes );
    swap( a._vals, b._vals );
    swap( a._spline, b._spline );
    swap( a._axis_limits, b._axis_limits );
}
