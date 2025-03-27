#pragma once

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/exceptions.hpp"
#include "NCPA/interpolation.hpp"

#include <array>
#include <stdexcept>
#include <unordered_map>

#define RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D \
    return static_cast<abstract_atmospheric_property_2d&>( *this );

namespace NCPA {
    namespace atmos {
        class abstract_atmospheric_property_2d;
        class grid_atmospheric_property_2d;
        class stratified_atmospheric_property_2d;
        class piecewise_stratified_atmospheric_property_2d;
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract_atmospheric_property_2d& a,
                  NCPA::atmos::abstract_atmospheric_property_2d& b ) noexcept;
static void swap( NCPA::atmos::grid_atmospheric_property_2d& a,
                  NCPA::atmos::grid_atmospheric_property_2d& b ) noexcept;
static void swap(
    NCPA::atmos::stratified_atmospheric_property_2d& a,
    NCPA::atmos::stratified_atmospheric_property_2d& b ) noexcept;
static void swap(
    NCPA::atmos::piecewise_stratified_atmospheric_property_2d& a,
    NCPA::atmos::piecewise_stratified_atmospheric_property_2d& b ) noexcept;
static void swap( NCPA::atmos::AtmosphericProperty2D& a,
                  NCPA::atmos::AtmosphericProperty2D& b ) noexcept;

namespace NCPA {
    namespace atmos {

        class abstract_atmospheric_property_2d
            : public abstract_atmospheric_property {
            public:
                virtual ~abstract_atmospheric_property_2d() {}

                friend void ::swap(
                    abstract_atmospheric_property_2d& a,
                    abstract_atmospheric_property_2d& b ) noexcept;

                virtual size_t dimensions() const override { return 2; }

                virtual std::unique_ptr<abstract_atmospheric_property_2d>
                    clone2d() const = 0;

                // Generate a 2-D profile from one or more 1-D profiles
                virtual abstract_atmospheric_property_2d& clear() = 0;
                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ranges,
                    const std::vector<
                        const abstract_atmospheric_property_1d *>& atmos1ds )
                    = 0;
                virtual abstract_atmospheric_property_2d& set(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    = 0;
                virtual abstract_atmospheric_property_2d& set(
                    const abstract_atmospheric_property_1d& atmos1d )
                    = 0;
                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const units_ptr_t units )
                    = 0;
                virtual abstract_atmospheric_property_2d& append(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    = 0;
                virtual _atm_prop_1d_ptr_t extract( double range ) = 0;

                virtual size_t size( size_t d ) const = 0;
                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    = 0;
                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;

                virtual abstract_atmospheric_property_2d& set(
                    size_t dim, const vector_u_t& ax )
                    = 0;
                virtual abstract_atmospheric_property_2d& set(
                    const vector2d_u_t& ax )
                    = 0;
                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector2d_u_t& vals )
                    = 0;
                virtual vector_u_t& axis( size_t n )             = 0;
                virtual const vector_u_t& axis( size_t n ) const = 0;
                virtual abstract_atmospheric_property_2d& set_limits(
                    size_t dim, double minax, double maxax )
                    = 0;
                virtual abstract_atmospheric_property_2d& reset_limits(
                    size_t dim )
                    = 0;
                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const
                    = 0;


                virtual vector2d_u_t& values()                 = 0;
                virtual const vector2d_u_t& values() const     = 0;
                virtual double get( double val1, double val2 ) = 0;

                virtual double f( double val1, double val2 ) {
                    return this->get( val1, val2 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     size_t rel )
                    = 0;

                virtual double df( double val1, double val2, size_t rel ) {
                    return this->get_first_derivative( val1, val2, rel );
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      size_t rel1,
                                                      size_t rel2 )
                    = 0;

                virtual double ddf( double val1, double val2, size_t rel1,
                                    size_t rel2 ) {
                    return this->get_second_derivative( val1, val2, rel1,
                                                        rel2 );
                }

                virtual const units_ptr_t get_units() const                = 0;
                virtual const units_ptr_t get_axis_units( size_t n ) const = 0;
                virtual abstract_atmospheric_property_2d& convert_units(
                    const NCPA::units::Unit& u )
                    = 0;
                virtual abstract_atmospheric_property_2d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u )
                    = 0;
                virtual abstract_atmospheric_property_2d& resample(
                    const vector_u_t& new_r, const vector_u_t& new_z )
                    = 0;

                virtual abstract_atmospheric_property_2d& resample(
                    size_t dim, const vector_u_t& new_ax ) {
                    if (dim == 0) {
                        return this->resample( new_ax, this->axis( 1 ) );
                    } else {
                        return this->resample( this->axis( 0 ), new_ax );
                    }
                }
        };

        typedef std::unique_ptr<abstract_atmospheric_property_2d>
            _atm_prop_2d_ptr_t;

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

                virtual abstract_atmospheric_property_2d& clear() override {
                    _axes[ 0 ].clear();
                    _axes[ 1 ].clear();
                    _vals.clear();
                    _axis_limits[ 0 ] = std::pair<double, double>( { 0, 0 } );
                    _axis_limits[ 1 ] = std::pair<double, double>( { 0, 0 } );
                    _spline.clear();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    return this->set(
                        vector_u_t( 1, range ),
                        std::vector<const abstract_atmospheric_property_1d *>(
                            1, &atmos1d ) );
                }

                virtual abstract_atmospheric_property_2d& set(
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    return this->set(
                        scalar_u_t( 0.0, _axes[ 0 ].get_units() ), atmos1d );
                }

                virtual abstract_atmospheric_property_2d& append(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                DECLARE_BOILERPLATE_METHODS( grid_atmospheric_property_2d,
                                             abstract_atmospheric_property )

                virtual abstract_atmospheric_property& copy(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual std::unique_ptr<abstract_atmospheric_property_2d>
                    clone2d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_2d>(
                        new grid_atmospheric_property_2d( *this ) );
                }

                virtual size_t size( size_t d ) const override {
                    return _axes[ d ].size();
                }

                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    this->_init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Grid property requires a 2-D interpolator" );
                }

                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const units_ptr_t units ) override {
                    _axes[ 0 ] = ax1;
                    _axes[ 1 ] = ax2;
                    _vals      = vector2d_u_t( ax1.size(), ax2.size(), units );
                    reset_limits( 0 );
                    reset_limits( 1 );
                    this->_init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    size_t dim, const vector_u_t& ax ) override {
                    _axes[ dim ] = ax;
                    if (_vals.dim( dim ) != ax.size()) {
                        auto shape   = _vals.shape();
                        shape[ dim ] = ax.size();
                        _vals.reshape( shape );
                    }
                    this->reset_limits( dim );
                    this->_init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                };

                virtual abstract_atmospheric_property_2d& set(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                };

                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector2d_u_t& vals ) override {
                    this->clear();
                    this->set( 0, ax1 );
                    this->set( 1, ax2 );
                    this->set( vals );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
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

                virtual abstract_atmospheric_property_2d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _vals.convert_units( u );
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D
                }

                virtual abstract_atmospheric_property_2d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) override {
                    _axes[ n ].convert_units( u );
                    this->reset_limits( n );
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D
                }

                virtual abstract_atmospheric_property_2d& resample(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D
                }

                virtual abstract_atmospheric_property_2d& set_limits(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D
                }

                virtual abstract_atmospheric_property_2d& reset_limits(
                    size_t dim ) override {
                    _axis_limits[ dim ].first  = _axes[ dim ].front();
                    _axis_limits[ dim ].second = _axes[ dim ].back();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
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

        class stratified_atmospheric_property_2d
            : public abstract_atmospheric_property_2d {
            public:
                using abstract_atmospheric_property_2d::resample;

                stratified_atmospheric_property_2d() {
                    _dummy = vector_u_t( { 0 }, NCPA::units::KILOMETERS );
                    set_interpolator(
                        NCPA::interpolation::interpolator_1d_type_t::
                            NEAREST_NEIGHBOR );
                }

                stratified_atmospheric_property_2d( const vector_u_t& z,
                                                    const vector_u_t& p ) :
                    stratified_atmospheric_property_2d() {
                    if (p.size() != z.size()) {
                        throw std::range_error( "Values dimensions do not "
                                                "match axes dimensions!" );
                    }
                    _z    = z;
                    _vals = p;
                    this->_init_spline();
                    _copy_to_2d();
                }

                stratified_atmospheric_property_2d( const vector_u_t& r,
                                                    const vector_u_t& z,
                                                    const vector2d_u_t& p ) :
                    stratified_atmospheric_property_2d() {
                    _z    = z;
                    _vals = vector_u_t( p[ 0 ], p.get_units() );
                    _copy_to_2d();
                }

                stratified_atmospheric_property_2d(
                    const stratified_atmospheric_property_2d& source ) :
                    stratified_atmospheric_property_2d() {
                    _z        = source._z;
                    _vals     = source._vals;
                    _z_limits = source._z_limits;
                    _spline   = NCPA::interpolation::InterpolatorFactory<
                          double, double>::build( source._spline.interptype() );
                    _as_2d = source._as_2d;
                    _dummy = source._dummy;
                }

                virtual abstract_atmospheric_property_2d& clear() override {
                    _z.clear();
                    _vals.clear();
                    _z_limits = std::pair<double, double>( { 0, 0 } );
                    _spline.clear();
                    _as_2d.clear();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ranges,
                    const std::vector<const abstract_atmospheric_property_1d
                                          *>& atmos1ds ) override {
                    if (atmos1ds.size() != 1) {
                        throw std::invalid_argument(
                            "Vector of component properties must be of size "
                            "1!" );
                    }
                    this->clear();
                    _z    = atmos1ds[ 0 ]->axis();
                    _vals = vector_u_t( atmos1ds.front()->values(),
                                        atmos1ds.front()->get_units() );
                    this->_init_spline();
                    _copy_to_2d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    return this->set(
                        vector_u_t( 1, range ),
                        std::vector<const abstract_atmospheric_property_1d *>(
                            1, &atmos1d ) );
                }

                virtual abstract_atmospheric_property_2d& set(
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    return this->set(
                        scalar_u_t( 0.0, NCPA::units::KILOMETERS ), atmos1d );
                }

                virtual abstract_atmospheric_property_2d& append(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    throw std::logic_error( "Cannot append more properties to "
                                            "a stratified 2-D property!" );
                }

                virtual _atm_prop_1d_ptr_t extract( double range ) override {
                    return _atm_prop_1d_ptr_t(
                        new tuple_atmospheric_property_1d( _z, _vals ) );
                }


                DECLARE_BOILERPLATE_METHODS(
                    stratified_atmospheric_property_2d,
                    abstract_atmospheric_property )

                virtual abstract_atmospheric_property& copy(
                    const abstract_atmospheric_property& source ) override {
                    if (source.dimensions() == 2) {
                        const abstract_atmospheric_property_2d *srcptr
                            = dynamic_cast<
                                const abstract_atmospheric_property_2d *>(
                                &source );
                        vector2d_u_t temp = srcptr->values();
                        temp.resize( 1 );
                        vector_u_t dummy_r = srcptr->axis( 0 );
                        dummy_r.resize( 1 );
                        this->set( dummy_r, srcptr->axis( 1 ), temp );
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
                            "stratified_atmospheric_property_2d.copy(): "
                            "Cannot "
                            "copy property with more than 2 dimension" );
                    }
                    _copy_to_2d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual std::unique_ptr<abstract_atmospheric_property_2d>
                    clone2d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_2d>(
                        new stratified_atmospheric_property_2d( *this ) );
                }

                virtual size_t size( size_t d ) const override {
                    return ( d == 0 ? 1 : _z.size() );
                }

                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Stratified property requires a 1-D interpolator" );
                }

                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    this->_init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D
                }

                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const units_ptr_t units ) override {
                    this->clear();
                    _z    = ax2;
                    _vals = vector_u_t( _z.size(), units );
                    reset_limits( 1 );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    size_t dim, const vector_u_t& ax ) override {
                    if (dim == 1) {
                        this->resample(
                            vector_u_t( { 0 }, NCPA::units::KILOMETERS ), ax );
                        this->reset_limits( dim );
                    }
                    _copy_to_2d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D
                };

                virtual abstract_atmospheric_property_2d& set(
                    const vector2d_u_t& ax ) override {
                    if (ax.dim( 1 ) != _z.size()) {
                        std::ostringstream oss;
                        oss << "2-D vector size Nx" << ax.dim( 1 )
                            << " not compatible with existing axes sizes Nx"
                            << _z.size();
                        throw std::logic_error( oss.str() );
                    }
                    _vals.assign( ax[ 0 ].cbegin(), ax[ 0 ].cend() );
                    _copy_to_2d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D
                };

                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector2d_u_t& vals ) override {
                    // this->set( 0, ax1 );
                    this->set( 1, ax2 );
                    this->set( vals );
                    _copy_to_2d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual vector_u_t& axis( size_t n ) override {
                    return ( n == 0 ? _dummy : _z );
                }

                virtual const vector_u_t& axis( size_t n ) const override {
                    return ( n == 0 ? _dummy : _z );
                }

                virtual vector2d_u_t& values() override { return _as_2d; }

                virtual const vector2d_u_t& values() const override {
                    return _as_2d;
                }

                virtual double get( double val1, double val2 ) override {
                    _check_limits( val1, val2 );
                    return _spline.eval_f( val2 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     size_t rel ) override {
                    _check_limits( val1, val2 );
                    if (rel != 0) {
                        return _spline.eval_df( val2 );
                    } else {
                        return 0.0;
                    }
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      size_t rel1,
                                                      size_t rel2 ) override {
                    _check_limits( val1, val2 );
                    if (rel1 != 0 && rel2 != 0) {
                        return _spline.eval_ddf( val2 );
                    } else {
                        return 0.0;
                    }
                }

                virtual const units_ptr_t get_units() const override {
                    return _vals.get_units();
                }

                virtual const units_ptr_t get_axis_units(
                    size_t n ) const override {
                    return _z.get_units();
                }

                virtual abstract_atmospheric_property_2d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _vals.convert_units( u );
                    _init_spline();
                    _copy_to_2d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) override {
                    _z.convert_units( u );
                    this->reset_limits( n );
                    _init_spline();
                    _copy_to_2d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& resample(
                    const vector_u_t& new_r,
                    const vector_u_t& new_z ) override {
                    if (!_spline) {
                        throw std::logic_error(
                            "Cannot resample, interpolator has not been "
                            "set." );
                    }
                    _z.convert_units( *new_z.get_units() );
                    vector_u_t new_vals( new_z.size() );
                    new_vals.set_units( *_vals.get_units() );
                    for (size_t j = 0; j < new_z.size(); j++) {
                        new_vals[ j ] = this->get( 0.0, new_z[ j ] );
                    }
                    _z    = new_z;
                    _vals = new_vals;
                    _init_spline();
                    _copy_to_2d();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set_limits(
                    size_t dim, double minax, double maxax ) override {
                    if (dim == 1) {
                        if (_z_limits.first < minax
                            || _z_limits.second > maxax) {
                            std::ostringstream oss;
                            oss << "Error setting axis " << dim
                                << " limits to [" << minax << "," << maxax
                                << "]: cannot contract existing axes.";
                            throw std::range_error( oss.str() );
                        }
                    }
                    _z_limits.first  = minax;
                    _z_limits.second = maxax;
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& reset_limits(
                    size_t dim ) override {
                    if (dim == 1) {
                        _z_limits.first  = _z.front();
                        _z_limits.second = _z.back();
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const override {
                    return ( dim == 0 ? std::pair<double, double> { 0.0, 0.0 }
                                      : _z_limits );
                }

            protected:
                void _init_spline() {
                    if (_z && _vals && _spline) {
                        _spline.init( this->size( 1 ) )
                            .fill( _z, _vals )
                            .ready();
                    }
                }

                void _copy_to_2d() {
                    _as_2d = vector2d_u_t( 1, _z.size(), _vals.get_units() );
                    _as_2d[ 0 ].assign( _vals.cbegin(), _vals.cend() );
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

                vector_u_t _z;
                vector_u_t _vals;
                std::pair<double, double> _z_limits;
                NCPA::interpolation::Interpolator1D<double, double> _spline;
                vector_u_t _dummy;
                vector2d_u_t _as_2d;
        };

        class piecewise_stratified_atmospheric_property_2d
            : public abstract_atmospheric_property_2d {
            public:
                using abstract_atmospheric_property_2d::resample;

                piecewise_stratified_atmospheric_property_2d() {
                    set_interpolator(
                        NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR );
                }

                piecewise_stratified_atmospheric_property_2d(
                    const vector_u_t& r, const vector_u_t& z,
                    const vector2d_u_t& p ) :
                    piecewise_stratified_atmospheric_property_2d() {
                    if (p.dim( 0 ) != r.size() || p.dim( 1 ) != z.size()) {
                        throw std::range_error( "Values dimensions do not "
                                                "match axes dimensions!" );
                    }
                    _ax1 = r;
                    _props.clear();
                    for (size_t i = 0; i < _ax1.size(); ++i) {
                        vector_u_t propvals( p[ i ], p.get_units() );
                        _props.push_back(
                            AtmosphericProperty1D( _atm_prop_1d_ptr_t(
                                new tuple_atmospheric_property_1d(
                                    z, propvals ) ) ) );
                    }
                    _make_tmpvals();
                    _calculate_change_points();
                }

                piecewise_stratified_atmospheric_property_2d(
                    const piecewise_stratified_atmospheric_property_2d&
                        source ) :
                    piecewise_stratified_atmospheric_property_2d() {
                    _ax1           = source._ax1;
                    _props         = source._props;
                    _change_points = source._change_points;
                    _tmpvals       = source._tmpvals;
                    set_interpolator( source._interptype );
                }

                virtual abstract_atmospheric_property_2d& clear() override {
                    _ax1.clear();
                    _props.clear();
                    _change_points.clear();
                    _tmpvals.clear();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ranges,
                    const std::vector<const abstract_atmospheric_property_1d
                                          *>& atmos1ds ) override {
                    if (ranges.size() != atmos1ds.size()) {
                        throw std::invalid_argument(
                            "Range and component vectors must be the same "
                            "size!" );
                    }
                    this->clear();
                    _ax1 = ranges;
                    for (auto it = atmos1ds.cbegin(); it != atmos1ds.cend();
                         ++it) {
                        _props.push_back( AtmosphericProperty1D( **it ) );
                    }
                    _calculate_change_points();
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    size_t new_ind = 0;
                    double r       = range.get_as( _ax1.get_units() );
                    while (new_ind < _ax1.size() && _ax1[ new_ind++ ] < r) {}
                    _ax1.insert( _ax1.begin() + new_ind, r );
                    _props.insert(
                        _props.begin() + new_ind,
                        AtmosphericProperty1D( atmos1d.clone1d() ) );
                    _calculate_change_points();
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    _change_points.clear();
                    _ax1.clear();
                    _props.clear();
                    _ax1.push_back( 0.0 );
                    _props.push_back(
                        AtmosphericProperty1D( atmos1d.clone1d() ) );
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& append(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    double r = range.get_as( _ax1.get_units() );
                    if (r <= _ax1.back()) {
                        throw std::range_error(
                            "Range to append is within existing boundaries, "
                            "use set() instead." );
                    }
                    _ax1.push_back( r );
                    _props.push_back(
                        AtmosphericProperty1D( atmos1d.clone1d() ) );
                    _calculate_change_points();
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                DECLARE_BOILERPLATE_METHODS(
                    piecewise_stratified_atmospheric_property_2d,
                    abstract_atmospheric_property )

                virtual abstract_atmospheric_property& copy(
                    const abstract_atmospheric_property& source ) override {
                    if (source.dimensions() == 2) {
                        const abstract_atmospheric_property_2d *srcptr
                            = dynamic_cast<
                                const abstract_atmospheric_property_2d *>(
                                &source );
                        return this->set( srcptr->axis( 0 ), srcptr->axis( 1 ),
                                          srcptr->values() );
                    } else if (source.dimensions() == 1) {
                        const abstract_atmospheric_property_1d *srcptr
                            = dynamic_cast<
                                const abstract_atmospheric_property_1d *>(
                                &source );
                        return this->set( *srcptr );
                    } else {
                        throw std::invalid_argument(
                            "piecewise_stratified_atmospheric_property_2d."
                            "copy(): Cannot "
                            "copy property with more than 2 dimension" );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual std::unique_ptr<abstract_atmospheric_property_2d>
                    clone2d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_2d>(
                        new piecewise_stratified_atmospheric_property_2d(
                            *this ) );
                }

                virtual size_t size( size_t d ) const override {
                    return ( d == 0 ? _ax1.size()
                                    : _props[ 0 ].axis().size() );
                }

                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Piecewise stratified property requires a 1-D "
                        "interpolator" );
                }

                virtual abstract_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _interptype = interp_type;
                    for (auto pit = _props.begin(); pit != _props.end();
                         ++pit) {
                        pit->set_interpolator( interp_type );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const units_ptr_t units ) override {
                    _change_points.clear();
                    _props.clear();
                    _ax1 = ax1;
                    for (auto it = _ax1.begin(); it != _ax1.end(); ++it) {
                        _props.push_back(
                            AtmosphericProperty1D( _atm_prop_1d_ptr_t(
                                new tuple_atmospheric_property_1d(
                                    ax2, vector_u_t( 0, units ) ) ) ) );
                    }
                    _calculate_change_points();
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& set(
                    size_t dim, const vector_u_t& ax ) override {
                    _validate_axis( dim );
                    if (dim == 0) {
                        auto newprops = _props;
                        newprops.clear();
                        _ax1.convert_units( *ax.get_units() );
                        for (auto i = 0; i < ax.size(); ++i) {
                            scalar_u_t r( ax[ i ], ax.get_units() );
                            newprops.push_back(
                                _props[ _profile_index( r ) ] );
                        }
                        _calculate_change_points();
                    } else {
                        for (auto it = _props.begin(); it != _props.end();
                             ++it) {
                            it->resample( ax );
                        }
                    }
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                };

                virtual abstract_atmospheric_property_2d& set(
                    const vector2d_u_t& ax ) override {
                    if (ax.dim( 0 ) != this->size( 0 )
                        || ax.dim( 1 ) != this->size( 1 )) {
                        std::ostringstream oss;
                        oss << "2-D vector size " << ax.dim( 0 ) << "x"
                            << ax.dim( 1 )
                            << " does not match existing axes sizes "
                            << this->size( 0 ) << "x" << this->size( 1 );
                        throw std::logic_error( oss.str() );
                    }
                    for (size_t i = 0; i < this->size( 0 ); ++i) {
                        vector_u_t propvals( ax[ i ], ax.get_units() );
                        _props[ i ].set( _props[ i ].axis(), propvals );
                    }
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                };

                virtual abstract_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector2d_u_t& vals ) override {
                    this->clear();
                    _ax1 = ax1;
                    for (size_t i = 0; i < _ax1.size(); ++i) {
                        _props.push_back(
                            AtmosphericProperty1D( _atm_prop_1d_ptr_t(
                                new tuple_atmospheric_property_1d(
                                    ax2,
                                    vector_u_t( vals[ i ],
                                                vals.get_units() ) ) ) ) );
                    }
                    _calculate_change_points();
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual vector_u_t& axis( size_t n ) override {
                    _validate_axis( n );
                    if (n == 0) {
                        return _ax1;
                    } else {
                        return _props[ 0 ].axis();
                    }
                }

                virtual const vector_u_t& axis( size_t n ) const override {
                    _validate_axis( n );
                    if (n == 0) {
                        return _ax1;
                    } else {
                        return _props[ 0 ].axis();
                    }
                }

                virtual vector2d_u_t& values() override {
                    vector_u_t z = _props.at( 0 ).axis();
                    return _tmpvals;
                }

                virtual const vector2d_u_t& values() const override {
                    vector_u_t z = _props.at( 0 ).axis();
                    return _tmpvals;
                }

                virtual double get( double val1, double val2 ) override {
                    return _props[ _profile_index( val1 ) ].get( val2 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     size_t rel ) override {
                    _validate_axis( rel );
                    return ( rel == 0 ? 0.0
                                      : _props[ _profile_index( val1 ) ]
                                            .get_first_derivative( val2 ) );
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      size_t rel1,
                                                      size_t rel2 ) override {
                    _validate_axis( rel1 );
                    _validate_axis( rel2 );
                    return ( rel1 == 0 || rel2 == 0
                                 ? 0.0
                                 : _props[ _profile_index( val1 ) ]
                                       .get_second_derivative( val2 ) );
                }

                virtual _atm_prop_1d_ptr_t extract( double range ) override {
                    return _props[ _profile_index( range ) ]
                        .internal()
                        ->clone1d();
                }

                virtual const units_ptr_t get_units() const override {
                    return _props[ 0 ].get_units();
                }

                virtual const units_ptr_t get_axis_units(
                    size_t n ) const override {
                    _validate_axis( n );
                    return ( n == 0 ? _ax1.get_units()
                                    : _props[ 0 ].get_axis_units() );
                }

                virtual abstract_atmospheric_property_2d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    for (auto it = _props.begin(); it != _props.end(); ++it) {
                        it->convert_units( u );
                    }
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D
                }

                virtual abstract_atmospheric_property_2d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) override {
                    _validate_axis( n );
                    if (n == 0) {
                        _ax1.convert_units( u );
                    } else {
                        for (auto it = _props.begin(); it != _props.end();
                             ++it) {
                            it->convert_axis_units( u );
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& resample(
                    const vector_u_t& new_r,
                    const vector_u_t& new_z ) override {
                    for (auto it = _props.begin(); it != _props.end(); ++it) {
                        it->resample( new_z );
                    }
                    _make_tmpvals();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D
                }

                virtual abstract_atmospheric_property_2d& set_limits(
                    size_t dim, double minax, double maxax ) override {
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual abstract_atmospheric_property_2d& reset_limits(
                    size_t dim ) override {
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_2D;
                }

                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const override {
                    _validate_axis( dim );
                    if (dim == 0) {
                        return std::pair<double, double> { -1.0, -1.0 };
                    } else {
                        auto z        = _props[ 0 ].axis();
                        double minlim = z.front();
                        double maxlim = z.back();
                        for (auto it = _props.cbegin(); it != _props.cend();
                             ++it) {
                            z      = it->axis();
                            minlim = std::max( minlim, z.front() );
                            maxlim = std::min( maxlim, z.back() );
                        }
                        return std::pair<double, double> { minlim, maxlim };
                    }
                }

            protected:
                size_t _profile_index( scalar_u_t r_s ) {
                    return _profile_index(
                        r_s.get_as( _change_points.get_units() ) );
                }

                size_t _profile_index( double r ) {
                    for (size_t i = 0; i < _change_points.size(); ++i) {
                        if (r < _change_points[ i ]) {
                            return i;
                        }
                    }
                    return _props.size() - 1;
                }

                void _calculate_change_points() {
                    _change_points.clear();
                    _change_points.set_units( *_ax1.get_units() );
                    for (size_t i = 1; i < _ax1.size(); ++i) {
                        _change_points.push_back(
                            0.5 * ( _ax1[ i ] - _ax1[ i - 1 ] ) );
                    }
                }

                void _make_tmpvals() {
                    vector_u_t z = _props.at( 0 ).axis();
                    _tmpvals     = vector2d_u_t( _ax1.size(), z.size(),
                                                 _props[ 0 ].get_units() );
                    for (size_t i = 0; i < _ax1.size(); ++i) {
                        AtmosphericProperty1D prop( _props[ i ] );
                        for (size_t j = 0; j < z.size(); ++j) {
                            _tmpvals[ i ][ j ] = prop.get( z[ j ] );
                        }
                    }
                }

                vector_u_t _ax1;
                vector_u_t _change_points;
                std::vector<AtmosphericProperty1D> _props;
                NCPA::interpolation::interpolator_1d_type_t _interptype;
                vector2d_u_t _tmpvals;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract_atmospheric_property_2d& a,
                  NCPA::atmos::abstract_atmospheric_property_2d& b ) noexcept {
}

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

static void swap(
    NCPA::atmos::stratified_atmospheric_property_2d& a,
    NCPA::atmos::stratified_atmospheric_property_2d& b ) noexcept {
    using std::swap;
    ::swap(
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_2d&>( a ),
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_2d&>( a ) );
    swap( a._z, b._z );
    swap( a._vals, b._vals );
    swap( a._spline, b._spline );
    swap( a._z_limits, b._z_limits );
    swap( a._dummy, b._dummy );
    swap( a._as_2d, b._as_2d );
}

static void swap(
    NCPA::atmos::piecewise_stratified_atmospheric_property_2d& a,
    NCPA::atmos::piecewise_stratified_atmospheric_property_2d& b ) noexcept {
    using std::swap;
    ::swap(
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_2d&>( a ),
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_2d&>( a ) );
    // swap( a._axis_limits, b._axis_limits );
    swap( a._ax1, b._ax1 );
    swap( a._change_points, b._change_points );
    swap( a._props, b._props );
    swap( a._interptype, b._interptype );
    swap( a._tmpvals, b._tmpvals );
}

namespace NCPA {
    namespace atmos {

        class AtmosphericProperty2D {
            public:
                AtmosphericProperty2D() {}

                AtmosphericProperty2D( _atm_prop_2d_ptr_t engine ) {
                    _ptr = std::move( engine->clone2d() );
                }

                AtmosphericProperty2D( const AtmosphericProperty2D& other ) :
                    AtmosphericProperty2D() {
                    _ptr = std::move( other._ptr->clone2d() );
                }

                AtmosphericProperty2D(
                    const abstract_atmospheric_property& prop ) {
                    if (prop.dimensions() == 2) {
                        _ptr = std::move(
                            dynamic_cast<
                                const abstract_atmospheric_property_2d&>(
                                prop )
                                .clone2d() );
                    } else {
                        throw std::invalid_argument(
                            "AtmosphericProperty2D: Can only ingest 2-D "
                            "atmospheric properties" );
                    }
                }

                AtmosphericProperty2D& set(
                    const vector_u_t& ranges,
                    const std::vector<const AtmosphericProperty1D *>& props ) {
                    _check_pointer();
                    std::vector<const abstract_atmospheric_property_1d *> ptrs;
                    for (auto it = props.cbegin(); it != props.cend(); ++it) {
                        ptrs.push_back( ( *it )->internal() );
                    }
                    _ptr->set( ranges, ptrs );
                    return *this;
                }

                AtmosphericProperty2D& set(
                    const scalar_u_t& range,
                    const AtmosphericProperty1D& prop ) {
                    _check_pointer();
                    return this->set(
                        vector_u_t( 1, range ),
                        std::vector<const AtmosphericProperty1D *>( 1,
                                                                    &prop ) );
                }

                AtmosphericProperty2D& set(
                    const AtmosphericProperty1D& prop ) {
                    _check_pointer();
                    return this->set(
                        scalar_u_t( 0.0, this->get_axis_units( 0 ) ), prop );
                }

                AtmosphericProperty2D& append(
                    const scalar_u_t& range,
                    const AtmosphericProperty1D& prop ) {
                    _check_pointer();
                    _ptr->append( range, *prop.internal() );
                    return *this;
                }

                virtual AtmosphericProperty1D extract( double range ) const {
                    _check_pointer();
                    auto ptr = _ptr->extract( range );
                    return AtmosphericProperty1D( _ptr->extract( range ) );
                }

                DECLARE_WRAPPER_BOILERPLATE_METHODS( AtmosphericProperty2D )

                virtual size_t size( size_t d ) const {
                    _check_pointer();
                    return _ptr->size( d );
                }

                virtual AtmosphericProperty2D& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type ) {
                    _check_pointer();
                    _ptr->set_interpolator( interp_type );
                    return *this;
                }

                virtual AtmosphericProperty2D& set( size_t dim,
                                                    const vector_u_t& ax ) {
                    _check_pointer();
                    _ptr->set( dim, ax );
                    return *this;
                }

                virtual AtmosphericProperty2D& set( const vector2d_u_t& ax ) {
                    _check_pointer();
                    _ptr->set( ax );
                    return *this;
                }

                virtual AtmosphericProperty2D& set( const vector_u_t& ax1,
                                                    const vector_u_t& ax2,
                                                    const vector2d_u_t vals ) {
                    _check_pointer();
                    _ptr->set( ax1, ax2, vals );
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

                virtual vector2d_u_t& values() {
                    _check_pointer();
                    return _ptr->values();
                }

                virtual const vector2d_u_t& values() const {
                    _check_pointer();
                    return _ptr->values();
                }

                virtual double get( double val1, double val2 ) {
                    _check_pointer();
                    return _ptr->get( val1, val2 );
                }

                virtual double f( double val1, double val2 ) {
                    _check_pointer();
                    return _ptr->f( val1, val2 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     size_t rel ) {
                    _check_pointer();
                    return _ptr->get_first_derivative( val1, val2, rel );
                }

                virtual double df( double val1, double val2, size_t rel ) {
                    _check_pointer();
                    return _ptr->df( val1, val2, rel );
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      size_t rel1,
                                                      size_t rel2 ) {
                    _check_pointer();
                    return _ptr->get_second_derivative( val1, val2, rel1,
                                                        rel2 );
                }

                virtual double ddf( double val1, double val2, size_t rel1,
                                    size_t rel2 ) {
                    _check_pointer();
                    return _ptr->ddf( val1, val2, rel1, rel2 );
                }

                virtual const units_ptr_t get_units() const {
                    _check_pointer();
                    return _ptr->get_units();
                }

                virtual const units_ptr_t get_axis_units( size_t n ) const {
                    _check_pointer();
                    return _ptr->get_axis_units( n );
                }

                virtual AtmosphericProperty2D& convert_units(
                    const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_units( u );
                    return *this;
                }

                virtual AtmosphericProperty2D& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_axis_units( n, u );
                    return *this;
                }

                virtual AtmosphericProperty2D& resample(
                    const vector_u_t& new_r, const vector_u_t& new_z ) {
                    _check_pointer();
                    _ptr->resample( new_r, new_z );
                    return *this;
                }

                virtual AtmosphericProperty2D& resample(
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

                virtual AtmosphericProperty2D& set_limits( size_t dim,
                                                           double minax,
                                                           double maxax ) {
                    _check_pointer();
                    _ptr->set_limits( dim, minax, maxax );
                    return *this;
                }

                virtual _atm_prop_2d_ptr_t& internal() { return _ptr; }

                virtual const _atm_prop_2d_ptr_t& internal() const {
                    return _ptr;
                }


            protected:
                _atm_prop_2d_ptr_t _ptr;

                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphericProperty2D: no pointer has been set" );
                    }
                }
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphericProperty2D& a,
                  NCPA::atmos::AtmosphericProperty2D& b ) noexcept {
    using std::swap;
    swap( a._ptr, b._ptr );
}
