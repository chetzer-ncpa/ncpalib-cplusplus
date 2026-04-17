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

static void swap(
    NCPA::atmos::stratified_atmospheric_property_2d& a,
    NCPA::atmos::stratified_atmospheric_property_2d& b ) noexcept;

namespace NCPA {
    namespace atmos {
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

                virtual stratified_atmospheric_property_2d& clear() override {
                    _z.clear();
                    _vals.clear();
                    _z_limits = std::pair<double, double>( { 0, 0 } );
                    _spline.clear();
                    _as_2d.clear();
                    return *this;
                }

                virtual stratified_atmospheric_property_2d& set(
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
                    return *this;
                }

                virtual stratified_atmospheric_property_2d& set(
                    const scalar_u_t& range,
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    return this->set(
                        vector_u_t( 1, range ),
                        std::vector<const abstract_atmospheric_property_1d *>(
                            1, &atmos1d ) );
                }

                virtual stratified_atmospheric_property_2d& set(
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    return this->set(
                        scalar_u_t( 0.0, NCPA::units::KILOMETERS ), atmos1d );
                }

                virtual stratified_atmospheric_property_2d& append(
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

                stratified_atmospheric_property_2d(
                    stratified_atmospheric_property_2d&& source ) noexcept :
                    stratified_atmospheric_property_2d() {
                    ::swap( *this, source );
                }

                virtual ~stratified_atmospheric_property_2d() {}

                stratified_atmospheric_property_2d& operator=(
                    stratified_atmospheric_property_2d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap(
                    stratified_atmospheric_property_2d& a,
                    stratified_atmospheric_property_2d& b ) noexcept;

                virtual std::unique_ptr<abstract_atmospheric_property> clone()
                    const override {
                    return std::unique_ptr<abstract_atmospheric_property>(
                        new stratified_atmospheric_property_2d( *this ) );
                }

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
                    return *this;
                }

                virtual std::unique_ptr<abstract_atmospheric_property_2d>
                    clone2d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_2d>(
                        new stratified_atmospheric_property_2d( *this ) );
                }

                virtual size_t size( size_t d ) const override {
                    return ( d == 0 ? 1 : _z.size() );
                }

                virtual stratified_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Stratified property requires a 1-D interpolator" );
                }

                virtual stratified_atmospheric_property_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    this->_init_spline();
                    return *this;
                }

                virtual stratified_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const units_ptr_t units ) override {
                    this->clear();
                    _z    = ax2;
                    _vals = vector_u_t( _z.size(), units );
                    reset_limits( 1 );
                    return *this;
                }

                virtual stratified_atmospheric_property_2d& set(
                    size_t dim, const vector_u_t& ax ) override {
                    if (dim == 1) {
                        this->resample(
                            vector_u_t( { 0 }, NCPA::units::KILOMETERS ), ax );
                        this->reset_limits( dim );
                    }
                    _copy_to_2d();
                    return *this;
                };

                virtual stratified_atmospheric_property_2d& set(
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
                    return *this;
                };

                virtual stratified_atmospheric_property_2d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    const vector2d_u_t& vals ) override {
                    // this->set( 0, ax1 );
                    this->set( 1, ax2 );
                    this->set( vals );
                    _copy_to_2d();
                    return *this;
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

                virtual stratified_atmospheric_property_2d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _vals.convert_units( u );
                    _init_spline();
                    _copy_to_2d();
                    return *this;
                }

                virtual stratified_atmospheric_property_2d& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) override {
                    _z.convert_units( u );
                    this->reset_limits( n );
                    _init_spline();
                    _copy_to_2d();
                    return *this;
                }

                virtual stratified_atmospheric_property_2d& resample(
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
                    return *this;
                }

                virtual stratified_atmospheric_property_2d& set_limits(
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
                    return *this;
                }

                virtual stratified_atmospheric_property_2d& reset_limits(
                    size_t dim ) override {
                    if (dim == 1) {
                        _z_limits.first  = _z.front();
                        _z_limits.second = _z.back();
                    }
                    return *this;
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
    }  // namespace atmos
}  // namespace NCPA

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
