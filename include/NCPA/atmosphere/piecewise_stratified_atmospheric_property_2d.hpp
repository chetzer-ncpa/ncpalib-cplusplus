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
    NCPA::atmos::piecewise_stratified_atmospheric_property_2d& a,
    NCPA::atmos::piecewise_stratified_atmospheric_property_2d& b ) noexcept;

namespace NCPA {
    namespace atmos {
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

                virtual piecewise_stratified_atmospheric_property_2d& clear()
                    override {
                    _ax1.clear();
                    _props.clear();
                    _change_points.clear();
                    _tmpvals.clear();
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d& set(
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
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d& set(
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
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d& set(
                    const abstract_atmospheric_property_1d& atmos1d )
                    override {
                    _change_points.clear();
                    _ax1.clear();
                    _props.clear();
                    _ax1.push_back( 0.0 );
                    _props.push_back(
                        AtmosphericProperty1D( atmos1d.clone1d() ) );
                    _make_tmpvals();
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d& append(
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
                    return *this;
                }

                piecewise_stratified_atmospheric_property_2d(
                    piecewise_stratified_atmospheric_property_2d&&
                        source ) noexcept :
                    piecewise_stratified_atmospheric_property_2d() {
                    ::swap( *this, source );
                }

                virtual ~piecewise_stratified_atmospheric_property_2d() {}

                piecewise_stratified_atmospheric_property_2d& operator=(
                    piecewise_stratified_atmospheric_property_2d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap(
                    piecewise_stratified_atmospheric_property_2d& a,
                    piecewise_stratified_atmospheric_property_2d& b ) noexcept;

                virtual std::unique_ptr<abstract_atmospheric_property> clone()
                    const override {
                    return std::unique_ptr<abstract_atmospheric_property>(
                        new piecewise_stratified_atmospheric_property_2d(
                            *this ) );
                }

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
                    return *this;
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

                virtual piecewise_stratified_atmospheric_property_2d&
                    set_interpolator(
                        NCPA::interpolation::interpolator_2d_type_t
                            interp_type ) override {
                    throw NCPA::NotImplementedError(
                        "Piecewise stratified property requires a 1-D "
                        "interpolator" );
                }

                virtual piecewise_stratified_atmospheric_property_2d&
                    set_interpolator(
                        NCPA::interpolation::interpolator_1d_type_t
                            interp_type ) override {
                    _interptype = interp_type;
                    for (auto pit = _props.begin(); pit != _props.end();
                         ++pit) {
                        pit->set_interpolator( interp_type );
                    }
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d& set(
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
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d& set(
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
                    return *this;
                };

                virtual piecewise_stratified_atmospheric_property_2d& set(
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
                    return *this;
                };

                virtual piecewise_stratified_atmospheric_property_2d& set(
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
                    return *this;
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

                virtual piecewise_stratified_atmospheric_property_2d&
                    convert_units( const NCPA::units::Unit& u ) override {
                    for (auto it = _props.begin(); it != _props.end(); ++it) {
                        it->convert_units( u );
                    }
                    _make_tmpvals();
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d&
                    convert_axis_units( size_t n,
                                        const NCPA::units::Unit& u ) override {
                    _validate_axis( n );
                    if (n == 0) {
                        _ax1.convert_units( u );
                    } else {
                        for (auto it = _props.begin(); it != _props.end();
                             ++it) {
                            it->convert_axis_units( u );
                        }
                    }
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d& resample(
                    const vector_u_t& new_r,
                    const vector_u_t& new_z ) override {
                    for (auto it = _props.begin(); it != _props.end(); ++it) {
                        it->resample( new_z );
                    }
                    _make_tmpvals();
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d&
                    set_limits( size_t dim, double minax,
                                double maxax ) override {
                    return *this;
                }

                virtual piecewise_stratified_atmospheric_property_2d&
                    reset_limits( size_t dim ) override {
                    return *this;
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
