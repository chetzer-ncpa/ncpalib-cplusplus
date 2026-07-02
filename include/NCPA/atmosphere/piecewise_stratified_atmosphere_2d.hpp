#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_2d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/AtmosphericModel.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/tuple_atmospheric_property_1d.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/units.hpp"

#include <memory>
#include <string>

static void swap( NCPA::atmos::piecewise_stratified_atmosphere_2d&,
                  NCPA::atmos::piecewise_stratified_atmosphere_2d& ) noexcept;

namespace NCPA {
    namespace atmos {
        class piecewise_stratified_atmosphere_2d
            : public abstract_atmosphere_2d {
            public:
                piecewise_stratified_atmosphere_2d() :
                    abstract_atmosphere_2d() {
                    _axes[ 0 ].set_units(
                        NCPA_ATMOSPHERE_DEFAULT_DISTANCE_UNITS );
                    _axes[ 1 ].set_units(
                        NCPA_ATMOSPHERE_DEFAULT_DISTANCE_UNITS );
                }

                piecewise_stratified_atmosphere_2d(
                    const piecewise_stratified_atmosphere_2d& other ) :
                    piecewise_stratified_atmosphere_2d() {
                    _components           = other._components;
                    _scalar_properties    = other._scalar_properties;
                    _1d_interpolator_type = other._1d_interpolator_type;
                    _axes                 = other._axes;
                    _change_points        = other._change_points;
                    _build_scalar_splines();
                }

                piecewise_stratified_atmosphere_2d(
                    piecewise_stratified_atmosphere_2d&& source ) noexcept :
                    piecewise_stratified_atmosphere_2d() {
                    ::swap( *this, source );
                }

                virtual ~piecewise_stratified_atmosphere_2d() {}

                piecewise_stratified_atmosphere_2d& operator=(
                    piecewise_stratified_atmosphere_2d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap(
                    piecewise_stratified_atmosphere_2d& a,
                    piecewise_stratified_atmosphere_2d& b ) noexcept;

                virtual std::unique_ptr<abstract_atmosphere_2d> clone()
                    const override {
                    return std::unique_ptr<abstract_atmosphere_2d>(
                        new piecewise_stratified_atmosphere_2d( *this ) );
                }

                virtual size_t size( size_t dim ) const override {
                    this->validate_axis( dim );
                    return _axes[ dim ].size();
                }

                virtual bool is_stratified() const override { return false; }

                virtual abstract_atmosphere_2d& set(
                    abstract_atmosphere_1d& atmos1d ) override {
                    _components.clear();
                    _scalar_properties.clear();
                    vector_u_t tmp_r = _axes[ 0 ];
                    _axes[ 0 ].clear();
                    _axes[ 1 ]       = atmos1d.get_axis_vector();
                    auto scalar_keys = atmos1d.get_scalar_keys();
                    for (auto it = scalar_keys.cbegin();
                         it != scalar_keys.cend(); ++it) {
                        _scalar_properties[ *it ]
                            = std::make_pair<vector_u_t, vector_u_t>(
                                vector_u_t( 0, _axes[ 0 ].get_units() ),
                                vector_u_t( 0, atmos1d.get_units( *it ) ) );
                    }
                    for (size_t i = 0; i < tmp_r.size(); ++i) {
                        this->append(
                            scalar_u_t( tmp_r[ i ], tmp_r.get_units() ),
                            atmos1d );
                    }
                    _build_scalar_splines();
                    _calculate_change_points();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& set(
                    const vector_u_t ranges,
                    std::vector<abstract_atmosphere_1d *> components )
                    override {
                    _axes[ 0 ].clear();
                    _components.clear();
                    _scalar_properties.clear();
                    _axes[ 1 ]       = components[ 0 ]->get_axis_vector();
                    auto scalar_keys = components[ 0 ]->get_scalar_keys();
                    for (auto it = scalar_keys.cbegin();
                         it != scalar_keys.cend(); ++it) {
                        _scalar_properties[ *it ]
                            = std::pair<vector_u_t, vector_u_t> {
                                  vector_u_t( 0, _axes[ 0 ].get_units() ),
                                  vector_u_t(
                                      0, components[ 0 ]->get_units( *it ) )
                              };
                    }
                    for (size_t i = 0; i < ranges.size(); ++i) {
                        this->append(
                            scalar_u_t( ranges[ i ], ranges.get_units() ),
                            *components[ i ] );
                    }
                    _axes[ 0 ] = ranges;
                    _build_scalar_splines();
                    _calculate_change_points();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& append(
                    scalar_u_t range,
                    abstract_atmosphere_1d& atmos1d ) override {
                    double r = range.get_as( _axes[ 0 ].get_units() );
                    if (_axes[ 0 ].size() > 0 && r <= _axes[ 0 ].back()) {
                        std::ostringstream oss;
                        oss << "Supplied range " << r
                            << " is inside existing limits";
                        throw std::range_error( oss.str() );
                    }
                    _axes[ 0 ].push_back( r );
                    if (_axes[ 1 ].size() == 0) {
                        _axes[ 1 ] = atmos1d.get_axis_vector();
                    }
                    _components.push_back( Atmosphere1D( atmos1d.clone() )
                                               .resample( _axes[ 1 ] ) );
                    auto scalar_keys = atmos1d.get_scalar_keys();
                    for (auto it = scalar_keys.cbegin();
                         it != scalar_keys.cend(); ++it) {
                        _scalar_properties[ *it ].first.push_back( r );
                        _scalar_properties[ *it ].second.push_back(
                            atmos1d.get( *it ) );
                    }
                    _build_scalar_splines();
                    _calculate_change_points();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual std::unique_ptr<abstract_atmosphere_1d> extract(
                    double range ) override {
                    return _at( range ).internal()->clone();
                }

                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "2-D interpolator not applicable to locally "
                        "stratified "
                        "atmospheres" );
                }

                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _1d_interpolator_type = interp_type;
                    for (auto it = _components.begin();
                         it != _components.end(); ++it) {
                        it->set_interpolator( interp_type );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& set_axis(
                    size_t axis, vector_u_t vals ) override {
                    this->validate_axis( axis );
                    if (axis == 1) {
                        this->resample( 1, vals );
                    } else {
                        if (vals.size() == _axes[ 0 ].size()) {
                            _axes[ 0 ] = vals;
                            _build_scalar_splines();
                            _calculate_change_points();
                        } else {
                            throw std::logic_error(
                                "Cannot resample axis 0 of locally stratified "
                                "atmosphere" );
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const AtmosphericProperty2D& property ) override {
                    _assert_does_not_contain( key );
                    if (_axes[ 1 ].empty()) {
                        _axes[ 1 ] = property.axis( 1 );
                    }
                    if (_components.size() == 0) {
                        throw std::logic_error(
                            "Cannot add property, no existing atmospheres "
                            "to add to" );
                    }
                    for (size_t i = 0; i < _axes[ 0 ].size(); ++i) {
                        AtmosphericProperty1D prop
                            = property.extract( _axes[ 0 ][ i ] )
                                  .resample( _axes[ 1 ] );
                        _components[ i ].add_property( key, prop );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const vector2d_u_t& property ) override {
                    if (_axes[ 0 ].size() != property.dim( 0 )) {
                        throw std::range_error(
                            "Size mismatch between internal axis 1 vector and "
                            "first dimension of provided property" );
                    }
                    if (_axes[ 1 ].size() != property.dim( 1 )) {
                        throw std::range_error(
                            "Size mismatch between internal Z vector and "
                            "second dimension of provided property" );
                    }
                    for (size_t i = 0; i < _components.size(); ++i) {
                        vector_u_t p( property[ i ], property.get_units() );
                        tuple_atmospheric_property_1d prop( _axes[ 1 ], p );
                        _components[ i ].add_property( key, prop );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const vector_u_t& property ) override {
                    if (property.size() != _axes[ 0 ].size()) {
                        throw std::range_error(
                            "Size mismatch between internal axis 1 vector and "
                            "new property vector" );
                    }
                    _assert_does_not_contain( key );
                    _scalar_properties[ key ]
                        = std::make_pair( _axes[ 0 ], property );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key, const vector_u_t& property,
                    const vector_u_t& index ) override {
                    _assert_does_not_contain( key );
                    if (property.size() != index.size()) {
                        throw std::range_error(
                            "add_property(): Size mismatch between supplied "
                            "vectors" );
                    }
                    _scalar_properties[ key ] = { index, property };
                    _build_scalar_splines();
                    _calculate_change_points();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const scalar_u_t& property ) override {
                    _assert_does_not_contain( key );
                    _scalar_properties[ key ].second
                        = vector_u_t( _axes[ 0 ].size(), property );
                    _build_scalar_splines();
                    _calculate_change_points();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& remove_property(
                    const std::string& key ) override {
                    _assert_contains( key );
                    if (this->contains_scalar( key )) {
                        _scalar_properties.erase( key );
                        _build_scalar_splines();
                        _calculate_change_points();
                    } else {
                        for (auto it = _components.begin();
                             it != _components.end(); ++it) {
                            it->remove_property( key );
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& copy_property(
                    const std::string& old_key,
                    const std::string& new_key ) override {
                    _assert_contains( old_key );
                    _assert_does_not_contain( new_key );
                    if (this->contains_scalar( old_key )) {
                        _scalar_properties[ new_key ]
                            = _scalar_properties[ old_key ];
                        _build_scalar_splines();
                        _calculate_change_points();
                    } else {
                        for (auto it = _components.begin();
                             it != _components.end(); ++it) {
                            it->copy_property( old_key, new_key );
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual vector_u_t get_axis_vector( size_t n ) override {
                    this->validate_axis( n );
                    return _axes[ n ];
                }

                virtual vector_u_t get_axis_vector( size_t n ) const override {
                    this->validate_axis( n );
                    return _axes[ n ];
                }

                virtual double get( const std::string& key,
                                    double r ) override {
                    if (this->contains_scalar( key )) {
                        return _scalar_splines[ key ].eval_f( r );
                    } else {
                        return _at( r ).get( key );
                    }
                }

                virtual double get( const std::string& key, double r,
                                    double z ) override {
                    return _at( r ).get( key, z );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double r, double z,
                                                     size_t wrt ) override {
                    if (_scalar_splines.count( key ) != 0) {
                        return ( wrt == 0 ? _scalar_splines[ key ].eval_df( r )
                                          : 0.0 );
                    } else {
                        return ( wrt == 0 ? 0.0
                                          : _at( r ).get_first_derivative(
                                                key, z ) );
                    }
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double r, double z,
                                                      size_t wrt1,
                                                      size_t wrt2 ) override {
                    if (_scalar_splines.count( key ) != 0) {
                        return ( wrt1 == 0 && wrt2 == 0
                                     ? _scalar_splines[ key ].eval_ddf( r )
                                     : 0.0 );
                    } else {
                        return (
                            wrt1 == 0 || wrt2 == 0
                                ? 0.0
                                : _at( r ).get_second_derivative( key, z ) );
                    }
                }

                virtual units_ptr_t get_axis_units( size_t n ) const override {
                    this->validate_axis( n );
                    return _axes[ n ].get_units();
                }

                virtual units_ptr_t get_units(
                    const std::string& key ) const override {
                    if (this->contains_scalar( key )) {
                        return _scalar_properties.at( key ).second.get_units();
                    } else {
                        return _components.at( 0 ).get_units( key );
                    }
                }

                virtual double get_minimum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    return _axes[ n ].front();
                }

                virtual double get_maximum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    return _axes[ n ].back();
                }

                virtual abstract_atmosphere_2d& convert_axis_units(
                    size_t n, units_ptr_t new_units ) override {
                    this->validate_axis( n );
                    if (n == 0) {
                        _axes[ 0 ].convert_units( *new_units );
                        for (auto it = _scalar_properties.begin();
                             it != _scalar_properties.end(); ++it) {
                            it->second.first.convert_units( *new_units );
                        }
                        _build_scalar_splines();
                        _calculate_change_points();

                    } else {
                        for (auto it = _components.begin();
                             it != _components.end(); ++it) {
                            it->convert_axis_units( new_units );
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& convert_units(
                    const std::string& key, units_ptr_t new_units ) override {
                    if (this->contains_scalar( key )) {
                        _scalar_properties[ key ].second.convert_units(
                            *new_units );
                        _build_scalar_splines();
                        _calculate_change_points();
                    } else {
                        for (auto it = _components.begin();
                             it != _components.end(); ++it) {
                            it->convert_units( key, new_units );
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    size_t axis, double new_d ) override {
                    this->validate_axis( axis );
                    vector_u_t new_axis = _axes[ axis ];
                    new_axis.resize( 0 );
                    double d = _axes[ axis ].front();
                    while (d <= _axes[ axis ].back()) {
                        new_axis.push_back( d );
                        d += new_d;
                    }
                    return this->resample( axis, new_axis );
                }

                virtual abstract_atmosphere_2d& resample(
                    size_t axis, vector_u_t new_z ) override {
                    this->validate_axis( axis );
                    if (axis == 1) {
                        for (auto it = _components.begin();
                             it != _components.end(); ++it) {
                            it->resample( new_z );
                        }
                    } else {
                        for (auto it = _scalar_properties.begin();
                             it != _scalar_properties.end(); ++it) {
                            it->second.first.convert_units(
                                *new_z.get_units() );
                            NCPA::interpolation::Interpolator1D<double, double>
                                interp = NCPA::interpolation::
                                    InterpolatorFactory<double, double>::build(
                                        _1d_interpolator_type );
                            interp.init( it->second.first.size() )
                                .fill( it->second.first, it->second.second )
                                .ready();
                            vector_u_t newv = it->second.second;
                            newv.resize( 0 );
                            for (auto it = new_z.cbegin(); it != new_z.cend();
                                 ++it) {
                                newv.push_back( interp.eval_f( *it ) );
                            }
                            it->second.first  = new_z;
                            it->second.second = newv;
                        }
                        _build_scalar_splines();
                        _calculate_change_points();
                    }
                    _axes[ axis ] = new_z;
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    vector_u_t new_ax1, vector_u_t new_ax2 ) override {
                    for (auto it = _scalar_properties.begin();
                         it != _scalar_properties.end(); ++it) {
                        it->second.first.convert_units( *new_ax1.get_units() );
                        NCPA::interpolation::Interpolator1D<double, double>
                            interp = NCPA::interpolation::
                                InterpolatorFactory<double, double>::build(
                                    _1d_interpolator_type );
                        interp.init( it->second.first.size() )
                            .fill( it->second.first, it->second.second )
                            .ready();
                        vector_u_t newv = it->second.second;
                        newv.resize( 0 );
                        for (auto it = new_ax1.cbegin(); it != new_ax1.cend();
                             ++it) {
                            newv.push_back( interp.eval_f( *it ) );
                        }
                        it->second.first  = new_ax1;
                        it->second.second = newv;
                    }
                    _axes[ 0 ] = new_ax1;
                    for (auto it = _components.begin();
                         it != _components.end(); ++it) {
                        it->resample( new_ax2 );
                    }
                    _axes[ 1 ] = new_ax2;
                    _build_scalar_splines();
                    _calculate_change_points();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual std::vector<std::string> get_keys() const override {
                    std::vector<std::string> vkeys = this->get_vector_keys(),
                                             skeys = this->get_scalar_keys(),
                                             finalkeys;
                    std::set_intersection( vkeys.cbegin(), vkeys.cend(),
                                           skeys.cbegin(), skeys.cend(),
                                           std::back_inserter( finalkeys ) );
                    return finalkeys;
                }

                virtual std::vector<std::string> get_vector_keys()
                    const override {
                    if (_components.empty()) {
                        return std::vector<std::string> {};
                    } else {
                        return _components.at( 0 ).get_vector_keys();
                    }
                }

                virtual std::vector<std::string> get_scalar_keys()
                    const override {
                    std::vector<std::string> spkeys;
                    for (auto it = _scalar_properties.cbegin();
                         it != _scalar_properties.cend(); ++it) {
                        spkeys.push_back( it->first );
                    }
                    return spkeys;
                }

                virtual bool contains_scalar(
                    const std::string& key ) const override {
                    return ( _scalar_properties.count( key ) > 0 );
                }

                virtual bool contains_vector(
                    const std::string& key ) const override {
                    return ( ( !_components.empty() )
                             && _components.at( 0 ).contains_vector( key ) );
                }

                virtual bool contains_key(
                    const std::string& key ) const override {
                    return this->contains_vector( key )
                        || this->contains_scalar( key );
                }

                virtual void print( std::ostream& os ) override {
                    throw NCPA::NotImplementedError(
                        "Print function not yet implemented." );
                }

                virtual bool same( scalar_u_t r1,
                                   scalar_u_t r2 ) const override {
                    return ( _profile_index( r1 ) == _profile_index( r2 ) );
                }

                virtual bool same( double r1, double r2 ) const override {
                    return ( _profile_index( r1 ) == _profile_index( r2 ) );
                }

                size_t _profile_index( scalar_u_t r_s ) const {
                    return _profile_index(
                        r_s.get_as( _change_points.get_units() ) );
                }

                size_t _profile_index( double r ) const {
                    for (size_t i = 0; i < _change_points.size(); ++i) {
                        if (r <= _change_points[ i ]) {
                            return i - 1;
                        }
                    }
                    return _components.size() - 1;
                }

                Atmosphere1D& _at( double r ) {
                    return _components[ _profile_index( r ) ];
                }

                Atmosphere1D& _at( scalar_u_t r ) {
                    return _components[ _profile_index( r ) ];
                }

            protected:
                void _assert_contains_vector( const std::string& key ) const {
                    if (!contains_vector( key )) {
                        throw std::range_error( "Key " + key
                                                + " does not exist!" );
                    }
                }

                void _assert_contains_scalar( const std::string& key ) const {
                    if (!contains_scalar( key )) {
                        throw std::range_error( "Key " + key
                                                + " does not exist!" );
                    }
                }

                void _assert_contains( const std::string& key ) const {
                    _assert_contains_vector( key );
                    _assert_contains_scalar( key );
                }

                void _assert_does_not_contain( const std::string& key ) const {
                    if (contains_key( key )) {
                        throw std::range_error( "Key " + key
                                                + " already exists!" );
                    }
                }

                void _build_scalar_splines() {
                    _scalar_splines.clear();
                    for (auto propit = _scalar_properties.cbegin();
                         propit != _scalar_properties.cend(); ++propit) {
                        _scalar_splines[ propit->first ]
                            = NCPA::interpolation::
                                InterpolatorFactory<double, double>::build(
                                    _1d_interpolator_type );
                        _scalar_splines[ propit->first ]
                            .init( propit->second.first.size() )
                            .fill( propit->second.first,
                                   propit->second.second )
                            .ready();
                    }
                }

                void _calculate_change_points() {
                    _change_points.clear();
                    _change_points.set_units( *_axes[ 0 ].get_units() );
                    _change_points.push_back( 0.0 );
                    for (size_t i = 1; i < _axes[ 0 ].size(); ++i) {
                        _change_points.push_back(
                            _change_points[ i - 1 ]
                            + 0.5
                                  * ( _axes[ 0 ][ i ]
                                      - _axes[ 0 ][ i - 1 ] ) );
                    }
                }

                // vector_u_t _ax1, _ax2;
                std::array<vector_u_t, 2> _axes;
                std::vector<Atmosphere1D> _components;
                std::unordered_map<std::string,
                                   std::pair<vector_u_t, vector_u_t>>
                    _scalar_properties;
                std::unordered_map<
                    std::string,
                    NCPA::interpolation::Interpolator1D<double, double>>
                    _scalar_splines;
                NCPA::interpolation::interpolator_1d_type_t
                    _1d_interpolator_type
                    = NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR;
                vector_u_t _change_points;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap(
    NCPA::atmos::piecewise_stratified_atmosphere_2d& a,
    NCPA::atmos::piecewise_stratified_atmosphere_2d& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::atmos::abstract_atmosphere_2d&>( a ),
            dynamic_cast<NCPA::atmos::abstract_atmosphere_2d&>( b ) );
    swap( a._components, b._components );
    swap( a._scalar_properties, b._scalar_properties );
    swap( a._scalar_splines, b._scalar_splines );
    swap( a._axes, b._axes );
    swap( a._change_points, b._change_points );
    swap( a._1d_interpolator_type, b._1d_interpolator_type );
}
