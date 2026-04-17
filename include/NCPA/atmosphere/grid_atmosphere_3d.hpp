#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_3d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/AtmosphericModel.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/grid_atmospheric_property_3d.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/units.hpp"

#include <cfloat>
#include <memory>
#include <string>

static void swap( NCPA::atmos::grid_atmosphere_3d&,
                  NCPA::atmos::grid_atmosphere_3d& ) noexcept;

namespace NCPA {
    namespace atmos {
        class grid_atmosphere_3d : public abstract_atmosphere_3d {
            public:
                grid_atmosphere_3d() : abstract_atmosphere_3d() {
                    for (auto it = _axes.begin(); it != _axes.end(); ++it) {
                        it->set_units(
                            NCPA_ATMOSPHERE_DEFAULT_DISTANCE_UNITS );
                    }
                }

                grid_atmosphere_3d( const grid_atmosphere_3d& other ) :
                    grid_atmosphere_3d() {
                    _properties           = other._properties;
                    _scalar_properties    = other._scalar_properties;
                    _interpolator_type    = other._interpolator_type;
                    _2d_interpolator_type = other._2d_interpolator_type;
                    _axes                 = other._axes;
                    _build_scalar_splines();
                }

                grid_atmosphere_3d( grid_atmosphere_3d&& source ) noexcept : grid_atmosphere_3d() { ::swap( *this, source ); } virtual ~grid_atmosphere_3d() {} grid_atmosphere_3d& operator=( grid_atmosphere_3d other ) { ::swap( *this, other ); return *this; } friend void ::swap( grid_atmosphere_3d& a, grid_atmosphere_3d& b ) noexcept; virtual std::unique_ptr<abstract_atmosphere_3d> clone() const override { return std::unique_ptr<abstract_atmosphere_3d>( new grid_atmosphere_3d( *this ) ); }

                virtual grid_atmosphere_3d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    NCPA::arrays::ndvector<2, abstract_atmosphere_1d *>&
                        components ) override {
                    _axes[ 0 ] = ax1;
                    _axes[ 1 ] = ax2;
                    return this->set( components );
                }

                virtual grid_atmosphere_3d& set(
                    NCPA::arrays::ndvector<2, abstract_atmosphere_1d *>&
                        components ) override {
                    auto dims = components.shape();
                    if (dims[ 0 ] != _axes[ 0 ].size()
                        || dims[ 1 ] != _axes[ 1 ].size()) {
                        throw std::range_error(
                            "Size mismatch between axis dimensions and "
                            "supplied components vector" );
                    }
                    if (_axes[ 2 ].size() == 0) {
                        _axes[ 2 ] = components[ 0 ][ 0 ]->get_axis_vector();
                    } else {
                        _axes[ 2 ].convert_units(
                            *components[ 0 ][ 0 ]->get_axis_units() );
                    }
                    _properties.clear();
                    _scalar_properties.clear();
                    _scalar_splines.clear();
                    auto vkeys = components[ 0 ][ 0 ]->get_vector_keys();
                    auto skeys = components[ 0 ][ 0 ]->get_scalar_keys();
                    for (auto vit = vkeys.cbegin(); vit != vkeys.cend();
                         ++vit) {
                        vector3d_u_t propvec(
                            _axes[ 0 ].size(), _axes[ 1 ].size(),
                            _axes[ 2 ].size(),
                            components[ 0 ][ 0 ]->get_units( *vit ) );
                        for (size_t i = 0; i < dims[ 0 ]; ++i) {
                            for (size_t j = 0; j < dims[ 1 ]; ++j) {
                                for (size_t k = 0; k < _axes[ 2 ].size();
                                     ++k) {
                                    propvec[ i ][ j ][ k ]
                                        = components[ i ][ j ]->get(
                                            *vit, _axes[ 2 ][ k ] );
                                }
                            }
                        }
                        _properties.emplace( std::make_pair(
                            *vit, grid_atmospheric_property_3d(
                                      _axes[ 0 ], _axes[ 1 ], _axes[ 2 ],
                                      propvec ) ) );
                        // _properties.emplace( *vit, _axes[0], _axes[1],
                        // _axes[2], propvec );
                    }
                    for (auto sit = skeys.cbegin(); sit != skeys.cend();
                         ++sit) {
                        vector2d_u_t propvec(
                            dims[ 0 ], dims[ 1 ],
                            components[ 0 ][ 0 ]->get_units( *sit ) );
                        for (size_t i = 0; i < dims[ 0 ]; ++i) {
                            for (size_t j = 0; j < dims[ 1 ]; ++j) {
                                propvec[ i ][ j ]
                                    = components[ i ][ j ]->get( *sit );
                            }
                        }
                        _scalar_properties.emplace(
                            std::make_pair( *sit, propvec ) );
                    }
                    this->set_interpolator( _interpolator_type );
                    _build_scalar_splines();
                    return *this;;
                }

                virtual size_t size( size_t dim ) const override {
                    this->validate_axis( dim );
                    return _axes[ dim ].size();
                }

                // set all horizontal points to the same 1-D atmosphere
                virtual grid_atmosphere_3d& set(
                    abstract_atmosphere_1d& atmos1d ) override {
                    auto keys = atmos1d.get_vector_keys();
                    _properties.clear();
                    std::unordered_map<std::string, vector3d_u_t> newprops;
                    for (auto it = keys.begin(); it != keys.end(); ++it) {
                        auto prop = atmos1d.get_property( *it );
                        prop.convert_axis_units( *_axes[ 2 ].get_units() );
                        vector3d_u_t v3d( _axes[ 0 ].size(), _axes[ 1 ].size(),
                                          _axes[ 2 ].size(),
                                          prop.get_units() );
                        for (size_t i = 0; i < _axes[ 0 ].size(); ++i) {
                            for (size_t j = 0; j < _axes[ 1 ].size(); ++j) {
                                v3d[ i ][ j ] = prop.values();
                            }
                        }
                        _properties[ *it ] = grid_atmospheric_property_3d(
                            _axes[ 0 ].size(), _axes[ 1 ].size(),
                            _axes[ 2 ].size(), v3d );
                    }
                    _scalar_properties.clear();
                    keys = atmos1d.get_scalar_keys();
                    for (auto it = keys.begin(); it != keys.end(); ++it) {
                        _scalar_properties[ *it ] = vector2d_u_t(
                            _axes[ 0 ].size(), _axes[ 1 ].size(),
                            scalar_u_t( atmos1d.get( *it ),
                                        atmos1d.get_units( *it ) ) );
                    }
                    this->set_interpolator( _interpolator_type );
                    _build_scalar_splines();
                    return *this;;
                }

                virtual grid_atmosphere_3d& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type )
                    override {
                    _interpolator_type = interp_type;
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.set_interpolator( interp_type );
                    }
                    return *this;;
                }

                virtual grid_atmosphere_3d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    _2d_interpolator_type = interp_type;
                    _build_scalar_splines();
                    return *this;;
                }

                virtual grid_atmosphere_3d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Cannot set a 1-D interpolator on a 3-D atmosphere" );
                }

                virtual grid_atmosphere_3d& set_axis(
                    size_t axis, vector_u_t vals ) override {
                    this->validate_axis( axis );
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.resample( axis, vals );
                    }
                    _axes[ axis ] = vals;
                    return *this;;
                }

                virtual grid_atmosphere_3d& add_property(
                    const std::string& key,
                    const AtmosphericProperty3D& property ) override {
                    _assert_does_not_contain( key );
                    _properties[ key ] = grid_atmospheric_property_3d();
                    _properties[ key ].copy( *property.internal() );
                    if (!_axes[ 0 ].empty() && !_axes[ 1 ].empty()
                        && !_axes[ 2 ].empty()) {
                        _properties[ key ].resample( _axes[ 0 ], _axes[ 1 ],
                                                     _axes[ 2 ] );
                    } else {
                        for (size_t i = 0; i < 3; ++i) {
                            _axes[ i ] = _properties[ key ].axis( i );
                        }
                    }
                    _properties[ key ].set_interpolator( _interpolator_type );
                    return *this;;
                }

                virtual grid_atmosphere_3d& add_property(
                    const std::string& key,
                    const vector3d_u_t& property ) override {
                    _assert_does_not_contain( key );
                    for (size_t ax = 0; ax < 3; ++ax) {
                        _assert_axis_is_set( ax );
                        if (_axes[ ax ].size() != property.dim( ax )) {
                            std::ostringstream oss;
                            oss << "Size mismatch between internal size "
                                   "vectors and "
                                   "provided property: axis "
                                << ax;
                            throw std::range_error( oss.str() );
                        }
                    }
                    _properties[ key ] = grid_atmospheric_property_3d();
                    _properties[ key ].set( _axes[ 0 ], _axes[ 1 ], _axes[ 2 ],
                                            property );
                    _properties[ key ].set_interpolator( _interpolator_type );
                    return *this;;
                }

                virtual grid_atmosphere_3d& add_property(
                    const std::string& key,
                    const vector2d_u_t& property ) override {
                    _assert_does_not_contain( key );
                    for (size_t ax = 0; ax < 2; ++ax) {
                        _assert_axis_is_set( ax );
                        if (_axes[ ax ].size() != property.dim( ax )) {
                            std::ostringstream oss;
                            oss << "Size mismatch between internal size "
                                   "vectors and "
                                   "provided property: axis "
                                << ax;
                            throw std::range_error( oss.str() );
                        }
                    }
                    _scalar_properties[ key ] = property;
                    _build_scalar_splines();
                    return *this;;
                }

                virtual grid_atmosphere_3d& add_property(
                    const std::string& key,
                    const scalar_u_t& property ) override {
                    _assert_does_not_contain( key );
                    if (_axes[ 0 ].empty() || _axes[ 1 ].empty()) {
                        throw std::range_error(
                            "No axes set to apply scalar property!" );
                    }
                    _scalar_properties.emplace( std::make_pair(
                        key, vector2d_u_t( _axes[ 0 ].size(),
                                           _axes[ 1 ].size(), property ) ) );
                    _build_scalar_splines();
                    return *this;;
                }

                virtual grid_atmosphere_3d& remove_property(
                    const std::string& key ) override {
                    _properties.erase( key );
                    _scalar_properties.erase( key );
                    _scalar_splines.erase( key );
                    return *this;;
                }

                virtual grid_atmosphere_3d& copy_property(
                    const std::string& old_key,
                    const std::string& new_key ) override {
                    _assert_does_not_contain( new_key );
                    if (contains_vector( old_key )) {
                        _properties[ new_key ] = _properties[ old_key ];
                    } else if (contains_scalar( old_key )) {
                        _scalar_properties[ new_key ]
                            = _scalar_properties[ old_key ];
                        _build_scalar_splines();
                    } else {
                        throw std::invalid_argument(
                            "Key " + old_key
                            + " does not exist in atmosphere!" );
                    }
                    return *this;;
                }

                virtual vector_u_t get_axis_vector( size_t n ) override {
                    this->validate_axis( n );
                    return _axes[ n ];
                }

                virtual vector_u_t get_axis_vector( size_t n ) const override {
                    this->validate_axis( n );
                    return _axes[ n ];
                }

                virtual vector3d_u_t get_values(
                    const std::string& key ) override {
                    return _properties.at( key ).values();
                }

                virtual vector3d_u_t get_values(
                    const std::string& key ) const override {
                    return _properties.at( key ).values();
                }

                virtual double get( const std::string& key, double x1,
                                    double x2 ) override {
                    return _scalar_splines.at( key ).eval_f( x1, x2 );
                }

                virtual double get( const std::string& key, double x1,
                                    double x2, double x3 ) override {
                    return _properties.at( key ).get( x1, x2, x3 );
                }

                virtual vector3d_u_t get(
                    const std::string& key, const std::vector<double>& v1,
                    const std::vector<double>& v2,
                    const std::vector<double>& v3 ) override {
                    return _properties.at( key ).get( v1, v2, v3 );
                }

                virtual vector2d_u_t get(
                    const std::string& key, const std::vector<double>& x1,
                    const std::vector<double>& x2 ) override {
                    vector2d_u_t vout( x1.size(), x2.size(),
                                       this->get_units( key ) );
                    for (size_t i = 0; i < x1.size(); ++i) {
                        for (size_t j = 0; j < x2.size(); ++j) {
                            vout[ i ][ j ] = _scalar_splines.at( key ).eval_f(
                                x1[ i ], x2[ j ] );
                        }
                    }
                    return vout;
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double x1, double x2,
                                                     size_t wrt ) override {
                    if (wrt != 2) {
                        return _scalar_splines.at( key ).eval_df( x1, x2,
                                                                  wrt );
                    } else {
                        return 0.0;
                    }
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double x1, double x2,
                                                     double x3,
                                                     size_t wrt ) override {
                    return _properties.at( key ).get_first_derivative(
                        x1, x2, x3, wrt );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double x1, double x2,
                                                      double x3, size_t wrt1,
                                                      size_t wrt2 ) override {
                    return _properties.at( key ).get_second_derivative(
                        x1, x2, x3, wrt1, wrt2 );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double x1, double x2,
                                                      size_t wrt1,
                                                      size_t wrt2 ) override {
                    if (wrt1 != 2 && wrt2 != 2) {
                        return _scalar_splines.at( key ).eval_ddf(
                            x1, x2, wrt1, wrt2 );
                    } else {
                        return 0.0;
                    }
                }

                virtual units_ptr_t get_axis_units( size_t n ) const override {
                    this->validate_axis( n );
                    return _axes[ n ].get_units();
                }

                virtual units_ptr_t get_units(
                    const std::string& key ) const override {
                    if (this->contains_vector( key )) {
                        return _properties.at( key ).get_units();
                    } else {
                        _assert_contains_scalar( key );
                        return _scalar_properties.at( key ).get_units();
                    }
                }

                virtual double get_minimum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if (n == 0) {
                        return 0.0;
                    }
                    double minax = _axes[ n ].back();
                    for (auto it = _properties.cbegin();
                         it != _properties.cend(); ++it) {
                        minax = std::min( minax,
                                          it->second.get_limits( n ).first );
                    }
                    return minax;
                }

                virtual double get_maximum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if (n == 0) {
                        return 0.0;
                    }
                    double maxax = _axes[ n ].front();
                    for (auto it = _properties.cbegin();
                         it != _properties.cend(); ++it) {
                        maxax = std::max( maxax,
                                          it->second.get_limits( n ).second );
                    }
                    return maxax;
                }

                virtual grid_atmosphere_3d& convert_axis_units(
                    size_t n, units_ptr_t new_units ) override {
                    this->validate_axis( n );
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.convert_axis_units( n, *new_units );
                    }
                    _axes[ n ].convert_units( *new_units );
                    _build_scalar_splines();
                    return *this;;
                }

                virtual grid_atmosphere_3d& convert_units(
                    const std::string& key, units_ptr_t new_units ) override {
                    if (this->contains_vector( key )) {
                        _properties[ key ].convert_units( *new_units );
                    } else if (this->contains_scalar( key )) {
                        _scalar_properties[ key ].convert_units( *new_units );
                    } else {
                        throw std::range_error( "Unknown key: " + key );
                    }
                    return *this;;
                }

                virtual grid_atmosphere_3d& resample(
                    size_t axis, double new_d ) override {
                    this->validate_axis( axis );
                    vector_u_t new_x;
                    new_x.set_units( *_axes[ axis ].get_units() );
                    double x = _axes[ axis ].front();
                    while (x <= _axes[ axis ].back()) {
                        new_x.push_back( x );
                        x += new_d;
                    }
                    return this->resample( axis, new_x );
                }

                virtual grid_atmosphere_3d& resample(
                    size_t axis, vector_u_t new_z ) override {
                    this->validate_axis( axis );
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.resample( axis, new_z );
                    }

                    if (axis <= 1) {
                        _axes[ axis ].convert_units( *new_z.get_units() );
                        auto scalar_copy = _scalar_properties;
                        std::array<std::vector<double>, 2> newdims {
                            _axes[ 0 ], _axes[ 1 ]
                        };
                        newdims[ axis ] = new_z;
                        for (auto it = scalar_copy.begin();
                             it != scalar_copy.end(); ++it) {
                            it->second.reshape(
                                { newdims[ 0 ].size(), newdims[ 1 ].size() } );
                            for (size_t i = 0; i < newdims[ 0 ].size(); ++i) {
                                for (size_t j = 0; j < newdims[ 1 ].size();
                                     ++j) {
                                    it->second[ i ][ j ]
                                        = _scalar_splines[ it->first ].eval_f(
                                            newdims[ 0 ][ i ],
                                            newdims[ 1 ][ j ] );
                                }
                            }
                        }
                        _scalar_properties = scalar_copy;
                    }
                    _axes[ axis ] = new_z;
                    if (axis <= 1) {
                        _build_scalar_splines();
                    }
                    return *this;;
                }

                virtual grid_atmosphere_3d& resample(
                    vector_u_t new_ax1, vector_u_t new_ax2,
                    vector_u_t new_ax3 ) override {
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.resample( new_ax1, new_ax2, new_ax3 );
                    }
                    _axes[ 0 ].convert_units( *new_ax1.get_units() );
                    _axes[ 1 ].convert_units( *new_ax2.get_units() );

                    auto scalar_copy = _scalar_properties;
                    std::array<std::vector<double>, 2> newdims { new_ax1,
                                                                 new_ax2 };
                    for (auto it = scalar_copy.begin();
                         it != scalar_copy.end(); ++it) {
                        it->second.reshape(
                            { newdims[ 0 ].size(), newdims[ 1 ].size() } );
                        for (size_t i = 0; i < newdims[ 0 ].size(); ++i) {
                            for (size_t j = 0; j < newdims[ 1 ].size(); ++j) {
                                it->second[ i ][ j ]
                                    = _scalar_splines[ it->first ].eval_f(
                                        newdims[ 0 ][ i ], newdims[ 1 ][ j ] );
                            }
                        }
                    }
                    _scalar_properties = scalar_copy;
                    _axes              = { new_ax1, new_ax2, new_ax3 };
                    _build_scalar_splines();
                    return *this;;
                }

                virtual std::vector<std::string> get_keys() const override {
                    std::vector<std::string> vkeys = get_vector_keys(),
                                             skeys = get_scalar_keys(), jkeys;
                    std::set_intersection( vkeys.cbegin(), vkeys.cend(),
                                           skeys.cbegin(), skeys.cend(),
                                           std::back_inserter( jkeys ) );
                    return jkeys;
                }

                virtual std::vector<std::string> get_vector_keys()
                    const override {
                    std::vector<std::string> vkeys;
                    for (auto it = _properties.cbegin();
                         it != _properties.cend(); ++it) {
                        vkeys.push_back( it->first );
                    }
                    return vkeys;
                }

                virtual std::vector<std::string> get_scalar_keys()
                    const override {
                    std::vector<std::string> skeys;
                    for (auto it = _scalar_properties.cbegin();
                         it != _scalar_properties.cend(); ++it) {
                        skeys.push_back( it->first );
                    }
                    return skeys;
                }

                virtual bool contains_scalar(
                    const std::string& key ) const override {
                    return _scalar_properties.count( key ) == 1;
                }

                virtual bool contains_vector(
                    const std::string& key ) const override {
                    return _properties.count( key ) == 1;
                }

                virtual bool contains_key(
                    const std::string& key ) const override {
                    return contains_vector( key ) || contains_scalar( key );
                }

                virtual bool is_stratified() const override { return false; }

                virtual void print( std::ostream& os ) override {
                    throw std::runtime_error(
                        "Print function not yet implemented." );
                }

                virtual bool same( scalar_u_t x1_1, scalar_u_t x2_1,
                                   scalar_u_t x1_2,
                                   scalar_u_t x2_2 ) const override {
                    return false;
                }

                virtual bool same( double x1_1, double x2_1, double x1_2,
                                   double x2_2 ) const override {
                    return false;
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

                void _assert_axis_is_set( size_t i ) const {
                    if (_axes[ i ].size() == 0) {
                        std::ostringstream oss;
                        oss << "Axis " << i << " not set!";
                        throw std::logic_error( oss.str() );
                    }
                }

                void _build_scalar_splines() {
                    _scalar_splines.clear();
                    for (auto propit = _scalar_properties.cbegin();
                         propit != _scalar_properties.cend(); ++propit) {
                        _scalar_splines[ propit->first ]
                            = NCPA::interpolation::
                                InterpolatorFactory<double, double>::build(
                                    _2d_interpolator_type );
                        _scalar_splines[ propit->first ]
                            .init( _axes[ 0 ].size(), _axes[ 1 ].size() )
                            .fill( _axes[ 0 ], _axes[ 1 ], propit->second )
                            .ready();
                    }
                }

                virtual abstract_atmospheric_property_3d& _property(
                    const std::string& key ) {
                    return _properties.at( key );
                }

                std::unordered_map<std::string, grid_atmospheric_property_3d>
                    _properties;
                std::unordered_map<std::string, vector2d_u_t>
                    _scalar_properties;
                std::unordered_map<
                    std::string,
                    NCPA::interpolation::Interpolator2D<double, double>>
                    _scalar_splines;
                std::array<vector_u_t, 3> _axes;
                NCPA::interpolation::interpolator_3d_type_t _interpolator_type
                    = NCPA_ATMOSPHERE_DEFAULT_3D_INTERPOLATOR;
                NCPA::interpolation::interpolator_2d_type_t
                    _2d_interpolator_type
                    = NCPA_ATMOSPHERE_DEFAULT_2D_INTERPOLATOR;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::grid_atmosphere_3d& a,
                  NCPA::atmos::grid_atmosphere_3d& b ) noexcept {
    using std::swap;
    swap( a._properties, b._properties );
    swap( a._scalar_properties, b._scalar_properties );
    swap( a._scalar_splines, b._scalar_splines );
    swap( a._axes, b._axes );
    swap( a._interpolator_type, b._interpolator_type );
    swap( a._2d_interpolator_type, b._2d_interpolator_type );
}
