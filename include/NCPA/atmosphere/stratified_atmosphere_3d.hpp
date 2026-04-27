#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_3d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/AtmosphericModel.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/tuple_atmospheric_property_1d.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/units.hpp"

#include <cfloat>
#include <memory>
#include <string>

static void swap( NCPA::atmos::stratified_atmosphere_3d&,
                  NCPA::atmos::stratified_atmosphere_3d& ) noexcept;

namespace NCPA {
    namespace atmos {
        class stratified_atmosphere_3d : public abstract_atmosphere_3d {
            public:
                stratified_atmosphere_3d() : abstract_atmosphere_3d() {
                    // for (auto it = _axes.begin(); it != _axes.end(); ++it) {
                    //     it->set_units(
                    //         NCPA_ATMOSPHERE_DEFAULT_DISTANCE_UNITS );
                    // }
                    this->set_axis( 0, NCPA::atmos::vector_u_t( 
                        std::vector<double>{ 0.0, 1.0e15 }, NCPA_ATMOSPHERE_DEFAULT_DISTANCE_UNITS )
                    );
                    this->set_axis( 1, NCPA::atmos::vector_u_t( 
                        std::vector<double>{ 0.0, 1.0e15 }, NCPA_ATMOSPHERE_DEFAULT_DISTANCE_UNITS )
                    );
                    _axes[2].set_units( NCPA_ATMOSPHERE_DEFAULT_DISTANCE_UNITS );
                }

                stratified_atmosphere_3d(
                    const stratified_atmosphere_3d& other ) :
                    stratified_atmosphere_3d() {
                    _1d_atmos             = other._1d_atmos;
                    _2d_scalar_properties = other._2d_scalar_properties;
                    _interpolator_type    = other._interpolator_type;
                    _axes                 = other._axes;
                    _build_scalar_splines();
                }

                stratified_atmosphere_3d( stratified_atmosphere_3d&& source ) noexcept : stratified_atmosphere_3d() { ::swap( *this, source ); } virtual ~stratified_atmosphere_3d() {} stratified_atmosphere_3d& operator=( stratified_atmosphere_3d other ) { ::swap( *this, other ); return *this; } friend void ::swap( stratified_atmosphere_3d& a, stratified_atmosphere_3d& b ) noexcept; virtual std::unique_ptr<abstract_atmosphere_3d> clone() const override { return std::unique_ptr<abstract_atmosphere_3d>( new stratified_atmosphere_3d( *this ) ); }

                virtual size_t size( size_t dim ) const override {
                    this->validate_axis( dim );
                    return _axes[ dim ].size();
                }

                // set all horizontal points to the same 1-D atmosphere
                virtual stratified_atmosphere_3d& set(
                    abstract_atmosphere_1d& atmos1d ) override {
                    for (auto it = _2d_scalar_properties.cbegin();
                         it != _2d_scalar_properties.cend(); ++it) {
                        if (atmos1d.contains_key( it->first )) {
                            throw std::logic_error(
                                "Profile already contains range-dependent "
                                "scalar property "
                                + it->first );
                        }
                    }
                    _1d_atmos = Atmosphere1D( atmos1d.clone() );
                    _1d_atmos.set_interpolator( _interpolator_type );
                    this->set_axis( 2, _1d_atmos.get_axis_vector() );
                    return *this;;
                }

                virtual stratified_atmosphere_3d& set(
                    const vector_u_t& ax1, const vector_u_t& ax2,
                    NCPA::arrays::ndvector<2, abstract_atmosphere_1d *>&
                        components ) override {
                    _axes[ 0 ] = ax1;
                    _axes[ 1 ] = ax2;
                    return this->set( components );
                }

                virtual stratified_atmosphere_3d& set(
                    NCPA::arrays::ndvector<2, abstract_atmosphere_1d *>&
                        components ) override {
                    auto dims = components.shape();
                    if (dims[ 0 ] != _axes[ 0 ].size()
                        || dims[ 1 ] != _axes[ 1 ].size()) {
                        throw std::range_error(
                            "Size mismatch between axis dimensions and "
                            "supplied components vector" );
                    }
                    _1d_atmos  = Atmosphere1D( components[ 0 ][ 0 ]->clone() );
                    _axes[ 2 ] = components[ 0 ][ 0 ]->get_axis_vector();
                    _2d_scalar_properties.clear();
                    _scalar_splines.clear();
                    auto skeys = components[ 0 ][ 0 ]->get_scalar_keys();
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
                        _2d_scalar_properties.emplace(
                            std::make_pair( *sit, propvec ) );
                    }
                    _build_scalar_splines();
                    return *this;;
                }

                virtual stratified_atmosphere_3d& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "Cannot set a 3-D interpolator on a stratified "
                        "atmosphere" );
                }

                virtual stratified_atmosphere_3d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    _2d_interpolator_type = interp_type;
                    _build_scalar_splines();
                    return *this;;
                }

                virtual stratified_atmosphere_3d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _interpolator_type = interp_type;
                    _1d_atmos.set_interpolator( interp_type );
                    return *this;;
                }

                virtual stratified_atmosphere_3d& set_axis(
                    size_t axis, vector_u_t vals ) override {
                    this->validate_axis( axis );
                    if (axis == 2) {
                        _1d_atmos.resample( vals );
                        _axes[ 2 ] = vals;
                    } else {
                        _axes[ axis ].set_units( *vals.get_units() );
                    }
                    return *this;;
                }

                virtual stratified_atmosphere_3d& add_property(
                    const std::string& key,
                    const AtmosphericProperty3D& property ) override {
                    _assert_does_not_contain( key );
                    tuple_atmospheric_property_1d prop1d;
                    AtmosphericProperty3D propcopy = property;
                    propcopy.convert_axis_units( 2, *_axes[ 2 ].get_units() );
                    if (_axes[ 2 ].empty()) {
                        _axes[ 2 ] = propcopy.axis( 2 );
                    }
                    vector_u_t propvec( _axes[ 2 ].size(),
                                        property.get_units() );
                    for (size_t i = 0; i < _axes[ 2 ].size(); ++i) {
                        propvec[ i ]
                            = propcopy.get( 0.0, 0.0, _axes[ 2 ][ i ] );
                    }
                    prop1d.set( _axes[ 2 ], propvec );
                    return *this;;
                }

                virtual stratified_atmosphere_3d& add_property(
                    const std::string& key,
                    const vector3d_u_t& property ) override {
                    _assert_does_not_contain( key );
                    _assert_axis_is_set( 2 );
                    if (_axes[ 2 ].size() != property.dim( 2 )) {
                        std::ostringstream oss;
                        oss << "Size mismatch between internal size "
                               "vectors and "
                               "provided property: axis 2";
                        throw std::range_error( oss.str() );
                    }
                    vector_u_t v1d( property[ 0 ][ 0 ], property.get_units() );
                    tuple_atmospheric_property_1d prop1d( _axes[ 2 ], v1d );
                    _1d_atmos.add_property( key, prop1d );
                    return *this;;
                }

                virtual stratified_atmosphere_3d& add_property(
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
                    _2d_scalar_properties[ key ] = property;
                    _build_scalar_splines();
                    return *this;;
                }

                virtual stratified_atmosphere_3d& add_property(
                    const std::string& key,
                    const scalar_u_t& property ) override {
                    _assert_does_not_contain( key );
                    if (_axes[ 0 ].empty() || _axes[ 1 ].empty()) {
                        throw std::range_error(
                            "No axes set to apply scalar property!" );
                    }
                    _2d_scalar_properties.emplace( std::make_pair(
                        key, vector2d_u_t( _axes[ 0 ].size(),
                                           _axes[ 1 ].size(), property ) ) );
                    _build_scalar_splines();
                    return *this;;
                }

                virtual stratified_atmosphere_3d& remove_property(
                    const std::string& key ) override {
                    _1d_atmos.remove_property( key );
                    _2d_scalar_properties.erase( key );
                    _scalar_splines.erase( key );
                    return *this;;
                }

                virtual stratified_atmosphere_3d& copy_property(
                    const std::string& old_key,
                    const std::string& new_key ) override {
                    _assert_does_not_contain( new_key );
                    if (contains_vector( old_key )) {
                        _1d_atmos.copy_property( old_key, new_key );
                    } else if (contains_scalar( old_key )) {
                        _2d_scalar_properties[ new_key ]
                            = _2d_scalar_properties[ old_key ];
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
                    _assert_contains_vector( key );
                    vector3d_u_t vals( _axes[ 0 ].size(), _axes[ 1 ].size(),
                                       _axes[ 2 ].size(),
                                       this->get_units( key ) );
                    for (size_t i = 0; i < _axes[ 0 ].size(); ++i) {
                        for (size_t j = 0; j < _axes[ 1 ].size(); ++j) {
                            vals[ i ][ j ]
                                = _1d_atmos.get_property( key ).values();
                        }
                    }
                    return vals;
                }

                virtual vector3d_u_t get_values(
                    const std::string& key ) const override {
                    _assert_contains_vector( key );
                    vector3d_u_t vals( _axes[ 0 ].size(), _axes[ 1 ].size(),
                                       _axes[ 2 ].size(),
                                       this->get_units( key ) );
                    for (size_t i = 0; i < _axes[ 0 ].size(); ++i) {
                        for (size_t j = 0; j < _axes[ 1 ].size(); ++j) {
                            vals[ i ][ j ]
                                = _1d_atmos.get_property( key ).values();
                        }
                    }
                    return vals;
                }

                virtual double get( const std::string& key, double x1,
                                    double x2 ) override {
                    // return _scalar_splines.at( key ).eval_f( x1, x2 );
                    return _get_scalar( key, x1, x2 );
                }

                virtual double get( const std::string& key, double x1,
                                    double x2, double x3 ) override {
                    return _1d_atmos.get( key, x3 );
                }

                virtual vector3d_u_t get(
                    const std::string& key, const std::vector<double>& v1,
                    const std::vector<double>& v2,
                    const std::vector<double>& v3 ) override {
                    vector3d_u_t out3d( v1.size(), v2.size(), v3.size(),
                                        this->get_units( key ) );
                    for (size_t k = 0; k < v3.size(); ++k) {
                        double v = _1d_atmos.get( key, v3[ k ] );
                        for (size_t i = 0; i < v1.size(); ++i) {
                            for (size_t j = 0; j < v2.size(); ++j) {
                                out3d[ i ][ j ][ k ] = v;
                            }
                        }
                    }
                    return out3d;
                }

                virtual vector2d_u_t get(
                    const std::string& key, const std::vector<double>& x1,
                    const std::vector<double>& x2 ) override {
                    vector2d_u_t vout( x1.size(), x2.size(),
                                       this->get_units( key ) );
                    for (size_t i = 0; i < x1.size(); ++i) {
                        for (size_t j = 0; j < x2.size(); ++j) {
                            // vout[ i ][ j ] = _scalar_splines.at( key
                            // ).eval_f(
                            //     x1[ i ], x2[ j ] );
                            vout[ i ][ j ]
                                = _get_scalar( key, x1[ i ], x2[ j ] );
                        }
                    }
                    return vout;
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double x1, double x2,
                                                     double x3,
                                                     size_t wrt ) override {
                    if (wrt == 2) {
                        return _1d_atmos.get_first_derivative( key, x3 );
                    } else {
                        return 0.0;
                    }
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double x1, double x2,
                                                     size_t wrt ) override {
                    if (wrt != 2 && contains_2d_scalar( key )) {
                        return _scalar_splines.at( key ).eval_df( x1, x2,
                                                                  wrt );
                    } else {
                        return 0.0;
                    }
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double x1, double x2,
                                                      double x3, size_t wrt1,
                                                      size_t wrt2 ) override {
                    if (wrt1 == 2 && wrt2 == 2) {
                        return _1d_atmos.get_second_derivative( key, x3 );
                    } else {
                        return 0.0;
                    }
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double x1, double x2,
                                                      size_t wrt1,
                                                      size_t wrt2 ) override {
                    if (wrt1 != 2 && wrt2 != 2 && contains_2d_scalar( key )) {
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
                        return _1d_atmos.get_units( key );
                    } else {
                        _assert_contains_scalar( key );
                        if (contains_1d_scalar( key )) {
                            return _1d_atmos.get_units( key );
                        } else {
                            return _2d_scalar_properties.at( key ).get_units();
                        }
                    }
                }

                virtual double get_minimum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if (n != 2) {
                        return 0.0;
                    } else {
                        return _1d_atmos.get_minimum_axis();
                    }
                }

                virtual double get_maximum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if (n != 2) {
                        return DBL_MAX;
                    } else {
                        return _1d_atmos.get_maximum_axis();
                    }
                }

                virtual bool is_stratified() const override { return true; }

                virtual stratified_atmosphere_3d& convert_axis_units(
                    size_t n, units_ptr_t new_units ) override {
                    this->validate_axis( n );
                    if (n == 2) {
                        _1d_atmos.convert_axis_units( new_units );
                    }
                    _axes[ n ].convert_units( *new_units );
                    _build_scalar_splines();
                    return *this;;
                }

                virtual stratified_atmosphere_3d& convert_units(
                    const std::string& key, units_ptr_t new_units ) override {
                    if (this->contains_vector( key )) {
                        _1d_atmos.convert_units( key, new_units );
                    } else if (this->contains_2d_scalar( key )) {
                        _2d_scalar_properties[ key ].convert_units(
                            *new_units );
                        _build_scalar_splines();
                    } else if (contains_1d_scalar( key )) {
                        _1d_atmos.convert_units( key, new_units );
                    } else {
                        throw std::range_error( "Unknown key: " + key );
                    }
                    return *this;;
                }

                virtual stratified_atmosphere_3d& resample(
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

                virtual stratified_atmosphere_3d& resample(
                    size_t axis, vector_u_t new_z ) override {
                    this->validate_axis( axis );
                    if (axis == 2) {
                        _1d_atmos.resample( new_z );
                        _axes[ 2 ] = new_z;
                    } else {
                        _axes[ axis ].convert_units( *new_z.get_units() );
                        auto scalar_copy = _2d_scalar_properties;
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
                        _2d_scalar_properties = scalar_copy;
                        _build_scalar_splines();
                    }
                    return *this;;
                }

                virtual stratified_atmosphere_3d& resample(
                    vector_u_t new_ax1, vector_u_t new_ax2,
                    vector_u_t new_ax3 ) override {
                    _1d_atmos.resample( new_ax3 );
                    _axes[ 2 ] = new_ax3;
                    _axes[ 0 ].convert_units( *new_ax1.get_units() );
                    _axes[ 1 ].convert_units( *new_ax2.get_units() );
                    auto scalar_copy = _2d_scalar_properties;
                    std::array<std::vector<double>, 2> newdims { _axes[ 0 ],
                                                                 _axes[ 1 ] };
                    newdims[ 0 ] = new_ax1;
                    newdims[ 1 ] = new_ax2;
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
                    _2d_scalar_properties = scalar_copy;
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
                    return _1d_atmos.get_vector_keys();
                }

                virtual std::vector<std::string> get_scalar_keys()
                    const override {
                    std::vector<std::string> skeys;
                    for (auto it = _2d_scalar_properties.cbegin();
                         it != _2d_scalar_properties.cend(); ++it) {
                        skeys.push_back( it->first );
                    }
                    std::vector<std::string> skeys1d
                        = _1d_atmos.get_scalar_keys();
                    skeys.insert( skeys.begin(), skeys1d.begin(),
                                  skeys1d.end() );
                    return skeys;
                }

                virtual bool contains_scalar(
                    const std::string& key ) const override {
                    return contains_1d_scalar( key )
                        || contains_2d_scalar( key );
                    // return _1d_atmos.contains_scalar( key )
                    //     || _2d_scalar_properties.count( key ) == 1;
                }

                virtual bool contains_1d_scalar(
                    const std::string& key ) const {
                    return _1d_atmos.contains_scalar( key );
                }

                virtual bool contains_2d_scalar(
                    const std::string& key ) const {
                    return _2d_scalar_properties.count( key ) == 1;
                }

                virtual bool contains_vector(
                    const std::string& key ) const override {
                    return _1d_atmos.contains_vector( key );
                }

                virtual bool contains_key(
                    const std::string& key ) const override {
                    return contains_vector( key ) || contains_scalar( key );
                }

                virtual void print( std::ostream& os ) override {
                    // throw std::runtime_error(
                    //     "Print function not yet implemented." );
                    _1d_atmos.print( os );
                }

                virtual bool same( scalar_u_t x1_1, scalar_u_t x2_1,
                                   scalar_u_t x1_2,
                                   scalar_u_t x2_2 ) const override {
                    return true;
                }

                virtual bool same( double x1_1, double x2_1, double x1_2,
                                   double x2_2 ) const override {
                    return true;
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
                    for (auto propit = _2d_scalar_properties.cbegin();
                         propit != _2d_scalar_properties.cend(); ++propit) {
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

                double _get_scalar( const std::string& key, double x1,
                                    double x2 ) {
                    if (contains_2d_scalar( key )) {
                        return _scalar_splines.at( key ).eval_f( x1, x2 );
                    } else {
                        return _1d_atmos.get( key );
                    }
                }

                Atmosphere1D _1d_atmos;
                std::unordered_map<std::string, vector2d_u_t>
                    _2d_scalar_properties;
                std::unordered_map<
                    std::string,
                    NCPA::interpolation::Interpolator2D<double, double>>
                    _scalar_splines;
                std::array<vector_u_t, 3> _axes;
                NCPA::interpolation::interpolator_1d_type_t _interpolator_type
                    = NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR;
                NCPA::interpolation::interpolator_2d_type_t
                    _2d_interpolator_type
                    = NCPA_ATMOSPHERE_DEFAULT_2D_INTERPOLATOR;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::stratified_atmosphere_3d& a,
                  NCPA::atmos::stratified_atmosphere_3d& b ) noexcept {
    using std::swap;
    swap( a._1d_atmos, b._1d_atmos );
    swap( a._2d_scalar_properties, b._2d_scalar_properties );
    swap( a._scalar_splines, b._scalar_splines );
    swap( a._axes, b._axes );
    swap( a._interpolator_type, b._interpolator_type );
    swap( a._2d_interpolator_type, b._2d_interpolator_type );
}
