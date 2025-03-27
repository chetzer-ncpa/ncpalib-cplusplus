#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_2d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/AtmosphericModel.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/units.hpp"

#include <memory>
#include <string>

static void swap( NCPA::atmos::stratified_atmosphere_2d&,
                  NCPA::atmos::stratified_atmosphere_2d& ) noexcept;
static void swap( NCPA::atmos::piecewise_stratified_atmosphere_2d&,
                  NCPA::atmos::piecewise_stratified_atmosphere_2d& ) noexcept;
static void swap( NCPA::atmos::grid_atmosphere_2d&,
                  NCPA::atmos::grid_atmosphere_2d& ) noexcept;
static void swap( NCPA::atmos::Atmosphere2D&,
                  NCPA::atmos::Atmosphere2D& ) noexcept;

namespace NCPA {
    namespace atmos {

        class grid_atmosphere_2d : public abstract_atmosphere_2d {
            public:
                grid_atmosphere_2d() : abstract_atmosphere_2d() {
                    _r.set_units( NCPA::units::KILOMETERS );
                    _z.set_units( NCPA::units::KILOMETERS );
                }

                grid_atmosphere_2d( const grid_atmosphere_2d& other ) :
                    grid_atmosphere_2d() {
                    _properties        = other._properties;
                    _scalar_properties = other._scalar_properties;
                    _z                 = other._z;
                    _interpolator_type = other._interpolator_type;
                    _r                 = other._r;
                    _internals         = other._internals;
                }

                DECLARE_BOILERPLATE_METHODS( grid_atmosphere_2d,
                                             abstract_atmosphere_2d )

                virtual std::unique_ptr<abstract_atmosphere_1d> extract(
                    double range ) override {
                    std::unique_ptr<abstract_atmosphere_1d> atm(
                        new tuple_atmosphere_1d() );
                    for (auto vit = _properties.begin();
                         vit != _properties.end(); ++vit) {
                        atm->add_property( vit->first,
                                           vit->second.extract( range ) );
                    }
                    for (auto sit = _scalar_splines.begin();
                         sit != _scalar_splines.end(); ++sit) {
                        atm->add_property(
                            sit->first,
                            scalar_u_t( sit->second.eval_f( range ),
                                        _scalar_properties[ sit->first ]
                                            .get_units() ) );
                    }
                    return atm;
                }

                virtual size_t size( size_t dim ) const override {
                    this->validate_axis( dim );
                    if (dim == 0) {
                        return _r.size();
                    } else {
                        return _z.size();
                    }
                }

                virtual abstract_atmosphere_2d& set(
                    abstract_atmosphere_1d& atmos1d ) override {
                    _r.resize( 1 );
                    _r[ 0 ] = 0.0;
                    _z      = atmos1d.get_axis_vector();
                    _internals.clear();
                    _internals.push_back( Atmosphere1D( atmos1d.clone() ) );
                    _properties.clear();
                    _scalar_properties.clear();
                    _translate_to_grid();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& set(
                    const vector_u_t ranges,
                    std::vector<abstract_atmosphere_1d *> components )
                    override {
                    _r = ranges;
                    _z = components[ 0 ]->get_axis_vector();
                    _internals.clear();
                    _properties.clear();
                    _scalar_properties.clear();
                    for (auto it = components.cbegin();
                         it != components.cend(); ++it) {
                        _internals.push_back( ( *it )->clone() );
                    }
                    _translate_to_grid();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& append(
                    scalar_u_t range,
                    abstract_atmosphere_1d& atmos1d ) override {
                    if (range.get_as( _r.get_units() ) <= _r.back()) {
                        throw std::range_error(
                            "Can't append inside existing range vector!" );
                    }
                    _properties.clear();
                    _scalar_properties.clear();
                    _r.push_back( range.get_as( _r.get_units() ) );
                    _internals.push_back( Atmosphere1D( atmos1d.clone() ) );
                    _internals.back().resample( _z );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    _interpolator_type = interp_type;
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.set_interpolator( interp_type );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    this->_translate_to_grid();
                    _1d_interpolator_type = interp_type;
                    _build_scalar_splines();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& set_axis(
                    size_t axis, vector_u_t vals ) override {
                    this->validate_axis( axis );
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.resample( axis, vals );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const AtmosphericProperty2D& property ) override {
                    _assert_does_not_contain( key );
                    _properties[ key ] = grid_atmospheric_property_2d();
                    _properties[ key ].copy( *property.internal() );
                    if (!_r.empty() && !_z.empty()) {
                        _properties[ key ].resample( _r, _z );
                    } else {
                        _r = _properties[ key ].axis( 0 );
                        _z = _properties[ key ].axis( 1 );
                    }
                    _properties[ key ].set_interpolator( _interpolator_type );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const vector2d_u_t& property ) override {
                    _assert_does_not_contain( key );
                    _assert_z_is_set();
                    if (_z.size() != property.dim( 1 )
                        || _r.size() != property.dim( 0 )) {
                        throw std::range_error(
                            "Size mismatch between internal size vectors and "
                            "provided property" );
                    }
                    _properties[ key ] = grid_atmospheric_property_2d();
                    _properties[ key ].set( _r, _z, property );
                    _properties[ key ].set_interpolator( _interpolator_type );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const vector_u_t& property ) override {
                    return this->add_property( key, property, _r );
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
                    if (index == _r) {
                        _scalar_properties[ key ] = property;
                    } else {
                        NCPA::interpolation::Interpolator1D<double, double>
                            interp = NCPA::interpolation::
                                InterpolatorFactory<double, double>::build(
                                    NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR );
                        interp.init( property.size() )
                            .fill( index.as( _r.get_units() ), property )
                            .ready();
                        vector_u_t newvec = property;
                        newvec.resize( _r.size() );
                        for (size_t i = 0; i < _r.size(); ++i) {
                            newvec[ i ] = interp.eval_f( _r[ i ] );
                        }
                        _scalar_properties[ key ] = newvec;
                    }
                    _build_scalar_splines();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const scalar_u_t& property ) override {
                    _assert_does_not_contain( key );
                    _scalar_properties[ key ]
                        = vector_u_t( _r.size(), property );
                    _build_scalar_splines();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& remove_property(
                    const std::string& key ) override {
                    _properties.erase( key );
                    _scalar_properties.erase( key );
                    _scalar_splines.erase( key );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& copy_property(
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual vector_u_t get_axis_vector( size_t n ) override {
                    this->validate_axis( n );
                    return ( n == 0 ? _r : _z );
                }

                virtual vector_u_t get_axis_vector( size_t n ) const override {
                    this->validate_axis( n );
                    return ( n == 0 ? _r : _z );
                }

                virtual double get( const std::string& key,
                                    double r ) override {
                    return _scalar_splines.at( key ).eval_f( r );
                }

                virtual double get( const std::string& key, double r,
                                    double z ) override {
                    return _properties.at( key ).get( r, z );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double r, double z,
                                                     size_t wrt ) override {
                    return _properties.at( key ).get_first_derivative( r, z,
                                                                       wrt );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double r, double z,
                                                      size_t wrt1,
                                                      size_t wrt2 ) override {
                    return _properties.at( key ).get_second_derivative(
                        r, z, wrt1, wrt2 );
                }

                virtual units_ptr_t get_axis_units( size_t n ) const override {
                    this->validate_axis( n );
                    return ( n == 0 ? _r.get_units() : _z.get_units() );
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
                    double minax = _z.back();
                    for (auto it = _properties.cbegin();
                         it != _properties.cend(); ++it) {
                        minax = std::min( minax,
                                          it->second.get_limits( 1 ).first );
                    }
                    return minax;
                }

                virtual double get_maximum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if (n == 0) {
                        return 0.0;
                    }
                    double maxax = _z.front();
                    for (auto it = _properties.cbegin();
                         it != _properties.cend(); ++it) {
                        maxax = std::max( maxax,
                                          it->second.get_limits( 1 ).second );
                    }
                    return maxax;
                }

                virtual abstract_atmosphere_2d& convert_axis_units(
                    size_t n, units_ptr_t new_units ) override {
                    this->validate_axis( n );
                    if (n == 0) {
                        _r.convert_units( *new_units );
                    } else {
                        for (auto it = _properties.begin();
                             it != _properties.end(); ++it) {
                            it->second.convert_axis_units( n, *new_units );
                        }
                        _z.convert_units( *new_units );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& convert_units(
                    const std::string& key, units_ptr_t new_units ) override {
                    if (this->contains_vector( key )) {
                        _properties[ key ].convert_units( *new_units );
                    } else if (this->contains_scalar( key )) {
                        _scalar_properties[ key ].convert_units( *new_units );
                    } else {
                        throw std::range_error( "Unknown key: " + key );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    size_t axis, double new_d ) override {
                    this->validate_axis( axis );
                    if (axis == 1) {
                        vector_u_t new_z;
                        new_z.set_units( *_z.get_units() );
                        double z = _z.front();
                        while (z <= _z.back()) {
                            new_z.push_back( z );
                            z += new_d;
                        }
                        return this->resample( axis, new_z );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    size_t axis, vector_u_t new_z ) override {
                    this->validate_axis( axis );
                    if (axis == 1) {
                        for (auto it = _properties.begin();
                             it != _properties.end(); ++it) {
                            it->second.resample( 1, new_z );
                        }
                        _z = new_z;
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    vector_u_t new_ax1, vector_u_t new_ax2 ) override {
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.resample( new_ax1, new_ax2 );
                    }
                    _z = new_ax2;
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
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

                virtual void print( std::ostream& os ) override {
                    throw std::runtime_error(
                        "Print function not yet implemented." );
                }

                virtual bool same( scalar_u_t r1,
                                   scalar_u_t r2 ) const override {
                    return ( r1 == r2 );
                }

                virtual bool same( double r1, double r2 ) const override {
                    return ( r1 == r2 );
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

                void _assert_z_is_set() const {
                    if (!_z) {
                        throw std::logic_error( "Altitude vector not set!" );
                    }
                }

                void _translate_to_grid() {
                    if (!( _properties.empty()
                           || _properties.begin()->second.size( 0 )
                                  != _r.size() )) {
                        return;
                    }
                    _properties.clear();
                    _scalar_properties.clear();

                    auto vector_keys = _internals[ 0 ].get_vector_keys();
                    for (auto vit = vector_keys.cbegin();
                         vit != vector_keys.cend(); ++vit) {
                        grid_atmospheric_property_2d temp;
                        temp.set(
                            _r.get_scalar( 0 ),
                            *_internals[ 0 ].get_property( *vit ).internal() );
                        for (size_t i = 1; i < _r.size(); ++i) {
                            temp.append( _r.get_scalar( i ),
                                         *_internals[ i ]
                                              .get_property( *vit )
                                              .internal() );
                        }
                        this->add_property( *vit, temp );
                    }

                    _scalar_properties.clear();
                    auto scalar_keys = _internals[ 0 ].get_scalar_keys();
                    for (auto sit = scalar_keys.cbegin();
                         sit != scalar_keys.cend(); ++sit) {
                        _scalar_properties[ *sit ].resize( _r.size() );
                        _scalar_properties[ *sit ].set_units(
                            *_internals[ 0 ].get_units( *sit ) );
                        for (size_t i = 0; i < _r.size(); ++i) {
                            scalar_u_t tempscalar(
                                _internals[ i ].get( *sit ),
                                _internals[ i ].get_units( *sit ) );
                            tempscalar.convert(
                                _scalar_properties[ *sit ].get_units() );
                            _scalar_properties[ *sit ][ i ] = tempscalar.get();
                        }
                    }
                    _build_scalar_splines();
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
                            .init( _r.size() )
                            .fill( _r, propit->second )
                            .ready();
                    }
                }

                std::unordered_map<std::string, grid_atmospheric_property_2d>
                    _properties;
                std::unordered_map<std::string, vector_u_t> _scalar_properties;
                std::unordered_map<
                    std::string,
                    NCPA::interpolation::Interpolator1D<double, double>>
                    _scalar_splines;
                vector_u_t _z, _r;
                std::vector<Atmosphere1D> _internals;
                NCPA::interpolation::interpolator_2d_type_t _interpolator_type
                    = NCPA_ATMOSPHERE_DEFAULT_2D_INTERPOLATOR;
                NCPA::interpolation::interpolator_1d_type_t
                    _1d_interpolator_type
                    = NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR;
        };

        class stratified_atmosphere_2d : public abstract_atmosphere_2d {
            public:
                stratified_atmosphere_2d() : abstract_atmosphere_2d() {
                    _dummy = vector_u_t(
                        1, scalar_u_t( 0.0, NCPA::units::KILOMETERS ) );
                }

                stratified_atmosphere_2d(
                    const stratified_atmosphere_2d& other ) :
                    stratified_atmosphere_2d() {
                    _1d                   = other._1d;
                    _dummy                = other._dummy;
                    _scalar_properties    = other._scalar_properties;
                    _1d_interpolator_type = other._1d_interpolator_type;
                    _build_scalar_splines();
                }

                DECLARE_BOILERPLATE_METHODS( stratified_atmosphere_2d,
                                             abstract_atmosphere_2d )

                virtual size_t size( size_t dim ) const override {
                    this->validate_axis( dim );
                    if (dim == 0) {
                        return 1;
                    } else {
                        return _1d.size();
                    }
                }

                virtual abstract_atmosphere_2d& set(
                    abstract_atmosphere_1d& atmos1d ) override {
                    _1d              = tuple_atmosphere_1d();
                    auto vector_keys = atmos1d.get_vector_keys();
                    vector_u_t z     = atmos1d.get_axis_vector();

                    for (auto vit = vector_keys.cbegin();
                         vit != vector_keys.cend(); ++vit) {
                        _1d.add_property( *vit,
                                          tuple_atmospheric_property_1d(
                                              atmos1d.get_axis_vector(),
                                              atmos1d.get_vector( *vit ) ) );
                    }
                    auto scalar_keys = atmos1d.get_scalar_keys();
                    for (auto sit = scalar_keys.cbegin();
                         sit != scalar_keys.cend(); ++sit) {
                        _1d.add_property(
                            *sit, scalar_u_t( atmos1d.get( *sit ),
                                              atmos1d.get_units( *sit ) ) );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& set(
                    const vector_u_t ranges,
                    std::vector<abstract_atmosphere_1d *> components )
                    override {
                    if (components.size() > 1) {
                        throw std::range_error(
                            "Cannot set a 2-D stratified atmosphere with "
                            "multiple 1-D profiles" );
                    }
                    return this->set( *components[ 0 ] );
                }

                virtual abstract_atmosphere_2d& append(
                    scalar_u_t range,
                    abstract_atmosphere_1d& atmos1d ) override {
                    throw std::logic_error(
                        "Cannot append to a stratified atmosphere!" );
                }

                virtual std::unique_ptr<abstract_atmosphere_1d> extract(
                    double range ) override {
                    return std::unique_ptr<abstract_atmosphere_1d>(
                        _1d.clone() );
                }

                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    throw NCPA::NotImplementedError(
                        "2-D interpolator not applicable to stratified "
                        "profiles" );
                }

                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _1d.set_interpolator( interp_type );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& set_axis(
                    size_t axis, vector_u_t vals ) override {
                    this->validate_axis( axis );
                    if (axis == 0) {
                        _dummy.convert_units( *vals.get_units() );
                    } else {
                        this->resample( 1, vals );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const AtmosphericProperty2D& property ) override {
                    AtmosphericProperty1D prop = property.extract( 0.0 );
                    _1d.add_property( key, *prop.internal() );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const vector2d_u_t& property ) override {
                    vector_u_t z = _1d.get_axis_vector();
                    if (z.size() != property.dim( 1 )) {
                        throw std::range_error(
                            "Size mismatch between internal Z vector and "
                            "second dimension of provided property" );
                    }
                    vector_u_t p( property[ 0 ], property.get_units() );
                    tuple_atmospheric_property_1d prop( z, p );
                    _1d.add_property( key, prop );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const vector_u_t& property ) override {
                    throw std::logic_error(
                        "stratified_atmosphere_2d.add_property(): Adding a "
                        "1-d vector property requires an index vector" );
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
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& add_property(
                    const std::string& key,
                    const scalar_u_t& property ) override {
                    _1d.add_property( key, property );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& remove_property(
                    const std::string& key ) override {
                    _1d.remove_property( key );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& copy_property(
                    const std::string& old_key,
                    const std::string& new_key ) override {
                    _1d.copy_property( old_key, new_key );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual vector_u_t get_axis_vector( size_t n ) override {
                    this->validate_axis( n );
                    return ( n == 0 ? _dummy : _1d.get_axis_vector() );
                }

                virtual vector_u_t get_axis_vector( size_t n ) const override {
                    return ( n == 0 ? _dummy : _1d.get_axis_vector() );
                }

                virtual double get( const std::string& key,
                                    double r ) override {
                    if (_scalar_splines.count( key ) != 0) {
                        return _scalar_splines[ key ].eval_f( r );
                    } else {
                        return _1d.get( key );
                    }
                }

                virtual double get( const std::string& key, double r,
                                    double z ) override {
                    return _1d.get( key, z );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double r, double z,
                                                     size_t wrt ) override {
                    if (_scalar_splines.count( key ) != 0) {
                        return ( wrt == 0 ? _scalar_splines[ key ].eval_df( r )
                                          : 0.0 );
                    } else {
                        return ( wrt == 0
                                     ? 0.0
                                     : _1d.get_first_derivative( key, z ) );
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
                        return ( wrt1 == 0 || wrt2 == 0
                                     ? 0.0
                                     : _1d.get_second_derivative( key, z ) );
                    }
                }

                virtual units_ptr_t get_axis_units( size_t n ) const override {
                    this->validate_axis( n );
                    return ( n == 0 ? _dummy.get_units()
                                    : _1d.get_axis_units() );
                }

                virtual units_ptr_t get_units(
                    const std::string& key ) const override {
                    if (_scalar_properties.count( key ) != 0) {
                        return _scalar_properties.at( key ).second.get_units();
                    } else {
                        return _1d.get_units( key );
                    }
                }

                virtual double get_minimum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if (n == 0) {
                        if (_scalar_properties.empty()) {
                            return 0.0;
                        } else {
                            double minax = _scalar_properties.cbegin()
                                               ->second.second.back();
                            for (auto it = _scalar_properties.cbegin();
                                 it != _scalar_properties.cend(); ++it) {
                                minax = std::min( minax,
                                                  it->second.second.front() );
                            }
                            return minax;
                        }
                    } else {
                        return _1d.get_minimum_axis();
                    }
                }

                virtual double get_maximum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if (n == 0) {
                        if (_scalar_properties.empty()) {
                            return 0.0;
                        } else {
                            double maxax = _scalar_properties.cbegin()
                                               ->second.second.front();
                            for (auto it = _scalar_properties.cbegin();
                                 it != _scalar_properties.cend(); ++it) {
                                maxax = std::min( maxax,
                                                  it->second.second.back() );
                            }
                            return maxax;
                        }
                    } else {
                        return _1d.get_maximum_axis();
                    }
                }

                virtual abstract_atmosphere_2d& convert_axis_units(
                    size_t n, units_ptr_t new_units ) override {
                    this->validate_axis( n );
                    if (n == 0) {
                        if (_scalar_properties.empty()) {
                            _dummy.convert_units( *new_units );
                        } else {
                            for (auto it = _scalar_properties.begin();
                                 it != _scalar_properties.end(); ++it) {
                                it->second.first.convert_units( *new_units );
                            }
                            _build_scalar_splines();
                        }
                    } else {
                        _1d.convert_axis_units( new_units );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& convert_units(
                    const std::string& key, units_ptr_t new_units ) override {
                    if (_scalar_properties.count( key ) != 0) {
                        _scalar_properties[ key ].second.convert_units(
                            *new_units );
                        _build_scalar_splines();
                    } else {
                        _1d.convert_units( key, new_units );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    size_t axis, double new_d ) override {
                    this->validate_axis( axis );
                    if (axis == 1) {
                        _1d.resample( new_d );
                    } else {
                        for (auto it = _scalar_properties.begin();
                             it != _scalar_properties.end(); ++it) {
                            NCPA::interpolation::Interpolator1D<double, double>
                                interp = NCPA::interpolation::
                                    InterpolatorFactory<double, double>::build(
                                        _1d_interpolator_type );
                            interp.init( it->second.first.size() )
                                .fill( it->second.first, it->second.second )
                                .ready();
                            vector_u_t newv = it->second.second;
                            vector_u_t newr = it->second.first;
                            newv.resize( 0 );
                            newr.resize( 0 );
                            for (double d = it->second.first.front();
                                 d <= it->second.first.back(); d += new_d) {
                                newr.push_back( d );
                                newv.push_back( interp.eval_f( d ) );
                            }
                            it->second.first  = newr;
                            it->second.second = newv;
                        }
                        _build_scalar_splines();
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    size_t axis, vector_u_t new_z ) override {
                    this->validate_axis( axis );
                    if (axis == 1) {
                        _1d.resample( new_z );
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
                    }
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
                    _build_scalar_splines();
                    _1d.resample( new_ax2 );
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
                    return _1d.get_vector_keys();
                }

                virtual std::vector<std::string> get_scalar_keys()
                    const override {
                    std::vector<std::string> skeys = _1d.get_scalar_keys(),
                                             spkeys, finalkeys;
                    for (auto it = _scalar_properties.cbegin();
                         it != _scalar_properties.cend(); ++it) {
                        spkeys.push_back( it->first );
                    }
                    std::set_intersection( spkeys.cbegin(), spkeys.cend(),
                                           skeys.cbegin(), skeys.cend(),
                                           std::back_inserter( finalkeys ) );
                    return finalkeys;
                }

                virtual bool contains_scalar(
                    const std::string& key ) const override {
                    return _1d.contains_scalar( key )
                        || ( _scalar_properties.count( key ) > 0 );
                }

                virtual bool contains_vector(
                    const std::string& key ) const override {
                    return _1d.contains_vector( key );
                }

                virtual bool contains_key(
                    const std::string& key ) const override {
                    return contains_vector( key ) || contains_scalar( key );
                }

                virtual void print( std::ostream& os ) override {
                    throw NCPA::NotImplementedError(
                        "Print function not yet implemented." );
                }

                virtual bool same( scalar_u_t r1,
                                   scalar_u_t r2 ) const override {
                    return true;
                }

                virtual bool same( double r1, double r2 ) const override {
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

                tuple_atmosphere_1d _1d;
                vector_u_t _dummy;
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
        };

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

                DECLARE_BOILERPLATE_METHODS(
                    piecewise_stratified_atmosphere_2d,
                    abstract_atmosphere_2d )

                virtual size_t size( size_t dim ) const override {
                    this->validate_axis( dim );
                    return _axes[ dim ].size();
                }

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
                            return i-1;
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

        class Atmosphere2D : public AtmosphericModel {
            public:
                Atmosphere2D() : AtmosphericModel() {}

                Atmosphere2D( std::unique_ptr<abstract_atmosphere_2d> ptr ) :
                    Atmosphere2D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                Atmosphere2D( const Atmosphere2D& other ) : Atmosphere2D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                Atmosphere2D( Atmosphere2D&& source ) noexcept :
                    Atmosphere2D() {
                    ::swap( *this, source );
                }

                virtual ~Atmosphere2D() {}

                friend void ::swap( Atmosphere2D& a,
                                    Atmosphere2D& b ) noexcept;

                Atmosphere2D& operator=( Atmosphere2D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                Atmosphere2D& set( Atmosphere1D& atmos1d ) {
                    check_pointer();
                    _ptr->set( *atmos1d.internal() );
                    return *this;
                }

                Atmosphere2D& set( const vector_u_t ranges,
                                   std::vector<Atmosphere1D>& atms ) {
                    check_pointer();
                    std::vector<abstract_atmosphere_1d *> ptrs;
                    for (auto it = atms.begin(); it != atms.end(); ++it) {
                        ptrs.push_back( it->internal() );
                    }
                    _ptr->set( ranges, ptrs );
                    return *this;
                }

                Atmosphere2D& append( const scalar_u_t& range,
                                      const Atmosphere1D& atmos1d ) {
                    check_pointer();
                    _ptr->append( range, *atmos1d.internal()->clone().get() );
                    return *this;
                }

                virtual Atmosphere2D set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type ) {
                    check_pointer();
                    if (NCPA::interpolation::InterpolatorFactory<
                            double, double>::can_build( interp_type )) {
                        _ptr->set_interpolator( interp_type );
                    } else {
                        throw std::logic_error(
                            "Selected interpolator type not available" );
                    }
                    return *this;
                }

                virtual Atmosphere2D set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type ) {
                    check_pointer();
                    if (NCPA::interpolation::InterpolatorFactory<
                            double, double>::can_build( interp_type )) {
                        _ptr->set_interpolator( interp_type );
                    } else {
                        throw std::logic_error(
                            "Selected interpolator type not available" );
                    }
                    return *this;
                }

                virtual Atmosphere2D& add_property(
                    const std::string& key,
                    const AtmosphericProperty2D& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere2D& add_property(
                    const std::string& key, const vector2d_u_t& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere2D& add_property(
                    const std::string& key, const scalar_u_t& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere2D remove_property(
                    const std::string& key ) {
                    check_pointer();
                    _ptr->remove_property( key );
                    return *this;
                }

                virtual Atmosphere2D& copy_property(
                    const std::string& old_key, const std::string& new_key ) {
                    check_pointer();
                    _ptr->copy_property( old_key, new_key );
                    return *this;
                }

                // virtual AtmosphericProperty2D& get_property(
                //     const std::string& key ) const {
                //     check_pointer();
                //     return _ptr->get_property( key );
                // }

                virtual vector_u_t get_axis_vector( size_t dim ) {
                    check_pointer();
                    return _ptr->get_axis_vector( dim );
                }

                virtual double get( const std::string& key ) {
                    check_pointer();
                    return _ptr->get( key, 0.0 );
                }

                virtual double get( const std::string& key, double range,
                                    double altitude ) {
                    check_pointer();
                    return _ptr->get( key, range, altitude );
                }

                virtual double get(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& range,
                    const NCPA::units::ScalarWithUnits<double>& altitude ) {
                    check_pointer();
                    return _ptr->get(
                        key, range.get_as( this->get_axis_units( 0 ) ),
                        altitude.get_as( this->get_axis_units( 1 ) ) );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double range,
                                                     double altitude,
                                                     size_t wrt1 ) {
                    check_pointer();
                    return _ptr->get_first_derivative( key, range, altitude,
                                                       wrt1 );
                }

                virtual double get_first_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& range,
                    const NCPA::units::ScalarWithUnits<double>& altitude,
                    size_t wrt1 ) {
                    check_pointer();
                    return _ptr->get_first_derivative(
                        key, range.get_as( this->get_axis_units( 0 ) ),
                        altitude.get_as( this->get_axis_units( 1 ) ), wrt1 );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double range,
                                                      double altitude,
                                                      size_t wrt1,
                                                      size_t wrt2 ) {
                    check_pointer();
                    return _ptr->get_second_derivative( key, range, altitude,
                                                        wrt1, wrt2 );
                }

                virtual double get_second_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& range,
                    const NCPA::units::ScalarWithUnits<double>& altitude,
                    size_t wrt1, size_t wrt2 ) {
                    check_pointer();
                    return _ptr->get_second_derivative(
                        key, range.get_as( this->get_axis_units( 0 ) ),
                        altitude.get_as( this->get_axis_units( 1 ) ), wrt1,
                        wrt2 );
                }

                virtual units_ptr_t get_axis_units( size_t dim ) {
                    check_pointer();
                    return _ptr->get_axis_units( dim );
                }

                virtual units_ptr_t get_units( const std::string& key ) {
                    check_pointer();
                    return _ptr->get_units( key );
                }

                virtual double get_minimum_axis( size_t dim ) const {
                    check_pointer();
                    return _ptr->get_minimum_axis( dim );
                }

                virtual double get_maximum_axis( size_t dim ) const {
                    check_pointer();
                    return _ptr->get_maximum_axis( dim );
                }

                virtual Atmosphere2D& convert_axis_units(
                    size_t dim, units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_axis_units( dim, new_units );
                    return *this;
                }

                virtual Atmosphere2D& convert_units( const std::string& key,
                                                     units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_units( key, new_units );
                    return *this;
                }

                virtual Atmosphere2D& resample( size_t dim, double new_dz ) {
                    check_pointer();
                    _ptr->resample( dim, new_dz );
                    return *this;
                }

                virtual Atmosphere2D& resample( size_t dim,
                                                vector_u_t new_z ) {
                    check_pointer();
                    _ptr->resample( dim, new_z );
                    return *this;
                }

                virtual Atmosphere2D& resample( vector_u_t new_ax1,
                                                vector_u_t new_ax2 ) {
                    check_pointer();
                    _ptr->resample( new_ax1, new_ax2 );
                    return *this;
                }

                virtual std::vector<std::string> get_keys() const {
                    check_pointer();
                    return _ptr->get_keys();
                }

                virtual std::vector<std::string> get_vector_keys() const {
                    check_pointer();
                    return _ptr->get_vector_keys();
                }

                virtual std::vector<std::string> get_scalar_keys() const {
                    check_pointer();
                    return _ptr->get_scalar_keys();
                }

                virtual bool contains_scalar( const std::string& key ) const {
                    check_pointer();
                    return _ptr->contains_scalar( key );
                }

                virtual bool contains_vector( const std::string& key ) const {
                    check_pointer();
                    return _ptr->contains_vector( key );
                }

                virtual bool contains_key( const std::string& key ) const {
                    check_pointer();
                    return _ptr->contains_key( key );
                }

                virtual bool same( double r1, double r2 ) const {
                    check_pointer();
                    return _ptr->same( r1, r2 );
                }

                virtual bool same( scalar_u_t r1, scalar_u_t r2 ) const {
                    check_pointer();
                    return _ptr->same( r1, r2 );
                }

                virtual void check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "Atmosphere2D: Internal pointer not set!" );
                    }
                }

                // friend binary operators
                friend std::ostream& operator<<( std::ostream& os,
                                                 const Atmosphere2D& atm ) {
                    if (atm) {
                        atm._ptr->print( os );
                    }
                    return os;
                }

                explicit operator bool() const {
                    return ( _ptr ? true : false );
                }

            private:
                std::unique_ptr<abstract_atmosphere_2d> _ptr;
        };


    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::stratified_atmosphere_2d& a,
                  NCPA::atmos::stratified_atmosphere_2d& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::atmos::abstract_atmosphere_2d&>( a ),
            dynamic_cast<NCPA::atmos::abstract_atmosphere_2d&>( b ) );
    swap( a._1d, b._1d );
    swap( a._dummy, b._dummy );
    swap( a._scalar_properties, b._scalar_properties );
    swap( a._scalar_splines, b._scalar_splines );
    swap( a._1d_interpolator_type, b._1d_interpolator_type );
}

static void swap( NCPA::atmos::grid_atmosphere_2d& a,
                  NCPA::atmos::grid_atmosphere_2d& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::atmos::abstract_atmosphere_2d&>( a ),
            dynamic_cast<NCPA::atmos::abstract_atmosphere_2d&>( b ) );
    swap( a._properties, b._properties );
    swap( a._scalar_properties, b._scalar_properties );
    swap( a._scalar_splines, b._scalar_splines );
    swap( a._z, b._z );
    swap( a._r, b._r );
    swap( a._internals, b._internals );
    swap( a._interpolator_type, b._interpolator_type );
    swap( a._1d_interpolator_type, b._1d_interpolator_type );
}

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

static void swap( NCPA::atmos::Atmosphere2D& a,
                  NCPA::atmos::Atmosphere2D& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::atmos::AtmosphericModel&>( a ),
          static_cast<NCPA::atmos::AtmosphericModel&>( b ) );
    a._ptr.swap( b._ptr );
}
