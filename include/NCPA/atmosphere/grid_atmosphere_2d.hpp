#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_1d.hpp"
#include "NCPA/atmosphere/abstract_atmosphere_2d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
// #include "NCPA/atmosphere/builders.hpp"
#include "NCPA/atmosphere/calculations.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/tuple_atmosphere_1d.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <string>
#include <vector>


static void swap( NCPA::atmos::grid_atmosphere_2d&,
                  NCPA::atmos::grid_atmosphere_2d& ) noexcept;

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
                    for ( auto vit = _properties.begin();
                          vit != _properties.end(); ++vit ) {
                        atm->add_property( vit->first,
                                           vit->second.extract( range ) );
                    }
                    for ( auto sit = _scalar_splines.begin();
                          sit != _scalar_splines.end(); ++sit ) {
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
                    if ( dim == 0 ) {
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
                }

                virtual abstract_atmosphere_2d& append(
                    scalar_u_t range,
                    abstract_atmosphere_1d& atmos1d ) override {
                    if ( range.get_as( _r.get_units() ) <= _r.back() ) {
                        throw std::range_error(
                            "Can't append inside existing range vector!" );
                    }
                    _properties.clear();
                    _scalar_properties.clear();
                    _r.push_back( range.get_as( _r.get_units() ) );
                    _internals.push_back( Atmosphere1D( atmos1d.clone() ) );
                    _internals.back().resample( _z );
                }

                virtual abstract_atmosphere_2d& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type )
                    override {
                    _interpolator_type = interp_type;
                    for ( auto it = _properties.begin();
                          it != _properties.end(); ++it ) {
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
                }

                virtual abstract_atmosphere_2d& set_axis(
                    size_t axis, vector_u_t vals ) override {
                    this->validate_axis( axis );
                    for ( auto it = _properties.begin();
                          it != _properties.end(); ++it ) {
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
                    if ( !_r.empty() && !_z.empty() ) {
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
                    if ( _z.size() != property.dim( 1 )
                         || _r.size() != property.dim( 0 ) ) {
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
                    if ( contains_vector( old_key ) ) {
                        _properties[ new_key ] = _properties[ old_key ];
                    } else if ( contains_scalar( old_key ) ) {
                        _scalar_properties[ new_key ]
                            = _scalar_properties[ old_key ];
                        _scalar_splines[ new_key ]
                            = _scalar_splines[ old_key ];
                    } else {
                        throw std::invalid_argument(
                            "Key " + old_key
                            + " does not exist in atmosphere!" );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                // virtual AtmosphericProperty2D& get_property(
                //     const std::string& key ) override {
                //     _assert_contains_vector( key );
                //     return AtmosphericProperty2D( _properties.at( key
                //     ).clone2d() );
                // }

                // virtual const AtmosphericProperty2D& get_property(
                //     const std::string& key ) const override {
                //     _assert_contains_vector( key );
                //     return _properties.at( key );
                // }

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
                    if ( this->contains_vector( key ) ) {
                        return _properties.at( key ).get_units();
                    } else {
                        _assert_contains_scalar( key );
                        return _scalar_properties.at( key ).get_units();
                    }
                }

                virtual double get_minimum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if ( n == 0 ) {
                        return 0.0;
                    }
                    double minax = _z.back();
                    for ( auto it = _properties.cbegin();
                          it != _properties.cend(); ++it ) {
                        minax = std::min( minax,
                                          it->second.get_limits( 1 ).first );
                    }
                    return minax;
                }

                virtual double get_maximum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if ( n == 0 ) {
                        return 0.0;
                    }
                    double maxax = _z.front();
                    for ( auto it = _properties.cbegin();
                          it != _properties.cend(); ++it ) {
                        maxax = std::max( maxax,
                                          it->second.get_limits( 1 ).second );
                    }
                    return maxax;
                }

                virtual abstract_atmosphere_2d& convert_axis_units(
                    size_t n, units_ptr_t new_units ) override {
                    this->validate_axis( n );
                    if ( n == 0 ) {
                        _r.convert_units( *new_units );
                    } else {
                        for ( auto it = _properties.begin();
                              it != _properties.end(); ++it ) {
                            it->second.convert_axis_units( n, *new_units );
                        }
                        _z.convert_units( *new_units );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& convert_units(
                    const std::string& key, units_ptr_t new_units ) override {
                    if ( this->contains_vector( key ) ) {
                        _properties[ key ].convert_units( *new_units );
                    } else if ( this->contains_scalar( key ) ) {
                        _scalar_properties[ key ].convert_units( *new_units );
                    } else {
                        throw std::range_error( "Unknown key: " + key );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    size_t axis, double new_d ) override {
                    this->validate_axis( axis );
                    if ( axis == 1 ) {
                        vector_u_t new_z;
                        new_z.set_units( *_z.get_units() );
                        double z = _z.front();
                        while ( z <= _z.back() ) {
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
                    if ( axis == 1 ) {
                        for ( auto it = _properties.begin();
                              it != _properties.end(); ++it ) {
                            it->second.resample( 1, new_z );
                        }
                        _z = new_z;
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    vector_u_t new_ax1, vector_u_t new_ax2 ) override {
                    for ( auto it = _properties.begin();
                          it != _properties.end(); ++it ) {
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
                    for ( auto it = _properties.cbegin();
                          it != _properties.cend(); ++it ) {
                        vkeys.push_back( it->first );
                    }
                    return vkeys;
                }

                virtual std::vector<std::string> get_scalar_keys()
                    const override {
                    std::vector<std::string> skeys;
                    for ( auto it = _scalar_properties.cbegin();
                          it != _scalar_properties.cend(); ++it ) {
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

            protected:
                void _assert_contains_vector( const std::string& key ) const {
                    if ( !contains_vector( key ) ) {
                        throw std::range_error( "Key " + key
                                                + " does not exist!" );
                    }
                }

                void _assert_contains_scalar( const std::string& key ) const {
                    if ( !contains_scalar( key ) ) {
                        throw std::range_error( "Key " + key
                                                + " does not exist!" );
                    }
                }

                void _assert_contains( const std::string& key ) const {
                    _assert_contains_vector( key );
                    _assert_contains_scalar( key );
                }

                void _assert_does_not_contain( const std::string& key ) const {
                    if ( contains_key( key ) ) {
                        throw std::range_error( "Key " + key
                                                + " already exists!" );
                    }
                }

                void _assert_z_is_set() const {
                    if ( !_z ) {
                        throw std::logic_error( "Altitude vector not set!" );
                    }
                }

                void _translate_to_grid() {
                    if ( !( _properties.empty()
                            || _properties.begin()->second.size( 0 )
                                   != _r.size() ) ) {
                        return;
                    }
                    _properties.clear();
                    _scalar_properties.clear();

                    auto vector_keys = _internals[ 0 ].get_vector_keys();
                    for ( auto vit = vector_keys.cbegin();
                          vit != vector_keys.cend(); ++vit ) {
                        grid_atmospheric_property_2d temp;
                        temp.set(
                            _r.get_scalar( 0 ),
                            *_internals[ 0 ].get_property( *vit ).internal() );
                        for ( size_t i = 1; i < _r.size(); ++i ) {
                            temp.append( _r.get_scalar( i ),
                                         *_internals[ i ]
                                              .get_property( *vit )
                                              .internal() );
                        }
                        this->add_property( *vit, temp );
                    }

                    _scalar_properties.clear();
                    auto scalar_keys = _internals[ 0 ].get_scalar_keys();
                    for ( auto sit = scalar_keys.cbegin();
                          sit != scalar_keys.cend(); ++sit ) {
                        _scalar_properties[ *sit ].resize( _r.size() );
                        _scalar_properties[ *sit ].set_units(
                            *_internals[ 0 ].get_units( *sit ) );
                        for ( size_t i = 0; i < _r.size(); ++i ) {
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
                    for ( auto propit = _scalar_properties.cbegin();
                          propit != _scalar_properties.cend(); ++propit ) {
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
                    = NCPA::interpolation::interpolator_2d_type_t::
                        LANL_BICUBIC;
                NCPA::interpolation::interpolator_1d_type_t
                    _1d_interpolator_type
                    = NCPA::interpolation::interpolator_1d_type_t::LANL_CUBIC;
        };


    }  // namespace atmos
}  // namespace NCPA

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
