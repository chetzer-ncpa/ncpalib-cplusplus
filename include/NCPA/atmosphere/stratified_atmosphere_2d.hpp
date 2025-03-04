#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_1d.hpp"
#include "NCPA/atmosphere/abstract_atmosphere_2d.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
// #include "NCPA/atmosphere/builders.hpp"
#include "NCPA/atmosphere/calculations.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/tuple_atmosphere_1d.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/exceptions.hpp"
#include "NCPA/interpolation.hpp"

#include <string>
#include <vector>


static void swap( NCPA::atmos::stratified_atmosphere_2d&,
                  NCPA::atmos::stratified_atmosphere_2d& ) noexcept;

namespace NCPA {
    namespace atmos {
        class stratified_atmosphere_2d : public abstract_atmosphere_2d {
            public:
                stratified_atmosphere_2d() : abstract_atmosphere_2d() {
                    _dummy = vector_u_t(
                        1, scalar_u_t( 0.0, NCPA::units::KILOMETERS ) );
                    // _z.set_units( NCPA::units::KILOMETERS );
                }

                stratified_atmosphere_2d(
                    const stratified_atmosphere_2d& other ) :
                    stratified_atmosphere_2d() {
                    _1d = other._1d;
                    // _properties        = other._properties;
                    // _scalar_properties = other._scalar_properties;
                    // _z                 = other._z;
                    // _interpolator_type = other._interpolator_type;
                    // _dummy             = other._dummy;
                }

                DECLARE_BOILERPLATE_METHODS( stratified_atmosphere_2d,
                                             abstract_atmosphere_2d )

                virtual size_t size( size_t dim ) const override {
                    this->validate_axis( dim );
                    if ( dim == 0 ) {
                        return 1;
                    } else {
                        return _1d.size();
                    }
                }

                virtual abstract_atmosphere_2d& set(
                    abstract_atmosphere_1d& atmos1d ) override {
                    _1d = tuple_atmosphere_1d();
                    auto vector_keys = atmos1d.get_vector_keys();
                    vector_u_t z = atmos1d.get_axis_vector();
                    
                    for ( auto vit = vector_keys.cbegin();
                          vit != vector_keys.cend(); ++vit ) {
                        _1d.add_property( *vit,
                                          tuple_atmospheric_property_1d(
                                              atmos1d.get_axis_vector(),
                                              atmos1d.get_vector( *vit ) ) );
                    }
                    auto scalar_keys = atmos1d.get_scalar_keys();
                    for ( auto sit = scalar_keys.cbegin();
                          sit != scalar_keys.cend(); ++sit ) {
                        _1d.add_property( *sit, scalar_u_t( atmos1d.get( *sit ), atmos1d.get_units( *sit ) ));
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
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
                    if ( axis == 0 ) {
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
                    if ( z.size() != property.dim( 1 ) ) {
                        throw std::range_error(
                            "Size mismatch between internal Z vector and "
                            "second dimension of provided property" );
                    }
                    // _properties[ key ] = AtmosphereFactory::build(
                    //     atmospheric_property_2d_t::STRATIFIED );
                    vector_u_t p( property[ 0 ], property.get_units() );
                    tuple_atmospheric_property_1d prop( z, p );
                    _1d.add_property( key, prop );
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

                // virtual AtmosphericProperty2D& get_property(
                //     const std::string& key ) override {
                //     _assert_contains_vector( key );
                //     return _properties.at( key );
                // }

                // virtual const AtmosphericProperty2D& get_property(
                //     const std::string& key ) const override {
                //     return _1d.get( key );
                // }

                virtual vector_u_t get_axis_vector( size_t n ) override {
                    this->validate_axis( n );
                    return ( n == 0 ? _dummy : _1d.get_axis_vector() );
                }

                virtual vector_u_t get_axis_vector( size_t n ) const override {
                    return ( n == 0 ? _dummy : _1d.get_axis_vector() );
                }

                virtual double get( const std::string& key,
                                    double r ) override {
                    return _1d.get( key );
                }

                virtual double get( const std::string& key, double r,
                                    double z ) override {
                    return _1d.get( key, z );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double r, double z,
                                                     size_t wrt ) override {
                    return ( wrt == 0 ? 0.0
                                      : _1d.get_first_derivative( key, z ) );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double r, double z,
                                                      size_t wrt1,
                                                      size_t wrt2 ) override {
                    return ( wrt1 == 0 || wrt2 == 0
                                 ? 0.0
                                 : _1d.get_second_derivative( key, z ) );
                }

                virtual units_ptr_t get_axis_units( size_t n ) const override {
                    this->validate_axis( n );
                    return ( n == 0 ? _dummy.get_units()
                                    : _1d.get_axis_units() );
                }

                virtual units_ptr_t get_units(
                    const std::string& key ) const override {
                    return _1d.get_units( key );
                }

                virtual double get_minimum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if ( n == 0 ) {
                        return 0.0;
                    }
                    return _1d.get_minimum_axis();
                }

                virtual double get_maximum_axis( size_t n ) const override {
                    this->validate_axis( n );
                    if ( n == 0 ) {
                        return 0.0;
                    }
                    return _1d.get_maximum_axis();
                }

                virtual abstract_atmosphere_2d& convert_axis_units(
                    size_t n, units_ptr_t new_units ) override {
                    this->validate_axis( n );
                    if ( n == 0 ) {
                        _dummy.convert_units( *new_units );

                    } else {
                        _1d.convert_axis_units( new_units );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& convert_units(
                    const std::string& key, units_ptr_t new_units ) override {
                    _1d.convert_units( key, new_units );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    size_t axis, double new_d ) override {
                    this->validate_axis( axis );
                    if ( axis == 1 ) {
                        _1d.resample( new_d );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    size_t axis, vector_u_t new_z ) override {
                    this->validate_axis( axis );
                    if ( axis == 1 ) {
                        _1d.resample( new_z );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual abstract_atmosphere_2d& resample(
                    vector_u_t new_ax1, vector_u_t new_ax2 ) override {
                    _1d.resample( new_ax2 );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_2D;
                }

                virtual std::vector<std::string> get_keys() const override {
                    std::vector<std::string> vkeys = this->get_vector_keys(),
                                             skeys = this->get_scalar_keys(),
                                             jkeys;
                    std::set_intersection( vkeys.cbegin(), vkeys.cend(),
                                           skeys.cbegin(), skeys.cend(),
                                           std::back_inserter( jkeys ) );
                    return jkeys;
                }

                virtual std::vector<std::string> get_vector_keys()
                    const override {
                    return _1d.get_vector_keys();
                }

                virtual std::vector<std::string> get_scalar_keys()
                    const override {
                    return _1d.get_scalar_keys();
                }

                virtual bool contains_scalar(
                    const std::string& key ) const override {
                    return _1d.contains_scalar( key );
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

                // void _assert_z_is_set() const {
                //     if ( !_z ) {
                //         throw std::logic_error( "Altitude vector not set!"
                //         );
                //     }
                // }

                // std::unordered_map<std::string, AtmosphericProperty2D>
                //     _properties;
                // std::unordered_map<std::string, scalar_u_t>
                // _scalar_properties; vector_u_t _z, _dummy;
                tuple_atmosphere_1d _1d;
                vector_u_t _dummy;

                // NCPA::interpolation::interpolator_1d_type_t
                // _interpolator_type
                //     = ( NCPA_INTERPOLATION_GSL_STEFFEN_SPLINE_AVAILABLE
                //             ? NCPA::interpolation::interpolator_1d_type_t::
                //                   GSL_STEFFEN
                //             : NCPA::interpolation::interpolator_1d_type_t::
                //                   LANL_CUBIC );
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
}
