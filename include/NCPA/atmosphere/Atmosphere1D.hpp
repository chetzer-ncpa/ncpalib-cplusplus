#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_1d.hpp"
#include "NCPA/atmosphere/AtmosphericModel.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/units.hpp"

#include <memory>
#include <string>

static void swap( NCPA::atmos::tuple_atmosphere_1d&,
                  NCPA::atmos::tuple_atmosphere_1d& ) noexcept;

static void swap( NCPA::atmos::Atmosphere1D&,
                  NCPA::atmos::Atmosphere1D& ) noexcept;

namespace NCPA {
    namespace atmos {

        class tuple_atmosphere_1d : public abstract_atmosphere_1d {
            public:
                tuple_atmosphere_1d() : abstract_atmosphere_1d() {}

                tuple_atmosphere_1d( const tuple_atmosphere_1d& other ) :
                    tuple_atmosphere_1d() {
                    _properties        = other._properties;
                    _scalar_properties = other._scalar_properties;
                    _z                 = other._z;
                    _interpolator_type = other._interpolator_type;
                }

                tuple_atmosphere_1d( tuple_atmosphere_1d&& source ) noexcept :
                    tuple_atmosphere_1d() {
                    ::swap( *this, source );
                }

                virtual ~tuple_atmosphere_1d() {}

                tuple_atmosphere_1d& operator=( tuple_atmosphere_1d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( tuple_atmosphere_1d& a,
                                    tuple_atmosphere_1d& b ) noexcept;

                virtual std::unique_ptr<abstract_atmosphere_1d> clone()
                    const override {
                    return std::unique_ptr<abstract_atmosphere_1d>(
                        new tuple_atmosphere_1d( *this ) );
                }

                virtual size_t size() const override { return _z.size(); }

                virtual abstract_atmosphere_1d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _interpolator_type = interp_type;
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.set_interpolator( _interpolator_type );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D;
                }

                virtual abstract_atmosphere_1d& set_axis(
                    vector_u_t z ) override {
                    _z = z;
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D;
                }

                virtual abstract_atmosphere_1d& add_property(
                    const std::string& key,
                    const AtmosphericProperty1D& property ) override {
                    _assert_does_not_contain( key );
                    _properties[ key ] = property;
                    _properties[ key ].set_interpolator( _interpolator_type );
                    if (_z.empty()) {
                        _z = property.axis();
                    } else {
                        _properties[ key ].resample( _z );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D
                }

                virtual abstract_atmosphere_1d& add_property(
                    const std::string& key,
                    const scalar_u_t& property ) override {
                    _assert_does_not_contain( key );
                    _scalar_properties[ key ] = property;
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D;
                }

                virtual abstract_atmosphere_1d& remove_property(
                    const std::string& key ) override {
                    _properties.erase( key );
                    _scalar_properties.erase( key );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D;
                }

                virtual abstract_atmosphere_1d& copy_property(
                    const std::string& old_key,
                    const std::string& new_key ) override {
                    _assert_does_not_contain( new_key );
                    if (contains_vector( old_key )) {
                        _properties[ new_key ] = _properties[ old_key ];
                    } else if (contains_scalar( old_key )) {
                        _scalar_properties[ new_key ]
                            = _scalar_properties[ old_key ];
                    } else {
                        throw std::invalid_argument(
                            "Key " + old_key
                            + " does not exist in atmosphere!" );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D;
                }

                virtual AtmosphericProperty1D& get_property(
                    const std::string& key ) override {
                    _assert_contains_vector( key );
                    return _properties.at( key );
                }

                virtual const AtmosphericProperty1D& get_property(
                    const std::string& key ) const override {
                    _assert_contains_vector( key );
                    return _properties.at( key );
                }

                virtual vector_u_t get_axis_vector() const override {
                    return _z;
                }

                virtual vector_u_t get_vector(
                    const std::string& key ) override {
                    return _properties.at( key ).values();
                }

                virtual double get( const std::string& key ) const override {
                    return _scalar_properties.at( key ).get();
                }

                virtual double get( const std::string& key,
                                    double altitude ) override {
                    // return _properties.at( key ).get( altitude );
                    return get_property( key ).get( altitude );
                }

                virtual double get_first_derivative(
                    const std::string& key, double altitude ) override {
                    return get_property( key ).get_first_derivative(
                        altitude );
                }

                virtual double get_second_derivative(
                    const std::string& key, double altitude ) override {
                    return get_property( key ).get_second_derivative(
                        altitude );
                }

                virtual units_ptr_t get_axis_units() const override {
                    return _z.get_units();
                }

                virtual units_ptr_t get_units(
                    const std::string& key ) const override {
                    if (contains_vector( key )) {
                        return get_property( key ).get_units();
                    } else if (contains_scalar( key )) {
                        scalar_u_t s = _scalar_properties.at( key );
                        return _scalar_properties.at( key ).get_units();
                    } else {
                        throw std::out_of_range( "Key " + key
                                                 + " not found!" );
                    }
                }

                virtual double get_minimum_axis() const override {
                    return _z.front();
                }

                virtual double get_maximum_axis() const override {
                    return _z.back();
                }

                virtual abstract_atmosphere_1d& convert_axis_units(
                    units_ptr_t new_units ) override {
                    _z.convert_units( *new_units );
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.convert_axis_units( *new_units );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D
                }

                virtual abstract_atmosphere_1d& convert_units(
                    const std::string& key, units_ptr_t new_units ) override {
                    if (contains_vector( key )) {
                        get_property( key ).convert_units( *new_units );
                    } else if (contains_scalar( key )) {
                        _scalar_properties[ key ].convert_units( *new_units );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D
                }

                virtual abstract_atmosphere_1d& resample(
                    vector_u_t new_z ) override {
                    for (auto it = _properties.begin();
                         it != _properties.end(); ++it) {
                        it->second.resample( new_z );
                    }
                    _z = new_z;
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERE_1D;
                }

                virtual abstract_atmosphere_1d& resample(
                    double new_dz ) override {
                    vector_u_t new_z;
                    new_z.set_units( *_z.get_units() );
                    double z = _z.front();
                    while (z <= _z.back()) {
                        new_z.push_back( z );
                        z += new_dz;
                    }
                    return resample( new_z );
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
                    this->print( os, this->get_vector_keys(), "Z" );
                }

                virtual void print(
                    std::ostream& os,
                    const std::vector<std::string>& columnorder,
                    const std::string& altitude_key ) const {
                    // check columnorder variable for key validity
                    std::vector<std::string>::const_iterator vit;
                    for (vit = columnorder.cbegin(); vit != columnorder.cend();
                         ++vit) {
                        if (!contains_vector( *vit )) {
                            throw std::invalid_argument(
                                "No vector quantity exists with key " + *vit );
                        }
                    }

                    // first we do the header.  That contains column
                    // descriptions as well as scalar values scalars first
                    for (auto mit = _scalar_properties.cbegin();
                         mit != _scalar_properties.cend(); ++mit) {
                        os << "#% 0, "
                           << NCPA::strings::deblank( mit->first, "_" ) << ", "
                           << *( mit->second.get_units() ) << ", "
                           << mit->second.get() << std::endl;
                    }

                    // Now column descriptors.  Altitude first
                    os << "#% 1, " << altitude_key << ", "
                       << *( _z.get_units() ) << std::endl;
                    unsigned int column = 2;
                    for (vit = columnorder.cbegin(); vit != columnorder.cend();
                         ++vit) {
                        os << "#% " << column << ", "
                           << NCPA::strings::deblank( *vit, "_" ) << ", "
                           << *( this->get_units( *vit ) ) << std::endl;
                        column++;
                    }

                    // Now columns
                    size_t nz_ = _z.size();
                    os.setf( std::ios::scientific, std::ios::floatfield );
                    os.setf( std::ios::right, std::ios::adjustfield );
                    os.precision( 6 );
                    os.width( 9 );
                    os.fill( ' ' );
                    for (size_t i = 0; i < nz_; i++) {
                        os << _z[ i ];
                        for (vit = columnorder.cbegin();
                             vit != columnorder.cend(); ++vit) {
                            os << " " << get_property( *vit ).values()[ i ];
                        }
                        os << std::endl;
                    }
                    os.flush();
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

                std::unordered_map<std::string, AtmosphericProperty1D>
                    _properties;
                std::unordered_map<std::string, scalar_u_t> _scalar_properties;
                vector_u_t _z;
                NCPA::interpolation::interpolator_1d_type_t _interpolator_type
                    = ( NCPA_INTERPOLATION_GSL_STEFFEN_SPLINE_AVAILABLE
                            ? NCPA::interpolation::interpolator_1d_type_t::
                                  GSL_STEFFEN
                            : NCPA::interpolation::interpolator_1d_type_t::
                                  LANL_CUBIC );
        };

        class Atmosphere1D : public AtmosphericModel {
            public:
                Atmosphere1D() : AtmosphericModel() {}

                Atmosphere1D( std::unique_ptr<abstract_atmosphere_1d> ptr ) :
                    Atmosphere1D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                Atmosphere1D( const Atmosphere1D& other ) : Atmosphere1D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                Atmosphere1D( Atmosphere1D&& source ) noexcept :
                    Atmosphere1D() {
                    ::swap( *this, source );
                }

                virtual ~Atmosphere1D() {}

                friend void ::swap( Atmosphere1D& a,
                                    Atmosphere1D& b ) noexcept;

                Atmosphere1D& operator=( Atmosphere1D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual Atmosphere1D set_interpolator(
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

                virtual Atmosphere1D& add_property(
                    const std::string& key,
                    const AtmosphericProperty1D& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                // virtual Atmosphere1D& add_property(
                //     const std::string& key, const vector_u_t& property ) {
                //     check_pointer();
                //     _ptr->add_property( key, property );
                //     return *this;
                // }

                virtual Atmosphere1D& add_property(
                    const std::string& key, const scalar_u_t& property ) {
                    check_pointer();
                    _ptr->add_property( key, property );
                    return *this;
                }

                virtual Atmosphere1D remove_property(
                    const std::string& key ) {
                    check_pointer();
                    _ptr->remove_property( key );
                    return *this;
                }

                virtual Atmosphere1D& copy_property(
                    const std::string& old_key, const std::string& new_key ) {
                    check_pointer();
                    _ptr->copy_property( old_key, new_key );
                    return *this;
                }

                virtual AtmosphericProperty1D& get_property(
                    const std::string& key ) const {
                    check_pointer();
                    return _ptr->get_property( key );
                }

                virtual size_t size() const {
                    check_pointer();
                    return _ptr->get_axis_vector().size(); 
                }

                virtual vector_u_t get_axis_vector() {
                    check_pointer();
                    return _ptr->get_axis_vector();
                }

                virtual double get( const std::string& key ) {
                    check_pointer();
                    return _ptr->get( key );
                }

                virtual double get( const std::string& key, double altitude ) {
                    check_pointer();
                    return _ptr->get( key, altitude );
                }

                virtual double get(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& altitude ) {
                    check_pointer();
                    return _ptr->get(
                        key, altitude.get_as( this->get_axis_units() ) );
                }

                virtual double get_first_derivative( const std::string& key,
                                                     double altitude ) {
                    check_pointer();
                    return _ptr->get_first_derivative( key, altitude );
                }

                virtual double get_first_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& altitude ) {
                    check_pointer();
                    return _ptr->get_first_derivative(
                        key, altitude.get_as( this->get_axis_units() ) );
                }

                virtual double get_second_derivative( const std::string& key,
                                                      double altitude ) {
                    check_pointer();
                    return _ptr->get_second_derivative( key, altitude );
                }

                virtual double get_second_derivative(
                    const std::string& key,
                    const NCPA::units::ScalarWithUnits<double>& altitude ) {
                    check_pointer();
                    return _ptr->get_second_derivative(
                        key, altitude.get_as( this->get_axis_units() ) );
                }

                virtual units_ptr_t get_axis_units() {
                    check_pointer();
                    return _ptr->get_axis_units();
                }

                virtual units_ptr_t get_units( const std::string& key ) {
                    check_pointer();
                    return _ptr->get_units( key );
                }

                virtual units_ptr_t get_units( const std::string& key ) const {
                    check_pointer();
                    return _ptr->get_units( key );
                }

                virtual double get_minimum_axis() const {
                    check_pointer();
                    return _ptr->get_minimum_axis();
                }

                virtual double get_maximum_axis() const {
                    check_pointer();
                    return _ptr->get_maximum_axis();
                }

                virtual Atmosphere1D& convert_axis_units(
                    units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_axis_units( new_units );
                    return *this;
                }

                virtual Atmosphere1D& convert_units( const std::string& key,
                                                     units_ptr_t new_units ) {
                    check_pointer();
                    _ptr->convert_units( key, new_units );
                    return *this;
                }

                virtual Atmosphere1D& resample( double new_dz ) {
                    check_pointer();
                    _ptr->resample( new_dz );
                    return *this;
                }

                virtual Atmosphere1D& resample( vector_u_t new_z ) {
                    check_pointer();
                    _ptr->resample( new_z );
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

                virtual void check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "Atmosphere1D: Internal pointer not set!" );
                    }
                }

                // friend binary operators
                friend std::ostream& operator<<( std::ostream& os,
                                                 const Atmosphere1D& atm ) {
                    if (atm) {
                        atm._ptr->print( os );
                    }
                    return os;
                }

                explicit operator bool() const {
                    return ( _ptr ? true : false );
                }

                abstract_atmosphere_1d *internal() { return _ptr.get(); }

                const abstract_atmosphere_1d *internal() const {
                    return _ptr.get();
                }

            private:
                std::unique_ptr<abstract_atmosphere_1d> _ptr;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::tuple_atmosphere_1d& a,
                  NCPA::atmos::tuple_atmosphere_1d& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::atmos::abstract_atmosphere_1d&>( a ),
            static_cast<NCPA::atmos::abstract_atmosphere_1d&>( b ) );
    swap( a._properties, b._properties );
    swap( a._scalar_properties, b._scalar_properties );
    swap( a._z, b._z );
    swap( a._interpolator_type, b._interpolator_type );
}

static void swap( NCPA::atmos::Atmosphere1D& a,
                  NCPA::atmos::Atmosphere1D& b ) noexcept {
    using std::swap;
    swap( static_cast<NCPA::atmos::AtmosphericModel&>( a ),
          static_cast<NCPA::atmos::AtmosphericModel&>( b ) );
    a._ptr.swap( b._ptr );
}
