#pragma once

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace atmos {
        class abstract_atmospheric_property_1d;
        class tuple_atmospheric_property_1d;
    }  // namespace atmos
}  // namespace NCPA

#define RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D \
    return static_cast<abstract_atmospheric_property_1d&>( *this );

static void swap( NCPA::atmos::abstract_atmospheric_property_1d& a,
                  NCPA::atmos::abstract_atmospheric_property_1d& b ) noexcept;
static void swap( NCPA::atmos::tuple_atmospheric_property_1d& a,
                  NCPA::atmos::tuple_atmospheric_property_1d& b ) noexcept;
static void swap( NCPA::atmos::AtmosphericProperty1D& a,
                  NCPA::atmos::AtmosphericProperty1D& b ) noexcept;

namespace NCPA {
    namespace atmos {
        class abstract_atmospheric_property_1d
            : public abstract_atmospheric_property {
            public:
                abstract_atmospheric_property_1d() {}

                virtual ~abstract_atmospheric_property_1d() {}

                friend void ::swap(
                    abstract_atmospheric_property_1d& a,
                    abstract_atmospheric_property_1d& b ) noexcept;

                virtual size_t size() const = 0;
                virtual std::unique_ptr<abstract_atmospheric_property_1d>
                    clone1d() const = 0;
                virtual abstract_atmospheric_property_1d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    = 0;
                virtual abstract_atmospheric_property_1d& set_extrapolator(
                    NCPA::interpolation::extrapolator_1d_type_t extrap_type )
                    = 0;
                virtual NCPA::interpolation::Interpolator1D<double, double>&
                    interpolator()
                    = 0;
                virtual abstract_atmospheric_property_1d& set(
                    const vector_u_t& z, const vector_u_t& p )
                    = 0;
                virtual vector_u_t& values()                            = 0;
                virtual const vector_u_t& values() const                = 0;
                virtual vector_u_t& axis()                              = 0;
                virtual const vector_u_t& axis() const                  = 0;
                virtual double get( double altitude )                   = 0;
                virtual double get_first_derivative( double altitude )  = 0;
                virtual double get_second_derivative( double altitude ) = 0;
                virtual const units_ptr_t get_units() const             = 0;
                virtual const units_ptr_t get_axis_units() const        = 0;
                virtual abstract_atmospheric_property_1d& convert_units(
                    const NCPA::units::Unit& u )
                    = 0;
                virtual abstract_atmospheric_property_1d& convert_axis_units(
                    const NCPA::units::Unit& u )
                    = 0;
                virtual abstract_atmospheric_property_1d& resample(
                    vector_u_t new_z )
                    = 0;

                virtual size_t dimensions() const override { return 1; }
        };

        typedef std::unique_ptr<abstract_atmospheric_property_1d>
            _atm_prop_1d_ptr_t;

        class tuple_atmospheric_property_1d
            : public abstract_atmospheric_property_1d {
            public:
                tuple_atmospheric_property_1d() {
                    set_interpolator(
                        NCPA::interpolation::interpolator_1d_type_t::
                            NEAREST_NEIGHBOR );
                    set_extrapolator( NCPA::interpolation::extrapolator_1d_type_t::CONSTANT );
                }

                tuple_atmospheric_property_1d( const vector_u_t& z,
                                               const vector_u_t& p ) :
                    tuple_atmospheric_property_1d() {
                    set( z, p );
                }

                tuple_atmospheric_property_1d(
                    const tuple_atmospheric_property_1d& source ) :
                    tuple_atmospheric_property_1d() {
                    set( source._z, source._x );
                }

                DECLARE_BOILERPLATE_METHODS( tuple_atmospheric_property_1d,
                                             abstract_atmospheric_property )

                virtual std::unique_ptr<abstract_atmospheric_property_1d>
                    clone1d() const override {
                    return std::unique_ptr<abstract_atmospheric_property_1d>(
                        new tuple_atmospheric_property_1d( *this ) );
                }

                virtual size_t size() const override { return _z.size(); }

                virtual abstract_atmospheric_property& copy(
                    const abstract_atmospheric_property& source ) override {
                    if ( source.dimensions() == 1 ) {
                        const abstract_atmospheric_property_1d *srcptr
                            = dynamic_cast<
                                const abstract_atmospheric_property_1d *>(
                                &source );
                        set( srcptr->axis(), srcptr->values() );
                    } else {
                        throw std::invalid_argument(
                            "tuple_atmospheric_property_1d.copy(): Cannot "
                            "copy property with more than 1 dimension" );
                    }
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D;
                }

                virtual abstract_atmospheric_property_1d& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type )
                    override {
                    _spline = NCPA::interpolation::InterpolatorFactory<
                        double, double>::build( interp_type );
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D;
                }

                virtual abstract_atmospheric_property_1d& set_extrapolator(
                    NCPA::interpolation::extrapolator_1d_type_t extrap_type )
                    override {
                    this->interpolator().set_extrapolation( extrap_type );
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D;
                }

                virtual NCPA::interpolation::Interpolator1D<double, double>&
                    interpolator() override {
                    return _spline;
                }

                virtual abstract_atmospheric_property_1d& set(
                    const vector_u_t& z, const vector_u_t& p ) override {
                    if ( z.size() != p.size() ) {
                        throw std::invalid_argument(
                            "Altitude and property vectors are not the same "
                            "size!" );
                    }
                    _z = z;
                    _x = p;
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D;
                }

                virtual vector_u_t& values() override { return _x; }

                virtual const vector_u_t& values() const override {
                    return _x;
                }

                virtual vector_u_t& axis() override { return _z; }

                virtual const vector_u_t& axis() const override { return _z; }

                virtual double get( double altitude ) override {
                    return _spline.eval_f( altitude );
                }

                virtual double get_first_derivative(
                    double altitude ) override {
                    return _spline.eval_df( altitude );
                }

                virtual double get_second_derivative(
                    double altitude ) override {
                    return _spline.eval_ddf( altitude );
                }

                virtual const units_ptr_t get_units() const override {
                    return _x.get_units();
                }

                virtual const units_ptr_t get_axis_units() const override {
                    return _z.get_units();
                }

                virtual abstract_atmospheric_property_1d& convert_units(
                    const NCPA::units::Unit& u ) override {
                    _x.convert_units( u );
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D;
                }

                virtual abstract_atmospheric_property_1d& convert_axis_units(
                    const NCPA::units::Unit& u ) override {
                    _z.convert_units( u );
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D;
                }

                virtual abstract_atmospheric_property_1d& resample(
                    vector_u_t new_z ) override {
                    if ( !_spline ) {
                        throw std::logic_error(
                            "Cannot resample, interpolator has not been "
                            "set." );
                    }
                    _z.convert_units( *new_z.get_units() );
                    vector_u_t new_x( new_z.size() );
                    new_x.set_units( *_x.get_units() );
                    for ( size_t i = 0; i < new_z.size(); i++ ) {
                        new_x[ i ] = this->get( new_z[ i ] );
                    }
                    _x = new_x;
                    _z = new_z;
                    _init_spline();
                    RETURN_THIS_AS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D;
                }

            protected:
                void _init_spline() {
                    if ( _z && _x && _spline ) {
                        _spline.init( size() ).fill( _z, _x ).ready();
                    }
                }

            protected:
                vector_u_t _z, _x;
                NCPA::interpolation::Interpolator1D<double, double> _spline;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract_atmospheric_property_1d& a,
                  NCPA::atmos::abstract_atmospheric_property_1d& b ) noexcept {
}

static void swap( NCPA::atmos::tuple_atmospheric_property_1d& a,
                  NCPA::atmos::tuple_atmospheric_property_1d& b ) noexcept {
    using std::swap;
    ::swap(
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_1d&>( a ),
        dynamic_cast<NCPA::atmos::abstract_atmospheric_property_1d&>( b ) );
    swap( a._z, b._z );
    swap( a._x, b._x );
    swap( a._spline, b._spline );
}

namespace NCPA {
    namespace atmos {
        class AtmosphericProperty1D {
            public:
                AtmosphericProperty1D() {}

                AtmosphericProperty1D( _atm_prop_1d_ptr_t engine ) {
                    _ptr = std::move( engine->clone1d() );
                }

                AtmosphericProperty1D( const AtmosphericProperty1D& other ) :
                    AtmosphericProperty1D() {
                    _ptr = std::move( other._ptr->clone1d() );
                }

                AtmosphericProperty1D(
                    const abstract_atmospheric_property& prop ) {
                    if ( prop.dimensions() == 1 ) {
                        _ptr = std::move(
                            dynamic_cast<
                                const abstract_atmospheric_property_1d&>(
                                prop )
                                .clone1d() );
                    } else {
                        throw std::invalid_argument(
                            "AtmosphericProperty1D: Can only ingest 1-D "
                            "atmospheric properties" );
                    }
                }

                DECLARE_WRAPPER_BOILERPLATE_METHODS( AtmosphericProperty1D )

                virtual size_t size() const {
                    this->_check_pointer();
                    return _ptr->size();
                }

                virtual const vector_u_t& axis() const {
                    this->_check_pointer();
                    return _ptr->axis();
                }

                virtual vector_u_t& values() {
                    this->_check_pointer();
                    return _ptr->values();
                }

                virtual const vector_u_t& values() const {
                    this->_check_pointer();
                    return _ptr->values();
                }

                virtual vector_u_t& axis() {
                    this->_check_pointer();
                    return _ptr->axis();
                }

                virtual AtmosphericProperty1D& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type ) {
                    _check_pointer();
                    _ptr->set_interpolator( interp_type );
                    return *this;
                }

                virtual AtmosphericProperty1D& set( const vector_u_t& z,
                                                    const vector_u_t& p ) {
                    _check_pointer();
                    _ptr->set( z, p );
                    return *this;
                }

                // virtual vector_u_t& values() {
                //     _check_pointer();
                //     return _ptr->values();
                // }

                // virtual const vector_u_t& values() const {
                //     _check_pointer();
                //     return _ptr->values();
                // }

                virtual double get( double x ) {
                    _check_pointer();
                    return _ptr->get( x );
                }

                virtual double get_first_derivative( double x ) {
                    _check_pointer();
                    return _ptr->get_first_derivative( x );
                }

                virtual double get_second_derivative( double x ) {
                    _check_pointer();
                    return _ptr->get_second_derivative( x );
                }

                virtual const units_ptr_t get_units() const {
                    _check_pointer();
                    return _ptr->get_units();
                }

                virtual const units_ptr_t get_axis_units() const {
                    _check_pointer();
                    return _ptr->get_axis_units();
                }

                virtual AtmosphericProperty1D& convert_units(
                    const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_units( u );
                    return *this;
                }

                virtual AtmosphericProperty1D& convert_axis_units(
                    const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_axis_units( u );
                    return *this;
                }

                virtual AtmosphericProperty1D& resample( vector_u_t new_z ) {
                    _check_pointer();
                    _ptr->resample( new_z );
                    return *this;
                }

                explicit operator bool() const {
                    return ( _ptr ? true : false );
                }

                void set_pointer( _atm_prop_1d_ptr_t newptr ) {
                    _ptr = std::move( newptr );
                }

                abstract_atmospheric_property_1d* internal() {
                    return _ptr.get();
                }

                const abstract_atmospheric_property_1d* internal() const {
                    return _ptr.get();
                }

            protected:
                _atm_prop_1d_ptr_t _ptr;

                void _check_pointer() const {
                    if ( !_ptr ) {
                        throw std::logic_error(
                            "AtmosphericProperty1D: no pointer has been set" );
                    }
                }
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap(
    NCPA::atmos::AtmosphericProperty1D& a,
    NCPA::atmos::AtmosphericProperty1D& b ) noexcept {
    using std::swap;
    swap( a._ptr, b._ptr );
}