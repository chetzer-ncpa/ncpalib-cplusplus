#pragma once

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/abstract_atmospheric_property_1d.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation.hpp"

#include <memory>
#include <stdexcept>

static void swap( NCPA::atmos::AtmosphericProperty1D& a,
                  NCPA::atmos::AtmosphericProperty1D& b ) noexcept;

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
                    if (prop.dimensions() == 1) {
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

                AtmosphericProperty1D(
                    AtmosphericProperty1D&& source ) noexcept :
                    AtmosphericProperty1D() {
                    ::swap( *this, source );
                }

                virtual ~AtmosphericProperty1D() {}

                AtmosphericProperty1D& operator=(
                    AtmosphericProperty1D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( AtmosphericProperty1D& a,
                                    AtmosphericProperty1D& b ) noexcept;

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

                abstract_atmospheric_property_1d *internal() {
                    return _ptr.get();
                }

                const abstract_atmospheric_property_1d *internal() const {
                    return _ptr.get();
                }

            protected:
                _atm_prop_1d_ptr_t _ptr;

                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphericProperty1D: no pointer has been set" );
                    }
                }
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphericProperty1D& a,
                  NCPA::atmos::AtmosphericProperty1D& b ) noexcept {
    using std::swap;
    swap( a._ptr, b._ptr );
}
