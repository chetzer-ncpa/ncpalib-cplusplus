#pragma once

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/abstract_atmospheric_property_2d.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/exceptions.hpp"
#include "NCPA/interpolation.hpp"

#include <array>
#include <stdexcept>
#include <unordered_map>

static void swap( NCPA::atmos::AtmosphericProperty2D& a,
                  NCPA::atmos::AtmosphericProperty2D& b ) noexcept;

namespace NCPA {
    namespace atmos {

        class AtmosphericProperty2D {
            public:
                AtmosphericProperty2D() {}

                AtmosphericProperty2D( _atm_prop_2d_ptr_t engine ) {
                    _ptr = std::move( engine->clone2d() );
                }

                AtmosphericProperty2D( const AtmosphericProperty2D& other ) :
                    AtmosphericProperty2D() {
                    _ptr = std::move( other._ptr->clone2d() );
                }

                AtmosphericProperty2D(
                    const abstract_atmospheric_property& prop ) {
                    if (prop.dimensions() == 2) {
                        _ptr = std::move(
                            dynamic_cast<
                                const abstract_atmospheric_property_2d&>(
                                prop )
                                .clone2d() );
                    } else {
                        throw std::invalid_argument(
                            "AtmosphericProperty2D: Can only ingest 2-D "
                            "atmospheric properties" );
                    }
                }

                AtmosphericProperty2D& set(
                    const vector_u_t& ranges,
                    const std::vector<const AtmosphericProperty1D *>& props ) {
                    _check_pointer();
                    std::vector<const abstract_atmospheric_property_1d *> ptrs;
                    for (auto it = props.cbegin(); it != props.cend(); ++it) {
                        ptrs.push_back( ( *it )->internal() );
                    }
                    _ptr->set( ranges, ptrs );
                    return *this;
                }

                AtmosphericProperty2D& set(
                    const scalar_u_t& range,
                    const AtmosphericProperty1D& prop ) {
                    _check_pointer();
                    return this->set(
                        vector_u_t( 1, range ),
                        std::vector<const AtmosphericProperty1D *>( 1,
                                                                    &prop ) );
                }

                AtmosphericProperty2D& set(
                    const AtmosphericProperty1D& prop ) {
                    _check_pointer();
                    return this->set(
                        scalar_u_t( 0.0, this->get_axis_units( 0 ) ), prop );
                }

                AtmosphericProperty2D& append(
                    const scalar_u_t& range,
                    const AtmosphericProperty1D& prop ) {
                    _check_pointer();
                    _ptr->append( range, *prop.internal() );
                    return *this;
                }

                virtual AtmosphericProperty1D extract( double range ) const {
                    _check_pointer();
                    auto ptr = _ptr->extract( range );
                    return AtmosphericProperty1D( _ptr->extract( range ) );
                }

                AtmosphericProperty2D(
                    AtmosphericProperty2D&& source ) noexcept :
                    AtmosphericProperty2D() {
                    ::swap( *this, source );
                }

                virtual ~AtmosphericProperty2D() {}

                AtmosphericProperty2D& operator=(
                    AtmosphericProperty2D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( AtmosphericProperty2D& a,
                                    AtmosphericProperty2D& b ) noexcept;

                virtual size_t size( size_t d ) const {
                    _check_pointer();
                    return _ptr->size( d );
                }

                virtual AtmosphericProperty2D& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type ) {
                    _check_pointer();
                    _ptr->set_interpolator( interp_type );
                    return *this;
                }

                virtual AtmosphericProperty2D& set( size_t dim,
                                                    const vector_u_t& ax ) {
                    _check_pointer();
                    _ptr->set( dim, ax );
                    return *this;
                }

                virtual AtmosphericProperty2D& set( const vector2d_u_t& ax ) {
                    _check_pointer();
                    _ptr->set( ax );
                    return *this;
                }

                virtual AtmosphericProperty2D& set( const vector_u_t& ax1,
                                                    const vector_u_t& ax2,
                                                    const vector2d_u_t vals ) {
                    _check_pointer();
                    _ptr->set( ax1, ax2, vals );
                    return *this;
                }

                virtual vector_u_t& axis( size_t n ) {
                    _check_pointer();
                    return _ptr->axis( n );
                }

                virtual const vector_u_t& axis( size_t n ) const {
                    _check_pointer();
                    return _ptr->axis( n );
                }

                virtual vector2d_u_t& values() {
                    _check_pointer();
                    return _ptr->values();
                }

                virtual const vector2d_u_t& values() const {
                    _check_pointer();
                    return _ptr->values();
                }

                virtual double get( double val1, double val2 ) {
                    _check_pointer();
                    return _ptr->get( val1, val2 );
                }

                virtual double f( double val1, double val2 ) {
                    _check_pointer();
                    return _ptr->f( val1, val2 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     size_t rel ) {
                    _check_pointer();
                    return _ptr->get_first_derivative( val1, val2, rel );
                }

                virtual double df( double val1, double val2, size_t rel ) {
                    _check_pointer();
                    return _ptr->df( val1, val2, rel );
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      size_t rel1,
                                                      size_t rel2 ) {
                    _check_pointer();
                    return _ptr->get_second_derivative( val1, val2, rel1,
                                                        rel2 );
                }

                virtual double ddf( double val1, double val2, size_t rel1,
                                    size_t rel2 ) {
                    _check_pointer();
                    return _ptr->ddf( val1, val2, rel1, rel2 );
                }

                virtual const units_ptr_t get_units() const {
                    _check_pointer();
                    return _ptr->get_units();
                }

                virtual const units_ptr_t get_axis_units( size_t n ) const {
                    _check_pointer();
                    return _ptr->get_axis_units( n );
                }

                virtual AtmosphericProperty2D& convert_units(
                    const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_units( u );
                    return *this;
                }

                virtual AtmosphericProperty2D& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_axis_units( n, u );
                    return *this;
                }

                virtual AtmosphericProperty2D& resample(
                    const vector_u_t& new_r, const vector_u_t& new_z ) {
                    _check_pointer();
                    _ptr->resample( new_r, new_z );
                    return *this;
                }

                virtual AtmosphericProperty2D& resample(
                    size_t axis, const vector_u_t& new_ax ) {
                    _check_pointer();
                    _ptr->resample( axis, new_ax );
                    return *this;
                }

                virtual const std::pair<double, double> get_limits(
                    size_t dim ) const {
                    _check_pointer();
                    return _ptr->get_limits( dim );
                }

                virtual AtmosphericProperty2D& set_limits( size_t dim,
                                                           double minax,
                                                           double maxax ) {
                    _check_pointer();
                    _ptr->set_limits( dim, minax, maxax );
                    return *this;
                }

                virtual _atm_prop_2d_ptr_t& internal() { return _ptr; }

                virtual const _atm_prop_2d_ptr_t& internal() const {
                    return _ptr;
                }


            protected:
                _atm_prop_2d_ptr_t _ptr;

                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphericProperty2D: no pointer has been set" );
                    }
                }
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphericProperty2D& a,
                  NCPA::atmos::AtmosphericProperty2D& b ) noexcept {
    using std::swap;
    swap( a._ptr, b._ptr );
}
