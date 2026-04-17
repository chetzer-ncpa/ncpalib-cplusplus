#pragma once

#include "NCPA/atmosphere/abstract_atmospheric_property.hpp"
#include "NCPA/atmosphere/abstract_atmospheric_property_3d.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/exceptions.hpp"
#include "NCPA/interpolation.hpp"

#include <array>
#include <stdexcept>
#include <unordered_map>

static void swap( NCPA::atmos::AtmosphericProperty3D& a,
                  NCPA::atmos::AtmosphericProperty3D& b ) noexcept;

namespace NCPA {
    namespace atmos {

        class AtmosphericProperty3D {
            public:
                AtmosphericProperty3D() {}

                AtmosphericProperty3D( _atm_prop_3d_ptr_t engine ) {
                    _ptr = std::move( engine->clone3d() );
                }

                AtmosphericProperty3D( const AtmosphericProperty3D& other ) :
                    AtmosphericProperty3D() {
                    _ptr = std::move( other._ptr->clone3d() );
                }

                AtmosphericProperty3D(
                    const abstract_atmospheric_property& prop ) {
                    if (prop.dimensions() == 3) {
                        _ptr = std::move(
                            dynamic_cast<
                                const abstract_atmospheric_property_3d&>(
                                prop )
                                .clone3d() );
                    } else {
                        throw std::invalid_argument(
                            "AtmosphericProperty3D: Can only ingest 3-D "
                            "atmospheric properties" );
                    }
                }

                AtmosphericProperty3D& set(
                    const vector_u_t& ax1, const vector_u_t ax2,
                    const vector_u_t ax3,
                    const NCPA::arrays::ndvector<
                        2, const AtmosphericProperty1D *>& props ) {
                    _check_pointer();
                    if (props.shape()[ 0 ] != ax1.size()
                        || props.shape()[ 1 ] != ax2.size()) {
                        throw std::range_error(
                            "AtmosphericProperty3D.set(): Size mismatch "
                            "between supplied axes and supplied 2-D vector of "
                            "properties" );
                    }
                    vector3d_u_t grid( ax1.size(), ax2.size(), ax3.size(),
                                       props[ 0 ][ 0 ]->get_units() );
                    for (size_t i = 0; i < ax1.size(); ++i) {
                        for (size_t j = 0; j < ax2.size(); ++j) {
                            AtmosphericProperty1D prop( *props[ i ][ j ] );
                            prop.convert_units( *grid.get_units() );
                            prop.convert_axis_units( *ax3.get_units() );
                            // prop.resample( ax3 );
                            for (size_t k = 0; k < ax3.size(); ++k) {
                                grid[ i ][ j ][ k ] = prop.get( ax3[ k ] );
                            }
                        }
                    }
                    _ptr->set( ax1, ax2, ax3, grid );
                    return *this;
                }

                AtmosphericProperty3D& append(
                    size_t axis, const scalar_u_t& range,
                    const AtmosphericProperty2D& prop ) {
                    _check_pointer();
                    _ptr->append( axis,
                                  range.get_as( _ptr->get_axis_units( axis ) ),
                                  prop.internal()->values() );
                    return *this;
                }

                AtmosphericProperty3D(
                    AtmosphericProperty3D&& source ) noexcept :
                    AtmosphericProperty3D() {
                    ::swap( *this, source );
                }

                virtual ~AtmosphericProperty3D() {}

                AtmosphericProperty3D& operator=(
                    AtmosphericProperty3D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( AtmosphericProperty3D& a,
                                    AtmosphericProperty3D& b ) noexcept;

                virtual size_t size( size_t d ) const {
                    _check_pointer();
                    return _ptr->dim( d );
                }

                virtual AtmosphericProperty3D& set_interpolator(
                    NCPA::interpolation::interpolator_3d_type_t interp_type ) {
                    _check_pointer();
                    _ptr->set_interpolator( interp_type );
                    return *this;
                }

                virtual AtmosphericProperty3D& set_interpolator(
                    NCPA::interpolation::interpolator_2d_type_t interp_type ) {
                    _check_pointer();
                    _ptr->set_interpolator( interp_type );
                    return *this;
                }

                virtual AtmosphericProperty3D& set_interpolator(
                    NCPA::interpolation::interpolator_1d_type_t interp_type ) {
                    _check_pointer();
                    _ptr->set_interpolator( interp_type );
                    return *this;
                }

                virtual AtmosphericProperty3D& set( size_t dim,
                                                    const vector_u_t& ax ) {
                    _check_pointer();
                    _ptr->set( dim, ax );
                    return *this;
                }

                virtual AtmosphericProperty3D& set( const vector3d_u_t& ax ) {
                    _check_pointer();
                    _ptr->set( ax );
                    return *this;
                }

                virtual AtmosphericProperty3D& set( const vector_u_t& ax1,
                                                    const vector_u_t& ax2,
                                                    const vector_u_t& ax3,
                                                    const vector3d_u_t vals ) {
                    _check_pointer();
                    _ptr->set( ax1, ax2, ax3, vals );
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

                virtual vector3d_u_t& values() {
                    _check_pointer();
                    return _ptr->values();
                }

                virtual const vector3d_u_t& values() const {
                    _check_pointer();
                    return _ptr->values();
                }

                virtual double get( double val1, double val2, double val3 ) {
                    _check_pointer();
                    return _ptr->get( val1, val2, val3 );
                }

                virtual double f( double val1, double val2, double val3 ) {
                    _check_pointer();
                    return _ptr->f( val1, val2, val3 );
                }

                virtual vector3d_u_t get( const vector_u_t& v1,
                                          const vector_u_t& v2,
                                          const vector_u_t& v3 ) {
                    _check_pointer();
                    return _ptr->get( v1, v2, v3 );
                }

                virtual vector3d_u_t get( const std::vector<double>& v1,
                                          const std::vector<double>& v2,
                                          const std::vector<double>& v3 ) {
                    _check_pointer();
                    return _ptr->get( v1, v2, v3 );
                }

                virtual double get_first_derivative( double val1, double val2,
                                                     double val3,
                                                     size_t rel ) {
                    _check_pointer();
                    return _ptr->get_first_derivative( val1, val2, val3, rel );
                }

                virtual double df( double val1, double val2, double val3,
                                   size_t rel ) {
                    _check_pointer();
                    return _ptr->df( val1, val2, val3, rel );
                }

                virtual double get_second_derivative( double val1, double val2,
                                                      double val3, size_t rel1,
                                                      size_t rel2 ) {
                    _check_pointer();
                    return _ptr->get_second_derivative( val1, val2, val3, rel1,
                                                        rel2 );
                }

                virtual double ddf( double val1, double val2, double val3,
                                    size_t rel1, size_t rel2 ) {
                    _check_pointer();
                    return _ptr->ddf( val1, val2, val3, rel1, rel2 );
                }

                virtual const units_ptr_t get_units() const {
                    _check_pointer();
                    return _ptr->get_units();
                }

                virtual const units_ptr_t get_axis_units( size_t n ) const {
                    _check_pointer();
                    return _ptr->get_axis_units( n );
                }

                virtual AtmosphericProperty3D& convert_units(
                    const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_units( u );
                    return *this;
                }

                virtual AtmosphericProperty3D& convert_axis_units(
                    size_t n, const NCPA::units::Unit& u ) {
                    _check_pointer();
                    _ptr->convert_axis_units( n, u );
                    return *this;
                }

                virtual AtmosphericProperty3D& resample(
                    const vector_u_t& new_ax1, const vector_u_t& new_ax2,
                    const vector_u_t& new_ax3 ) {
                    _check_pointer();
                    _ptr->resample( new_ax1, new_ax2, new_ax3 );
                    return *this;
                }

                virtual AtmosphericProperty3D& resample(
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

                virtual AtmosphericProperty3D& set_limits( size_t dim,
                                                           double minax,
                                                           double maxax ) {
                    _check_pointer();
                    _ptr->set_limits( dim, minax, maxax );
                    return *this;
                }

                virtual _atm_prop_3d_ptr_t& internal() { return _ptr; }

                virtual const _atm_prop_3d_ptr_t& internal() const {
                    return _ptr;
                }


            protected:
                _atm_prop_3d_ptr_t _ptr;

                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphericProperty3D: no pointer has been set" );
                    }
                }
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphericProperty3D& a,
                  NCPA::atmos::AtmosphericProperty3D& b ) noexcept {
    using std::swap;
    swap( a._ptr, b._ptr );
}
