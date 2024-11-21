#pragma once

#include "NCPA/atmosphere/types.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/defines.hpp"

#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace atmos {
        class AtmosphericProperty1D {
            public:
                AtmosphericProperty1D() {}

                AtmosphericProperty1D( const vector_t& z, const vector_t& p ) :
                    AtmosphericProperty1D() {
                    set( z, p );
                    // if ( z.size() != p.size() ) {
                    //     throw std::invalid_argument(
                    //         "Altitude and property vectors are not the same "
                    //         "size!" );
                    // }
                    // _z = z;
                    // _x = p;
                    // _init_spline();
                }

                AtmosphericProperty1D( const AtmosphericProperty1D& source ) :
                    AtmosphericProperty1D() {
                    _z = source._z;
                    _x = source._x;
                    _init_spline();
                }

                virtual ~AtmosphericProperty1D() {}

                virtual size_t size() const {
                    return _z.size();
                }

                virtual AtmosphericProperty1D& set_interpolator(
                    NCPA::interpolation::interpolator_type_t interp_type ) {
                    _spline = NCPA::interpolation::InterpolatorFactory::build<
                        double, double>( interp_type );
                    _init_spline();
                    return *this;
                }

                virtual AtmosphericProperty1D& set(const vector_t& z, const vector_t& p ) {
                    if ( z.size() != p.size() ) {
                        throw std::invalid_argument(
                            "Altitude and property vectors are not the same "
                            "size!" );
                    }
                    _z = z;
                    _x = p;
                    _init_spline();
                }

                scalar_t get( scalar_t altitude ) {
                    return scalar_t(
                        _spline.eval_f( altitude.get_as(* _z.get_units()) ),
                        *_z.get_units() );
                }
		        scalar_t get_first_derivative( scalar_t altitude ) {
                    return scalar_t(
                        _spline.eval_df( altitude.get_as(* _z.get_units()) ),
                        *_z.get_units() 
                    );
                }
		        scalar_t get_second_derivative( scalar_t altitude ) {
                    return scalar_t(
                        _spline.eval_ddf( altitude.get_as(* _z.get_units()) ),
                        *_z.get_units() 
                    );
                }



            protected:
                void _init_spline() {
                    if (_z && _x && _spline) {
                        _spline.init(size());
                        _spline.fill( _z.as_std_vector(), _x.as_std_vector() );
                    }
                }

            private:
                NCPA::units::VectorWithUnits<double> _z, _x;
                NCPA::interpolation::Interpolator1D<double, double> _spline;
        };
    }  // namespace atmos
}  // namespace NCPA
