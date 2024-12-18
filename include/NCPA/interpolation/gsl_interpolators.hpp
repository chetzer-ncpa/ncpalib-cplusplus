#pragma once

#ifndef NCPA_INTERPOLATION_IGNORE_GSL

#include "NCPA/interpolation/defines.hpp"

#ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_version.h"
#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/abstract_spline_1d.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

namespace NCPA {
    namespace interpolation {
        namespace GSL {

            DECLARE_GENERIC_INTERPOLATOR_TEMPLATE_WITH_PARAM(
                gsl_spline_1d, NCPA::interpolation::_spline_1d,
                const gsl_interp_type * );

            _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION_WITH_PARAM  //
                class gsl_spline_1d<INDEPTYPE, DEPTYPE, PARAMTYPE,
                                    ENABLE_IF_REAL( INDEPTYPE ),
                                    ENABLE_IF_REAL( DEPTYPE )>
                : public NCPA::interpolation::_spline_1d<INDEPTYPE, DEPTYPE> {
                public:
                    gsl_spline_1d() :
                        NCPA::interpolation::_spline_1d<INDEPTYPE, DEPTYPE>(),
                        _interptype { nullptr },
                        _spline { nullptr },
                        _accel { nullptr },
                        _ready { false },
                        _xmin { 0.0 },
                        _xmax { 0.0 } {}

                    gsl_spline_1d( PARAMTYPE interptype ) :
                        gsl_spline_1d<INDEPTYPE, DEPTYPE, PARAMTYPE>() {
                        set_interpolation_type( interptype );
                    }

                    virtual ~gsl_spline_1d() { clear(); }

                    virtual void set_interpolation_type(
                        PARAMTYPE interptype ) {
                        _interptype = interptype;
                    }

                    virtual void fill( size_t N, const INDEPTYPE *x,
                                       const DEPTYPE *f ) override {
                        if ( N != _size ) {
                            init( N );
                        }
                        _set_spline( N, x, f );
                        _ready = true;
                    }

                    virtual void fill(
                        const std::vector<INDEPTYPE>& x,
                        const std::vector<DEPTYPE>& f ) override {
                        size_t N = x.size();
                        if ( N != f.size() ) {
                            throw std::invalid_argument(
                                "Vectors must be same size" );
                        }
                        fill( N, &x[ 0 ], &f[ 0 ] );
                    }

                    virtual void init( size_t N ) override {
                        clear();
                        _check_interptype();
                        _allocate_spline( N );
                    }

                    virtual void clear() override {
                        _free_spline();
                        _free_accel();
                    }

                    virtual void ready() override {
                        _check_spline();
                        _check_accel();
                    }

                    DEPTYPE eval_f( INDEPTYPE x ) override {
                        _check_ready();
                        return (DEPTYPE)gsl_spline_eval( _spline, (double)x,
                                                         _accel );
                    }

                    DEPTYPE eval_df( INDEPTYPE x ) override {
                        _check_ready();
                        return (DEPTYPE)gsl_spline_eval_deriv(
                            _spline, (double)x, _accel );
                    }

                    DEPTYPE eval_ddf( INDEPTYPE x ) override {
                        _check_ready();
                        return (DEPTYPE)gsl_spline_eval_deriv2(
                            _spline, (double)x, _accel );
                    }

                    DEPTYPE eval_dddf( INDEPTYPE x ) override {
                        throw std::out_of_range( "GSL splines have no third "
                                                 "derivative capability." );
                    }

                protected:
                    void _allocate_spline( size_t n ) {
                        if ( n < _interptype->min_size ) {
                            std::ostringstream oss;
                            oss << "GSL interpolation type "
                                << _interptype->name << " requires at least "
                                << _interptype->min_size << " points.";
                            throw std::logic_error( oss.str() );
                        }
                        _spline = gsl_spline_alloc( _interptype, n );
                        _accel  = gsl_interp_accel_alloc();
                        _size   = n;
                    }

                    void _set_spline( size_t n, const INDEPTYPE *x,
                                      const DEPTYPE *y ) {
                        double *xd = NCPA::arrays::zeros<double>( n );
                        double *yd = NCPA::arrays::zeros<double>( n );
                        _to_double<INDEPTYPE>( n, x, xd );
                        _to_double<DEPTYPE>( n, y, yd );
                        gsl_spline_init( _spline, xd, yd, n );
                        _xmin = (INDEPTYPE)x[ 0 ];
                        _xmax = (INDEPTYPE)x[ n - 1 ];
                        delete[] xd;
                        delete[] yd;
                    }

                    template<typename T>
                    void _to_double( size_t n, const T *input,
                                     double *& output ) {
                        for ( auto i = 0; i < n; i++ ) {
                            output[ i ] = (double)( input[ i ] );
                        }
                    }

                    // template<>
                    // void _to_double( size_t n, const double *input,
                    //                  double *& output ) {
                    //     std::memcpy( output, input, n * sizeof( double ) );
                    // }

                    void _free_interptype() { _interptype = nullptr; }

                    void _free_spline() {
                        if ( _spline != nullptr ) {
                            gsl_spline_free( _spline );
                            _spline = nullptr;
                        }
                        _free_accel();
                        _ready = false;
                    }

                    void _free_accel() {
                        if ( _accel != nullptr ) {
                            gsl_interp_accel_free( _accel );
                            _accel = nullptr;
                        }
                        _ready = false;
                    }

                    void _check_ready() const {
                        if ( !_ready ) {
                            throw std::logic_error( "Spline not ready" );
                        }
                    }

                    void _check_interptype() const {
                        if ( _interptype == nullptr ) {
                            throw std::logic_error(
                                "GSL Spline type has not been set!" );
                        }
                    }

                    void _check_spline() const {
                        if ( _spline == nullptr ) {
                            throw std::logic_error(
                                "Spline has not been initialized!" );
                        }
                    }

                    void _check_accel() const {
                        if ( _accel == nullptr ) {
                            throw std::logic_error(
                                "Spline accelerator has not been "
                                "initialized" );
                        }
                    }


                private:
                    PARAMTYPE _interptype    = nullptr;
                    gsl_spline *_spline      = nullptr;
                    gsl_interp_accel *_accel = nullptr;
                    bool _ready              = false;
                    INDEPTYPE _xmin = 0.0, _xmax = 0.0;
                    size_t _size = 0;
            };

            DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAM(
                gsl_spline_1d, NCPA::interpolation::_spline_1d,
                const gsl_interp_type * );

        }  // namespace GSL
    }
}

#endif

#endif