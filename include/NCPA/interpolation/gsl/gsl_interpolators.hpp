#pragma once

#ifndef NCPA_INTERPOLATION_IGNORE_GSL

#  include "NCPA/interpolation/defines.hpp"

#  ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE

#    include "gsl/gsl_interp.h"
#    include "gsl/gsl_spline.h"
#    include "gsl/gsl_version.h"
#    include "NCPA/arrays.hpp"
#    include "NCPA/defines.hpp"
#    include "NCPA/interpolation/abstract_spline_1d.hpp"
#    include "NCPA/interpolation/gsl/gsl_declarations.hpp"
#    include "NCPA/interpolation/types.hpp"
#    include "NCPA/math.hpp"
#    include "NCPA/types.hpp"

template<typename T, typename U>
static void swap( NCPA::interpolation::GSL::gsl_spline_1d<T, U>& a,
                  NCPA::interpolation::GSL::gsl_spline_1d<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {
        namespace GSL {

            // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE_WITH_PARAM(
            //     gsl_spline_1d, NCPA::interpolation::_spline_1d,
            //     const gsl_interp_type * );

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
                        _spline_size { 0 } {}

                    gsl_spline_1d( PARAMTYPE interptype ) :
                        gsl_spline_1d<INDEPTYPE, DEPTYPE, PARAMTYPE>() {
                        set_interpolation_type( interptype );
                    }

                    gsl_spline_1d(
                        const gsl_spline_1d<INDEPTYPE, DEPTYPE>& other ) :
                        gsl_spline_1d<INDEPTYPE, DEPTYPE>() {
                        _interptype = other._interptype;
                        if ( _interptype != nullptr && other._x.size() > 0
                             && other._f.size() > 0 ) {
                            this->init( other._x.size() );
                            this->fill( other._x, other._f );
                        }
                    }

                    gsl_spline_1d(
                        gsl_spline_1d<INDEPTYPE, DEPTYPE>&& source ) noexcept :
                        gsl_spline_1d<INDEPTYPE, DEPTYPE>() {
                        ::swap( *this, source );
                    }

                    friend void ::swap<INDEPTYPE, DEPTYPE>(
                        gsl_spline_1d<INDEPTYPE, DEPTYPE>& a,
                        gsl_spline_1d<INDEPTYPE, DEPTYPE>& b ) noexcept;

                    gsl_spline_1d<INDEPTYPE, DEPTYPE>& operator=(
                        gsl_spline_1d<INDEPTYPE, DEPTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual ~gsl_spline_1d() { clear(); }

                    virtual void set_interpolation_type(
                        PARAMTYPE interptype ) {
                        _interptype = interptype;
                    }

                    virtual void fill(
                        const std::vector<INDEPTYPE>& x,
                        const std::vector<DEPTYPE>& f ) override {
                        if ( x.size() != f.size() ) {
                            throw std::range_error(
                                "gsl_spline_1d.fill(): Input vectors not the "
                                "same size!" );
                        }
                        this->fill( x.size(), &x[ 0 ], &f[ 0 ] );
                    }

                    virtual void fill( size_t N, const INDEPTYPE *x,
                                       const DEPTYPE *f ) override {
                        if ( N != _spline_size ) {
                            init( N );
                        }
                        _set_spline( N, x, f );
                        _ready = true;
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

                    virtual std::pair<INDEPTYPE, INDEPTYPE> limits()
                        const override {
                        _check_ready();
                        return std::pair<INDEPTYPE, INDEPTYPE>(
                            { _x.front(), _x.back() } );
                    }

                    virtual interpolator_1d_type_t interptype() const override {
                        if (_interptype == gsl_interp_linear) {
                            return interpolator_1d_type_t::GSL_LINEAR;
                        } else if (_interptype == gsl_interp_cspline) {
                            return interpolator_1d_type_t::GSL_CUBIC;
                        } else if (_interptype == gsl_interp_polynomial) {
                            return interpolator_1d_type_t::GSL_POLYNOMIAL;
                        } else if (_interptype == gsl_interp_cspline_periodic) {
                            return interpolator_1d_type_t::GSL_CUBIC_PERIODIC;
                        } else if (_interptype == gsl_interp_akima) {
                            return interpolator_1d_type_t::GSL_AKIMA;
                        } else if (_interptype == gsl_interp_akima_periodic) {
                            return interpolator_1d_type_t::GSL_AKIMA_PERIODIC;
                        } 
#ifdef NCPA_INTERPOLATION_GSL_STEFFEN_SPLINE_AVAILABLE
#if NCPA_INTERPOLATION_GSL_STEFFEN_SPLINE_AVAILABLE
                        else if (_interptype == gsl_interp_steffen) {
                            return interpolator_1d_type_t::GSL_STEFFEN;
                        }
#endif
#endif
                        throw std::range_error( "Unknown or unimplemented interpolator" );
                    } 

                    virtual  std::vector<INDEPTYPE> source_x()
                        const override {
                        return _x;
                    }

                    virtual  std::vector<DEPTYPE> source_f()
                        const override {
                        return _f;
                    }

                protected:
                    void _allocate_spline( size_t n ) {
                        _check_interptype();
                        if ( n < _interptype->min_size ) {
                            std::ostringstream oss;
                            oss << "GSL interpolation type "
                                << _interptype->name << " requires at least "
                                << _interptype->min_size << " points.";
                            throw std::logic_error( oss.str() );
                        }
                        _spline      = gsl_spline_alloc( _interptype, n );
                        _accel       = gsl_interp_accel_alloc();
                        _spline_size = n;
                    }

                    void _set_spline( size_t n, const INDEPTYPE *x,
                                      const DEPTYPE *y ) {
                        _check_interptype();
                        double *xd = NCPA::arrays::zeros<double>( n );
                        double *yd = NCPA::arrays::zeros<double>( n );
                        _to_double<INDEPTYPE>( n, x, xd );
                        _to_double<DEPTYPE>( n, y, yd );
                        gsl_spline_init( _spline, xd, yd, n );
                        _x = std::vector<INDEPTYPE>( x, x + n );
                        _f = std::vector<DEPTYPE>( y, y + n );
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

                    PARAMTYPE _interptype    = nullptr;
                    gsl_spline *_spline      = nullptr;
                    gsl_interp_accel *_accel = nullptr;
                    bool _ready              = false;
                    size_t _spline_size      = 0;
                    std::vector<INDEPTYPE> _x;
                    std::vector<DEPTYPE> _f;
            };

            DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR_WITH_PARAM(
                gsl_spline_1d, NCPA::interpolation::_spline_1d,
                const gsl_interp_type * );

        }  // namespace GSL
    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap( NCPA::interpolation::GSL::gsl_spline_1d<T, U>& a,
                  NCPA::interpolation::GSL::gsl_spline_1d<T, U>& b ) noexcept {
    using std::swap;
    swap( a._interptype, b._interptype );
    swap( a._spline, b._spline );
    swap( a._accel, b._accel );
    swap( a._ready, b._ready );
    swap( a._spline_size, b._spline_size );
    swap( a._x, b._x );
    swap( a._f, b._f );
}

#  endif

#endif
