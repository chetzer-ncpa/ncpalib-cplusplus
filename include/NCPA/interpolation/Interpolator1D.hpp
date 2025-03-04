#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/builders.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/Extrapolator1D.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <memory>
#include <stdexcept>
#include <utility>

#ifndef NCPA_INTERPOLATION_DEFAULT_1D_EXTRAPOLATOR
#define NCPA_INTERPOLATION_DEFAULT_1D_EXTRAPOLATOR FORBIDDEN
#endif

template<typename T, typename U>
static void swap( NCPA::interpolation::Interpolator1D<T, U>& a,
                  NCPA::interpolation::Interpolator1D<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {

        template<typename INDEPTYPE, typename DEPTYPE>
        class Interpolator1D {
            public:
                Interpolator1D() {
                    this->set_extrapolation(
                        extrapolator_1d_type_t::NCPA_INTERPOLATION_DEFAULT_1D_EXTRAPOLATOR );
                }

                Interpolator1D(
                    spline_engine_1d_t<INDEPTYPE, DEPTYPE> engine ) :
                    Interpolator1D<INDEPTYPE, DEPTYPE>() {
                    set_engine( engine );
                }

                Interpolator1D(
                    const Interpolator1D<INDEPTYPE, DEPTYPE>& other ) :
                    Interpolator1D<INDEPTYPE, DEPTYPE>() {
                    _engine = InterpolatorFactory<INDEPTYPE,DEPTYPE>::build( other.interptype() );
                    std::vector<INDEPTYPE> x = other._engine->source_x();
                    std::vector<DEPTYPE> f = other._engine->source_f();
                    this->init(x.size()).fill( x, f ).ready();
                }

                Interpolator1D(
                    Interpolator1D<INDEPTYPE, DEPTYPE>&& source ) noexcept :
                    Interpolator1D<INDEPTYPE, DEPTYPE>() {
                    ::swap( *this, source );
                }

                friend void ::swap<INDEPTYPE, DEPTYPE>(
                    Interpolator1D<INDEPTYPE, DEPTYPE>& a,
                    Interpolator1D<INDEPTYPE, DEPTYPE>& b ) noexcept;

                Interpolator1D<INDEPTYPE, DEPTYPE>& operator=(
                    Interpolator1D<INDEPTYPE, DEPTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual ~Interpolator1D() {}

                virtual Interpolator1D<INDEPTYPE, DEPTYPE>& set_engine(
                    spline_engine_1d_t<INDEPTYPE, DEPTYPE> engine ) {
                    _engine = std::move( engine );
                    return *this;
                }

                virtual Interpolator1D<INDEPTYPE, DEPTYPE>& set_extrapolation(
                    extrapolator_1d_type_t extraptype ) {
                    _extrap = InterpolatorFactory<INDEPTYPE, DEPTYPE>::build(
                        extraptype );
                    return *this;
                }

                virtual Interpolator1D<INDEPTYPE, DEPTYPE>& fill(
                    size_t N, const INDEPTYPE *x, const DEPTYPE *f ) {
                    check_engine();
                    _engine->fill( N, x, f );
                    return *this;
                }

                virtual Interpolator1D<INDEPTYPE, DEPTYPE>& fill(
                    const std::vector<INDEPTYPE>& x,
                    const std::vector<DEPTYPE>& f ) {
                    check_engine();
                    _engine->fill( x, f );
                    return *this;
                }

                virtual Interpolator1D<INDEPTYPE, DEPTYPE>& init( size_t N ) {
                    check_engine();
                    _engine->init( N );
                    return *this;
                }

                virtual Interpolator1D<INDEPTYPE, DEPTYPE>& clear() {
                    _engine.reset();
                    return *this;
                }

                virtual Interpolator1D<INDEPTYPE, DEPTYPE>& ready() {
                    check_engine();
                    _engine->ready();
                    return *this;
                }

                virtual DEPTYPE eval_f( INDEPTYPE x ) {
                    check_engine();
                    return ( this->_in_limits( x )
                                 ? _engine->eval_f( x )
                                 : _extrap.extrapolate( x, *this ) );
                }

                virtual DEPTYPE eval_df( INDEPTYPE x ) {
                    check_engine();
                    // return _engine->eval_df( x );
                    return ( this->_in_limits( x )
                                 ? _engine->eval_df( x )
                                 : _extrap.extrapolate_df( x, *this ) );
                }

                virtual DEPTYPE eval_ddf( INDEPTYPE x ) {
                    check_engine();
                    // return _engine->eval_ddf( x );
                    return ( this->_in_limits( x )
                                 ? _engine->eval_ddf( x )
                                 : _extrap.extrapolate_ddf( x, *this ) );
                }

                virtual DEPTYPE eval_dddf( INDEPTYPE x ) {
                    check_engine();
                    // return _engine->eval_dddf( x );
                    return ( this->_in_limits( x )
                                 ? _engine->eval_dddf( x )
                                 : _extrap.extrapolate_dddf( x, *this ) );
                }

                virtual std::pair<INDEPTYPE, INDEPTYPE> limits() const {
                    check_engine();
                    return _engine->limits();
                }

                virtual void check_engine() const {
                    if ( !_engine ) {
                        throw std::logic_error(
                            "Interpolator1D: no engine has been set" );
                    }
                }

                explicit operator bool() const {
                    return ( _engine ? true : false );
                }

                virtual interpolator_1d_type_t interptype() const {
                    check_engine();
                    return _engine->interptype();
                }

            protected:
                spline_engine_1d_t<INDEPTYPE, DEPTYPE> _engine;
                Extrapolator1D<INDEPTYPE, DEPTYPE> _extrap;

                bool _in_limits( INDEPTYPE x ) const {
                    std::pair<INDEPTYPE, INDEPTYPE> lim = this->limits();
                    return ( x >= lim.first && x <= lim.second );
                }
        };

    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap( NCPA::interpolation::Interpolator1D<T, U>& a,
                  NCPA::interpolation::Interpolator1D<T, U>& b ) noexcept {
    using std::swap;
    swap( a._engine, b._engine );
}
