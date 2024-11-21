#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/abstract_spline_1d.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace interpolation {

        template<typename INDEPTYPE, typename DEPTYPE>
        class Interpolator1D {
            public:
                Interpolator1D() {}

                Interpolator1D(
                    std::unique_ptr<_spline_1d<INDEPTYPE, DEPTYPE>> engine ) {
                    set_engine( engine );
                }

                void set_engine(
                    std::unique_ptr<_spline_1d<INDEPTYPE, DEPTYPE>> engine ) {
                    _engine = std::move( engine );
                }

                virtual void fill( size_t N, const INDEPTYPE *x,
                                   const DEPTYPE *f ) {
                    check_engine();
                    _engine->fill( N, x, f );
                }

                virtual void fill( const std::vector<INDEPTYPE>& x,
                                   const std::vector<DEPTYPE>& f ) {
                    check_engine();
                    _engine->fill( x, f );
                }

                virtual void init( size_t N ) {
                    check_engine();
                    _engine->init( N );
                }

                virtual void clear() { _engine.reset(); }

                virtual void ready() {
                    check_engine();
                    _engine->ready();
                }

                virtual DEPTYPE eval_f( INDEPTYPE x ) {
                    check_engine();
                    return _engine->eval_f( x );
                }

                virtual DEPTYPE eval_df( INDEPTYPE x ) {
                    check_engine();
                    return _engine->eval_df( x );
                }

                virtual DEPTYPE eval_ddf( INDEPTYPE x ) {
                    check_engine();
                    return _engine->eval_ddf( x );
                }

                virtual DEPTYPE eval_dddf( INDEPTYPE x ) {
                    check_engine();
                    return _engine->eval_dddf( x );
                }

                virtual void check_engine() {
                    if ( !_engine ) {
                        throw std::logic_error(
                            "Interpolator1D: no engine has been set" );
                    }
                }

                explicit operator bool() const {
                    return ( _engine ? true : false );
                }

            private:
                std::unique_ptr<_spline_1d<INDEPTYPE, DEPTYPE>> _engine;
        };

    }  // namespace interpolation
}  // namespace NCPA
