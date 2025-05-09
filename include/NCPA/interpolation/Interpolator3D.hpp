#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/abstract_spline_3d.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <memory>
#include <stdexcept>

template<typename T, typename U>
static void swap( NCPA::interpolation::Interpolator3D<T, U>& a,
                  NCPA::interpolation::Interpolator3D<T, U>& b ) noexcept;

namespace NCPA {
    namespace interpolation {

        template<typename INDEPTYPE, typename DEPTYPE>
        class Interpolator3D {
            public:
                Interpolator3D() {}

                Interpolator3D( spline_engine_3d_t<INDEPTYPE,DEPTYPE> engine ) {
                    set_engine( engine );
                }

                virtual Interpolator3D& set_engine(
                    spline_engine_3d_t<INDEPTYPE,DEPTYPE> engine ) {
                    _engine = std::move( engine );
                    return *this;
                }

                virtual Interpolator3D& fill( size_t N1, size_t N2, size_t N3,
                                              const INDEPTYPE *x1,
                                              const INDEPTYPE *x2,
                                              const INDEPTYPE *x3,
                                              const DEPTYPE ***f ) {
                    check_engine();
                    _engine->fill( N1, N2, N3, x1, x2, x3, f );
                    return *this;
                }

                virtual Interpolator3D& fill(
                    const std::vector<INDEPTYPE>& x1,
                    const std::vector<INDEPTYPE>& x2,
                    const std::vector<INDEPTYPE>& x3,
                    const NCPA::arrays::vector3d_t<DEPTYPE>& f ) {
                    check_engine();
                    _engine->fill( x1, x2, x3, f );
                    return *this;
                }

                virtual Interpolator3D& init( size_t N1, size_t N2,
                                              size_t N3 ) {
                    check_engine();
                    _engine->init( N1, N2, N3 );
                    return *this;
                }

                virtual Interpolator3D& clear() {
                    _engine->clear();
                    return *this;
                }

                virtual Interpolator3D& ready() {
                    check_engine();
                    _engine->ready();
                    return *this;
                }

                virtual DEPTYPE eval_f( INDEPTYPE x1, INDEPTYPE x2,
                                        INDEPTYPE x3 ) {
                    check_engine();
                    return _engine->eval_f( x1, x2, x3 );
                }

                virtual DEPTYPE eval_df( INDEPTYPE x1, INDEPTYPE x2,
                                         INDEPTYPE x3, size_t dim1 ) {
                    check_engine();
                    return _engine->eval_df( x1, x2, x3, dim1 );
                }

                virtual DEPTYPE eval_ddf( INDEPTYPE x1, INDEPTYPE x2,
                                          INDEPTYPE x3, size_t dim1,
                                          size_t dim2 ) {
                    check_engine();
                    return _engine->eval_ddf( x1, x2, x3, dim1, dim2 );
                }

                virtual DEPTYPE eval_dddf( INDEPTYPE x1, INDEPTYPE x2,
                                           INDEPTYPE x3, size_t dim1,
                                           size_t dim2, size_t dim3 ) {
                    check_engine();
                    return _engine->eval_dddf( x1, x2, x3, dim1, dim2, dim3 );
                }

                virtual interpolator_3d_type_t interptype() const {
                    check_engine();
                    return _engine->interptype();
                }

                virtual void check_engine() const {
                    if ( !_engine ) {
                        throw std::logic_error(
                            "Interpolator3D: no engine has been set" );
                    }
                }

                explicit operator bool() const {
                    return ( _engine ? true : false );
                }

            private:
                spline_engine_3d_t<INDEPTYPE,DEPTYPE> _engine;
        };

    }  // namespace interpolation
}  // namespace NCPA

template<typename T, typename U>
static void swap( NCPA::interpolation::Interpolator3D<T, U>& a,
                  NCPA::interpolation::Interpolator3D<T, U>& b ) noexcept {
    using std::swap;
    swap( a._engine, b._engine );
}
