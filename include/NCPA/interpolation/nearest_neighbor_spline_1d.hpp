#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/types.hpp"
#include "NCPA/interpolation/abstract_spline_1d.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace NCPA {
    namespace interpolation {
        // Simple nearest neighbor interpolator
        DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( nearest_neighbor_spline_1d,
                                               details::_abstract_spline_1d );

        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
            class nearest_neighbor_spline_1d<INDEPTYPE, DEPTYPE, void,
                                             ENABLE_IF_REAL( INDEPTYPE ),
                                             ENABLE_IF_REAL( DEPTYPE )>
            : public details::_abstract_spline_1d<INDEPTYPE, DEPTYPE> {
            public:
                virtual ~nearest_neighbor_spline_1d() {}

                virtual void fill( size_t N, const INDEPTYPE *x,
                                   const DEPTYPE *f ) override {
                    if ( N != _size ) {
                        init( N );
                    }
                    std::memcpy( _x.get(), x, N * sizeof( INDEPTYPE ) );
                    std::memcpy( _f.get(), f, N * sizeof( DEPTYPE ) );
                    _size = N;
                }

                virtual void fill( const std::vector<INDEPTYPE>& x,
                                   const std::vector<DEPTYPE>& f ) override {
                    if ( x.size() != f.size() ) {
                        throw std::invalid_argument(
                            "Vector sizes must be equal!" );
                    }
                    fill( x.size(), &x[ 0 ], &f[ 0 ] );
                }

                virtual void init( size_t N ) override {
                    _x = std::unique_ptr<INDEPTYPE>(
                        NCPA::arrays::zeros<INDEPTYPE>( N ) );
                    _f = std::unique_ptr<DEPTYPE>(
                        NCPA::arrays::zeros<DEPTYPE>( N ) );
                    _size = N;
                }

                virtual void clear() override {
                    _x.reset();
                    _f.reset();
                    _size = 0;
                }

                virtual void ready() override {
                    if ( ( !_x ) || ( !_f ) ) {
                        throw std::logic_error(
                            "Interpolator has not been set up!" );
                    }
                }

                virtual DEPTYPE eval_f( INDEPTYPE x ) override {
                    return _f[ NCPA::math::find_closest_index<INDEPTYPE>(
                        _x.get(), _size, x ) ];
                }

                virtual DEPTYPE eval_df( INDEPTYPE x ) override { return 0; }

                virtual DEPTYPE eval_ddf( INDEPTYPE x ) override { return 0; }

                virtual DEPTYPE eval_dddf( INDEPTYPE x ) override { return 0; }

            private:
                std::unique_ptr<INDEPTYPE> _x;
                std::unique_ptr<DEPTYPE> _f;
                size_t _size = 0;
        };
        DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR( nearest_neighbor_spline_1d,
                                                details::_abstract_spline_1d )
    }
}

