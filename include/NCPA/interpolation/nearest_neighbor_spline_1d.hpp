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
        // DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( nearest_neighbor_spline_1d,
        //                                        _spline_1d );

        _INTERPOLATOR_SPECIALIZED_TEMPLATE_DECLARATION  //
            class nearest_neighbor_spline_1d<INDEPTYPE, DEPTYPE, void,
                                             ENABLE_IF_REAL( INDEPTYPE ),
                                             ENABLE_IF_REAL( DEPTYPE )>
            : public NCPA::interpolation::_spline_1d<INDEPTYPE, DEPTYPE> {
            public:
                nearest_neighbor_spline_1d() : _size{ 0 } {}

                virtual ~nearest_neighbor_spline_1d() {}

                virtual void fill( size_t N, const INDEPTYPE *x,
                                   const DEPTYPE *f ) override {
                    if ( N != _size ) {
                        init( N );
                    }
                    _x.assign( x, x+N );
                    _f.assign( f, f+N );
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
                    clear();
                    _size = N;
                }

                virtual void clear() override {
                    _x.clear();
                    _f.clear();
                    _size = 0;
                }

                virtual void ready() override {
                    if ( _x.empty() || _f.empty() ) {
                        throw std::logic_error(
                            "Interpolator has not been set up!" );
                    }
                }

                virtual DEPTYPE eval_f( INDEPTYPE x ) override {
                    return _f[ NCPA::math::find_closest_index<INDEPTYPE>(
                        _x, x ) ];
                }

                virtual DEPTYPE eval_df( INDEPTYPE x ) override { return 0; }

                virtual DEPTYPE eval_ddf( INDEPTYPE x ) override { return 0; }

                virtual DEPTYPE eval_dddf( INDEPTYPE x ) override { return 0; }

            private:
                std::vector<INDEPTYPE> _x;
                std::vector<DEPTYPE> _f;
                size_t _size = 0;
        };
        DEFINE_COMPLEX_VERSION_OF_INTERPOLATOR( nearest_neighbor_spline_1d,
                                                _abstract_spline_1d )
    }
}

