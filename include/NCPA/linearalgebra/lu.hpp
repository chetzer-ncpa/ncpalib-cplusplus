#pragma once

#include "NCPA/linearalgebra/builders.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/vector.hpp"
#include "NCPA/math.hpp"

#include <stdexcept>

#ifndef LU_DECOMPOSITION_TOLERANCE
#  define LU_DECOMPOSITION_TOLERANCE 1.0e-20
#endif

namespace NCPA {
    namespace linear {

        template<typename ELEMENTTYPE>
        class LUDecomposition;

    }  // namespace linear
}  // namespace NCPA

template<typename ELEMENTTYPE>
static void swap( NCPA::linear::LUDecomposition<ELEMENTTYPE>& a,
                  NCPA::linear::LUDecomposition<ELEMENTTYPE>& b ) noexcept;

namespace NCPA {
    namespace linear {
        template<typename ELEMENTTYPE>
        class LUDecomposition {
            public:
                LUDecomposition() :
                    _tolerance { NCPA::math::one<ELEMENTTYPE>()
                                 * LU_DECOMPOSITION_TOLERANCE },
                    _family { family_t::INVALID } {}

                LUDecomposition( family_t family ) :
                    LUDecomposition(), _family { family } {
                    this->init();
                }

                LUDecomposition( const LUDecomposition<ELEMENTTYPE>& other ) :
                    LUDecomposition<ELEMENTTYPE>() {
                    _tolerance   = other._tolerance;
                    _family      = other._family;
                    _upper       = other._upper;
                    _lower       = other._lower;
                    _permutation = other._permutation;
                }

                LUDecomposition<ELEMENTTYPE>& operator=(
                    LUDecomposition<ELEMENTTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual LUDecomposition<ELEMENTTYPE>& init( family_t family ) {
                    _family = family;
                    return init();
                }

                virtual LUDecomposition<ELEMENTTYPE>& init() {
                    if ( _family != family_t::INVALID ) {
                        _lower = MatrixFactory<ELEMENTTYPE>::build( _family );
                        _upper = MatrixFactory<ELEMENTTYPE>::build( _family );
                        _permutation
                            = MatrixFactory<ELEMENTTYPE>::build( _family );
                    }
                    return *this;
                }

                virtual LUDecomposition<ELEMENTTYPE>& clear() {
                    _lower.clear();
                    _upper.clear();
                    _permutation.clear();
                    return *this;
                }

                explicit operator bool() const {
                    return ( _upper && _lower && _permutation
                             && !_upper.is_empty() && !_lower.is_empty()
                             && !_permutation.is_empty() );
                }

                // virtual LUDecomposition<ELEMENTTYPE>& set(
                //     const Matrix<ELEMENTTYPE>& mat ) {
                //     if ( !mat.is_square() ) {
                //         throw std::invalid_argument(
                //             "LUDecomposition.set(): Matrix must be square"
                //             );
                //     }
                //     _base = mat;
                //     return *this;
                // }

                virtual LUDecomposition<ELEMENTTYPE>& set_tolerance(
                    ELEMENTTYPE tol ) {
                    _tolerance = tol;
                    return *this;
                }

                virtual ELEMENTTYPE tolerance() const { return _tolerance; }

                virtual LUDecomposition<ELEMENTTYPE>& decompose(
                    const Matrix<ELEMENTTYPE>& base, bool pivot = true ) {
                    if ( !base ) {
                        throw std::logic_error(
                            "LUDecomposition.compute(): base matrix has not "
                            "been set up!" );
                    }
                    if ( !base.is_square() ) {
                        throw std::logic_error(
                            "LUDecomposition.set(): Base matrix must be "
                            "square" );
                    }
                    // if (_upper.is_empty() || _lower.is_empty() ||
                    // _permutation.is_empty()) {
                    if ( !_upper || !_lower || !_permutation ) {
                        _upper       = base;
                        _lower       = base;
                        _permutation = base;
                    }
                    size_t N = base.rows();
                    _upper.copy( base );
                    _lower.clear().resize( N, N );
                    _permutation.identity( N, N );

                    size_t i, j, k;

                    ELEMENTTYPE maxA, absA;
                    // _upper *= base;

                    for ( k = 0; k < N; k++ ) {
                        if ( pivot ) {
                            size_t pivotRow    = k;
                            ELEMENTTYPE maxVal = _upper.get( k, k );

                            // find pivot row and swap
                            for ( i = k + 1; i < N; i++ ) {
                                ELEMENTTYPE x = _upper.get( i, k );
                                if ( std::abs( x ) > std::abs( maxVal ) ) {
                                    maxVal   = x;
                                    pivotRow = i;
                                }
                            }
                            if ( std::abs( maxVal )
                                 < std::abs( tolerance() ) ) {
                                clear().init();
                                throw std::invalid_argument(
                                    "LUDecomposition.compute(): Matrix is "
                                    "degenerate" );
                            }
                            _upper.swap_rows( k, pivotRow );
                            _lower.swap_rows( k, pivotRow );
                            _permutation.swap_rows( k, pivotRow );
                        }

                        // perform elimination
                        for ( i = k + 1; i < N; i++ ) {
                            ELEMENTTYPE ratio
                                = _upper.get( i, k ) / _upper.get( k, k );
                            _lower.set( i, k, ratio );
                            for ( j = 0; j < N; j++ ) {
                                ELEMENTTYPE diff = _upper.get( i, j )
                                                 - ratio * _upper.get( k, j );
                                _upper.set( i, j, diff );
                            }
                        }
                    }
                    for ( k = 0; k < N; k++ ) {
                        _lower.set( k, k, NCPA::math::one<ELEMENTTYPE>() );
                    }
                    return *this;
                }

                virtual const Matrix<ELEMENTTYPE>& lower() const {
                    return _lower;
                }

                virtual const Matrix<ELEMENTTYPE>& upper() const {
                    return _upper;
                }

                virtual const Matrix<ELEMENTTYPE>& permutation() const {
                    return _permutation;
                }

            private:
                ELEMENTTYPE _tolerance;
                family_t _family;
                Matrix<ELEMENTTYPE> _lower, _upper, _permutation;
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::LUDecomposition<T>& a,
                  NCPA::linear::LUDecomposition<T>& b ) noexcept {
    using std::swap;
    swap( a._lower, b._lower );
    swap( a._upper, b._upper );
    swap( a._permutation, b._permutation );
    swap( a._tolerance, b._tolerance );
    swap( a._family, b._family );
}
