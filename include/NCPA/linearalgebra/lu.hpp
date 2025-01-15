#pragma once

#include "NCPA/linearalgebra/builders.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/vector.hpp"
#include "NCPA/math.hpp"

#include <stdexcept>

#ifndef LU_DECOMPOSITION_TOLERANCE
#  define LU_DECOMPOSITION_TOLERANCE 1.0e-20
#endif

template<typename ELEMENTTYPE>
static void swap( NCPA::linear::LUDecomposition<ELEMENTTYPE>& a,
                  NCPA::linear::LUDecomposition<ELEMENTTYPE>& b ) noexcept;

template<typename ELEMENTTYPE>
static void swap(
    NCPA::linear::BandDiagonalLUDecomposition<ELEMENTTYPE>& a,
    NCPA::linear::BandDiagonalLUDecomposition<ELEMENTTYPE>& b ) noexcept;

namespace NCPA {
    namespace linear {
        template<typename ELEMENTTYPE>
        class LUDecomposition {
            public:
                LUDecomposition() :
                    _tolerance { NCPA::math::one<ELEMENTTYPE>()
                                 * LU_DECOMPOSITION_TOLERANCE } {}

                LUDecomposition( const LUDecomposition<ELEMENTTYPE>& other ) :
                    LUDecomposition<ELEMENTTYPE>() {
                    _tolerance   = other._tolerance;
                    _upper       = other._upper;
                    _lower       = other._lower;
                    _permutation = other._permutation;
                }

                LUDecomposition<ELEMENTTYPE>& operator=(
                    LUDecomposition<ELEMENTTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual LUDecomposition<ELEMENTTYPE>& clear() {
                    _lower.clear();
                    _upper.clear();
                    _permutation.clear();
                    return *this;
                }

                friend void ::swap<ELEMENTTYPE>(
                    LUDecomposition<ELEMENTTYPE>& a,
                    LUDecomposition<ELEMENTTYPE>& b ) noexcept;

                explicit operator bool() const {
                    return ( _upper && _lower && _permutation
                             && !_upper.is_empty() && !_lower.is_empty()
                             && !_permutation.is_empty() );
                }

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
                    // _upper       = base;
                    // _lower       = base;
                    // _permutation = base;
                    size_t N = base.rows();
                    _upper
                        = MatrixFactory<ELEMENTTYPE>::build( matrix_t::DENSE );
                    _lower
                        = MatrixFactory<ELEMENTTYPE>::build( matrix_t::DENSE );
                    _permutation
                        = MatrixFactory<ELEMENTTYPE>::build( matrix_t::DENSE );
                    _lower.clear().resize( N, N );
                    _upper.clear().resize( N, N ).copy( base );
                    _permutation.identity( N, N );

                    size_t i, j, k;

                    ELEMENTTYPE maxA, absA;

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
                                clear();
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

                // protected:
                virtual Matrix<ELEMENTTYPE>& lower() { return _lower; }

                virtual Matrix<ELEMENTTYPE>& upper() { return _upper; }

                virtual Matrix<ELEMENTTYPE>& permutation() {
                    return _permutation;
                }

            protected:
                ELEMENTTYPE _tolerance;
                Matrix<ELEMENTTYPE> _lower, _upper, _permutation;
        };

        template<typename ELEMENTTYPE>
        class BandDiagonalLUDecomposition
            : public LUDecomposition<ELEMENTTYPE> {
            public:
                friend class details::basic_band_diagonal_linear_system_solver<
                    ELEMENTTYPE>;

                BandDiagonalLUDecomposition() :
                    LUDecomposition<ELEMENTTYPE>() {}

                BandDiagonalLUDecomposition(
                    const BandDiagonalLUDecomposition<ELEMENTTYPE>& other ) :
                    LUDecomposition<ELEMENTTYPE>( other ) {}

                BandDiagonalLUDecomposition<ELEMENTTYPE>& operator=(
                    BandDiagonalLUDecomposition<ELEMENTTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual LUDecomposition<ELEMENTTYPE>& clear() override {
                    _A.clear();
                    return static_cast<LUDecomposition<ELEMENTTYPE>&>( *this );
                }

                virtual LUDecomposition<ELEMENTTYPE>& decompose(
                    const Matrix<ELEMENTTYPE>& Mbase,
                    bool pivot = true ) override {
                    if ( !Mbase ) {
                        throw std::logic_error(
                            "BandDiagonalLUDecomposition.compute(): base "
                            "matrix has not "
                            "been set up!" );
                    }

                    if ( !Mbase.is_band_diagonal() ) {
                        throw std::logic_error(
                            "BandDiagonalLUDecomposition.decompose(): Base "
                            "matrix must be band-diagonal" );
                    }

                    if ( !Mbase.is_square() ) {
                        throw std::logic_error(
                            "BandDiagonalLUDecomposition.decompose(): Base "
                            "matrix must be "
                            "square" );
                    }

                    _A.clear();
                    _A.copy(dynamic_cast<const band_diagonal_matrix<ELEMENTTYPE>&>(
                        Mbase.internal() ));

                    size_t n = _A.rows(), p = _A.lower_bandwidth(),
                           q = _A.upper_bandwidth();

                    for ( auto k = 1; k <= n - 1; k++ ) {
                        size_t ni = std::min( k + p, n );
                        for ( auto i = k + 1; i <= ni; i++ ) {
                            _A.set( i - 1, k - 1,
                                    _A.get( i - 1, k - 1 )
                                        / _A.get( k - 1, k - 1 ) );
                        }
                        size_t nj = std::min( k + q, n );
                        for ( auto j = k + 1; j <= nj; j++ ) {
                            ni = std::min( k + p, n );
                            for ( auto i = k + 1; i <= ni; i++ ) {
                                _A.set( i - 1, j - 1,
                                        _A.get( i - 1, j - 1 )
                                            - _A.get( i - 1, k - 1 )
                                                  * _A.get( k - 1, j - 1 ) );
                            }
                        }
                    }
                    return static_cast<LUDecomposition<ELEMENTTYPE>&>( *this );
                }

            protected:
                band_diagonal_matrix<ELEMENTTYPE> _A;
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
}

template<typename T>
static void swap( NCPA::linear::BandDiagonalLUDecomposition<T>& a,
                  NCPA::linear::BandDiagonalLUDecomposition<T>& b )
                  noexcept {
    using std::swap;
    std::swap( static_cast<NCPA::linear::LUDecomposition<T>&>( a ),
               static_cast<NCPA::linear::LUDecomposition<T>&>( b ) );
    swap( a._A, b._A );
}
