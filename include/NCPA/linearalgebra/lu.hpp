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

// template<typename ELEMENTTYPE>
// static void swap(
//     NCPA::linear::BandDiagonalLUDecomposition<ELEMENTTYPE>& a,
//     NCPA::linear::BandDiagonalLUDecomposition<ELEMENTTYPE>& b ) noexcept;

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
                    _upper       = base;
                    _lower       = base;
                    _permutation = base;
                    size_t N     = base.rows();
                    _lower.clear().resize( N, N );
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

            private:
                ELEMENTTYPE _tolerance;
                Matrix<ELEMENTTYPE> _lower, _upper, _permutation;
        };

        // template<typename ELEMENTTYPE>
        // class BandDiagonalLUDecomposition
        //     : public LUDecomposition<ELEMENTTYPE> {
        //     public:
        //         friend class details::basic_band_diagonal_linear_system_solver<ELEMENTTYPE>;

        //         BandDiagonalLUDecomposition() :
        //             LUDecomposition<ELEMENTTYPE>() {}

        //         BandDiagonalLUDecomposition(
        //             const BandDiagonalLUDecomposition<ELEMENTTYPE>& other ) :
        //             LUDecomposition<ELEMENTTYPE>( other ) {}

        //         BandDiagonalLUDecomposition<ELEMENTTYPE>& operator=(
        //             BandDiagonalLUDecomposition<ELEMENTTYPE> other ) {
        //             swap( *this, other );
        //             return *this;
        //         }

        //         virtual LUDecomposition<ELEMENTTYPE>& decompose(
        //             const Matrix<ELEMENTTYPE>& Mbase,
        //             bool pivot = true ) override {
        //             if ( !Mbase ) {
        //                 throw std::logic_error(
        //                     "BandDiagonalLUDecomposition.compute(): base "
        //                     "matrix has not "
        //                     "been set up!" );
        //             }

        //             // if ( !Mbase.is_band_diagonal() ) {
        //             //     throw std::logic_error(
        //             //         "BandDiagonalLUDecomposition.decompose(): Base "
        //             //         "matrix must be band-diagonal" );
        //             // }

        //             if ( !Mbase.is_square() ) {
        //                 throw std::logic_error(
        //                     "BandDiagonalLUDecomposition.decompose(): Base "
        //                     "matrix must be "
        //                     "square" );
        //             }

        //             Matrix<ELEMENTTYPE> base
        //                 = MatrixFactory<ELEMENTTYPE>::build(
        //                     family_t::NCPA_BAND_DIAGONAL );
        //             base.copy( Mbase );
        //             const band_diagonal_matrix<ELEMENTTYPE> *amat = dynamic_cast<
        //                       const band_diagonal_matrix<ELEMENTTYPE> *>(
        //                       base.internal() );
                    
        //             int i, j, k, l, mm;
        //             double dum;

        //             n = base.rows();
        //             _indx.resize( n );
        //             std::vector<std::vector<ELEMENTTYPE>> a = amat->_contents;
        //             _au = a;
        //             m1 = amat->_n_lower;
        //             m2 = amat->_n_upper;
        //             _al = std::vector<std::vector<ELEMENTTYPE>>>( m1 );
        //             for (i = 0; i < m1; i++) {
        //                 _al[ i ] = std::vector<ELEMENTTYPE>( n );
        //             }

        //             mm = m1 + m2;
        //             l  = m1;
        //             for ( i = 0; i < m1; i++ ) {
        //                 for ( j = m1 - i; j < mm; j++ ) {
        //                     _au[ j - 1 ][ i ]
        //                         = _au[ j ][ i ];
        //                 }
        //                 l--;
        //                 for ( j = mm - l - 1; j < mm; j++ ) {
        //                     _au[ j ][ i ] = _zero;
        //                 }
        //             }
        //             d = NCPA::math::one<double>();
        //             l = m1;
        //             for ( k = 0; k < n; k++ ) {
        //                 dum = _au[ 0 ][ k ];
        //                 i   = k;
        //                 if ( l < n ) {
        //                     l++;
        //                 }
        //                 for ( j = k + 1; j < l; j++ ) {
        //                     if ( std::abs( _au[ 0 ][ j ] )
        //                          > std::abs( dum ) ) {
        //                         dum = _au[ 0 ][ j ];
        //                         i   = j;
        //                     }
        //                 }
        //                 _indx[ k ] = i + 1;
        //                 if ( NCPA::math::is_zero( dum ) ) {
        //                     // matrix is algorithmically singular but keep
        //                     // going with tiny pivot
        //                     _au[ 0 ][ k ] = tolerance();
        //                 }
        //                 if ( i != k ) {
        //                     // interchange rows
        //                     d = -d;
        //                     for ( j = 0; j < mm; j++ ) {
        //                         std::swap( _au[ j ][ k ],
        //                                    _au[ j ][ i ] );
        //                     }
        //                 }
        //                 for ( i = k + 1; i < l; i++ ) {
        //                     dum = _au[ 0 ][ i ]
        //                         / _au[ 0 ][ k ];
        //                     _al[ i - k - 1 ][ k ] = dum;
        //                     for ( j = 1; j < mm; j++ ) {
        //                         _au[ j - 1 ][ i ]
        //                             = _au[ j ][ i ]
        //                             - dum * _au[ j ][ k ];
        //                         _au[ mm - 1 ][ i ] = _zero;
        //                     }
        //                 }
        //             }
        //             _bdlower.clear();
        //             _bdupper.clear();
        //         }

        //         virtual const Matrix<ELEMENTTYPE>& lower() const override {
        //             if (!_bdlower) {
        //                 band_diagonal_matrix<ELEMENTTYPE> templower( n, n, m1, 0 );
        //                 templower._contents = _al;
        //                 _bdlower = Matrix<ELEMENTTYPE>( templower.clone() );
        //             }
        //             return _bdlower;
        //         }

        //         virtual const Matrix<ELEMENTTYPE>& upper() const override {
        //             if (!_bdlower) {
        //                 band_diagonal_matrix<ELEMENTTYPE> tempupper( n, n, m1, m2 );
        //                 tempupper._contents = _au;
        //                 _bdupper = Matrix<ELEMENTTYPE>( tempupper.clone() );
        //             }
        //             return _bdupper;
        //         }

        //     private:
        //         // using nomenclature from Numerical Recipes for now, make
        //         // better later
        //         int n, m1, m2;
        //         double d;
        //         Vector<int> _indx;
        //         std::vector<std::vector<ELEMENTTYPE>> _au, _al;
        //         const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();
        //         Matrix<ELEMENTTYPE> _bdlower, _bdupper;
        // };
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

// template<typename T>
// static void swap( NCPA::linear::BandDiagonalLUDecomposition<T>& a,
//                   NCPA::linear::BandDiagonalLUDecomposition<T>& b ) noexcept {
//     using std::swap;
//     std::swap( static_cast<NCPA::linear::LUDecomposition<ELEMENTTYPE>&>( a ),
//                static_cast<NCPA::linear::LUDecomposition<ELEMENTTYPE>&>( b ) );
//     swap( a.n, b.n );
//     swap( a.m1, b.m1 );
//     swap( a.m2, b.m2 );
//     swap( a.d, b.d );
//     swap( a._indx, b._indx );
//     swap( a._au, b._au );
//     swap( a._al, b._al );
//     swap( a._bdlower, b._bdlower );
//     swap( a._bdupper, b._bdupper );
// }
