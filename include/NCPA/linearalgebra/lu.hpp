#pragma once

#include "NCPA/linearalgebra/builders.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/Matrix.hpp"
#include "NCPA/linearalgebra/Vector.hpp"
#include "NCPA/logging.hpp"
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
                                 * (ELEMENTTYPE)LU_DECOMPOSITION_TOLERANCE } {}

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
                    const Matrix<ELEMENTTYPE>& base, bool pivot = false ) {
                    if (!base) {
                        throw std::logic_error(
                            "LUDecomposition.compute(): base matrix has not "
                            "been set up!" );
                    }
                    if (!base.is_square()) {
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
                    // NCPA_DEBUG << "Setting LU lower matrix to size " << N << " x " << N << std::endl;
                    _lower.clear().resize( N, N );
                    // NCPA_DEBUG << "Setting LU upper matrix to size " << N << " x " << N << " and copying input matrix" << std::endl;
                    _upper.clear().resize( N, N ).copy( base );
                    // NCPA_DEBUG << "Setting permutation matrix to size " << N << " x " << N << " identity matrix" << std::endl;
                    _permutation.identity( N, N );

                    size_t i, j, k;

                    ELEMENTTYPE maxA, absA;

                    for (k = 0; k < N; k++) {
                        if (pivot) {
                            size_t pivotRow    = k;
                            ELEMENTTYPE maxVal = _upper.get( k, k );

                            // find pivot row and swap
                            for (i = k + 1; i < N; i++) {
                                ELEMENTTYPE x = _upper.get( i, k );
                                if (std::abs( x ) > std::abs( maxVal )) {
                                    maxVal   = x;
                                    pivotRow = i;
                                }
                            }
                            if (std::abs( maxVal ) < std::abs( tolerance() )) {
                                clear();
                                throw std::invalid_argument(
                                    "LUDecomposition.compute(): Matrix is "
                                    "degenerate" );
                            }
                            NCPA_DEBUG << "Pivoting rows " << k << " and " << pivotRow << std::endl;
                            _upper.swap_rows( k, pivotRow );
                            _lower.swap_rows( k, pivotRow );
                            _permutation.swap_rows( k, pivotRow );
                        }

                        // perform elimination
                        for (i = k + 1; i < N; i++) {
                            ELEMENTTYPE ratio
                                = _upper.get( i, k ) / _upper.get( k, k );
                            _lower.set( i, k, ratio );
                            for (j = 0; j < N; j++) {
                                ELEMENTTYPE diff = _upper.get( i, j )
                                                 - ratio * _upper.get( k, j );
                                _upper.set( i, j, diff );
                            }
                        }
                    }
                    for (k = 0; k < N; k++) {
                        _lower.set( k, k, NCPA::math::one<ELEMENTTYPE>() );
                    }
                    return *this;
                }

                // virtual const Matrix<ELEMENTTYPE>& lower() const {
                //     return _lower;
                // }

                // virtual const Matrix<ELEMENTTYPE>& upper() const {
                //     return _upper;
                // }

                // virtual const Matrix<ELEMENTTYPE>& permutation() const {
                //     return _permutation;
                // }

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
                friend class basic_band_diagonal_linear_system_solver<
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

                virtual BandDiagonalLUDecomposition<ELEMENTTYPE>& clear() override {
                    _A.clear();
                    return *this;
                    // return static_cast<LUDecomposition<ELEMENTTYPE>&>( *this );
                }

                virtual BandDiagonalLUDecomposition<ELEMENTTYPE>& decompose(
                    const Matrix<ELEMENTTYPE>& Mbase,
                    bool pivot = true ) override {
                    if (!Mbase) {
                        throw std::logic_error(
                            "BandDiagonalLUDecomposition.compute(): base "
                            "matrix has not "
                            "been set up!" );
                    }

                    if (!Mbase.is_band_diagonal()) {
                        throw std::logic_error(
                            "BandDiagonalLUDecomposition.decompose(): Base "
                            "matrix must be band-diagonal" );
                    }

                    if (!Mbase.is_square()) {
                        throw std::logic_error(
                            "BandDiagonalLUDecomposition.decompose(): Base "
                            "matrix must be "
                            "square" );
                    }

                    _A.clear();
                    _lower.clear();
                    _upper.clear();
                    if (auto bdm_ptr = dynamic_cast<
                            const band_diagonal_matrix<ELEMENTTYPE> *>(
                            Mbase.internal() )) {
                        NCPA_DEBUG << "Copying internal band-diagonal matrix pointer" << std::endl;
                        _A.copy( *bdm_ptr );
                    } else {
                        NCPA_DEBUG << "Copying band-diagonal matrix by diagonals" << std::endl;
                        _A.resize( Mbase.rows(), Mbase.columns() );
                        std::vector<int> diags = Mbase.diagonals();
                        for (auto it = diags.begin(); it != diags.end(); ++it) {
                            NCPA_DEBUG << "Setting diagonal " << *it << std::endl;
                            _A.set_diagonal( Mbase.get_diagonal( *it )->as_std(), (size_t)(*it) );
                        }
                    }

                    size_t nrows = _A.rows();
                    NCPA_DEBUG << "Original vector has " << nrows << " rows" << std::endl;
                    size_t lbw   = _A.lower_bandwidth();
                    NCPA_DEBUG << "Original vector has lower bandwidth " << lbw << std::endl;
                    size_t ubw   = _A.upper_bandwidth();
                    NCPA_DEBUG << "Original vector has upper bandwidth " << ubw << std::endl;

                    for (auto k = 1; k <= nrows - 1; k++) {
                        NCPA_DEBUG << "Decomposing row " << k << std::endl;
                        size_t ni = std::min( k + lbw, nrows );
                        for (auto i = k + 1; i <= ni; i++) {
                            _A.set( i - 1, k - 1,
                                    _A.get( i - 1, k - 1 )
                                        / _A.get( k - 1, k - 1 ) );
                        }
                        size_t nj = std::min( k + ubw, nrows );
                        for (auto j = k + 1; j <= nj; j++) {
                            ni = std::min( k + lbw, nrows );
                            for (auto i = k + 1; i <= ni; i++) {
                                _A.set( i - 1, j - 1,
                                        _A.get( i - 1, j - 1 )
                                            - _A.get( i - 1, k - 1 )
                                                  * _A.get( k - 1, j - 1 ) );
                            }
                        }
                    }
                    NCPA_DEBUG << "LU decomposition complete" << std::endl;
                    return *this;
                    // return static_cast<LUDecomposition<ELEMENTTYPE>&>( *this );
                }

                virtual  Matrix<ELEMENTTYPE>& lower() override {
                    if (_A.is_empty()) {
                        throw std::logic_error( "Decomposition has not yet been performed" );
                    }
                    if (_lower.is_empty()) {
                        _lower = MatrixFactory<ELEMENTTYPE>::build(matrix_t::BAND_DIAGONAL);
                        _lower.identity( _A.rows(), _A.columns() );
                    }
                    std::vector<int> diags = _A.diagonals();
                    for (auto it = diags.begin(); it != diags.end() && *it < 0; ++it) {
                        _lower.set_diagonal( *_A.get_diagonal( *it ), *it );
                    }
                    return _lower;
                }

                virtual  Matrix<ELEMENTTYPE>& upper() override {
                    if (_A.is_empty()) {
                        throw std::logic_error( "Decomposition has not yet been performed" );
                    }
                    if (_upper.is_empty()) {
                        _upper = MatrixFactory<ELEMENTTYPE>::build(matrix_t::BAND_DIAGONAL);
                        _upper.resize( _A.rows(), _A.columns() );
                    }
                    std::vector<int> diags = _A.diagonals();
                    for (auto it = diags.rbegin(); it != diags.rend() && *it >= 0; ++it) {
                        _upper.set_diagonal( *_A.get_diagonal( *it ), *it );
                    }
                    return _upper;
                }

                virtual  Matrix<ELEMENTTYPE>& permutation() override {
                    if (_A.is_empty()) {
                        throw std::logic_error( "Decomposition has not yet been performed" );
                    }
                    if (_permutation.is_empty()) {
                        _permutation.identity( _A.rows(), _A.columns() );
                    }
                    return _permutation;
                }

                

            protected:
                band_diagonal_matrix<ELEMENTTYPE> _A;
                Matrix<ELEMENTTYPE> _permutation;
                Matrix<ELEMENTTYPE> _lower, _upper;


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
                  NCPA::linear::BandDiagonalLUDecomposition<T>& b ) noexcept {
    using std::swap;
    std::swap( static_cast<NCPA::linear::LUDecomposition<T>&>( a ),
               static_cast<NCPA::linear::LUDecomposition<T>&>( b ) );
    swap( a._A, b._A );
}
