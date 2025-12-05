#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_linear_system_solver.hpp"
#include "NCPA/linearalgebra/BlockMatrix.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/Matrix.hpp"
#include "NCPA/linearalgebra/Vector.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <initializer_list>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::block_outrigger_linear_system_solver, ELEMENTTYPE );

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class block_outrigger_linear_system_solver<
                ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public abstract_linear_system_solver<ELEMENTTYPE> {
            public:
                block_outrigger_linear_system_solver() :
                    abstract_linear_system_solver<ELEMENTTYPE>() {}

                block_outrigger_linear_system_solver( matrix_t blocktype ) :
                    abstract_linear_system_solver<ELEMENTTYPE>() {
                    _blocktype = blocktype;
                    this->init();
                }

                // copy constructor
                block_outrigger_linear_system_solver(
                    const block_outrigger_linear_system_solver<ELEMENTTYPE>&
                        other ) :
                    abstract_linear_system_solver<ELEMENTTYPE>() {
                    _mat->copy( other._mat );
                }

                /**
                 * Move constructor.
                 * @param source The vector to assimilate.
                 */
                block_outrigger_linear_system_solver(
                    block_outrigger_linear_system_solver<ELEMENTTYPE>&&
                        source ) noexcept :
                    abstract_linear_system_solver<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                virtual ~block_outrigger_linear_system_solver() {}

                friend void ::swap<ELEMENTTYPE>(
                    block_outrigger_linear_system_solver<ELEMENTTYPE>& a,
                    block_outrigger_linear_system_solver<ELEMENTTYPE>&
                        b ) noexcept;

                /**
                 * Assignment operator.
                 * @param other The vector to assign to this.
                 */
                block_outrigger_linear_system_solver<ELEMENTTYPE>& operator=(
                    block_outrigger_linear_system_solver<ELEMENTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual block_outrigger_linear_system_solver<ELEMENTTYPE>&
                    clear() override {
                    _mat.reset();
                    this->init();
                    return *this;
                }

                virtual block_outrigger_linear_system_solver<ELEMENTTYPE>&
                    init() {
                    if (_blocktype != matrix_t::INVALID) {
                        _mat = std::unique_ptr<BlockMatrix<ELEMENTTYPE>>(
                            new BlockMatrix<ELEMENTTYPE>( _blocktype ) );
                    }
                    return *this;
                }

                virtual block_outrigger_linear_system_solver<ELEMENTTYPE>&
                    set_system_matrix( const Matrix<ELEMENTTYPE>& M,
                                       bool check = true ) override {
                    if (auto BM
                        = dynamic_cast<const BlockMatrix<ELEMENTTYPE>&>( M )) {
                        if (check) {
                            if (!BM.is_square()) {
                                throw std::invalid_argument(
                                    "System matrix must be square!" );
                            }
                            size_t n = BM.block_rows();
                            if (n != BM.block_columns()) {
                                throw std::invalid_argument(
                                    "System matrix must be blockwise square" );
                            }
                            if (BM.rows_per_block()
                                != BM.columns_per_block()) {
                                throw std::invalid_argument(
                                    "System matrix blocks must be square "
                                    "matrices" );
                            }
                            for (size_t r = 0; r < n; ++r) {
                                size_t cmin = ( r == 0 ? r : r - 1 );
                                size_t cmax = ( r == n - 1 ? r : r + 1 );
                                for (size_t c = cmin; c <= cmax; ++c) {
                                    _check_block_is_diagonal( BM, r, c );
                                }
                            }
                            _check_block_is_diagonal( BM, 0, n - 1 );
                            _check_block_is_diagonal( BM, n - 1, 0 );
                        }

                        _mat = std::unique_ptr<BlockMatrix<ELEMENTTYPE>>(
                            new BlockMatrix<ELEMENTTYPE>( BM ) );

                    } else {
                        throw std::invalid_argument(
                            "Block tridiagonal solver only works for block "
                            "matrices" );
                    }
                    return *this;
                }

                virtual Vector<ELEMENTTYPE> solve(
                    const Vector<ELEMENTTYPE>& b ) override {
                    if (b.size() != _mat->rows()) {
                        throw std::invalid_argument(
                            "Size mismatch between matrix and vector "
                            "arguments" );
                    }

                    size_t m = _mat->block_rows();
                    size_t n = _mat->get_block( 0, 0 ).rows();
                    Matrix<ELEMENTTYPE> blankcell
                        = MatrixFactory<ELEMENTTYPE>::build(
                            matrix_t::BAND_DIAGONAL );
                    blankcell.resize( n, n );
                    Matrix<ELEMENTTYPE> eye = blankcell;
                    eye.identity();
                    std::vector<Matrix<ELEMENTTYPE>> blankvec( m, blankcell );
                    BlockMatrix<ELEMENTTYPE> blankmat(
                        matrix_t::BAND_DIAGONAL );
                    blankmat.resize( m, m, n, n );

                    std::vector<Matrix<ELEMENTTYPE>> D = blankvec;
                    for (size_t i = 0; i < m; ++i) {
                        D.at( i ).set_diagonal( b.subvector( n * i, n ) );
                    }

                    Matrix<ELEMENTTYPE> gamma  = -_mat->get_block( 0, 0 );
                    Matrix<ELEMENTTYPE> igamma = gamma.inverse();

                    BlockMatrix<ELEMENTTYPE> B = *_mat;
                    B.get_block( 0, m - 1 ).zero();
                    B.get_block( m - 1, 0 ).zero();
                    B.get_block( 0, 0 )         -= gamma;
                    B.get_block( m - 1, m - 1 ) -= _mat->get_block( 0, m - 1 )
                                                 * _mat->get_block( m - 1, 0 )
                                                 * igamma;

                    std::vector<Matrix<ELEMENTTYPE>> U = blankvec;
                    U.front()                          = gamma;
                    U.back() = _mat->get_block( m - 1, 0 );
                    std::vector<Matrix<ELEMENTTYPE>> V = blankvec;
                    V.front()                          = eye;
                    V.back() = _mat->get_block( 0, m - 1 ) * igamma;

                    BlockMatrix<ELEMENTTYPE> UVT = blankmat;
                    UVT.set_block( 0, 0, U.front() * V.front() )
                        .set_block( 0, m - 1, U.front() * V.back() )
                        .set_block( m - 1, 0, U.back() * V.front() )
                        .set_block( m - 1, m - 1, U.back() * V.back() );
                    std::vector<Matrix<ELEMENTTYPE>> Y
                        = _block_vector_solve( B, D );
                    std::vector<Matrix<ELEMENTTYPE>> Q
                        = _block_vector_solve( B, U );
                    Matrix<ELEMENTTYPE> vTq = blankcell;
                    Matrix<ELEMENTTYPE> vTy = blankcell;
                    for (size_t i = 0; i < m; ++i) {
                        vTq += V[ i ] * Q[ i ];
                        vTy += V[ i ] * Y[ i ];
                    }
                    vTq                      += eye;
                    // ivTq.invert();
                    Matrix<ELEMENTTYPE> QVTY  = blankcell;
                    Vector<ELEMENTTYPE> X     = b;
                    X.zero();
                    for (size_t i = 0; i < m; ++i) {
                        QVTY.zero();
                        for (size_t j = 0; j < m; ++j) {
                            QVTY += Q[ i ] * V[ j ] * Y[ j ];
                        }
                        for (size_t k = 0; k < n; ++k) {
                            X.set( i * n + k,
                                   Y[ i ].get( k, k )
                                       - QVTY.get( k, k ) / vTq.get( k, k ) );
                        }
                    }
                    return X;
                }

                virtual Vector<ELEMENTTYPE> solve(
                    const Matrix<ELEMENTTYPE>& b ) override {
                    if (b.is_column_matrix()) {
                        return solve( *b.get_column( 0 ) );
                    } else if (b.is_row_matrix()) {
                        return solve( *b.get_row( 0 ) );
                    } else {
                        throw std::logic_error(
                            "solve(): input Matrix is neither a column "
                            "matrix nor a row matrix!" );
                    }
                }

            private:
                matrix_t _blocktype = matrix_t::INVALID;
                std::unique_ptr<BlockMatrix<ELEMENTTYPE>> _mat;

                void _check_block_is_diagonal(
                    const BlockMatrix<ELEMENTTYPE>& mat, size_t r, size_t c ) {
                    if (!mat.get_block( r, c ).is_diagonal()) {
                        std::ostringstream oss;
                        oss << "Block [ " << r << ", " << c
                            << " ] is not diagonal!";
                        throw std::invalid_argument( oss.str() );
                    }
                }

                std::vector<Matrix<ELEMENTTYPE>> _block_vector_solve(
                    const BlockMatrix<ELEMENTTYPE>& A,
                    const std::vector<Matrix<ELEMENTTYPE>>& D ) const {
                    if (!( A.is_square()
                           && A.block_rows() == A.block_columns() )) {
                        throw std::invalid_argument(
                            "Matrix must be square and block square!" );
                    }
                    if (A.block_rows() != D.size()) {
                        throw std::invalid_argument(
                            "Matrix-vector size mismatch" );
                    }
                    size_t m = D.size();
                    std::vector<Matrix<ELEMENTTYPE>> Cprime( m - 1 ),
                        Dprime( m );
                    Matrix<ELEMENTTYPE> A00inv = A.get_block( 0, 0 ).inverse();
                    Cprime[ 0 ]                = A.get_block( 0, 1 ) * A00inv;
                    Dprime[ 0 ]                = D.at( 0 ) * A00inv;
                    int i;
                    Matrix<ELEMENTTYPE> idenom
                        = ( A.get_block( 1, 1 )
                            - A.get_block( 1, 0 ) * Cprime[ 0 ] )
                              .inverse();
                    for (i = 1; i < m - 1; ++i) {
                        Cprime[ i ] = idenom * A.get_block( i, i + 1 );
                        Dprime[ i ]
                            = idenom
                            * ( D[ i ]
                                - A.get_block( i, i - 1 ) * Dprime[ i - 1 ] );
                        idenom = ( A.get_block( i + 1, i + 1 )
                                   - A.get_block( i + 1, i ) * Cprime[ i ] )
                                     .inverse();
                    }
                    i = m - 1;
                    std::vector<Matrix<ELEMENTTYPE>> X( m );
                    X[ i ] = idenom
                           * ( D[ i ]
                               - A.get_block( i, i - 1 ) * Dprime[ i - 1 ] );
                    for (i = m - 2; i >= 0; --i) {
                        X[ i ] = Dprime[ i ] - Cprime[ i ] * X[ i + 1 ];
                    }
                    return X;
                }
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap(
    NCPA::linear::block_outrigger_linear_system_solver<T>& a,
    NCPA::linear::block_outrigger_linear_system_solver<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::linear::abstract_linear_system_solver<T>&>( a ),
        static_cast<NCPA::linear::abstract_linear_system_solver<T>&>( b ) );
    swap( a._mat, b._mat );
}
