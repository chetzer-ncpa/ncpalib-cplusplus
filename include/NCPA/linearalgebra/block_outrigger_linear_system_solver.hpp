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
                    BlockMatrix<ELEMENTTYPE> A = *_mat;

                    // we're going to mess with the input vector too
                    Vector<ELEMENTTYPE> f = b;
                    size_t rows_per_block = A.rows_per_block();
                    size_t rows_of_blocks = A.block_rows();

                    // active overall row
                    size_t r;

                    // active overall column
                    size_t c;

                    // row to be scaled and subtracted from r
                    size_t rs;

                    // first we eliminate the lower-left corner
                    for (size_t row_in_block = 0;
                         row_in_block < rows_per_block; ++row_in_block) {
                        c = row_in_block;
                        r = rows_per_block * ( rows_of_blocks - 1 )
                          + row_in_block;

                        // repeat until we reach the lower off-diagonal
                        while (c < ( r - rows_per_block )) {
                            rs = c + rows_per_block;

                            ELEMENTTYPE factor
                                = A.get( r, c ) / A.get( rs, c );
                            A.zero( r, c );
                            A.set( r, rs,
                                   A.get( r, rs ) - factor * A.get( rs, rs ) );
                            A.set(
                                r, rs + rows_per_block,
                                A.get( r, rs + rows_per_block )
                                    - factor
                                          * A.get( rs, rs + rows_per_block ) );
                            f.set( r, f.get( r ) - factor * f.get( rs ) );
                            c += rows_per_block;
                        }
                    }

                    // explicitly remove the lower left corner
                    A.get_block( A.block_rows() - 1, 0 ).zero();

                    // now the upper right
                    for (r = 0; r < rows_per_block; ++r) {
                        c = ( rows_of_blocks - 1 ) * rows_per_block + r;
                        while (c > ( r + rows_per_block )) {
                            rs = c - rows_per_block;
                            ELEMENTTYPE factor
                                = A.get( r, c ) / A.get( rs, c );
                            A.set(
                                r, rs - rows_per_block,
                                A.get( r, rs - rows_per_block )
                                    - factor
                                          * A.get( rs, rs - rows_per_block ) );
                            A.set( r, rs,
                                   A.get( r, rs ) - factor * A.get( rs, rs ) );
                            A.zero( r, rs + rows_per_block );
                            f.set( r, f.get( r ) - factor * f.get( rs ) );
                            c -= rows_per_block;
                        }
                    }
                    A.get_block( 0, A.block_columns() - 1 ).zero();
                    Solver<ELEMENTTYPE> tdsolver
                        = SolverFactory<ELEMENTTYPE>::build(
                            solver_t::BLOCK_OUTRIGGER, _mat->block_type() );
                    tdsolver.set_system_matrix( A );

                    // internal is now block tridiagonal, solve as normal
                    return tdsolver.solve( f );
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
