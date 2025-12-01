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
    NCPA::linear::block_tridiagonal_linear_system_solver, ELEMENTTYPE );

namespace NCPA {
    namespace linear {


        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class block_tridiagonal_linear_system_solver<
                ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public abstract_linear_system_solver<ELEMENTTYPE> {
                friend class block_outrigger_linear_system_solver<ELEMENTTYPE>;

            public:
                block_tridiagonal_linear_system_solver() :
                    abstract_linear_system_solver<ELEMENTTYPE>() {}

                block_tridiagonal_linear_system_solver( matrix_t blocktype ) :
                    abstract_linear_system_solver<ELEMENTTYPE>() {
                    _blocktype = blocktype;
                    this->init();
                }

                // copy constructor
                block_tridiagonal_linear_system_solver(
                    const block_tridiagonal_linear_system_solver<ELEMENTTYPE>&
                        other ) :
                    abstract_linear_system_solver<ELEMENTTYPE>() {
                    _mat->copy( other._mat );
                    // _lu  = other._lu;
                }

                /**
                 * Move constructor.
                 * @param source The solver to assimilate.
                 */
                block_tridiagonal_linear_system_solver(
                    block_tridiagonal_linear_system_solver<ELEMENTTYPE>&&
                        source ) noexcept :
                    abstract_linear_system_solver<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                virtual ~block_tridiagonal_linear_system_solver() {}

                friend void ::swap<ELEMENTTYPE>(
                    block_tridiagonal_linear_system_solver<ELEMENTTYPE>& a,
                    block_tridiagonal_linear_system_solver<ELEMENTTYPE>&
                        b ) noexcept;

                /**
                 * Assignment operator.
                 * @param other The vector to assign to this.
                 */
                block_tridiagonal_linear_system_solver<ELEMENTTYPE>& operator=(
                    block_tridiagonal_linear_system_solver<ELEMENTTYPE>
                        other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual block_tridiagonal_linear_system_solver<ELEMENTTYPE>&
                    clear() override {
                    _mat.reset();
                    this->init();
                    // _lu.reset();
                    return *this;
                    // return *static_cast<
                    //     abstract_linear_system_solver<ELEMENTTYPE> *>( this
                    //     );
                }

                virtual block_tridiagonal_linear_system_solver<ELEMENTTYPE>&
                    init() {
                    if (_blocktype != matrix_t::INVALID) {
                        _mat = std::unique_ptr<BlockMatrix<ELEMENTTYPE>>(
                            new BlockMatrix<ELEMENTTYPE>( _blocktype ) );
                    }
                    return *this;
                }

                virtual block_tridiagonal_linear_system_solver<ELEMENTTYPE>&
                    set_system_matrix( const Matrix<ELEMENTTYPE>& M,
                                       bool check = true ) override {
                    if (auto BM
                        = dynamic_cast<const BlockMatrix<ELEMENTTYPE>&>( M )) {
                        if (check && !BM.is_square()) {
                            throw std::logic_error(
                                "System matrix must be square!" );
                        }
                        if (check && !BM.is_block_tridiagonal()) {
                            throw std::logic_error(
                                "System matrix must be block tridiagonal!" );
                        }

                        _mat = std::unique_ptr<BlockMatrix<ELEMENTTYPE>>(
                            new BlockMatrix<ELEMENTTYPE>( BM ) );

                        // return *static_cast<
                        //     abstract_linear_system_solver<ELEMENTTYPE> *>(
                        //     this );
                    } else {
                        throw std::logic_error(
                            "Block tridiagonal solver only works for block "
                            "matrices" );
                    }
                    return *this;
                }

                virtual NCPA::linear::Vector<ELEMENTTYPE> solve(
                    const NCPA::linear::Vector<ELEMENTTYPE>& b ) override {
                    size_t N = b.size();
                    if (N != _mat->rows()) {
                        std::ostringstream oss;
                        oss << "solver: size mismatch between system "
                               "matrix "
                               "size ["
                            << _mat->rows() << "x" << _mat->columns()
                            << "] and input vector size " << b.size();
                        throw std::logic_error( oss.str() );
                    }

                    Vector<ELEMENTTYPE> x = b;
                    x.zero();
                    size_t blocksize = _mat->rows_per_block();
                    size_t nblocks   = _mat->block_rows();

                    // hold the modified versions of b and C
                    std::vector<Vector<ELEMENTTYPE>> Riss( nblocks );
                    std::vector<Matrix<ELEMENTTYPE>> Cis( nblocks );
                    Vector<ELEMENTTYPE> *Ris;
                    Matrix<ELEMENTTYPE> Bi_inv
                        = _mat->get_block( 0, 0 ).inverse();
                    Matrix<ELEMENTTYPE> Ci = _mat->get_block(0,1);
                    Cis[ 0 ] = Bi_inv * Ci;
                    Vector<ELEMENTTYPE> Ris0
                        = Bi_inv * b.subvector( 0, blocksize );

                    for (size_t i = 1; i < nblocks; ++i) {
                        const Matrix<ELEMENTTYPE>& Ai
                            = _mat->get_block( i, i - 1 );
                        const Matrix<ELEMENTTYPE>& Bi
                            = _mat->get_block( i, i );
                        Bi_inv = ( Bi - Ai * Cis[ i - 1 ] ).inverse();
                        if (i < nblocks - 1) {
                            const Matrix<ELEMENTTYPE>& Ci
                                = _mat->get_block( i, i + 1 );
                            Cis[ i ] = Bi_inv * Ci;
                        }
                        Vector<ELEMENTTYPE> Ri
                            = b.subvector( i * blocksize, blocksize );
                        if (i == 1) {
                            Ris = &Ris0;
                        } else {
                            Ris = &( Riss.at( i - 1 ) );
                        }
                        Riss[ i ] = Bi_inv * ( Ri - ( Ai * ( *Ris ) ) );
                    }

                    x.splice( Riss.back(), ( nblocks - 1 ) * blocksize,
                              blocksize );
                    Vector<ELEMENTTYPE> xsub;
                    for (size_t i = nblocks - 2; i >= 1; --i) {
                        xsub = x.subvector( ( i + 1 ) * blocksize, blocksize );
                        Vector<ELEMENTTYPE> xseg
                            = Riss.at( i ) - ( Cis.at( i ) * xsub );
                        x.splice( xseg, i * blocksize, blocksize );
                    }

                    xsub = x.subvector( blocksize, blocksize );
                    x.splice( Ris0 - ( Cis.at( 0 ) * xsub ), 0, blocksize );
                    return x;
                }

                virtual NCPA::linear::Vector<ELEMENTTYPE> solve(
                    const NCPA::linear::Matrix<ELEMENTTYPE>& b ) override {
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

            protected:
                BlockMatrix<ELEMENTTYPE>& internal() { return *_mat; }

            private:
                matrix_t _blocktype = matrix_t::INVALID;
                std::unique_ptr<BlockMatrix<ELEMENTTYPE>> _mat;
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap(
    NCPA::linear::block_tridiagonal_linear_system_solver<T>& a,
    NCPA::linear::block_tridiagonal_linear_system_solver<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::linear::abstract_linear_system_solver<T>&>( a ),
        static_cast<NCPA::linear::abstract_linear_system_solver<T>&>( b ) );
    swap( a._mat, b._mat );
}
