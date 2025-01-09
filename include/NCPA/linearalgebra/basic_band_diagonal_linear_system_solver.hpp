#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_linear_system_solver.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/lu.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/vector.hpp"
#include "NCPA/logging.hpp"
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
    NCPA::linear::details::basic_band_diagonal_linear_system_solver,
    ELEMENTTYPE );

namespace NCPA {
    namespace linear {
        namespace details {

            NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
                class basic_band_diagonal_linear_system_solver<
                    ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
                : public abstract_linear_system_solver<ELEMENTTYPE> {
                public:
                    basic_band_diagonal_linear_system_solver() :
                        abstract_linear_system_solver<ELEMENTTYPE>() {}

                    // copy constructor
                    basic_band_diagonal_linear_system_solver(
                        const basic_band_diagonal_linear_system_solver<
                            ELEMENTTYPE>& other ) :
                        abstract_linear_system_solver<ELEMENTTYPE>() {
                        _mat = other._mat;
                    }

                    /**
                     * Move constructor.
                     * @param source The vector to assimilate.
                     */
                    basic_band_diagonal_linear_system_solver(
                        basic_band_diagonal_linear_system_solver<ELEMENTTYPE>&&
                            source ) noexcept :
                        abstract_linear_system_solver<ELEMENTTYPE>() {
                        ::swap( *this, source );
                    }

                    virtual ~basic_band_diagonal_linear_system_solver() {
                        this->clear();
                    }

                    friend void ::swap<ELEMENTTYPE>(
                        basic_band_diagonal_linear_system_solver<ELEMENTTYPE>&
                            a,
                        basic_band_diagonal_linear_system_solver<ELEMENTTYPE>&
                            b ) noexcept;

                    /**
                     * Assignment operator.
                     * @param other The vector to assign to this.
                     */
                    basic_band_diagonal_linear_system_solver<ELEMENTTYPE>&
                        operator=( basic_band_diagonal_linear_system_solver<
                                   ELEMENTTYPE>
                                       other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual abstract_linear_system_solver<ELEMENTTYPE>& clear()
                        override {
                        _mat.clear();
                        return *static_cast<
                            abstract_linear_system_solver<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_linear_system_solver<ELEMENTTYPE>&
                        set_system_matrix(
                            const Matrix<ELEMENTTYPE>& M ) override {
                        if ( !M.is_band_diagonal() ) {
                            throw std::logic_error(
                                "System matrix must be band-diagonal" );
                        }
                        if ( !M.is_square() ) {
                            throw std::logic_error(
                                "System matrix must be square!" );
                        }
                        this->clear();
                        _mat.copy( M.internal() );
                        return *static_cast<
                            abstract_linear_system_solver<ELEMENTTYPE> *>(
                            this );
                    }

                    // Implements an adaptation of the Thomas algorithm as shown in 
                    // https://github.com/aababaei/Thomas-Algorithm-Banded-Matrix/blob/main/THMS.f90
                    virtual NCPA::linear::Vector<ELEMENTTYPE> solve(
                        const NCPA::linear::Vector<ELEMENTTYPE>& b ) override {
                        size_t N = b.size();
                        if ( N != _mat.rows() ) {
                            std::ostringstream oss;
                            oss << "solver: size mismatch between system "
                                   "matrix "
                                   "size ["
                                << _mat.rows() << "x" << _mat.columns()
                                << "] and input vector size " << b.size();
                            throw std::logic_error( oss.str() );
                        }

                        NCPA::linear::Vector<ELEMENTTYPE> x = b;
                        std::vector<std::vector<ELEMENTTYPE>> T
                            = _mat._contents;
                        size_t nrows = _mat.rows();
                        size_t kl = _mat._n_lower, ku = _mat._n_upper;
                        size_t ni, nj;

                        for ( size_t k = 0; k < nrows - 1; k++ ) {
                            ni = std::min( k + kl + 1, nrows );
                            for ( size_t i = k + 1; i < ni; i++ ) {
                                ELEMENTTYPE u
                                    = T[ k + kl - i ][ i ] / T[ kl ][ k ];
                                nj = k + kl + ku - i + 1;
                                for ( size_t j = k + kl - i + 1; j < nj;
                                      j++ ) {
                                    T[ j ][ i ] -= T[ i + j - k ][ k ] * u;
                                }
                                x[ i ] -= x[ k ] * u;
                            }
                        }

                        for ( int i = nrows - 1; i >= 0; i-- ) {
                            int k    = i + 1;
                            ELEMENTTYPE S = 0.0;
                            int j    = kl + 1;
                            while ( j < kl + ku + 1 && k < nrows ) {
                                S += T[ j++ ][ i ] * x[ k++ ];
                            }
                            x[ i ] = ( x[ i ] - S ) / T[ kl ][ i ];
                        }
                        return x;
                    }

                    virtual NCPA::linear::Vector<ELEMENTTYPE> solve(
                        const NCPA::linear::Matrix<ELEMENTTYPE>& b ) override {
                        if ( b.is_column_matrix() ) {
                            return solve( *b.get_column( 0 ) );
                        } else if ( b.is_row_matrix() ) {
                            return solve( *b.get_row( 0 ) );
                        } else {
                            throw std::logic_error(
                                "solve(): input Matrix is neither a column "
                                "matrix nor a row matrix!" );
                        }
                    }

                private:
                    NCPA::linear::band_diagonal_matrix<ELEMENTTYPE> _mat;
                    const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap(
    NCPA::linear::details::basic_band_diagonal_linear_system_solver<T>& a,
    NCPA::linear::details::basic_band_diagonal_linear_system_solver<T>&
        b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::linear::details::abstract_linear_system_solver<T>&>(
            a ),
        static_cast<NCPA::linear::details::abstract_linear_system_solver<T>&>(
            b ) );
    swap( a._mat, b._mat );
    swap( a._lu, b._lu );
}
