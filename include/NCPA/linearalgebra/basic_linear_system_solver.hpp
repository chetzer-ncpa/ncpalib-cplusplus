#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_linear_system_solver.hpp"
#include "NCPA/linearalgebra/builders.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/lu.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/vector.hpp"
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

namespace NCPA {
    namespace linear {
        namespace details {
            NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE(
                basic_linear_system_solver, abstract_linear_system_solver );
        }
    }  // namespace linear
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::basic_linear_system_solver, ELEMENTTYPE );

namespace NCPA {
    namespace linear {
        namespace details {

            NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
                class basic_linear_system_solver<
                    ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
                : public abstract_linear_system_solver<ELEMENTTYPE> {
                public:
                    basic_linear_system_solver() :
                        abstract_linear_system_solver<ELEMENTTYPE>() {}

                    // copy constructor
                    basic_linear_system_solver(
                        const basic_linear_system_solver<ELEMENTTYPE>&
                            other ) :
                        abstract_linear_system_solver<ELEMENTTYPE>() {
                        _mat = other._mat;
                        _lu  = other._lu;
                    }

                    /**
                     * Move constructor.
                     * @param source The vector to assimilate.
                     */
                    basic_linear_system_solver(
                        basic_linear_system_solver<ELEMENTTYPE>&&
                            source ) noexcept :
                        abstract_linear_system_solver<ELEMENTTYPE>() {
                        ::swap( *this, source );
                    }

                    virtual ~basic_linear_system_solver() {}

                    friend void ::swap<ELEMENTTYPE>(
                        basic_linear_system_solver<ELEMENTTYPE>& a,
                        basic_linear_system_solver<ELEMENTTYPE>& b ) noexcept;

                    /**
                     * Assignment operator.
                     * @param other The vector to assign to this.
                     */
                    basic_linear_system_solver<ELEMENTTYPE>& operator=(
                        basic_linear_system_solver<ELEMENTTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual abstract_linear_system_solver<ELEMENTTYPE>& clear()
                        override {
                        _mat.reset();
                        _lu.reset();
                        return *static_cast<abstract_linear_system_solver<ELEMENTTYPE>*>( this );
                    }

                    virtual abstract_linear_system_solver<ELEMENTTYPE>&
                        set_system_matrix(
                            const Matrix<ELEMENTTYPE>& M ) override {
                        if ( !M.is_square() ) {
                            throw std::logic_error(
                                "System matrix must be square!" );
                        }
                        _mat = std::unique_ptr<Matrix<ELEMENTTYPE>>(
                            new Matrix<ELEMENTTYPE>( M ) );
                        return *static_cast<abstract_linear_system_solver<ELEMENTTYPE>*>( this );
                    }

                    virtual NCPA::linear::Matrix<ELEMENTTYPE> _solve_using_lu(
                        const NCPA::linear::Vector<ELEMENTTYPE>& b ) {
                        size_t N = b.size();
                        NCPA::linear::Matrix<ELEMENTTYPE> x
                            = NCPA::linear::MatrixFactory<ELEMENTTYPE>::build( family_t::NCPA_DENSE );
                        x.resize( N, 1 );
                        int i = 0, j = 0, iN = (int)N;

                        // temporary vectors
                        std::vector<ELEMENTTYPE> y(
                            N, NCPA::math::zero<ELEMENTTYPE>() ),
                            Pb( N, NCPA::math::zero<ELEMENTTYPE>() );

                        // forward substitution
                        for ( i = 0; i < iN; i++ ) {
                            // @todo is this loop necessary if no pivoting?
                            for ( j = 0; j < iN; j++ ) {
                                Pb[ i ] += ( _lu->permutation().get( i, j ) )
                                         * b.get( j );
                            }
                            for ( j = 0; j < i; j++ ) {
                                Pb[ i ] -= _lu->lower().get( i, j ) * y[ j ];
                            }
                            y[ i ] = Pb[ i ] / _lu->lower().get( i, i );
                        }

                        for ( i = iN - 1; i >= 0; i-- ) {
                            for ( j = i + 1; j < iN; j++ ) {
                                y[ i ]
                                    -= _lu->upper().get( i, j ) * x.get( j );
                            }
                            x.set( i, 0,  y[ i ] / _lu->upper().get( i, i ) );
                        }

                        return x;
                    }

                    virtual NCPA::linear::Matrix<ELEMENTTYPE> solve(
                        const NCPA::linear::Vector<ELEMENTTYPE>& b ) override {
                        size_t N = b.size();
                        if ( N != _mat->rows() ) {
                            std::ostringstream oss;
                            oss << "solver: size mismatch between system "
                                   "matrix "
                                   "size ["
                                << _mat->rows() << "x" << _mat->columns()
                                << "] and input vector size " << b.size();
                            throw std::logic_error( oss.str() );
                        }
                        if ( !_lu ) {
                            _build_lu();
                            _lu->init( family_t::NCPA_SPARSE );
                            _lu->decompose( *_mat, false );
                        }
                        try {
                            return _solve_using_lu( b );
                        } catch ( std::invalid_argument& e1 ) {
                            _lu->clear().init( family_t::NCPA_SPARSE );
                            _lu->decompose( *_mat, true );
                            return _solve_using_lu( b );
                        }
                    }

                    virtual NCPA::linear::Matrix<ELEMENTTYPE> solve(
                        const NCPA::linear::Matrix<ELEMENTTYPE>& b ) override {
                        if ( b.is_column_matrix() ) {
                            return solve( *b.get_column_vector( 0 ) );
                        } else if ( b.is_row_matrix() ) {
                            return solve( *b.get_row_vector( 0 ) );
                        } else {
                            throw std::logic_error(
                                "solve(): input Matrix is neither a column "
                                "matrix nor a row matrix!" );
                        }
                    }

                protected:
                    void _build_lu() {
                        _lu = std::unique_ptr<NCPA::linear::LUDecomposition<ELEMENTTYPE>>( new NCPA::linear::LUDecomposition<ELEMENTTYPE>() );
                    }

                private:
                    std::unique_ptr<NCPA::linear::Matrix<ELEMENTTYPE>> _mat;
                    std::unique_ptr<NCPA::linear::LUDecomposition<ELEMENTTYPE>> _lu;
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap(
    NCPA::linear::details::basic_linear_system_solver<T>& a,
    NCPA::linear::details::basic_linear_system_solver<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::linear::details::abstract_linear_system_solver<T>&>(
            a ),
        static_cast<NCPA::linear::details::abstract_linear_system_solver<T>&>(
            b ) );
    swap( a._mat, b._mat );
    swap( a._lu, b._lu );
}
