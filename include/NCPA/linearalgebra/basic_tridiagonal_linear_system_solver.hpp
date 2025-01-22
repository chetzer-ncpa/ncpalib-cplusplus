#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_linear_system_solver.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
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

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::basic_tridiagonal_linear_system_solver, ELEMENTTYPE );

namespace NCPA {
    namespace linear {


        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class basic_tridiagonal_linear_system_solver<
                ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public abstract_linear_system_solver<ELEMENTTYPE> {
            public:
                basic_tridiagonal_linear_system_solver() :
                    abstract_linear_system_solver<ELEMENTTYPE>() {}

                // copy constructor
                basic_tridiagonal_linear_system_solver(
                    const basic_tridiagonal_linear_system_solver<ELEMENTTYPE>&
                        other ) :
                    abstract_linear_system_solver<ELEMENTTYPE>() {
                    _mat = other._mat;
                    // _lu  = other._lu;
                }

                /**
                 * Move constructor.
                 * @param source The vector to assimilate.
                 */
                basic_tridiagonal_linear_system_solver(
                    basic_tridiagonal_linear_system_solver<ELEMENTTYPE>&&
                        source ) noexcept :
                    abstract_linear_system_solver<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                virtual ~basic_tridiagonal_linear_system_solver() {}

                friend void ::swap<ELEMENTTYPE>(
                    basic_tridiagonal_linear_system_solver<ELEMENTTYPE>& a,
                    basic_tridiagonal_linear_system_solver<ELEMENTTYPE>&
                        b ) noexcept;

                /**
                 * Assignment operator.
                 * @param other The vector to assign to this.
                 */
                basic_tridiagonal_linear_system_solver<ELEMENTTYPE>& operator=(
                    basic_tridiagonal_linear_system_solver<ELEMENTTYPE>
                        other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual abstract_linear_system_solver<ELEMENTTYPE>& clear()
                    override {
                    _mat.reset();
                    // _lu.reset();
                    return *static_cast<
                        abstract_linear_system_solver<ELEMENTTYPE> *>( this );
                }

                virtual abstract_linear_system_solver<ELEMENTTYPE>&
                    set_system_matrix(
                        const Matrix<ELEMENTTYPE>& M ) override {
                    if ( !M.is_square() ) {
                        throw std::logic_error(
                            "System matrix must be square!" );
                    }
                    if ( !M.is_tridiagonal() ) {
                        throw std::logic_error(
                            "System matrix must be tridiagonal!" );
                    }
                    _mat = std::unique_ptr<Matrix<ELEMENTTYPE>>(
                        new Matrix<ELEMENTTYPE>( M ) );
                    return *static_cast<
                        abstract_linear_system_solver<ELEMENTTYPE> *>( this );
                }

                virtual NCPA::linear::Vector<ELEMENTTYPE> solve(
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
                    ELEMENTTYPE cup, pot;
                    std::vector<ELEMENTTYPE> diag
                        = _mat->get_diagonal()->as_std(),
                        lower = _mat->get_diagonal( -1 )->as_std(),
                        upper = _mat->get_diagonal( 1 )->as_std();

                    size_t L = diag.size();
                    if ( b.size() != L ) {
                        throw std::logic_error(
                            "Matrix diagonal and RHS vector are not the "
                            "same size!" );
                    }
                    if ( b.size() != L ) {
                        throw std::logic_error(
                            "Matrix diagonal and solution vector are not "
                            "the same size!" );
                    }
                    std::vector<ELEMENTTYPE> vec( L );
                    NCPA::linear::Vector<ELEMENTTYPE> x = b;
                    x.resize( N ).zero();
                    int i;

                    if ( NCPA::math::is_zero<ELEMENTTYPE>( diag[ 0 ] ) ) {
                        throw std::range_error(
                            "First diagonal element is zero." );
                    }
                    cup      = diag[ 0 ];
                    vec[ 0 ] = cup;
                    x.set( 0, b[ 0 ] );
                    for ( i = 1; i < L; i++ ) {
                        pot = lower[ i - 1 ] / cup;
                        x.set( i, b[ i ] - ( x.get( i - 1 ) * pot ) );
                        cup = diag[ i ] - ( upper[ i - 1 ] * pot );
                        if ( NCPA::math::is_zero<ELEMENTTYPE>( cup ) ) {
                            throw std::invalid_argument(
                                "Matrix needs to be pivoted" );
                        }
                        vec[ i ] = cup;
                    }

                    x.set( L - 1, x.get( L - 1 ) / vec[ L - 1 ] );
                    for ( i = L - 2; i >= 0; i-- ) {
                        x.set( i,
                               ( x.get( i ) - ( upper[ i ] * x.get( i + 1 ) ) )
                                   / vec[ i ] );
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
                std::unique_ptr<NCPA::linear::Matrix<ELEMENTTYPE>> _mat;
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap(
    NCPA::linear::basic_tridiagonal_linear_system_solver<T>& a,
    NCPA::linear::basic_tridiagonal_linear_system_solver<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::linear::abstract_linear_system_solver<T>&>( a ),
        static_cast<NCPA::linear::abstract_linear_system_solver<T>&>( b ) );
    swap( a._mat, b._mat );
}
