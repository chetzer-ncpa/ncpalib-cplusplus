#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_linear_system_solver.hpp"
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

namespace NCPA {
    namespace linear {
        namespace details {
            NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE(
                basic_tridiagonal_linear_system_solver,
                abstract_linear_system_solver );
        }
    }  // namespace linear
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::basic_tridiagonal_linear_system_solver,
    ELEMENTTYPE );

namespace NCPA {
    namespace linear {
        namespace details {

            NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
                class basic_tridiagonal_linear_system_solver<
                    ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
                : public abstract_linear_system_solver<ELEMENTTYPE> {
                public:
                    basic_tridiagonal_linear_system_solver() :
                        abstract_linear_system_solver<ELEMENTTYPE>() {}

                    // copy constructor
                    basic_tridiagonal_linear_system_solver(
                        const basic_tridiagonal_linear_system_solver<
                            ELEMENTTYPE>& other ) :
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
                    basic_tridiagonal_linear_system_solver<ELEMENTTYPE>&
                        operator=(
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
                            abstract_linear_system_solver<ELEMENTTYPE> *>(
                            this );
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
                            abstract_linear_system_solver<ELEMENTTYPE> *>(
                            this );
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
                        ELEMENTTYPE cup, pot;
                        std::vector<ELEMENTTYPE> diag  = _mat->get_diagonal()->as_std(),
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
                        NCPA::linear::Matrix<ELEMENTTYPE> x = *_mat;
                        x.resize( N, 1 );
                        int i;

                        if ( NCPA::math::is_zero<ELEMENTTYPE>( diag[ 0 ] ) ) {
                            throw std::range_error(
                                "First diagonal element is zero." );
                        }
                        cup      = diag[ 0 ];
                        vec[ 0 ] = cup;
                        x[ 0 ][ 0 ]   = b[ 0 ];
                        for ( i = 1; i < L; i++ ) {
                            pot    = lower[ i - 1 ] / cup;
                            x[ i ][ 0 ] = b[ i ] - ( x[ i - 1 ][ 0 ] * pot );
                            cup    = diag[ i ] - ( upper[ i - 1 ] * pot );
                            if ( NCPA::math::is_zero<ELEMENTTYPE>(cup ) ) {
                                throw std::invalid_argument(
                                    "Matrix needs to be pivoted" );
                            }
                            vec[ i ] = cup;
                        }

                        x[ L - 1 ][ 0 ] = x[ L - 1 ][ 0 ] / vec[ L - 1 ];
                        for ( i = L - 2; i >= 0; i-- ) {
                            x[ i ][ 0 ] = ( x[ i ][ 0 ] - ( upper[ i ] * x[ i + 1 ][ 0 ] ) )
                                   / vec[ i ];
                        }
                        return x;
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

                // protected:
                    // void _build_lu() {
                    //     _lu = std::unique_ptr<
                    //         NCPA::linear::LUDecomposition<ELEMENTTYPE>>(
                    //         new NCPA::linear::LUDecomposition<ELEMENTTYPE>() );
                    // }

                private:
                    std::unique_ptr<NCPA::linear::Matrix<ELEMENTTYPE>> _mat;
                    // std::unique_ptr<NCPA::linear::LUDecomposition<ELEMENTTYPE>>
                    //     _lu;
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap(
    NCPA::linear::details::basic_tridiagonal_linear_system_solver<T>& a,
    NCPA::linear::details::basic_tridiagonal_linear_system_solver<T>&
        b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::linear::details::abstract_linear_system_solver<T>&>(
            a ),
        static_cast<NCPA::linear::details::abstract_linear_system_solver<T>&>(
            b ) );
    swap( a._mat, b._mat );
    // swap( a._lu, b._lu );
}
