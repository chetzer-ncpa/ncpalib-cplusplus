#pragma once

#include <memory>
#include "NCPA/linearalgebra/defines.hpp"

namespace NCPA {
    namespace linear {

        // main recommended API
        enum class vector_t { INVALID, DENSE, SPARSE };
        enum class matrix_t {
            INVALID,
            DENSE,
            BAND_DIAGONAL,
            BLOCK,
            SYMMETRIC,
            SYMMETRIC_FINITE_DIFFERENCE
        };
        enum class solver_t { INVALID, BASIC, BAND_DIAGONAL, TRIDIAGONAL };
        enum class matrix_polynomial_algorithm_t {
            INVALID,
            MULTIPLY,
            SYMMETRIC,
            SYMMETRIC_REFLECTED,
            FINITE_DIFFERENCE,
            FINITE_DIFFERENCE_REFLECTED
        };

        // main public-facing classes
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Vector );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Matrix );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Solver );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( BlockMatrix );
        template<typename T>
        class MatrixPolynomial;
        template<typename T>
        class BlockMatrixPolynomial;

        // factories
        template<typename ELEMENTTYPE>
        class VectorFactory;
        template<typename ELEMENTTYPE>
        class MatrixFactory;
        template<typename ELEMENTTYPE>
        class SolverFactory;

        // vectors (behind the scenes)
        template<typename ELEMENTTYPE>
        class abstract_vector;
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( dense_vector,
                                                     abstract_vector );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( sparse_vector,
                                                     abstract_vector );


        // matrices (behind the scenes)
        template<typename ELEMENTTYPE>
        class abstract_matrix;
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( dense_matrix,
                                                     abstract_matrix );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( band_diagonal_matrix,
                                                     abstract_matrix );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( symmetric_matrix,
                                                     band_diagonal_matrix );

        // solvers
        template<typename ELEMENTTYPE>
        class abstract_linear_system_solver;
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE(
            basic_linear_system_solver, abstract_linear_system_solver );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE(
            basic_band_diagonal_linear_system_solver,
            abstract_linear_system_solver );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE(
            basic_tridiagonal_linear_system_solver,
            abstract_linear_system_solver );

        template<typename ELEMENTTYPE>
        class LUDecomposition;
        template<typename ELEMENTTYPE>
        class BandDiagonalLUDecomposition;

        static size_t rc2index( size_t r, size_t c, size_t rows, size_t cols ) {
            return r * cols + c;
        }

        template<typename ELEMENTTYPE>
        using matrix_ptr_t = std::unique_ptr<NCPA::linear::Matrix<ELEMENTTYPE>>;
        template<typename ELEMENTTYPE>
        using vector_ptr_t = std::unique_ptr<NCPA::linear::Vector<ELEMENTTYPE>>;

    }  // namespace linear
}  // namespace NCPA

// friend function declarations go here
NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::Matrix,
                                             ELEMENTTYPE );

template<typename ELEMENTTYPE>
std::ostream& operator<<( std::ostream& os,
                          const NCPA::linear::Matrix<ELEMENTTYPE>& obj );

NCPA_LINEARALGEBRA_DECLARE_FRIEND_BINARY_OPERATORS( NCPA::linear::Matrix,
                                                    ELEMENTTYPE )
