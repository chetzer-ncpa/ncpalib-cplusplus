#pragma once

#include "NCPA/linearalgebra/defines.hpp"

namespace NCPA {
    namespace linear {
        enum class vector_t { INVALID, DENSE, SPARSE };
        enum class matrix_t { INVALID, DENSE, BAND_DIAGONAL, TRIDIAGONAL };
        enum class family_t { INVALID, NCPA_DENSE, NCPA_SPARSE, NCPA_BAND_DIAGONAL };
        enum class solver_t { INVALID, BASIC, BAND_DIAGONAL, TRIDIAGONAL };

        namespace details {

            // vectors
            template<typename ELEMENTTYPE>
            class abstract_vector;
            NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( dense_vector,
                                                         abstract_vector );
            NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( sparse_vector,
                                                         abstract_vector );

            // matrices
            template<typename ELEMENTTYPE>
            class abstract_matrix;
            NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( dense_matrix,
                                                         abstract_matrix );
            
            // NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( transposed_matrix_view,
            //                                              abstract_matrix );

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

            

        }  // namespace details

        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( band_diagonal_matrix,
                                                         details::abstract_matrix );

        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Vector );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( WrapperVector, Vector );


        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Matrix );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Solver );

        template<typename ELEMENTTYPE>
        class LUDecomposition;
        // template<typename ELEMENTTYPE>
        // class BandDiagonalLUDecomposition;

        // factories
        template<typename ELEMENTTYPE>
        class VectorFactory;
        template<typename ELEMENTTYPE>
        class MatrixFactory;
        template<typename ELEMENTTYPE>
        class SolverFactory;

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