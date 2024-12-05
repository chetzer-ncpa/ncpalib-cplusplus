#pragma once

#include "NCPA/linearalgebra/defines.hpp"

namespace NCPA {
    namespace linear {
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
        template<typename ELEMENTTYPE>
        class BandDiagonalLUDecomposition;

    }  // namespace linear
}  // namespace NCPA
