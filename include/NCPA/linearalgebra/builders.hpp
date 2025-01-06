#pragma once

#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/dense_matrix.hpp"
#include "NCPA/linearalgebra/dense_vector.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/band_diagonal_matrix.hpp"
#include "NCPA/linearalgebra/sparse_vector.hpp"
#include "NCPA/linearalgebra/vector.hpp"
#include "NCPA/linearalgebra/basic_linear_system_solver.hpp"
#include "NCPA/linearalgebra/basic_tridiagonal_linear_system_solver.hpp"
#include "NCPA/linearalgebra/solver.hpp"

#include <stdexcept>

namespace NCPA {
    namespace linear {
        

        template<typename ELEMENTTYPE>
        class VectorFactory {
            public:
                static Vector<ELEMENTTYPE> build( vector_t vectype ) {
                    switch ( vectype ) {
                        case vector_t::DENSE:
                            return Vector<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_vector<ELEMENTTYPE>>(
                                    new details::dense_vector<
                                        ELEMENTTYPE>() ) );
                            break;
                        case vector_t::SPARSE:
                            return Vector<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_vector<ELEMENTTYPE>>(
                                    new details::sparse_vector<
                                        ELEMENTTYPE>() ) );
                            break;
                        default:
                            throw std::logic_error(
                                "Unknown or unsupported linear algebra family "
                                "requested" );
                    }
                }
        };

        template<typename ELEMENTTYPE>
        class MatrixFactory {
            public:
                static Matrix<ELEMENTTYPE> build( matrix_t mattype ) {
                    switch ( mattype ) {
                        case matrix_t::DENSE:
                            return Matrix<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_matrix<ELEMENTTYPE>>(
                                    new details::dense_matrix<
                                        ELEMENTTYPE>() ) );
                            break;
                        case matrix_t::BAND_DIAGONAL:
                            return Matrix<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_matrix<ELEMENTTYPE>>(
                                    new band_diagonal_matrix<
                                        ELEMENTTYPE>() ) );
                            break;
                        default:
                            throw std::logic_error(
                                "Unknown or unsupported linear algebra family "
                                "requested" );
                    }
                }
        };

        template<typename ELEMENTTYPE>
        class SolverFactory {
            public:
                static Solver<ELEMENTTYPE> build( solver_t family ) {
                    switch ( family ) {
                        case solver_t::BASIC:
                            return Solver<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_linear_system_solver<
                                        ELEMENTTYPE>>(
                                    new details::basic_linear_system_solver<
                                        ELEMENTTYPE>() ) );
                            break;
                        case solver_t::TRIDIAGONAL:
                            return Solver<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_linear_system_solver<
                                        ELEMENTTYPE>>(
                                    new details::
                                        basic_tridiagonal_linear_system_solver<
                                            ELEMENTTYPE>() ) );
                            break;
                        case solver_t::BAND_DIAGONAL:
                            return Solver<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_linear_system_solver<
                                        ELEMENTTYPE>>(
                                    new details::
                                        basic_band_diagonal_linear_system_solver<
                                            ELEMENTTYPE>() ) );
                            break;
                        default:
                            throw std::logic_error( "Unknown or unsupported "
                                                    "solver type requested" );
                    }
                }
        };
    }  // namespace linear
}  // namespace NCPA
