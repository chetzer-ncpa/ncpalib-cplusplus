#pragma once

#include "NCPA/linearalgebra/band_diagonal_matrix.hpp"
#include "NCPA/linearalgebra/basic_linear_system_solver.hpp"
#include "NCPA/linearalgebra/basic_tridiagonal_linear_system_solver.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/dense_matrix.hpp"
#include "NCPA/linearalgebra/dense_vector.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/solver.hpp"
#include "NCPA/linearalgebra/sparse_vector.hpp"
#include "NCPA/linearalgebra/vector.hpp"

#include <stdexcept>

namespace NCPA {
    namespace linear {


        template<typename ELEMENTTYPE>
        class VectorFactory {
            public:
                static Vector<ELEMENTTYPE> build( vector_t vectype,
                                                  size_t n = 0 ) {
                    switch ( vectype ) {
                        case vector_t::DENSE:
                            return Vector<ELEMENTTYPE>(
                                std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                                    new dense_vector<ELEMENTTYPE>( n ) ) );
                            break;
                        case vector_t::SPARSE:
                            return Vector<ELEMENTTYPE>(
                                std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                                    new sparse_vector<ELEMENTTYPE>( n ) ) );
                            break;
                        default:
                            throw std::logic_error( "Unknown or unsupported "
                                                    "vector type requested" );
                    }
                }
        };

        template<typename ELEMENTTYPE>
        class MatrixFactory {
            public:
                static Matrix<ELEMENTTYPE> build( matrix_t mattype,
                                                  size_t nrows = 0,
                                                  size_t ncols = 0 ) {
                    switch ( mattype ) {
                        case matrix_t::DENSE:
                            return Matrix<ELEMENTTYPE>(
                                std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                                    new dense_matrix<ELEMENTTYPE>( nrows,
                                                                   ncols ) ) );
                            break;
                        case matrix_t::BAND_DIAGONAL:
                            return Matrix<ELEMENTTYPE>(
                                std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                                    new band_diagonal_matrix<ELEMENTTYPE>(
                                        nrows, ncols ) ) );
                            break;
                        case matrix_t::SYMMETRIC:
                            return Matrix<ELEMENTTYPE>(
                                std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                                    new symmetric_matrix<ELEMENTTYPE>(
                                        nrows, ncols ) ) );
                            break;
                        case matrix_t::SYMMETRIC_FINITE_DIFFERENCE:
                            return Matrix<ELEMENTTYPE>(
                                std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                                    new symmetric_matrix<ELEMENTTYPE>(
                                        nrows, ncols, true ) ) );
                            break;
                        default:
                            throw std::logic_error( "Unknown or unsupported "
                                                    "matrix type requested" );
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
                                std::unique_ptr<abstract_linear_system_solver<
                                    ELEMENTTYPE>>(
                                    new basic_linear_system_solver<
                                        ELEMENTTYPE>() ) );
                            break;
                        case solver_t::TRIDIAGONAL:
                            return Solver<ELEMENTTYPE>(
                                std::unique_ptr<abstract_linear_system_solver<
                                    ELEMENTTYPE>>(
                                    new basic_tridiagonal_linear_system_solver<
                                        ELEMENTTYPE>() ) );
                            break;
                        case solver_t::BAND_DIAGONAL:
                            return Solver<
                                ELEMENTTYPE>( std::unique_ptr<
                                              abstract_linear_system_solver<
                                                  ELEMENTTYPE>>(
                                new basic_band_diagonal_linear_system_solver<
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
