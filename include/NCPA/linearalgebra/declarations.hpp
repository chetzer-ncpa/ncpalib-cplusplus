#pragma once

#include "NCPA/linearalgebra/defines.hpp"

#include <memory>

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
        enum class solver_t {
            INVALID,
            BASIC,
            BAND_DIAGONAL,
            TRIDIAGONAL,
            BLOCK_TRIDIAGONAL
        };
        enum class matrix_polynomial_algorithm_t {
            INVALID,
            MULTIPLY,
            SYMMETRIC,
            SYMMETRIC_REFLECTED,
            FINITE_DIFFERENCE,
            FINITE_DIFFERENCE_REFLECTED
        };

        typedef struct matrix_size_t {
                size_t rows    = 0;
                size_t columns = 0;

                matrix_size_t( size_t r, size_t c ) :
                    rows { r }, columns { c } {}

                matrix_size_t() : matrix_size_t( 0, 0 ) {}
        } matrix_size_t;

        typedef struct matrix_coordinate_t {
                size_t row    = 0;
                size_t column = 0;

                matrix_coordinate_t( size_t r, size_t c ) :
                    row { r }, column { c } {}

                matrix_coordinate_t() : matrix_coordinate_t( 0, 0 ) {}
        } matrix_coordinate_t;

        typedef struct {
                matrix_coordinate_t topleft;
                matrix_coordinate_t bottomright;
        } matrix_coordinate_span_t;

        typedef struct matrix_diagonal_coordinate_t {
                int diagonal = 0;
                size_t index = 0;

                matrix_diagonal_coordinate_t( int d, size_t i ) :
                    diagonal { d }, index { i } {}

                matrix_diagonal_coordinate_t() :
                    matrix_diagonal_coordinate_t( 0, 0 ) {}
        } matrix_diagonal_coordinate_t;

        typedef struct {
                matrix_coordinate_t block_coordinates;
                matrix_coordinate_t coordinates_in_block;
        } block_matrix_coordinate_t;

        typedef struct block_matrix_indexed_coordinate_t {
                size_t index = 0;
                matrix_coordinate_t coordinates_in_block;

                block_matrix_indexed_coordinate_t( size_t i ) : index { i } {}

                block_matrix_indexed_coordinate_t() :
                    block_matrix_indexed_coordinate_t( 0 ) {}
        } block_matrix_indexed_coordinate_t;

        typedef struct {
                matrix_coordinate_t block_coordinates;
                matrix_diagonal_coordinate_t diagonal_coordinates_in_block;
        } block_matrix_diagonal_coordinate_t;

        // main public-facing classes
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Vector );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Matrix );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Solver );
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS(
            BlockMatrix );
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
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE(
            block_tridiagonal_linear_system_solver,
            abstract_linear_system_solver );

        class TridiagonalEigensolver;

        template<typename ELEMENTTYPE>
        class LUDecomposition;
        template<typename ELEMENTTYPE>
        class BandDiagonalLUDecomposition;

        template<typename ELEMENTTYPE>
        using matrix_ptr_t
            = std::unique_ptr<NCPA::linear::Matrix<ELEMENTTYPE>>;
        template<typename ELEMENTTYPE>
        using vector_ptr_t
            = std::unique_ptr<NCPA::linear::Vector<ELEMENTTYPE>>;

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
