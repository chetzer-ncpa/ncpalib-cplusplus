#pragma once

// general definitions
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/debug.hpp"

// vector classes
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/dense_vector.hpp"
#include "NCPA/linearalgebra/sparse_vector.hpp"
#include "NCPA/linearalgebra/Vector.hpp"

// matrix classes
#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/dense_matrix.hpp"
#include "NCPA/linearalgebra/band_diagonal_matrix.hpp"
#include "NCPA/linearalgebra/symmetric_matrix.hpp"
#include "NCPA/linearalgebra/Matrix.hpp"
#include "NCPA/linearalgebra/BlockMatrix.hpp"

// matrix-adjacent classes
#include "NCPA/linearalgebra/MatrixPolynomial.hpp"
#include "NCPA/linearalgebra/BlockMatrixPolynomial.hpp"

// factories
#include "NCPA/linearalgebra/builders.hpp"

// LU decomposition
#include "NCPA/linearalgebra/lu.hpp"

// solvers
#include "NCPA/linearalgebra/abstract_linear_system_solver.hpp"
#include "NCPA/linearalgebra/basic_linear_system_solver.hpp"
#include "NCPA/linearalgebra/basic_band_diagonal_linear_system_solver.hpp"
#include "NCPA/linearalgebra/basic_tridiagonal_linear_system_solver.hpp"
#include "NCPA/linearalgebra/Solver.hpp"