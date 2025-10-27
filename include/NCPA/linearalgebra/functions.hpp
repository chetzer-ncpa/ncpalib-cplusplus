#pragma once

#include "NCPA/linearalgebra/declarations.hpp"

#include <cmath>

namespace NCPA {
    namespace linear {
        inline size_t matrix_diagonal_size( size_t rows, size_t columns,
                                            int offset = 0 ) {
            return std::min( rows, columns ) - (size_t)std::abs( offset );
        }

        inline matrix_coordinate_t matrix_diagonal_start( int offset = 0 ) {
            matrix_coordinate_t coord;
            if (offset <= 0) {
                coord.row    = (size_t)( -offset );
                coord.column = 0;
            } else {
                coord.row    = 0;
                coord.column = (size_t)offset;
            }
            return coord;
        }

        inline matrix_coordinate_t diag2rc( int diag, size_t element ) {
            matrix_coordinate_t coord = matrix_diagonal_start( diag );
            coord.row += element;
            coord.column += element;
            return coord;
        }

        inline void rc2diag( size_t r, size_t c, int& diag, size_t& element ) {
            diag = c - r;
            element = ( diag <= 0 ? c : r );
        }

        inline size_t rc2index( size_t r, size_t c, size_t rows,
                                size_t cols ) {
            return r * cols + c;
        }
    }  // namespace linear
}  // namespace NCPA
