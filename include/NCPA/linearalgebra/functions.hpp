#pragma once

#include "NCPA/linearalgebra/declarations.hpp"

#include <cmath>

namespace NCPA {
    namespace linear {
        inline size_t matrix_diagonal_size( size_t rows, size_t columns,
                                            int diag = 0 ) {
                                                int ic = (int)columns;
                                                int ir = (int)rows;
            int diff = ic - ir;
            int maxdiag = std::min( ic, ir );
            int dsize = 0;
            if (diff <= 0) {
                if ( diag < 0 ) {
                    dsize = maxdiag - std::max( diff - diag, 0 ); 
                } else {
                    dsize = maxdiag - diag;
                }
            } else {
                if (diag <= 0) {
                    dsize = maxdiag + diag;
                } else {
                    dsize = maxdiag - std::max( diag - diff, 0 );
                }
            }
            return (dsize > 0 ? (size_t)dsize : 0);
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
            matrix_coordinate_t coord  = matrix_diagonal_start( diag );
            coord.row                 += element;
            coord.column              += element;
            return coord;
        }

        inline void rc2diag( size_t r, size_t c, int& diag, size_t& element ) {
            diag    = c - r;
            element = ( diag <= 0 ? c : r );
        }

        inline size_t rc2index( size_t r, size_t c, size_t rows,
                                size_t cols ) {
            return r * cols + c;
        }
    }  // namespace linear
}  // namespace NCPA
