#pragma once

#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/Matrix.hpp"
#include "NCPA/linearalgebra/Vector.hpp"
#include "NCPA/logging.hpp"

#include <functional>

namespace NCPA {
    namespace linear {

        // template<typename arg1_t, typename arg2_t, typename return_t>
        // using algorithm_t
        //     = std::function<return_t( const arg1_t&, const arg2_t& )>;

        // template<typename T>
        // using mvv_algorithm_t = algorithm_t<Matrix<T>, Vector<T>,
        // Vector<T>>;

        class algorithms {
            public:
                template<typename element_t>
                static Vector<element_t> multiply(
                    const Matrix<element_t>& A, const Vector<element_t>& b ) {
                    if (A.is_diagonal()) {
                        NCPA_DEBUG
                            << "Returning diagonal matrix-vector product"
                            << std::endl;
                        return _mult_diagonal_matrix_vector<element_t>( A, b );
                    } else if (A.is_tridiagonal() || A.is_band_diagonal()) {
                        NCPA_DEBUG
                            << "Returning band-diagonal matrix-vector product"
                            << std::endl;
                        return _mult_band_diagonal_matrix_vector<element_t>(
                            A, b );
                    } else {
                        NCPA_DEBUG << "Returning general matrix-vector product"
                                   << std::endl;
                        return _mult_matrix_vector<element_t>( A, b );
                    }
                }

                template<typename element_t>
                static Vector<element_t> _mult_diagonal_matrix_vector(
                    const Matrix<element_t>& A, const Vector<element_t>& b ) {
                    Vector<element_t> c( b, true );
                    c.resize( b.size() );
                    for (auto nzpair : b.nonzero()) {
                        c.set( nzpair.first,
                               A.get( nzpair.first, nzpair.first )
                                   * nzpair.second );
                    }
                    return c;
                }

                template<typename element_t>
                static Vector<element_t> _mult_band_diagonal_matrix_vector(
                    const Matrix<element_t>& A, const Vector<element_t>& b ) {
                    Vector<element_t> product( b, true );
                    product.resize( b.size() );
                    for (int r = 0; r < A.rows(); ++r) {
                        int cmin = std::max( 0, r - (int)A.lower_bandwidth() );
                        int cmax = std::min( (int)A.rows() - 1,
                                             r + (int)A.upper_bandwidth() );
                        element_t element = 0;
                        for (int c = cmin; c <= cmax; ++c) {
                            element += A.get( r, c ) * b.get( c );
                        }
                        product.set( r, element );
                    }
                    return product;
                }

                template<typename element_t>
                static Vector<element_t> _mult_matrix_vector(
                    const Matrix<element_t>& A, const Vector<element_t>& b ) {
                    Vector<element_t> product( b, true );
                    product.resize( b.size() );
                    auto ncpairs = b.nonzero();
                    for (size_t r = 0; r < A.rows(); ++r) {
                        element_t element = 0;
                        for (auto nzpair : b.nonzero()) {
                            element
                                += A.get( r, nzpair.first ) * nzpair.second;
                        }
                        product.set( r, element );
                    }
                    return product;
                }
        };
    }  // namespace linear
}  // namespace NCPA
