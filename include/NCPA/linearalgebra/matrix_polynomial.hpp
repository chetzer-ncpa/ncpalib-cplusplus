#pragma once

#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/matrix.hpp"

#include <vector>

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class MatrixPolynomial<ELEMENTTYPE,
                                   _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC> {
            public:
                MatrixPolynomial() {}

                MatrixPolynomial( const Matrix<ELEMENTTYPE>& mat ) {
                    _powers[ 0 ] = mat;
                }

                MatrixPolynomial( const Matrix<ELEMENTTYPE>& mat,
                                  size_t order ) :
                    MatrixPolynomial<ELEMENTTYPE>( mat ) {
                    this->compute( order );
                }

                MatrixPolynomial(
                    const MatrixPolynomial<ELEMENTTYPE>& other ) :
                    MatrixPolynomial<ELEMENTTYPE>() {
                    for ( auto it = other._powers.cbegin();
                          it != other._powers.cend(); ++it ) {
                        _powers.emplace_back( *it );
                    }
                }

                MatrixPolynomial(
                    MatrixPolynomial<ELEMENTTYPE>&& source ) noexcept :
                    MatrixPolynomial<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                virtual ~MatrixPolynomial() {}

                friend void ::swap<ELEMENTTYPE>(
                    MatrixPolynomial<ELEMENTTYPE>& a,
                    MatrixPolynomial<ELEMENTTYPE>& b ) noexcept;

                MatrixPolynomial<ELEMENTTYPE>& operator=(
                    MatrixPolynomial<ELEMENTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual size_t order() const { return _powers.size(); }

                virtual MatrixPolynomial<ELEMENTTYPE>& compute(
                    const Matrix<ELEMENTTYPE>& mat, size_t order ) {
                    _powers[ 0 ] = mat;
                    return this->compute( order );
                }

                virtual MatrixPolynomial<ELEMENTTYPE>& compute(
                    size_t order ) {
                    _powers.resize( 1 );
                    _powers.reserve( order );
                    if ( _powers[ 0 ].is_square()
                         && _powers[ 0 ].is_symmetrical() ) {
                        return _compute_symmetric( order );
                    } else {
                        return _compute_with_multiply( order );
                    }
                }


            protected:
                std::vector<Matrix<ELEMENTTYPE>> _powers;

                virtual MatrixPolynomial<ELEMENTTYPE>& _compute_with_multiply(
                    size_t order ) {
                    for ( size_t i = 1; i < order; i++ ) {
                        _powers[ i ] = _powers[ i - 1 ] * _powers[ 0 ];
                    }
                    return *this;
                }

                virtual MatrixPolynomial<ELEMENTTYPE>& _compute_symmetric(
                    size_t order ) {
                    int bw  = _powers[ 0 ].lower_bandwidth();
                    int dim = _powers[ 0 ].rows();

                    // all higher powers will be band-diagonal by definition
                    for ( auto i = 1; i < order; i++ ) {
                        _powers[ i ]
                            = NCPA::linear::MatrixFactory<double>::build(
                                  NCPA::linear::matrix_t::SYMMETRIC )
                                  .resize( dim, dim );
                    }

                    // the NCPA::linear::symmetric_matrix template class automatically
                    // reflects the upper off-diagonals onto the lower ones, so we can
                    // save a few write operations by checking to see if that's what we
                    // have, and not doing the reflected cases
                    bool unreflected = true;
                    if ( auto *derived = dynamic_cast<
                             const symmetric_matrix<ELEMENTTYPE> *>(
                             &(powers[0]) ) ) {
                        unreflected = false;
                    }

                    // dummy variable loops through the length of the diagonal
                    // plus extra steps.  On the first loop, computes the first
                    // column of the second power below the diagonal (and, by
                    // symmetry, the first row above the diagonal).  On the
                    // second loop, computes the second column below and second
                    // row above the diagonal for the second power, and the
                    // first row and column of the third power (which depends
                    // on the second row and column of the lower power).  Third
                    // loop adds the fourth power, and so on until all powers
                    // are in the loop.  This continues until the end, where
                    // the lower powers are sequentially completed and those
                    // loops ignored.
                    for ( int dummy = 0; dummy < dim + order - 2; dummy++ ) {
                        int min_p = std::min( std::max( 1, dummy - dim + 2 ),
                                              order - 1 );
                        int max_p = std::min( dummy + 2, order );
                        for ( int p = min_p; p < max_p; p++ ) {
                            int i = dummy - p + 1;

                            // set most distal diagonals to 1
                            Matrix<ELEMENTTYPE> *thispower = &( _powers[ p ] );
                            if ( i + p + 1 < dim ) {
                                thispower->set( i, i + p + 1, 1.0 );
                                if (unreflected) {
                                    thispower->set( i + p + 1, i, 1.0 );
                                }
                            }
                            Matrix<ELEMENTTYPE> *lastpower
                                = &( _powers[ p - 1 ] );

                            // compute main diagonal
                            int diag        = 0;
                            ELEMENTTYPE val = 0.0;
                            for ( int k = 0; k <= std::min( p, bw ); k++ ) {
                                if ( k == 0 ) {
                                    val += lastpower->get( i, i )
                                         * _powers[ 0 ].get( i, i );
                                } else {
                                    if ( i + k < dim ) {
                                        val += lastpower->get( i, i + k )
                                             * _powers[ 0 ].get( i + k, i );
                                    }
                                    if ( i - k >= 0 ) {
                                        val += lastpower->get( i, i - k )
                                             * _powers[ 0 ].get( i, i - k );
                                    }
                                }
                            }
                            thispower->set( i, i, val );

                            // compute off-diagonals that aren't 1
                            for ( int diag = 1; diag <= p; diag++ ) {
                                int row_ind = i + diag, col_ind = i;
                                // index pairs are [i+diag][i] and [i][i+diag]
                                if ( row_ind < dim ) {
                                    val = 0.0;
                                    for ( int k = std::max(
                                              i - std::min( diag, bw ), 0 );
                                          k < std::min(
                                              i + std::min( diag, bw ) + 1,
                                              dim );
                                          k++ ) {
                                        val += lastpower->get( row_ind, k )
                                             * _powers[ 0 ].get( k, i );
                                    }
                                    thispower->set( row_ind, i, val );
                                    if (unreflected) {
                                        thispower->set( i, row_ind, val );
                                    }
                                }
                            }
                        }
                    }
                    return *this;
                }
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::MatrixPolynomial<T>& a,
                  NCPA::linear::MatrixPolynomial<T>& b ) noexcept {
    using std::swap;
    swap( a._powers, b._powers );
}
