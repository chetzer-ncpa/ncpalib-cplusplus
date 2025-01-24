#pragma once

#include "NCPA/linearalgebra/builders.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/matrix.hpp"

#include <vector>

DECLARE_SWAP_FUNCTION( NCPA::linear::MatrixPolynomial )

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class MatrixPolynomial<ELEMENTTYPE,
                                   _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC> {
            public:
                MatrixPolynomial() { _powers.resize( 1 ); }

                // MatrixPolynomial( const Matrix<ELEMENTTYPE>& mat, size_t
                // order,
                //                   matrix_polynomial_algorithm_t algo ) :
                //     MatrixPolynomial<ELEMENTTYPE>() {
                //     _powers[ 0 ] = mat;
                //     this->set_order( order );
                //     if ( algo != matrix_polynomial_algorithm_t::INVALID ) {
                //         this->set_algorithm( algo );
                //     }
                // }


                // MatrixPolynomial( const Matrix<ELEMENTTYPE>& mat,
                //                   size_t order ) :
                //     MatrixPolynomial<ELEMENTTYPE>(
                //         mat, order, matrix_polynomial_algorithm_t::INVALID )
                //         {}

                // MatrixPolynomial( const Matrix<ELEMENTTYPE>& mat ) :
                //     MatrixPolynomial(
                //         mat, 1, matrix_polynomial_algorithm_t::INVALID ) {}

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

                // API from here
                virtual MatrixPolynomial<ELEMENTTYPE>& compute(
                    const Matrix<ELEMENTTYPE>& mat, size_t order,
                    matrix_polynomial_algorithm_t algo
                    = matrix_polynomial_algorithm_t::INVALID ) {
                    _powers[ 0 ] = mat;
                    this->set_order( order );
                    return this->compute( algo );
                }

                virtual MatrixPolynomial<ELEMENTTYPE>& compute(
                    matrix_polynomial_algorithm_t algo ) {
                    this->set_algorithm( algo );
                    return this->compute();
                }

                virtual MatrixPolynomial<ELEMENTTYPE>& compute() {
                    if ( this->order() <= 1 ) {
                        return *this;
                    }
                    if ( _algorithm
                         == matrix_polynomial_algorithm_t::INVALID ) {
                        // std::cout << "Selecting best algorithm" <<
                        // std::endl;
                        _algorithm = _select_best_algorithm();
                    }
                    
                    switch ( _algorithm ) {
                        case matrix_polynomial_algorithm_t::MULTIPLY:
                            // std::cout << "MULTIPLY" << std::endl;
                            return _compute_with_multiply();
                            break;
                        case matrix_polynomial_algorithm_t::SYMMETRIC:
                            // std::cout << "SYMMETRIC" << std::endl;
                            return _compute_symmetric( false, false );
                            break;
                        case matrix_polynomial_algorithm_t::
                            SYMMETRIC_REFLECTED:
                            // std::cout << "SYMMETRIC_REFLECTED" << std::endl;
                            return _compute_symmetric( true, false );
                            break;
                        case matrix_polynomial_algorithm_t::FINITE_DIFFERENCE:
                            // std::cout << "FINITE_DIFFERENCE" << std::endl;
                            return _compute_symmetric( false, true );
                            break;
                        case matrix_polynomial_algorithm_t::
                            FINITE_DIFFERENCE_REFLECTED:
                            // std::cout << "FINITE_DIFFERENCE_REFLECTED" <<
                            // std::endl;
                            return _compute_symmetric( true, true );
                            break;
                        default:
                            throw std::logic_error(
                                "Unknown or unsupported polynomial algorithm "
                                "requested." );
                    }
                }

                virtual size_t order() const { return _powers.size(); }

                virtual MatrixPolynomial<ELEMENTTYPE>& scale(
                    const std::vector<ELEMENTTYPE>& factors ) {
                    if ( factors.size() != _powers.size() ) {
                        throw std::range_error(
                            "MatrixPolynomial.scale(): factors vector and "
                            "polynomial vector not the same size." );
                    }
                    for ( size_t i = 0; i < _powers.size(); i++ ) {
                        _powers[ i ].scale( factors[ i ] );
                    }
                    return *this;
                }

                virtual MatrixPolynomial<ELEMENTTYPE>& scale(
                    const ELEMENTTYPE& factor ) {
                    for ( size_t i = 0; i < _powers.size(); i++ ) {
                        _powers[ i ].scale( factor );
                    }
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE> scaled_sum(
                    const std::vector<ELEMENTTYPE>& factors ) const {
                    if ( factors.size() != _powers.size() ) {
                        throw std::range_error(
                            "MatrixPolynomial.scaled_sum(): factors vector "
                            "and "
                            "polynomial vector not the same size." );
                    }
                    Matrix<ELEMENTTYPE> sum = _powers[ 0 ] * factors[ 0 ];
                    for ( size_t i = 1; i < _powers.size(); i++ ) {
                        sum += _powers[ i ] * factors[ i ];
                    }
                    return sum;
                }

                virtual Matrix<ELEMENTTYPE> scaled_sum(
                    const ELEMENTTYPE& factor ) const {
                    Matrix<ELEMENTTYPE> sum = _powers[ 0 ] * factor;
                    for ( size_t i = 1; i < _powers.size(); i++ ) {
                        sum += _powers[ i ] * factor;
                    }
                    return sum;
                }

                virtual MatrixPolynomial<ELEMENTTYPE>& set_algorithm(
                    matrix_polynomial_algorithm_t algo
                    = matrix_polynomial_algorithm_t::INVALID ) {
                    if ( algo == matrix_polynomial_algorithm_t::INVALID ) {
                        _algorithm = _select_best_algorithm();
                    } else {
                        _algorithm = algo;
                    }
                    return *this;
                }

                virtual MatrixPolynomial<ELEMENTTYPE>& set_base(
                    const Matrix<ELEMENTTYPE>& mat ) {
                    if ( _powers.size() == 0 ) {
                        _powers.assign( 1, mat );
                    } else {
                        _powers[ 0 ] = mat;
                    }
                    return *this;
                }

                virtual MatrixPolynomial<ELEMENTTYPE>& set_order(
                    size_t order ) {
                    _powers.resize( 1 );
                    _powers.reserve( order );
                    for ( size_t i = 1; i < order; i++ ) {
                        _powers.push_back( _powers[ 0 ] );
                        _powers.back().zero();
                    }
                    return *this;
                }

                virtual std::vector<Matrix<ELEMENTTYPE>>& vector() {
                    return _powers;
                }

                virtual const std::vector<Matrix<ELEMENTTYPE>>& vector()
                    const {
                    return _powers;
                }

                // virtual MatrixPolynomial<ELEMENTTYPE>& compute(
                //     bool use_symmetric ) {
                //     if ( this->order() <= 1 ) {
                //         return *this;
                //     }
                //     if ( use_symmetric ) {
                //         return _compute_symmetric();
                //     } else {
                //         return _compute_with_multiply();
                //     }
                // }


                // virtual MatrixPolynomial<ELEMENTTYPE>& reset( size_t
                // neworder
                //                                               = 0 ) {
                //     size_t order = _powers.size();
                //     if ( neworder > 0 && neworder != order ) {
                //         Matrix<ELEMENTTYPE> base = _powers[ 0 ];
                //         _powers.assign( order, base );
                //     }
                //     for ( size_t i = 1; i < order; i++ ) {
                //         _powers[ i ].zero();
                //     }
                //     return *this;
                // }

            protected:
                std::vector<Matrix<ELEMENTTYPE>> _powers;
                matrix_polynomial_algorithm_t _algorithm
                    = matrix_polynomial_algorithm_t::INVALID;

                matrix_polynomial_algorithm_t _select_best_algorithm() {
                    return matrix_polynomial_algorithm_t::MULTIPLY;
                }


                //     if ( !_powers[ 0 ].is_square()
                //          || !_powers[ 0 ].is_symmetric() ) {
                //         return matrix_polynomial_algorithm_t::MULTIPLY;
                //     }
                //     bool reflected = false, finite_difference = false;

                //     // right now only the symmetric_matrix class is guaranteed
                //     // reflected
                //     if ( auto *derived
                //          = dynamic_cast<const symmetric_matrix<ELEMENTTYPE> *>(
                //              &( _powers[ 0 ] ) ) ) {
                //         reflected = true;
                //     }

                //     // to be a finite difference matrix, it must be
                //     // tridiagonal, and the off-diagonals must all be 1.
                //     if ( _powers[ 0 ].is_tridiagonal()
                //          && _powers[ 0 ].bandwidth() == 3 ) {
                //         auto super      = _powers[ 0 ].get_diagonal( 1 );
                //         ELEMENTTYPE one = NCPA::math::one<ELEMENTTYPE>();
                //         size_t counter = 0, limit = super->size() - 1;

                //         // short-circuit if we hit a non-one value
                //         while ( !finite_difference && counter < limit ) {
                //             if ( super->get( counter++ ) != one ) {
                //                 finite_difference = true;
                //             }
                //         }

                //         // don't bother if we already know
                //         if ( !finite_difference ) {
                //             auto sub = _powers[ 0 ].get_diagonal( -1 );
                //             counter  = 0;
                //             // short-circuit if we hit a non-one value
                //             while ( !finite_difference && counter < limit ) {
                //                 if ( sub->get( counter++ ) != one ) {
                //                     finite_difference = true;
                //                 }
                //             }
                //         }
                //     }

                //     if ( reflected ) {
                //         if ( finite_difference ) {
                //             return matrix_polynomial_algorithm_t::
                //                 FINITE_DIFFERENCE_REFLECTED;
                //         } else {
                //             return matrix_polynomial_algorithm_t::
                //                 SYMMETRIC_REFLECTED;
                //         }
                //     } else {
                //         if ( finite_difference ) {
                //             return matrix_polynomial_algorithm_t::
                //                 FINITE_DIFFERENCE;
                //         } else {
                //             return matrix_polynomial_algorithm_t::SYMMETRIC;
                //         }
                //     }
                // }

                virtual MatrixPolynomial<ELEMENTTYPE>&
                    _compute_with_multiply() {
                    for ( size_t i = 1; i < this->order(); i++ ) {
                        _powers[ i ] = _powers[ i - 1 ] * _powers[ 0 ];
                    }
                    return *this;
                }

                /**
                 * Compute the matrix polynomial assuming symmetry, saving
                 * some multiplications at the cost of a more complex
                 * algorithm.  Jury is still out on whether we can make this
                 * faster than just multiplying matrices recursively.
                 *
                 * The boolean flags indicate whether certain
                 * shortcuts can be taken.  reflected indicates that for each
                 * symmetric pair [i,j] and [j,i], we only have to set one of
                 * them and the other will be set automatically (or in this
                 * case, indices are flipped internally when requesting
                 * subdiagonal elements).  This saves a little bit of time. The
                 * finite_difference flag indicates that the outermost
                 * diagonals are 1, so when computing a power of the matrix,
                 * the outermost diagonals of that matrix will also be 1
                 * (although they won't be the same diagonals).  This saves a
                 * surprising amount of time, actually.
                 */
                virtual MatrixPolynomial<ELEMENTTYPE>& _compute_symmetric(
                    bool reflected, bool finite_difference ) {
                    int bw           = _powers[ 0 ].upper_bandwidth();
                    int dim          = _powers[ 0 ].rows();
                    int order        = (int)this->order();
                    ELEMENTTYPE one  = NCPA::math::one<ELEMENTTYPE>();
                    // the NCPA::linear::symmetric_matrix template class
                    // automatically reflects the upper off-diagonals onto the
                    // lower ones, so we can save a few write operations by
                    // checking to see if that's what we have, and not doing
                    // the reflected cases
                    bool unreflected = !reflected;

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
                        // On the first few and last few loops, we ignore
                        // progressively higher and lower-order terms,
                        // respectively
                        int min_p = std::min(
                            std::max( 1, dummy - dim + 2 ),
                            (int)order - 1 );  // lowest power to compute
                        int max_p = std::min(
                            dummy + 2,
                            (int)order );  // highest power to compute + 1

                        // loop over the powers to compute
                        for ( int p = min_p; p < max_p; p++ ) {
                            // the row/column we're working on for this power
                            // on this loop
                            int i = dummy - p + 1;

                            // get a pointer to the current and previous power,
                            // for convenience
                            Matrix<ELEMENTTYPE> *thispower = &( _powers[ p ] ),
                                                *lastpower
                                                = &( _powers[ p - 1 ] );

                            // If the matrix is a finite-difference operator,
                            // we can save a loop and set the most distal
                            // diagonals to 1
                            if ( finite_difference && i + p + 1 < dim ) {
                                thispower->set( i, i + p + 1, one );

                                // If we're not working with a matrix type
                                // that reflects internally, we have to set
                                // the other value to maintain symmetry
                                if ( unreflected ) {
                                    thispower->set( i + p + 1, i, one );
                                }
                            }

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

                            // compute off-diagonals
                            // If we're working with a finite-difference
                            // matrix, we don't have to do the last diagonal,
                            // cause we set it to one above
                            int diaglimit = ( finite_difference ? p : p + 1 );
                            for ( int diag = 1; diag <= diaglimit; diag++ ) {
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

                                    // If we're not working with a matrix type
                                    // that reflects internally, we have to set
                                    // the other value to maintain symmetry
                                    if ( unreflected ) {
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
