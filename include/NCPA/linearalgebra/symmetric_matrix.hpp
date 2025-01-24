#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/sparse_vector.hpp"
#include "NCPA/logging.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <initializer_list>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::symmetric_matrix,
                                             ELEMENTTYPE );

namespace NCPA {
    namespace linear {


        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class symmetric_matrix<ELEMENTTYPE,
                                   _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public band_diagonal_matrix<ELEMENTTYPE> {
            public:
                using band_diagonal_matrix<ELEMENTTYPE>::set;
                using band_diagonal_matrix<ELEMENTTYPE>::set_row;
                using band_diagonal_matrix<ELEMENTTYPE>::set_column;
                using band_diagonal_matrix<ELEMENTTYPE>::set_diagonal;
                using band_diagonal_matrix<ELEMENTTYPE>::contents;
                friend class BandDiagonalLUDecomposition<ELEMENTTYPE>;
                friend class basic_band_diagonal_linear_system_solver<
                    ELEMENTTYPE>;

                symmetric_matrix( size_t nrows, size_t ncols,
                                  size_t n_offdiags ) {
                    if ( nrows != ncols ) {
                        throw std::invalid_argument(
                            "Symmetric matrix must be square!" );
                    }
                    _nrows     = nrows;
                    _ncols     = ncols;
                    _n_offdiag = n_offdiags;
                    _prep_diagonals();
                }

                symmetric_matrix( size_t nrows, size_t ncols ) :
                    symmetric_matrix<ELEMENTTYPE>( nrows, ncols, 0 ) {}

                symmetric_matrix() :
                    symmetric_matrix<ELEMENTTYPE>( 0, 0, 0 ) {}

                symmetric_matrix(
                    const symmetric_matrix<ELEMENTTYPE>& other ) :
                    symmetric_matrix<ELEMENTTYPE>( other.rows(),
                                                   other.columns() ) {
                    _n_offdiag = other._n_offdiag;
                    contents()  = other.contents();
                }

                symmetric_matrix( const abstract_matrix<ELEMENTTYPE>& other ) :
                    symmetric_matrix<ELEMENTTYPE>() {
                    this->copy( other );
                }

                symmetric_matrix(
                    symmetric_matrix<ELEMENTTYPE>&& source ) noexcept :
                    symmetric_matrix<ELEMENTTYPE>() {
                    std::cout << "Moving" << std::endl;
                    ::swap( *this, source );
                }

                // constructors, destructors, copying, and assignment
                virtual ~symmetric_matrix() {}

                friend void ::swap<ELEMENTTYPE>(
                    symmetric_matrix<ELEMENTTYPE>& a,
                    symmetric_matrix<ELEMENTTYPE>& b ) noexcept;

                symmetric_matrix<ELEMENTTYPE>& operator=(
                    symmetric_matrix<ELEMENTTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> clone()
                    const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new symmetric_matrix<ELEMENTTYPE>( *this ) );
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                    fresh_clone() const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new symmetric_matrix<ELEMENTTYPE>() );
                }

                virtual abstract_matrix<ELEMENTTYPE>& copy(
                    const abstract_matrix<ELEMENTTYPE>& other ) override {
                    resize( other.rows(), other.columns() );
                    this->set_diagonal( *( other.get_diagonal() ) );
                    int countdown = other.max_off_diagonal();
                    int ndiag     = 1;
                    auto diag     = other.get_diagonal( ndiag );
                    while ( ndiag <= countdown ) {
                        if ( !diag->is_zero() ) {
                            set_diagonal( *diag, ndiag );
                        }
                        ndiag++;
                        diag = other.get_diagonal( ndiag );
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                const symmetric_matrix<ELEMENTTYPE> *downcast(
                    const abstract_matrix<ELEMENTTYPE> *in ) const {
                    return dynamic_cast<const symmetric_matrix<ELEMENTTYPE> *>(
                        in );
                }

                // virtual abstract_matrix<ELEMENTTYPE>& upcast()
                //     override {
                //     return *static_cast<abstract_matrix<ELEMENTTYPE> *>(
                //         this );
                // }

                virtual std::string id() const override {
                    return "NCPA symmetric band-diagonal matrix";
                }

                virtual size_t rows() const override { return _nrows; }

                virtual size_t columns() const override { return _ncols; }

                virtual abstract_matrix<ELEMENTTYPE>& clear() override {
                    _nrows     = 0;
                    _ncols     = 0;
                    _n_offdiag = 0;
                    contents().clear();
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const override {
                    // this->check_size( row, col );
                    int ind1, ind2;
                    if ( rowcol2internal( row, col, ind1, ind2 ) ) {
                        return contents()[ ind1 ][ ind2 ];
                    } else {
                        return _zero;
                    }
                }

                virtual bool rowcol2internal( const size_t& row,
                                              const size_t& col, int& ind1,
                                              int& ind2 ) const override {
                    if ( col < row ) {
                        return rowcol2internal( col, row, ind1, ind2 );
                    }
                    ind1 = (int)col - (int)row;
                    ind2 = (int)row;
                    // std::cout << "[" << row << "," << col << "] -> [" <<
                    // ind1 << "," << ind2 << "]" << std::endl;
                    return ( ind1 >= 0 && ind1 < (int)( contents().size() )
                             && row < this->rows() && col < this->columns() );
                }

                virtual bool internal2rowcol( const size_t& ind1,
                                              const size_t& ind2, int& row,
                                              int& col ) const override {
                    col = (int)ind1 + (int)ind2;
                    row = (int)ind2;
                    return ( col >= 0 && col < this->columns()
                             && row < this->rows() );
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_row(
                    size_t row, size_t nvals, const size_t *column_inds,
                    const ELEMENTTYPE *vals ) override {
                    for ( size_t i = 0; i < nvals; i++ ) {
                        if ( row <= column_inds[ i ] ) {
                            this->set_safe( row, column_inds[ i ], vals[ i ] );
                        } else {
                            throw std::range_error(
                                "symmetric_matrix.set_row(): cannot directly "
                                "set subdiagonals of symmetric_matrix" );
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    throw std::logic_error(
                        "symmetric_matrix.set_row(): connot set row without "
                        "indices, because subdiagonals of symmetric matrix "
                        "cannot be set." );
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_row(
                    size_t row, ELEMENTTYPE val ) override {
                    for (size_t col = row; col < std::min( this->columns(), row + _n_offdiag + 1); col++) {
                        this->set( row, col, val );
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    throw std::logic_error(
                        "symmetric_matrix.set_column(): connot set column "
                        "without indices, because subdiagonals of symmetric "
                        "matrix cannot be set." );
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_column(
                    size_t column, size_t nvals, const size_t *row_inds,
                    const ELEMENTTYPE *vals ) override {
                    for ( size_t i = 0; i < nvals; i++ ) {
                        // std::cout << "Trying to set [" << row_inds[i] << "," << column << "]" << std::endl;
                        if ( column >= row_inds[ i ] ) {
                            this->set_safe( row_inds[ i ], column, vals[ i ] );
                        } else {
                            throw std::range_error(
                                "symmetric_matrix.set_column(): cannot "
                                "directly set subdiagonals of "
                                "symmetric_matrix" );
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_column(
                    size_t col, ELEMENTTYPE val ) override {
                    int icol = (int)col;
                    for (int row = icol; row >= std::max( icol - (int)_n_offdiag, 0 ); row--) {
                        this->set( row, col, val );
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& as_array(
                    size_t& nrows, size_t& ncols,
                    ELEMENTTYPE **& vals ) override {
                    if ( vals == nullptr ) {
                        nrows = rows();
                        ncols = columns();
                        vals
                            = NCPA::arrays::zeros<ELEMENTTYPE>( nrows, ncols );
                    } else {
                        if ( nrows == 0 ) {
                            nrows = rows();
                        } else if ( nrows != rows() ) {
                            throw std::invalid_argument(
                                "Wrong number of rows requested" );
                        }
                        if ( ncols == 0 ) {
                            ncols = columns();
                        } else if ( ncols != columns() ) {
                            throw std::invalid_argument(
                                "Wrong number of columns requested" );
                        }
                        NCPA::arrays::fill( vals, nrows, ncols, _zero );
                    }
                    int row, col;
                    for ( size_t ind1 = 0; ind1 < contents().size(); ind1++ ) {
                        for ( size_t ind2 = 0; ind2 < _max_ind2( ind1 );
                              ind2++ ) {
                            if ( internal2rowcol( ind1, ind2, row, col ) ) {
                                vals[ row ][ col ] = contents()[ ind1 ][ ind2 ];
                                if ( row != col ) {
                                    vals[ col ][ row ]
                                        = contents()[ ind1 ][ ind2 ];
                                }
                            }
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                // @todo make read-write vector view for columns
                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>> get_row(
                    size_t row ) const override {
                    // std::cout << "row = " << row << std::endl;
                    std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                        new sparse_vector<ELEMENTTYPE>( this->columns() ) );
                    int ind1, ind2;
                    v->set( row, contents()[ 0 ][ row ] );
                    for ( int offd = 1; offd <= _n_offdiag; offd++ ) {
                        if ( rowcol2internal( row, row + offd, ind1, ind2 ) ) {
                            v->set( row + offd, contents()[ ind1 ][ ind2 ] );
                        }
                        if ( rowcol2internal( row, row - offd, ind1, ind2 ) ) {
                            v->set( row - offd, contents()[ ind1 ][ ind2 ] );
                        }
                    }
                    return v;
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    get_column( size_t col ) const override {
                    // std::cout << "row = " << row << std::endl;
                    std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                        new sparse_vector<ELEMENTTYPE>( this->rows() ) );
                    int ind1, ind2;
                    v->set( col, contents()[ 0 ][ col ] );
                    for ( int offd = 1; offd <= _n_offdiag; offd++ ) {
                        if ( rowcol2internal( col + offd, col, ind1, ind2 ) ) {
                            v->set( col + offd, contents()[ ind1 ][ ind2 ] );
                        }
                        if ( rowcol2internal( col - offd, col, ind1, ind2 ) ) {
                            v->set( col - offd, contents()[ ind1 ][ ind2 ] );
                        }
                    }
                    return v;
                }

                virtual abstract_matrix<ELEMENTTYPE>& resize(
                    size_t r, size_t c ) override {
                    if ( r != c ) {
                        throw std::invalid_argument(
                            "Symmetric matrix must be square!" );
                    }
                    if ( contents().size() == 0 ) {
                        // starting from scratch
                        _nrows     = r;
                        _ncols     = c;
                        _n_offdiag = 0;
                        _prep_diagonals();
                    } else {
                        // did the diagonal size increase or decrease?
                        int ddiag = (int)std::min( r, c )
                                  - (int)std::min( rows(), columns() );
                        if ( ddiag > 0 ) {
                            // increase each diagonal by ddiag elements
                            for ( auto it = contents().begin();
                                  it != contents().end(); ++it ) {
                                it->insert( it->cend(), (size_t)ddiag, _zero );
                            }
                        } else if ( ddiag < 0 ) {
                            for ( auto it = contents().begin();
                                  it != contents().end(); ++it ) {
                                it->erase( it->cend() - size_t( -ddiag ),
                                           it->cend() );
                            }
                        }
                        _nrows = r;
                        _ncols = c;
                    }

                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set(
                    size_t row, size_t col, ELEMENTTYPE val ) override {
                    this->check_size( row, col );
                    int ind1, ind2;
                    while ( !rowcol2internal( row, col, ind1, ind2 ) ) {
                        this->_add_offdiagonal();
                    }
                    contents().at( ind1 ).at( ind2 ) = val;
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_safe(
                    size_t row, size_t col, ELEMENTTYPE val ) override {
                    this->check_size( row, col );
                    int ind1, ind2;
                    if ( rowcol2internal( row, col, ind1, ind2 ) ) {
                        contents().at( ind1 ).at( ind2 ) = val;
                    } else {
                        std::ostringstream oss;
                        oss << "Element [" << row << ", " << col
                            << "] is out of band.";
                        throw std::out_of_range( oss.str() );
                    }

                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set(
                    ELEMENTTYPE val ) override {
                    std::vector<ELEMENTTYPE> diag( this->rows(), val );
                    for ( auto i = 0; i <= _n_offdiag; i++ ) {
                        contents()[ i ].assign( diag.begin(), diag.end() );
                        diag[ this->rows() - i - 1 ] = _zero;
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& transpose() override {
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                // virtual abstract_matrix<ELEMENTTYPE>& swap_rows(
                //     size_t ind1, size_t ind2 ) override {
                //     throw std::logic_error(
                //         "Cannot swap rows of a band-diagonal matrix!" );
                // }

                // virtual abstract_matrix<ELEMENTTYPE>& swap_columns(
                //     size_t ind1, size_t ind2 ) override {
                //     throw std::logic_error(
                //         "Cannot swap columns of a band-diagonal matrix!" );
                // }

                virtual abstract_matrix<ELEMENTTYPE>& zero(
                    size_t row, size_t col ) override {
                    return set( row, col, _zero );
                }

                virtual abstract_matrix<ELEMENTTYPE>& zero() override {
                    _n_offdiag = 0;
                    contents().clear();
                    contents().push_back(
                        std::vector<ELEMENTTYPE>( this->diagonal_size() ) );
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                // virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                //     build_vector( size_t n = 0 ) const override {
                //     return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                //         new sparse_vector<ELEMENTTYPE>( n ) );
                // }

                // overrides/specializations of non-pure virtual methods
                virtual bool equals( const abstract_matrix<ELEMENTTYPE>&
                                         other ) const override {
                    if ( !is_this_subclass( other ) ) {
                        return band_diagonal_matrix<ELEMENTTYPE>::equals(
                            other );
                    }
                    const symmetric_matrix<ELEMENTTYPE> *b
                        = downcast( &other );

                    return ( rows() == b->rows() && columns() == b->columns()
                             && contents() == b->contents() );
                }

                virtual bool is_identity() const override {
                    if ( !this->is_square() || !this->is_diagonal() ) {
                        return false;
                    }
                    for ( size_t i = 0; i < contents()[ 0 ].size(); i++ ) {
                        if ( contents()[ 0 ][ i ]
                             != NCPA::math::one<ELEMENTTYPE>() ) {
                            return false;
                        }
                    }
                    return true;
                }

                virtual bool is_diagonal() const override {
                    if ( this->is_empty() ) {
                        return true;
                    }
                    return _n_offdiag == 0;
                }

                virtual bool is_tridiagonal() const override {
                    if ( this->is_empty() ) {
                        return true;
                    }
                    return _n_offdiag <= 1;
                }

                virtual bool is_upper_triangular() const override {
                    return this->is_diagonal();
                }

                virtual bool is_lower_triangular() const override {
                    return this->is_diagonal();
                }

                virtual bool is_band_diagonal() const override { return true; }

                virtual bool is_symmetric() const override { return true; }

                virtual abstract_matrix<ELEMENTTYPE>& add(
                    const abstract_matrix<ELEMENTTYPE>& babs,
                    ELEMENTTYPE modifier = 1.0 ) override {
                    if ( !is_this_subclass( babs ) ) {
                        return abstract_matrix<ELEMENTTYPE>::add( babs,
                                                                  modifier );
                    }
                    const symmetric_matrix<ELEMENTTYPE> *b = downcast( &babs );

                    this->check_size( *b );
                    // upper half
                    size_t skip = 0, add_on = 0, both = 0;
                    if ( _n_offdiag > b->_n_offdiag ) {
                        both = b->_n_offdiag;
                    } else if ( _n_offdiag < b->_n_offdiag ) {
                        both   = _n_offdiag;
                        add_on = b->_n_offdiag - _n_offdiag;
                    } else {
                        both = _n_offdiag;
                    }
                    for ( size_t ind1 = 0; ind1 < both; ind1++ ) {
                        contents()[ ind1 + 1 ] = NCPA::math::add_vectors(
                            contents()[ ind1 + 1 ],
                            NCPA::math::scale_vector( b->contents()[ ind1 + 1 ],
                                                      modifier ) );
                    }
                    contents().insert( contents().end(),
                                      b->contents().end() - add_on,
                                      b->contents().end() );
                    _n_offdiag += add_on;
                    for ( size_t ind1 = _n_offdiag - add_on;
                          ind1 < contents().size(); ind1++ ) {
                        contents()[ ind1 ] = NCPA::math::scale_vector(
                            contents()[ ind1 ], modifier );
                    }
                    contents()[ 0 ] = NCPA::math::add_vectors(
                        contents()[ 0 ], NCPA::math::scale_vector(
                                            b->contents()[ 0 ], modifier ) );
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& add(
                    ELEMENTTYPE b ) override {
                    for ( size_t ind1 = 0; ind1 < contents().size(); ind1++ ) {
                        for ( size_t ind2 = 0; ind2 < _max_ind2( ind1 );
                              ind2++ ) {
                            contents()[ ind1 ][ ind2 ] += b;
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual bool is_this_subclass(
                    const abstract_matrix<ELEMENTTYPE>& b ) const override {
                    if ( auto *derived
                         = dynamic_cast<const symmetric_matrix<ELEMENTTYPE> *>(
                             &b ) ) {
                        return true;
                    } else {
                        return false;
                    }
                }

                using band_diagonal_matrix<ELEMENTTYPE>::right_multiply;
                // virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                //     right_multiply( const abstract_vector<ELEMENTTYPE>& x )
                //         const override {
                //     if ( columns() != x.size() ) {
                //         std::ostringstream oss;
                //         oss << "Size mismatch in matrix-vector "
                //                "multiplication: "
                //             << columns() << " columns in matrix vs "
                //             << x.size() << " elements in vector";
                //         throw std::invalid_argument( oss.str() );
                //     }
                //     std::unique_ptr<abstract_vector<ELEMENTTYPE>> b(
                //         new dense_vector<ELEMENTTYPE>( rows() ) );
                //     int n = (int)rows();
                //     for ( int i = 0; i < n; i++ ) {
                //         int minloop = std::max( 0, i - (int)this->lower_bandwidth() ),
                //             maxloop = std::min( i + (int)this->upper_bandwidth() + 1, (int)this->columns() );
                //         ELEMENTTYPE bval = _zero;
                //         for (int k = minloop; k < maxloop; k++) {
                //             bval += this->get( i, k ) * x.get( k );
                //         }
                //         b->set( i, bval );
                //     }
                //         // int k            = i;
                //         // int tmploop      = std::min( bw, n - k );
                //         // std::cout << "i = " << i << ", tmploop = " << tmploop << ", looping j in [" << std::max(0,-k)<< std::endl;
                //         // ELEMENTTYPE bval = _zero;

                //         // for ( int j = std::max( 0, -k ); j < tmploop; j++ ) {
                //         //     std::cout << "Getting contents()[" << j << "," << i << "] and x[" << j+k << "]" << std::endl;
                //         //     bval += contents()[ j ][ i ] * x.get( j + k );
                //         // }
                //         // b->set( i, bval );
                //     // }
                //     return b;
                // }

                using band_diagonal_matrix<ELEMENTTYPE>::left_multiply;
                // virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                //     left_multiply( const abstract_vector<ELEMENTTYPE>& x )
                //         const override {
                //     if ( this->rows() != x.size() ) {
                //         std::ostringstream oss;
                //         oss << "Size mismatch in vector-matrix "
                //                "multiplication: "
                //             << x.size() << " elements in vector vs "
                //             << columns() << " rows in matrix";
                //         throw std::invalid_argument( oss.str() );
                //     }
                //     std::unique_ptr<abstract_vector<ELEMENTTYPE>> b(
                //         new dense_vector<ELEMENTTYPE>( this->columns() ) );
                //     int n = (int)columns();
                //     for ( int i = 0; i < n; i++ ) {
                //         int minloop = std::max( i - (int)_n_offdiag, 0 ),
                //             maxloop = std::min( i + (int)_n_offdiag + 1, (int)this->rows() );
                //         ELEMENTTYPE bval = _zero;
                //         for ( int k = minloop; k < maxloop; k++) {
                //             bval += x.get( k ) * this->get( k, i );
                //         }
                //         b->set( i, bval );
                //     }
                //     return b;
                // }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> multiply(
                    const abstract_matrix<ELEMENTTYPE>& other )
                    const override {
                    if ( !is_this_subclass( other ) ) {
                        return abstract_matrix<ELEMENTTYPE>::multiply( other );
                    }
                    const symmetric_matrix<ELEMENTTYPE> *b
                        = downcast( &other );

                    // size_t new_n_lower = this->_n_lower + b->_n_lower,
                    //        new_n_upper = this->_n_upper + b->_n_upper;
                    size_t new_n_offdiag = this->_n_offdiag + b->_n_offdiag;
                    std::unique_ptr<abstract_matrix<ELEMENTTYPE>> product(
                        new symmetric_matrix<ELEMENTTYPE>(
                            this->rows(), b->columns(), new_n_offdiag ) );
                    int prows = (int)( product->rows() );
                    int pcols = (int)( product->columns() );
                    for ( int r = 0; r < prows; r++ ) {
                        for ( int c = r;
                              c < std::min( (int)( product->columns() ),
                                            r + (int)new_n_offdiag + 1 );
                              c++ ) {
                            ELEMENTTYPE val = _zero;
                            int kmin        = std::max(
                                std::max( 0, r - (int)_n_offdiag ),
                                std::max( 0, c - (int)( b->_n_offdiag ) ) );
                            int kmax = std::min(
                                std::min( (int)columns(),
                                          r + (int)_n_offdiag + 1 ),
                                std::min( (int)( b->rows() ),
                                          c + (int)( b->_n_offdiag ) + 1 ) );
                            for ( int k = kmin; k < kmax; k++ ) {
                                val += this->get( r, k ) * b->get( k, c );
                            }
                            product->set( (size_t)r, (size_t)c, val );
                        }
                    }
                    return product;
                }

                // virtual abstract_matrix<ELEMENTTYPE>& scale(
                //     const abstract_matrix<ELEMENTTYPE>& b ) override {
                //     this->check_size( b );
                //     int r, c;
                //     for ( size_t ind1 = 0; ind1 < contents().size(); ind1++ )
                //     {
                //         for ( size_t ind2 = _min_ind2( ind1 );
                //               ind2 < _max_ind2( ind1 ); ind2++ ) {
                //             if ( internal2rowcol( ind1, ind2, r, c ) ) {
                //                 contents()[ ind1 ][ ind2 ] *= b.get( r, c );
                //             }
                //         }
                //     }
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;;
                // }

                // virtual abstract_matrix<ELEMENTTYPE>& scale(
                //     ELEMENTTYPE val ) override {
                //     for ( auto it1 = contents().begin(); it1 !=
                //     contents().end();
                //           ++it1 ) {
                //         for ( auto it2 = it1->begin(); it2 != it1->end();
                //               ++it2 ) {
                //             *it2 *= val;
                //         }
                //     }
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;;
                // }

                virtual abstract_matrix<ELEMENTTYPE>& identity(
                    size_t nrows, size_t ncols ) override {
                    this->clear();
                    _nrows     = nrows;
                    _ncols     = ncols;
                    _n_offdiag = 0;
                    std::vector<ELEMENTTYPE> diag(
                        this->diagonal_size( 0 ),
                        NCPA::math::one<ELEMENTTYPE>() );
                    contents().push_back( diag );
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_diagonal(
                    size_t nvals, const ELEMENTTYPE *vals,
                    int offset = 0 ) override {
                    if ( nvals > this->diagonal_size( offset ) ) {
                        throw std::out_of_range(
                            "Too many values for requested diagonal" );
                    }
                    offset = std::abs( offset );

                    while ( offset > _n_offdiag ) {
                        this->_add_offdiagonal();
                    }

                    for ( size_t i = 0; i < nvals; i++ ) {
                        contents()[ offset ][ i ] = vals[ i ];
                    }

                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    get_diagonal( int offset = 0 ) const override {
                    std::vector<ELEMENTTYPE> diag(
                        this->diagonal_size( offset ), _zero );
                    int ind1 = std::abs( offset );
                    if ( ind1 < contents().size() ) {
                        diag.assign( contents()[ ind1 ].begin(),
                                     contents()[ ind1 ].begin()
                                         + _max_ind2( ind1 ) );
                    }
                    return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                        new dense_vector<ELEMENTTYPE>( diag ) );
                }

                virtual size_t bandwidth() const override {
                    return 1 + 2 * _n_offdiag;
                }

                virtual size_t lower_bandwidth() const override {
                    return _n_offdiag;
                }

                virtual size_t upper_bandwidth() const override {
                    return _n_offdiag;
                }

                virtual std::vector<size_t> band_column_indices(
                    size_t row ) const override {
                    std::vector<size_t> inds;
                    for ( int i = std::max( 0, (int)row - (int)_n_offdiag );
                          i < std::min( (int)columns(),
                                        (int)row + (int)_n_offdiag + 1 );
                          i++ ) {
                        inds.push_back( (size_t)i );
                    }
                    return inds;
                }

                virtual std::vector<size_t> band_row_indices(
                    size_t col ) override {
                    std::vector<size_t> inds;
                    for ( int i = std::max( 0, (int)col - (int)_n_offdiag );
                          i < std::min( (int)rows(),
                                        (int)col + (int)_n_offdiag + 1 );
                          i++ ) {
                        inds.push_back( (size_t)i );
                    }
                    return inds;
                }

            protected:
                size_t _nrows, _ncols, _n_offdiag;
                // std::vector<std::vector<ELEMENTTYPE>> _contents;

                const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();

                void _prep_diagonals() {
                    contents().clear();
                    size_t n_diag = _n_offdiag + 1;
                    std::vector<ELEMENTTYPE> diag( this->diagonal_size( 0 ) );
                    for ( size_t i = 0; i < n_diag; i++ ) {
                        contents().push_back( diag );
                    }
                }

                // void _add_subdiagonal() {
                //     std::vector<ELEMENTTYPE> diag( this->diagonal_size( 0 )
                //     ); contents().insert( contents().cbegin(), diag );
                //     _n_lower++;
                // }

                void _add_offdiagonal() {
                    std::vector<ELEMENTTYPE> diag( this->diagonal_size( 0 ) );
                    contents().push_back( diag );
                    _n_offdiag++;
                }

                // size_t _min_ind2( size_t ind1 ) const {
                //     return (size_t)std::max( (int)_n_offdiag - (int)ind1, 0
                //     );
                // }

                size_t _max_ind2( size_t ind1 ) const {
                    return (size_t)std::min(
                        this->rows(), this->rows() + _n_offdiag - ind1 );
                }

                virtual abstract_matrix<ELEMENTTYPE>& _add(
                    const symmetric_matrix<ELEMENTTYPE>& b,
                    ELEMENTTYPE modifier = 1.0 ) {
                    this->check_size( b );
                    size_t new_n_offdiag
                        = std::max( _n_offdiag, b._n_offdiag );
                    std::vector<std::vector<ELEMENTTYPE>> newcontents(
                        new_n_offdiag + 1 );

                    // diagonals first
                    newcontents[ 0 ] = NCPA::math::add_vectors(
                        contents()[ 0 ], NCPA::math::scale_vector(
                                            b.contents()[ 0 ], modifier ) );

                    // superdiagonals
                    for ( size_t n = 1; n < new_n_offdiag; n++ ) {
                        if ( n <= _n_offdiag && n <= b._n_offdiag ) {
                            newcontents[ n ] = NCPA::math::add_vectors(
                                contents()[ n ],
                                NCPA::math::scale_vector( b.contents()[ n ],
                                                          modifier ) );
                        } else if ( n <= _n_offdiag ) {
                            newcontents[ n ] = contents()[ n ];
                        } else if ( n <= b._n_offdiag ) {
                            newcontents[ n ] = NCPA::math::scale_vector(
                                b.contents()[ n ], modifier );
                        }
                    }
                    contents()  = newcontents;
                    _n_offdiag = new_n_offdiag;
                    RETURN_THIS_AS_ABSTRACT_MATRIX;;
                }
        };

    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::symmetric_matrix<T>& a,
                  NCPA::linear::symmetric_matrix<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::abstract_matrix<T>&>( a ),
            static_cast<NCPA::linear::abstract_matrix<T>&>( b ) );
    swap( a._nrows, b._nrows );
    swap( a._ncols, b._ncols );
    swap( a._n_offdiag, b._n_offdiag );
    swap( a.contents(), b.contents() );
}
