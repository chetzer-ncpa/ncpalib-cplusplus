#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/band_diagonal_matrix.hpp"
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

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::tridiagonal_matrix,
                                             ELEMENTTYPE );

namespace NCPA {
    namespace linear {


        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class tridiagonal_matrix<ELEMENTTYPE,
                                     _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public band_diagonal_matrix<ELEMENTTYPE> {
            public:
                friend class basic_tridiagonal_linear_system_solver<
                    ELEMENTTYPE>;

                // explicitly inherit these
                using band_diagonal_matrix<ELEMENTTYPE>::set;
                using band_diagonal_matrix<ELEMENTTYPE>::set_row;
                using band_diagonal_matrix<ELEMENTTYPE>::set_column;
                using band_diagonal_matrix<ELEMENTTYPE>::set_diagonal;
                using band_diagonal_matrix<ELEMENTTYPE>::zero;
                using band_diagonal_matrix<ELEMENTTYPE>::scale;

                tridiagonal_matrix( size_t nrows, size_t ncols,
                                    bool haslower = false,
                                    bool hasupper = false ) :
                    band_diagonal_matrix<ELEMENTTYPE>( nrows, ncols,
                                                       ( haslower ? 1 : 0 ),
                                                       ( hasupper ? 1 : 0 ) ) {
                }

                tridiagonal_matrix() :
                    band_diagonal_matrix<ELEMENTTYPE>( 0, 0, 0, 0 ) {}

                tridiagonal_matrix(
                    const tridiagonal_matrix<ELEMENTTYPE>& other ) :
                    band_diagonal_matrix<ELEMENTTYPE>( other ) {}

                tridiagonal_matrix(
                    const abstract_matrix<ELEMENTTYPE>& other ) :
                    tridiagonal_matrix<ELEMENTTYPE>() {
                    this->copy( other );
                }

                tridiagonal_matrix(
                    tridiagonal_matrix<ELEMENTTYPE>&& source ) noexcept :
                    tridiagonal_matrix<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                // constructors, destructors, copying, and assignment
                virtual ~tridiagonal_matrix() {}

                friend void ::swap<ELEMENTTYPE>(
                    tridiagonal_matrix<ELEMENTTYPE>& a,
                    tridiagonal_matrix<ELEMENTTYPE>& b ) noexcept;

                tridiagonal_matrix<ELEMENTTYPE>& operator=(
                    tridiagonal_matrix<ELEMENTTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> clone()
                    const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new tridiagonal_matrix( *this ) );
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                    fresh_clone() const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new tridiagonal_matrix() );
                }

                virtual abstract_matrix<ELEMENTTYPE>& copy(
                    const abstract_matrix<ELEMENTTYPE>& other ) override {
                    this->resize( other.rows(), other.columns() );
                    this->set_diagonal( *other.get_diagonal() );
                    if ( other.upper_bandwidth() > 0 ) {
                        this->set_diagonal( *other.get_diagonal( 1 ), 1 );
                    }
                    if ( other.lower_bandwidth() > 0 ) {
                        this->set_diagonal( *other.get_diagonal( -1 ), -1 );
                    }

                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                }

                const tridiagonal_matrix<ELEMENTTYPE> *downcast(
                    const abstract_matrix<ELEMENTTYPE> *in ) const {
                    return dynamic_cast<
                        const tridiagonal_matrix<ELEMENTTYPE> *>( in );
                }

                // virtual abstract_matrix<ELEMENTTYPE>& upcast()
                //     override {
                //     return *static_cast<abstract_matrix<ELEMENTTYPE> *>(
                //         this );
                // }

                virtual std::string id() const override {
                    return "NCPA tridiagonal matrix";
                }

                // virtual size_t rows() const override { return _nrows; }

                // virtual size_t columns() const override { return _ncols; }

                // virtual abstract_matrix<ELEMENTTYPE>& clear() override {
                //     _nrows   = 0;
                //     _ncols   = 0;
                //     _n_lower = 0;
                //     _n_upper = 0;
                //     this->contents().clear();
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }

                // virtual bool rowcol2internal( const size_t& row,
                //                               const size_t& col, int& ind1,
                //                               int& ind2 ) const {
                //     ind1 = (int)col - (int)row + (int)_n_lower;
                //     ind2 = (int)row;
                //     return ( ind1 >= 0
                //              && ind1 < (int)( this->contents().size() )
                //              && row < this->rows() && col < this->columns()
                //              );
                // }

                // virtual bool internal2rowcol( const size_t& ind1,
                //                               const size_t& ind2, int& row,
                //                               int& col ) const {
                //     col = (int)ind1 + (int)ind2 - (int)_n_lower;
                //     row = (int)ind2;
                //     return ( col >= 0 && col < this->columns()
                //              && row < this->rows() );
                //     // if ( icol < 0 || icol <= columns() ) {
                //     //     std::ostringstream oss;
                //     //     oss << "Attempted to retrieve invalid internal
                //     value
                //     //     "
                //     //            "at ind1="
                //     //         << ind1 << ", ind2=" << ind2
                //     //         << " (_n_lower = " << _n_lower << ")";
                //     //     throw std::out_of_range( oss.str() );
                //     // }
                //     // return (size_t)icol;
                // }

                // virtual bool row_in_range( size_t row ) const {
                //     return ( row >= 0 && row < this->rows() );
                // }


                virtual abstract_matrix<ELEMENTTYPE>& set(
                    size_t row, size_t col, ELEMENTTYPE val ) override {
                    this->check_size( row, col );
                    if ( std::abs( (int)row - (int)col ) > 1 ) {
                        std::ostringstream oss;
                        oss << "tridiagonal_matrix.set(): point [" << row
                            << "," << col << "] is out of band";
                        throw std::range_error( oss.str() );
                    }
                    int ind1, ind2;
                    while ( !this->rowcol2internal( row, col, ind1, ind2 ) ) {
                        if ( ind1 < 0 ) {
                            this->_add_subdiagonal();
                        } else {
                            this->_add_superdiagonal();
                        }
                    }
                    this->contents().at( ind1 ).at( row ) = val;
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                }

                // virtual abstract_matrix<ELEMENTTYPE>& set_safe(
                //     size_t row, size_t col, ELEMENTTYPE val ) {
                //     this->check_size( row, col );
                //     int ind1, ind2;
                //     if ( this->rowcol2internal( row, col, ind1, ind2 ) ) {
                //         this->contents().at( ind1 ).at( row ) = val;
                //     } else {
                //         std::ostringstream oss;
                //         oss << "Element [" << row << ", " << col
                //             << "] is out of band.";
                //         throw std::out_of_range( oss.str() );
                //     }
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }

                // virtual abstract_matrix<ELEMENTTYPE>& set_row(
                //     size_t row, size_t nvals, const size_t *column_inds,
                //     const ELEMENTTYPE *vals ) override {
                //     for ( size_t i = 0; i < nvals; i++ ) {
                //         this->set_safe( row, column_inds[ i ], vals[ i ] );
                //     }
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }

                virtual abstract_matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    size_t nelems;
                    if ( row == 0 || row == this->rows() - 1 ) {
                        nelems = 2;
                    } else {
                        nelems = 3;
                    }
                    if ( vals.size() == nelems ) {
                        return this->set_row(
                            row,
                            NCPA::arrays::index_vector<size_t>(
                                nelems, ( row == 0 ? 0 : row - 1 ) ),
                            vals );
                    }
                    throw std::invalid_argument(
                        "tridiagonal_matrix.set_row(): value vector size must "
                        "be 3 (2 for first or last row)" );
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_row(
                    size_t row, ELEMENTTYPE val ) override {
                    return set_row(
                        row,
                        std::vector<ELEMENTTYPE>(
                            ( row == 0 || row == this->rows() - 1 ? 2 : 3 ),
                            val ) );
                }

                // virtual abstract_matrix<ELEMENTTYPE>& set_column(
                //     size_t column, size_t nvals, const size_t *row_inds,
                //     const ELEMENTTYPE *vals ) override {
                //     for ( size_t i = 0; i < nvals; i++ ) {
                //         this->set_safe( row_inds[ i ], column, vals[ i ] );
                //     }
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }

                virtual abstract_matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    size_t nelems;
                    if ( col == 0 || col == this->columns() - 1 ) {
                        nelems = 2;
                    } else {
                        nelems = 3;
                    }
                    if ( vals.size() == nelems ) {
                        return this->set_column(
                            col,
                            NCPA::arrays::index_vector<size_t>(
                                nelems, ( col == 0 ? 0 : col - 1 ) ),
                            vals );
                    }
                    throw std::invalid_argument(
                        "tridiagonal_matrix.set_column(): value vector size "
                        "must be 3 (2 for first or last column)" );
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_column(
                    size_t col, ELEMENTTYPE val ) override {
                    return set_column(
                        col,
                        std::vector<ELEMENTTYPE>(
                            ( col == 0 || col == ( this->columns() - 1 ) ? 2
                                                                         : 3 ),
                            val ) );
                }

                // virtual abstract_matrix<ELEMENTTYPE>& as_array(
                //     size_t& nrows, size_t& ncols,
                //     ELEMENTTYPE **& vals ) override {
                //     if ( vals == nullptr ) {
                //         nrows = rows();
                //         ncols = columns();
                //         vals
                //             = NCPA::arrays::zeros<ELEMENTTYPE>( nrows, ncols
                //             );
                //     } else {
                //         if ( nrows == 0 ) {
                //             nrows = rows();
                //         } else if ( nrows != rows() ) {
                //             throw std::invalid_argument(
                //                 "Wrong number of rows requested" );
                //         }
                //         if ( ncols == 0 ) {
                //             ncols = columns();
                //         } else if ( ncols != columns() ) {
                //             throw std::invalid_argument(
                //                 "Wrong number of columns requested" );
                //         }
                //         NCPA::arrays::fill( vals, nrows, ncols, _zero );
                //     }
                //     int row, col;
                //     for ( size_t ind1 = 0; ind1 < this->contents().size();
                //           ind1++ ) {
                //         for ( size_t ind2 = _min_ind2( ind1 );
                //               ind2 < _max_ind2( ind1 ); ind2++ ) {
                //             if ( internal2rowcol( ind1, ind2, row, col ) ) {
                //                 vals[ row ][ col ]
                //                     = this->contents()[ ind1 ][ ind2 ];
                //             }
                //         }
                //     }
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }

                // @todo make read-write vector view for columns
                // virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                // get_row(
                //     size_t row ) const override {
                //     std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                //         new sparse_vector<ELEMENTTYPE>( this->columns() ) );
                //     int r, c;
                //     for ( size_t ind1 = 0; ind1 < this->contents().size();
                //           ind1++ ) {
                //         if ( internal2rowcol( ind1, row, r, c ) ) {
                //             v->set( c, this->contents()[ ind1 ][ row ] );
                //         }
                //     }
                //     return v;
                // }

                // virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                //     get_column( size_t column ) const override {
                //     std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                //         new sparse_vector<ELEMENTTYPE>( this->rows() ) );
                //     int ind1, ind2;
                //     for ( size_t r = 0; r < this->rows(); r++ ) {
                //         if ( rowcol2internal( r, column, ind1, ind2 ) ) {
                //             v->set( r, this->contents()[ ind1 ][ ind2 ] );
                //         }
                //     }
                //     // for ( size_t ind1 = 0; ind1 <
                //     this->contents().size();
                //     // ind1++ )
                //     // {
                //     //     size_t ind2 = (size_t)std::max(
                //     //         (int)( column + _n_lower ) - (int)ind1, 0 );
                //     //     ind2 = std::min( ind2, this->contents()[ ind1
                //     //     ].size() - 1
                //     //     ); v->set( ind2, this->contents()[ ind1 ][ ind2 ]
                //     );
                //     // }
                //     return v;
                // }

                // virtual const ELEMENTTYPE& get( size_t row,
                //                                 size_t col ) const override
                //                                 {
                //     // this->check_size( row, col );
                //     int ind1, ind2;
                //     if ( rowcol2internal( row, col, ind1, ind2 ) ) {
                //         return this->contents()[ ind1 ][ ind2 ];
                //     } else {
                //         return _zero;
                //     }
                // }

                // virtual abstract_matrix<ELEMENTTYPE>& resize(
                //     size_t r, size_t c ) override {
                //     if ( this->contents().size() == 0 ) {
                //         // starting from scratch
                //         _nrows   = r;
                //         _ncols   = c;
                //         _n_lower = 0;
                //         _n_upper = 0;
                //         _prep_diagonals();
                //     } else {
                //         // did the diagonal size increase or decrease?
                //         int ddiag = (int)std::min( r, c )
                //                   - (int)std::min( rows(), columns() );
                //         if ( ddiag > 0 ) {
                //             // increase each diagonal by ddiag elements
                //             for ( auto it = this->contents().begin();
                //                   it != this->contents().end(); ++it ) {
                //                 it->insert( it->cend(), (size_t)ddiag, _zero
                //                 );
                //             }
                //         } else if ( ddiag < 0 ) {
                //             for ( auto it = this->contents().begin();
                //                   it != this->contents().end(); ++it ) {
                //                 it->erase( it->cend() - size_t( -ddiag ),
                //                            it->cend() );
                //             }
                //         }
                //         _nrows = r;
                //         _ncols = c;
                //     }

                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }

                // virtual abstract_matrix<ELEMENTTYPE>& transpose() override {
                //     std::vector<std::vector<ELEMENTTYPE>> newcontents
                //         = this->contents();
                //     std::reverse( newcontents.begin(), newcontents.end() );
                //     std::swap( _n_lower, _n_upper );
                //     std::swap( _nrows, _ncols );
                //     if ( _n_lower > 0 ) {
                //         // take the zeros from the end and put them up front
                //         for ( size_t ind1 = 0; ind1 < _n_lower; ind1++ ) {
                //             size_t invalid = _n_lower - ind1;
                //             // NCPA_DEBUG << "For ind1 = " << ind1 << ":
                //             moving
                //             // "
                //             //            << invalid
                //             //            << " zeros from back to front"
                //             //            << std::endl;
                //             for ( size_t i = 0; i < invalid; i++ ) {
                //                 newcontents[ ind1 ].pop_back();
                //                 newcontents[ ind1 ].insert(
                //                     newcontents[ ind1 ].begin(), 1, _zero );
                //                 // NCPA_DEBUG << "Cycled one." << std::endl;
                //             }
                //         }
                //     }
                //     if ( _n_upper > 0 ) {
                //         // take the zeros from the front and put them at the
                //         // end
                //         for ( size_t ind1 = _n_lower + 1;
                //               ind1 < newcontents.size(); ind1++ ) {
                //             size_t invalid = ind1 - _n_lower;
                //             // NCPA_DEBUG << "For ind1 = " << ind1 << ":
                //             moving
                //             // "
                //             //            << invalid
                //             //            << " zeros from back to front"
                //             //            << std::endl;
                //             for ( size_t i = 0; i < invalid; i++ ) {
                //                 newcontents[ ind1 ].erase(
                //                     newcontents[ ind1 ].cbegin() );
                //                 newcontents[ ind1 ].push_back( _zero );
                //                 // NCPA_DEBUG << "Cycled one." << std::endl;
                //             }
                //         }
                //     }
                //     _contents = newcontents;
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }

                virtual abstract_matrix<ELEMENTTYPE>& swap_rows(
                    size_t ind1, size_t ind2 ) override {
                    throw std::logic_error(
                        "tridiagonal_matrix.swap_rows(): cannot swap rows of "
                        "a tridiagonal matrix!" );
                }

                virtual abstract_matrix<ELEMENTTYPE>& swap_columns(
                    size_t ind1, size_t ind2 ) override {
                    throw std::logic_error(
                        "tridiagonal_matrix.swap_columns(): cannot swap "
                        "columns of a tridiagonal matrix!" );
                }

                // virtual abstract_matrix<ELEMENTTYPE>& zero(
                //     size_t row, size_t col ) override {
                //     return set( row, col, _zero );
                // }

                virtual abstract_matrix<ELEMENTTYPE>& zero() override {
                    this->_n_lower = 0;
                    this->_n_upper = 0;
                    this->contents().assign(
                        1, std::vector<ELEMENTTYPE>( this->diagonal_size() ) );
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
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
                        return abstract_matrix<ELEMENTTYPE>::equals( other );
                    }
                    const tridiagonal_matrix<ELEMENTTYPE> *b
                        = downcast( &other );

                    return ( this->rows() == b->rows()
                             && this->columns() == b->columns()
                             && this->contents() == b->contents() );
                }

                virtual bool is_identity() const override {
                    if ( !this->is_square() || !this->is_diagonal() ) {
                        return false;
                    }
                    for ( size_t i = 0; i < this->contents()[ 0 ].size();
                          i++ ) {
                        if ( this->contents()[ 0 ][ i ]
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
                    return this->_n_lower == 0 && this->_n_upper == 0;
                }

                virtual bool is_tridiagonal() const override { return true; }

                virtual bool is_upper_triangular() const override {
                    return ( this->is_empty() || this->lower_bandwidth() == 0
                             || this->get_diagonal( -1 )->is_zero() );
                }

                virtual bool is_lower_triangular() const override {
                    // if ( this->is_empty() || this->lower_bandwidth() == 0) {
                    //     return true;
                    // }
                    return ( this->is_empty() || this->upper_bandwidth() == 0
                             || this->get_diagonal( 1 )->is_zero() );
                }

                virtual abstract_matrix<ELEMENTTYPE>& add(
                    const abstract_matrix<ELEMENTTYPE>& babs,
                    ELEMENTTYPE modifier = 1.0 ) override {
                    if ( !is_this_subclass( babs ) ) {
                        return band_diagonal_matrix<ELEMENTTYPE>::add(
                            babs, modifier );
                    }
                    const tridiagonal_matrix<ELEMENTTYPE> *b
                        = downcast( &babs );

                    this->check_size( *b );
                    // lower half
                    if ( b->lower_bandwidth() > 0 ) {
                        if ( this->lower_bandwidth() > 0 ) {
                            this->contents()[ 0 ] = NCPA::math::add_vectors(
                                this->contents()[ 0 ],
                                NCPA::math::scale_vector( b->contents()[ 0 ],
                                                          modifier ) );
                        } else {
                            this->contents().insert( this->contents().cbegin(),
                                                     b->contents()[ 0 ] );
                        }
                    }
                    // upper half
                    if ( b->upper_bandwidth() > 0 ) {
                        if ( this->lower_bandwidth() > 0 ) {
                            this->contents().back() = NCPA::math::add_vectors(
                                this->contents().back(),
                                NCPA::math::scale_vector( b->contents().back(),
                                                          modifier ) );
                        } else {
                            this->contents().push_back( b->contents().back() );
                        }
                    }
                    this->contents()[ this->_n_lower ]
                        = NCPA::math::add_vectors(
                            this->contents()[ this->_n_lower ],
                            NCPA::math::scale_vector(
                                b->contents()[ b->_n_lower ], modifier ) );
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                }

                virtual abstract_matrix<ELEMENTTYPE>& add(
                    ELEMENTTYPE b ) override {
                    for ( size_t ind1 = 0; ind1 < this->contents().size();
                          ind1++ ) {
                        for ( size_t ind2 = this->_min_ind2( ind1 );
                              ind2 < this->_max_ind2( ind1 ); ind2++ ) {
                            this->contents()[ ind1 ][ ind2 ] += b;
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                }

                virtual bool is_band_diagonal() const override { return true; }

                virtual bool is_this_subclass(
                    const abstract_matrix<ELEMENTTYPE>& b ) const override {
                    if ( auto *derived = dynamic_cast<
                             const tridiagonal_matrix<ELEMENTTYPE> *>( &b ) ) {
                        return true;
                    } else {
                        return false;
                    }
                }

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
                //     int n = (int)rows();  // , bw = (int)bandwidth();
                //     for ( int i = 0; i < n; i++ ) {
                //         int minloop
                //             = std::max( 0, i - (int)this->lower_bandwidth()
                //             ), maxloop = std::min( i +
                //             (int)this->upper_bandwidth() + 1,
                //                         (int)this->columns() );
                //         ELEMENTTYPE bval = _zero;
                //         for ( int k = minloop; k < maxloop; k++ ) {
                //             bval += this->get( i, k ) * x.get( k );
                //         }
                //         b->set( i, bval );
                //     }
                //     // for ( int i = 0; i < n; i++ ) {
                //     //     int k            = i - (int)_n_lower;
                //     //     int tmploop      = std::min( bw, n - k );
                //     //     ELEMENTTYPE bval = _zero;
                //     //     for ( int j = std::max( 0, -k ); j < tmploop; j++
                //     )
                //     //     {
                //     //         bval += this->contents()[ j ][ i ] * x.get( j
                //     +
                //     //         k );
                //     //     }
                //     //     b->set( i, bval );
                //     // }
                //     return b;
                // }

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
                //         int minloop
                //             = std::max( i - (int)this->upper_bandwidth(), 0
                //             ), maxloop = std::min( i +
                //             this->lower_bandwidth() + 1,
                //                         this->rows() );
                //         ELEMENTTYPE bval = _zero;
                //         for ( int k = minloop; k < maxloop; k++ ) {
                //             bval += x.get( k ) * this->get( k, i );
                //         }
                //         b->set( i, bval );
                //         // int k            = i - (int)_n_upper;
                //         // int tmploop      = std::min( bw, n - k );
                //         // ELEMENTTYPE bval = _zero;
                //         // for ( int j = std::max( 0, -k ); j < tmploop; j++
                //         )
                //         // {
                //         //     bval += this->contents()[ j ][ i ] * x.get( j
                //         +
                //         //     k );
                //         // }
                //         // b->set( i, bval );
                //     }
                //     return b;
                // }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> multiply(
                    const abstract_matrix<ELEMENTTYPE>& other )
                    const override {
                    if ( !other.is_diagonal() ) {
                        return band_diagonal_matrix<ELEMENTTYPE>::multiply(
                            other );
                    }
                    if ( auto *derived = dynamic_cast<
                             const band_diagonal_matrix<ELEMENTTYPE> *>(
                             &other ) ) {
                        std::unique_ptr<abstract_matrix<ELEMENTTYPE>> product(
                            new tridiagonal_matrix<ELEMENTTYPE>(
                                this->rows(), other.columns(),
                                ( this->lower_bandwidth() == 1 ),
                                ( this->upper_bandwidth() == 1 ) ) );
                        band_diagonal_matrix<ELEMENTTYPE> *bproduct
                            = dynamic_cast<
                                band_diagonal_matrix<ELEMENTTYPE> *>(
                                product.get() );

                        band_diagonal_matrix<ELEMENTTYPE>::
                            _do_band_diagonal_multiply( *this, *derived,
                                                        *bproduct );
                        return product;
                    } else {
                        return band_diagonal_matrix<ELEMENTTYPE>::multiply(
                            other );
                    }
                }

                // int prows = (int)( product->rows() );
                // int pcols = (int)( product->columns() );
                // for ( int r = 0; r < prows; r++ ) {
                //     for ( int c = std::max( 0, r - (int)new_n_lower );
                //           c < std::min( (int)( product->columns() ),
                //                         r + (int)new_n_upper + 1 );
                //           c++ ) {
                //         ELEMENTTYPE val = _zero;
                //         int kmin        = std::max(
                //             std::max( 0,
                //                              r -
                //                              (int)this->lower_bandwidth()
                //                              ),
                //             std::max(
                //                 0, c - (int)( b->upper_bandwidth() ) )
                //                 );
                //         int kmax = std::min(
                //             std::min( (int)columns(),
                //                       r + (int)this->upper_bandwidth()
                //                           + 1 ),
                //             std::min( (int)( b->rows() ),
                //                       c + (int)( b->lower_bandwidth() )
                //                           + 1 ) );
                //         for ( int k = kmin; k < kmax; k++ ) {
                //             val += this->get( r, k ) * b->get( k, c );
                //         }
                //         product->set( (size_t)r, (size_t)c, val );
                //     }
                // }


                virtual abstract_matrix<ELEMENTTYPE>& scale(
                    const abstract_matrix<ELEMENTTYPE>& b ) override {
                    this->check_size( b );
                    int r, c;
                    for ( size_t ind1 = 0; ind1 < this->contents().size();
                          ind1++ ) {
                        for ( size_t ind2 = this->_min_ind2( ind1 );
                              ind2 < this->_max_ind2( ind1 ); ind2++ ) {
                            if ( this->internal2rowcol( ind1, ind2, r, c ) ) {
                                this->contents()[ ind1 ][ ind2 ]
                                    *= b.get( r, c );
                            }
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                }

                // virtual abstract_matrix<ELEMENTTYPE>& scale(
                //     ELEMENTTYPE val ) override {
                //     for ( auto it1 = this->contents().begin();
                //           it1 != this->contents().end(); ++it1 ) {
                //         for ( auto it2 = it1->begin(); it2 != it1->end();
                //               ++it2 ) {
                //             *it2 *= val;
                //         }
                //     }
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }

                // virtual abstract_matrix<ELEMENTTYPE>& identity(
                //     size_t nrows, size_t ncols ) override {
                //     this->clear();
                //     this->_nrows   = nrows;
                //     _ncols   = ncols;
                //     _n_lower = 0;
                //     _n_upper = 0;
                //     std::vector<ELEMENTTYPE> diag(
                //         this->diagonal_size( 0 ),
                //         NCPA::math::one<ELEMENTTYPE>() );
                //     this->contents().push_back( diag );
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }

                virtual abstract_matrix<ELEMENTTYPE>& set_diagonal(
                    size_t nvals, const ELEMENTTYPE *vals,
                    int offset = 0 ) override {
                    if ( std::abs( offset ) <= 1 ) {
                        return band_diagonal_matrix<ELEMENTTYPE>::set_diagonal(
                            nvals, vals, offset );
                    } else {
                        throw std::range_error(
                            "tridiagonal_matrix.set_diagonal(): requested "
                            "diagonal is out of band" );
                    }
                }

                // virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                //     get_diagonal( int offset = 0 ) const override {
                //     // std::cout << "Calling get_diagonal( " << offset << "
                //     )"
                //     //           << std::endl;
                //     std::vector<ELEMENTTYPE> diag(
                //         this->diagonal_size( offset ), _zero );
                //     // if ( offset > _n_upper || -offset > _n_lower ) {
                //     //     diag.resize( this->diagonal_size( offset ), _zero
                //     );
                //     // } else {
                //     int ind1 = (int)_n_lower + offset;
                //     if ( ind1 >= 0 && ind1 < this->contents().size() ) {
                //         // std::cout
                //         //     << "ind1 = " << ind1
                //         //     << " is in range for _n_upper = " << _n_upper
                //         //     << " and _n_lower = " << _n_lower
                //         //     << ", this->contents().size() = " <<
                //         //     this->contents().size()
                //         //     << std::endl;
                //         // std::cout << "Fetching internal index " << ind1
                //         //           << std::endl;
                //         // std::cout << "Fetching ind2 = " << _min_ind2(
                //         ind1 )
                //         //           << " to " << _max_ind2( ind1 ) <<
                //         //           std::endl;
                //         diag.assign( this->contents()[ ind1 ].begin()
                //                          + _min_ind2( ind1 ),
                //                      this->contents()[ ind1 ].begin()
                //                          + _max_ind2( ind1 ) );
                //         // } else {
                //         //     std::cout << "ind1 = " << ind1
                //         //               << " out of rangefor _n_upper = "
                //         <<
                //         //               _n_upper
                //         //               << " and _n_lower = " << _n_lower
                //         //               << ", returning zeros" <<
                //         std::endl;
                //     }
                //     return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                //         new dense_vector<ELEMENTTYPE>( diag ) );
                //     // std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                //     //     vdiag = build_vector();
                //     // vdiag->set( diag );
                //     // return vdiag;
                // }

                // virtual size_t bandwidth() const override {
                //     return this->contents().size();
                // }

                // virtual size_t lower_bandwidth() const override {
                //     return _n_lower;
                // }

                // virtual size_t upper_bandwidth() const override {
                //     return _n_upper;
                // }

                // virtual std::vector<size_t> band_column_indices(
                //     size_t row ) const {
                //     std::vector<size_t> inds;
                //     for ( int i = std::max( 0, (int)row - (int)_n_lower );
                //           i < std::min( (int)columns(),
                //                         (int)row + (int)_n_upper + 1 );
                //           i++ ) {
                //         inds.push_back( (size_t)i );
                //     }
                //     return inds;
                // }

                // virtual std::vector<size_t> band_row_indices( size_t col ) {
                //     std::vector<size_t> inds;
                //     // NCPA_DEBUG << "Getting nonzero indices for column "
                //     <<
                //     // col << std::endl;
                //     for ( int i = std::max( 0, (int)col - (int)_n_upper );
                //           i < std::min( (int)rows(),
                //                         (int)col + (int)_n_lower + 1 );
                //           i++ ) {
                //         inds.push_back( (size_t)i );
                //         // NCPA_DEBUG << "Adding index " << i << std::endl;
                //     }
                //     return inds;
                // }

                // virtual std::vector<std::vector<ELEMENTTYPE>>& contents() {
                //     return _contents;
                // }

                // virtual const std::vector<std::vector<ELEMENTTYPE>>&
                // contents()
                //     const {
                //     return _contents;
                // }

            protected:
                // size_t _nrows, _ncols, _n_lower, _n_upper;

                // std::vector<std::vector<ELEMENTTYPE>> _contents;

                const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();

                // void _prep_diagonals() {
                //     this->contents().clear();
                //     size_t n_diag = _n_lower + _n_upper + 1;
                //     // std::cout << "Adding " << n_diag << " diagonals" <<
                //     // std::endl;
                //     std::vector<ELEMENTTYPE> diag( this->diagonal_size( 0 )
                //     ); for ( size_t i = 0; i < n_diag; i++ ) {
                //         this->contents().push_back( diag );
                //     }
                //     // std::cout << "Done" << std::endl;
                // }

                // void _add_subdiagonal() {
                //     std::vector<ELEMENTTYPE> diag( this->diagonal_size( 0 )
                //     ); this->contents().insert( this->contents().cbegin(),
                //     diag ); _n_lower++;
                // }

                // void _add_superdiagonal() {
                //     std::vector<ELEMENTTYPE> diag( this->diagonal_size( 0 )
                //     ); this->contents().push_back( diag ); _n_upper++;
                // }

                // size_t _min_ind2( size_t ind1 ) const {
                //     return (size_t)std::max( (int)_n_lower - (int)ind1, 0 );
                // }

                // size_t _max_ind2( size_t ind1 ) const {
                //     return (size_t)std::min( this->rows(),
                //                              this->rows() + _n_lower - ind1
                //                              );
                // }

                // virtual abstract_matrix<ELEMENTTYPE>& _add(
                //     const band_diagonal_matrix<ELEMENTTYPE>& b,
                //     ELEMENTTYPE modifier = 1.0 ) {
                //     this->check_size( b );
                //     size_t new_n_lower = std::max( _n_lower, b._n_lower );
                //     size_t new_n_upper = std::max( _n_upper, b._n_upper );
                //     std::vector<std::vector<ELEMENTTYPE>> newcontents(
                //         new_n_lower + new_n_upper + 1 );

                //     // diagonals first
                //     newcontents[ new_n_lower ] = NCPA::math::add_vectors(
                //         this->contents()[ _n_lower ],
                //         NCPA::math::scale_vector( b.contents()[ b._n_lower
                //         ],
                //                                   modifier ) );

                //     // subdiagonals
                //     for ( size_t n = 1; n < new_n_lower; n++ ) {
                //         if ( n <= _n_lower && n <= b._n_lower ) {
                //             newcontents[ new_n_lower - n ]
                //                 = NCPA::math::add_vectors(
                //                     this->contents()[ _n_lower - n ],
                //                     NCPA::math::scale_vector(
                //                         b.contents()[ b._n_lower - n ],
                //                         modifier ) );
                //         } else if ( n <= _n_lower ) {
                //             newcontents[ new_n_lower - n ]
                //                 = this->contents()[ _n_lower - n ];
                //         } else if ( n <= b._n_lower ) {
                //             newcontents[ new_n_lower - n ]
                //                 = NCPA::math::scale_vector(
                //                     b.contents()[ b._n_lower - n ], modifier
                //                     );
                //         }
                //     }

                //     // superdiagonals
                //     for ( size_t n = 1; n < new_n_upper; n++ ) {
                //         if ( n <= _n_upper && n <= b._n_upper ) {
                //             newcontents[ new_n_lower + n ]
                //                 = NCPA::math::add_vectors(
                //                     this->contents()[ _n_lower + n ],
                //                     NCPA::math::scale_vector(
                //                         b.contents()[ b._n_lower + n ],
                //                         modifier ) );
                //         } else if ( n <= _n_upper ) {
                //             newcontents[ new_n_lower + n ]
                //                 = this->contents()[ _n_lower + n ];
                //         } else if ( n <= b._n_upper ) {
                //             newcontents[ new_n_lower + n ]
                //                 = NCPA::math::scale_vector(
                //                     b.contents()[ b._n_lower + n ], modifier
                //                     );
                //         }
                //     }
                //     _contents = newcontents;
                //     _n_lower  = new_n_lower;
                //     _n_upper  = new_n_upper;
                //     RETURN_THIS_AS_ABSTRACT_MATRIX;
                // }
        };

    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::tridiagonal_matrix<T>& a,
                  NCPA::linear::tridiagonal_matrix<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::band_diagonal_matrix<T>&>( a ),
            static_cast<NCPA::linear::band_diagonal_matrix<T>&>( b ) );
    // swap( a._nrows, b._nrows );
    // swap( a._ncols, b._ncols );
    // swap( a._n_lower, b._n_lower );
    // swap( a._n_upper, b._n_upper );
    // swap( a._contents, b._contents );
}
