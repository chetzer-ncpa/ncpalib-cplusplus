#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/sparse_vector.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <cmath>
#include <complex>
#include <initializer_list>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace NCPA {
    namespace linear {
        namespace details {
            NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( sparse_matrix,
                                                         abstract_matrix );
        }
    }  // namespace linear
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::sparse_matrix, ELEMENTTYPE );

#define RETURN_SUPERCLASS_POINTER \
    return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>( this );

namespace NCPA {
    namespace linear {

        namespace details {
            NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
                class sparse_matrix<ELEMENTTYPE,
                                    _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
                : public abstract_matrix<ELEMENTTYPE> {
                public:
                    using abstract_matrix<ELEMENTTYPE>::set_row;
                    using abstract_matrix<ELEMENTTYPE>::set_column;
                    using abstract_matrix<ELEMENTTYPE>::is_empty;
                    using abstract_matrix<ELEMENTTYPE>::diagonal_size;

                    sparse_matrix() :
                        abstract_matrix<ELEMENTTYPE>(),
                        _nrows { 0 },
                        _ncols { 0 } {}

                    sparse_matrix( size_t nrows, size_t ncols ) :
                        sparse_matrix() {
                        resize( nrows, ncols );
                    }

                    sparse_matrix( const sparse_matrix<ELEMENTTYPE>& other ) :
                        sparse_matrix<ELEMENTTYPE>() {
                        _nrows = other._nrows;
                        _ncols = other._ncols;
                        _rows.clear();
                        for ( auto it = other._rows.cbegin();
                              it != other._rows.cend(); ++it ) {
                            _rows.emplace_back( *it );
                        }
                        _cols.clear();
                        _build_column_cache();
                    }

                    sparse_matrix( const abstract_matrix<ELEMENTTYPE>& other ) :
                        dense_matrix<ELEMENTTYPE>() {
                        resize( other.rows(), other.columns() );
                        for ( auto i = 0; i < other.rows(); i++ ) {
                            auto row = other.get_row();
                            auto inds = row->nonzero_indices();
                            for (auto it = inds.begin(); it != inds.end(); ++it) {
                                set( i, *it, row->get( *it ) );
                            }
                        }
                        _build_column_cache();
                    }

                    virtual ~sparse_matrix() {}

                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        clone() const override {
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new sparse_matrix( *this ) );
                    }

                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        fresh_clone() const override {
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new sparse_matrix() );
                    }

                    sparse_matrix<ELEMENTTYPE>& operator=(
                        sparse_matrix<ELEMENTTYPE> other ) {
                        swap( *this, other );
                        _rebuild_column_cache();
                        return *this;
                    }

                    sparse_matrix<ELEMENTTYPE>& operator=(
                        abstract_matrix<ELEMENTTYPE> other ) {
                        sparse_matrix<ELEMENTTYPE> copy( other );
                        swap( *this, copy );
                        _rebuild_column_cache();
                        return *this;
                    }

                    virtual std::string id() const override {
                        return "NCPA Basic Sparse Matrix";
                    }

                    friend void ::swap<ELEMENTTYPE>(
                        sparse_matrix<ELEMENTTYPE>& a,
                        sparse_matrix<ELEMENTTYPE>& b ) noexcept;

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        build_vector( size_t n = 0 ) const override {
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new sparse_vector<ELEMENTTYPE>( n ) );
                    }

                    virtual bool equals( const abstract_matrix<ELEMENTTYPE>&
                                             other ) const override {
                        if ( rows() != other.rows() ) {
                            return false;
                        }
                        if ( columns() != other.columns() ) {
                            return false;
                        }
                        for ( auto r = 0; r < rows(); r++ ) {
                            for ( auto c = 0; c < columns(); c++ ) {
                                if ( !NCPA::math::equals(
                                         get( r, c ), other.get( r, c ) ) ) {
                                    return false;
                                }
                            }
                        }
                        return true;
                    }

                    virtual bool is_diagonal() const override {
                        if ( is_empty() ) {
                            return true;
                        }
                        for ( size_t i = 0; i < rows(); i++ ) {
                            // if there's more than one nonzero value, it's not
                            // diagonal
                            std::vector<size_t> inds
                                = _rows[ i ].nonzero_indices();
                            if ( inds.size() > 1 ) {
                                return false;
                            }
                            // if the nonzero index is not on the diagonal,
                            // it's not diagonal
                            if ( inds.size() == 1 && inds[ 0 ] != i ) {
                                return false;
                            }
                        }

                        return true;
                    }

                    virtual bool is_tridiagonal() const override {
                        if ( is_empty() ) {
                            return true;
                        }
                        size_t nz_max = 0, min_ok = 0, max_ok = 0,
                               n_diag = diagonal_size( 0 );
                        for ( int row = 0; row < n_diag; row++ ) {
                            std::vector<size_t> inds
                                = _rows[ row ].nonzero_indices();
                            if ( row == 0 ) {
                                nz_max = 2;
                                min_ok = 0;
                                max_ok = 1;
                            } else if ( row == rows() - 1 ) {
                                nz_max = 2;
                                min_ok = n_diag - 2;
                                max_ok = std::max( n_diag, columns() - 1 );
                            } else {
                                nz_max = 3;
                                min_ok = row - 1;
                                max_ok = row + 1;
                            }
                            if ( inds.size() > nz_max ) {
                                // std::cout << "Row " << row
                                //           << ": Too many nonzero entries ("
                                //           << inds.size() << " vs " << nz_max
                                //           << ")" << std::endl;
                                return false;
                            } else {
                                for ( auto it = inds.cbegin();
                                      it != inds.cend(); ++it ) {
                                    if ( *it < min_ok || *it > max_ok ) {
                                        // std::cout << "Row " << row
                                        //           << ": Nonzero index " <<
                                        //           *it
                                        //           << " exceeds allowed range
                                        //           ["
                                        //           << min_ok << "," << max_ok
                                        //           << "]" << std::endl;
                                        return false;
                                    }
                                }
                            }
                        }
                        return true;
                    }

                    virtual bool is_upper_triangular() const override {
                        if ( is_empty() ) {
                            return true;
                        }
                        size_t ndiag = diagonal_size();
                        for ( size_t i = 1; i < ndiag; i++ ) {
                            std::vector<size_t> inds
                                = _rows[ i ].nonzero_indices();
                            // std::cout << "Row " << i << ": Nonzero = {";
                            // for (auto it = inds.begin(); it != inds.end(); ++it) {
                            //     std::cout << *it << "(" << _rows[i].get(*it) << "), ";
                            // } 
                            // std::cout << "}";
                            if ( inds.size() > 0 && inds.front() < i ) {
                                // std::cout << " - failed" << std::endl;
                                return false;
                            }
                            // std::cout << std::endl;
                        }
                        return true;
                    }

                    virtual bool is_lower_triangular() const override {
                        if ( is_empty() ) {
                            return true;
                        }
                        size_t ndiag = diagonal_size();
                        for ( size_t i = 0; i < ndiag; i++ ) {
                            std::vector<size_t> inds
                                = _rows[ i ].nonzero_indices();
                            if ( inds.size() > 0 && inds.back() > i ) {
                                return false;
                            }
                        }
                        return true;
                    }

                    virtual size_t rows() const override { return _nrows; }

                    virtual size_t columns() const override { return _ncols; }

                    virtual abstract_matrix<ELEMENTTYPE>& clear() override {
                        _clear_column_cache();
                        return resize( 0, 0 );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& zero(
                        size_t row, size_t col ) override {
                        _rows[ row ].zero( col );
                        _cols[ col ].zero( row );
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& zero() override {
                        for ( auto it = _rows.begin(); it != _rows.end();
                              ++it ) {
                            it->zero();
                        }
                        _rebuild_column_cache();
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& resize(
                        size_t newrows, size_t newcols ) override {
                        for ( size_t i = 0; i < std::min( newrows, rows() ); i++) {
                            _rows[ i ].resize( newcols );
                        }
                        _rows.resize( newrows );
                        if (newrows > rows()) {
                            for (size_t i = rows(); i < newrows; i++) {
                                _rows[ i ] = sparse_vector<ELEMENTTYPE>( newcols );
                            }
                        }
                        _nrows = newrows;
                        _ncols = newcols;
                        _rebuild_column_cache();

                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& transpose()
                        override {
                        std::swap( _rows, _cols );
                        std::swap( _nrows, _ncols );
                        // _clear_column_cache();
                        RETURN_SUPERCLASS_POINTER; 
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& swap_rows(size_t ind1, size_t ind2 ) override {
                        std::swap( _rows[ ind1 ], _rows[ ind2 ]); 
                        _rebuild_column_cache();
                        RETURN_SUPERCLASS_POINTER;
                    }
                    virtual abstract_matrix<ELEMENTTYPE>& swap_columns(size_t ind1, size_t ind2) override {
                        auto col1 = get_column( ind1 );
                        auto col2 = get_column( ind2 );
                        set_column( ind1, *col2 );
                        set_column( ind2, *col1 );
                        _rebuild_column_cache();
                        RETURN_SUPERCLASS_POINTER;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& add(
                        const abstract_matrix<ELEMENTTYPE>& b,
                        ELEMENTTYPE modifier = 1.0 ) override {
                        this->check_size( b );
                        for ( size_t r = 0; r < rows(); r++ ) {
                            _rows[ r ].add( *( b.get_row( r ) ), modifier );
                            // auto brow = b.get_row( r );
                            // std::vector<size_t> inds
                            //     = _rows[ r ].nonzero_index_intersection(
                            //         *brow );
                            // for ( auto cit = inds.cbegin(); cit !=
                            // inds.cend();
                            //       ++cit ) {
                            //     set( r, *cit,
                            //          get( r, *cit )
                            //              + brow.get( *cit ) * modifier );
                            // }
                        }
                        _rebuild_column_cache();
                        RETURN_SUPERCLASS_POINTER;
                    }

                    template<typename ANYTYPE,
                             ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                                           ELEMENTTYPE )>
                    abstract_matrix<ELEMENTTYPE>& add( ANYTYPE b ) {
                        for ( size_t row = 0; row < rows(); row++ ) {
                            for ( size_t col = 0; col < columns(); col++ ) {
                                set( row, col, get( row, col ) + b );
                                // _rows[ row ].get( col ) += b;
                            }
                        }
                        _rebuild_column_cache();
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    template<typename ANYTYPE,
                             ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                                           ELEMENTTYPE )>
                    abstract_matrix<ELEMENTTYPE>& scale( ANYTYPE val ) {
                        for ( auto rit = _rows.begin(); rit != _rows.end();
                              ++rit ) {
                            rit->scale( val );
                        }
                        // for ( auto cit = _cols.begin(); cit != _cols.end();
                        //       ++cit ) {
                        //     ( *cit ).scale( val );
                        // }
                        _rebuild_column_cache();
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& scale(
                        const abstract_matrix<ELEMENTTYPE>& b ) override {
                        // for ( auto cit = _cols.begin(); cit != _cols.end();
                        //       ++cit ) {
                        //     ( *cit ).clear();
                        // }
                        for ( size_t r = 0; r < rows(); r++ ) {
                            _rows[ r ].scale( *( b.get_row( r ) ) );
                            // std::vector<size_t> inds
                            //     = _rows[ r ].nonzero_index_union(
                            //         b.get_row( r ) );
                            // std::unique_ptr<sparse_vector<ELEMENTTYPE>>
                            // newrow(
                            //     new sparse_vector<ELEMENTTYPE>( columns() )
                            //     );

                            // for ( auto cit = inds.cbegin(); cit !=
                            // inds.cend();
                            //       ++cit ) {
                            //     newrow.set( *cit, _rows[ r ].get( *cit )
                            //                            * b.get( r, *cit ) );
                            //     _cols[ *cit ].set( r, newrow.get( *cit )
                            //     );
                            // }
                            // _rows[ r ] = newrow;
                        }
                        _rebuild_column_cache();
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& identity(
                        size_t nrows, size_t ncols ) override {
                        clear();
                        resize( nrows, ncols );
                        ELEMENTTYPE one = NCPA::math::one<ELEMENTTYPE>();
                        for ( auto i = 0; i < diagonal_size( 0 ); i++ ) {
                            set( i, i, one );
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        size_t row, size_t col, ELEMENTTYPE val ) override {
                        if ( NCPA::math::is_zero( val ) ) {
                            _rows[ row ].zero( col );
                            _cols[ col ].zero( row );
                        } else {
                            _rows[ row ].set( col, val );
                            _cols[ col ].set( row, val );
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        ELEMENTTYPE val ) override {
                        if ( NCPA::math::is_zero( val ) ) {
                            zero();
                        } else {
                            for ( size_t r = 0; r < rows(); r++ ) {
                                for ( size_t c = 0; c < columns(); c++ ) {
                                    set( r, c, val );
                                }
                            }
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, size_t nvals, const size_t *column_inds,
                        const ELEMENTTYPE *vals ) override {
                        for ( size_t i = 0; i < nvals; i++ ) {
                            set( row, column_inds[ i ], vals[ i ] );
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row,
                        const std::vector<ELEMENTTYPE>& vals ) override {
                        // _rows[ row ].zero();
                        for ( size_t i = 0; i < vals.size(); i++ ) {
                            // if ( !NCPA::math::is_zero( vals[ i ] ) ) {
                            set( row, i, vals[ i ] );
                            // }
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, size_t nvals, const size_t *row_inds,
                        const ELEMENTTYPE *vals ) override {
                        for ( size_t i = 0; i < nvals; i++ ) {
                            set( row_inds[ i ], column, vals[ i ] );
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column,
                        const std::vector<ELEMENTTYPE>& vals ) override {
                        // _cols[ column ].zero();
                        for ( size_t i = 0; i < vals.size(); i++ ) {
                            // if ( !NCPA::math::is_zero( vals[ i ] ) ) {
                            set( i, column, vals[ i ] );
                            // }
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    // virtual ELEMENTTYPE& get( size_t row, size_t col )
                    // override {
                    //
                    //     return _rows[ row ].get( col );
                    // }

                    virtual const ELEMENTTYPE& get(
                        size_t row, size_t col ) const override {
                        return _rows[ row ].get( col );
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_row( size_t row ) const override {
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new sparse_vector<ELEMENTTYPE>( _rows[ row ] ) );
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_column( size_t column ) const override {
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new sparse_vector<ELEMENTTYPE>(
                                _cols[ column ] ) );
                    }

                    // virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    //     get_diagonal( int offset ) const override {
                    //
                    //     size_t absoffset = (size_t)( std::abs( offset ) );
                    //     this->check_size( absoffset, absoffset );
                    //     size_t ndiag = diagonal_size( offset );
                    //     std::vector<ELEMENTTYPE> diag( ndiag );
                    //     for ( size_t i = 0; i < ndiag; i++ ) {
                    //         if ( offset >= 0 ) {
                    //             // diagonal and upper triangle
                    //             diag[ i ] = _elements[ i ][ i + absoffset ];
                    //         } else {
                    //             diag[ i ] = _elements[ i + absoffset ][ i ];
                    //         }
                    //     }
                    //     return
                    //     std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                    //         new sparse_vector<ELEMENTTYPE>( diag ) );
                    // }

                    virtual abstract_matrix<ELEMENTTYPE>& as_array(
                        size_t& nrows, size_t& ncols,
                        ELEMENTTYPE **& vals ) override {
                        if ( nrows == 0 && ncols == 0 ) {
                            nrows = rows();
                            ncols = columns();
                            vals  = NCPA::arrays::zeros<ELEMENTTYPE>( nrows,
                                                                      ncols );
                        } else if ( nrows != rows() || ncols != columns() ) {
                            throw std::invalid_argument(
                                "Size mismatch between vector and target "
                                "array" );
                        }
                        ELEMENTTYPE zero = NCPA::math::zero<ELEMENTTYPE>();
                        NCPA::arrays::fill( vals, nrows, ncols, zero );
                        for ( size_t row = 0; row < rows(); row++ ) {
                            auto nonzero = _rows[ row ].nonzero();
                            for ( auto it = nonzero.begin();
                                  it != nonzero.end(); ++it ) {
                                vals[ row ][ it->first ] = it->second;
                            }
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& operator[](
                        size_t i ) override {
                        this->check_size( i, 0 );
                        return _rows[ i ];
                    }

                    // // this is a self-consistency check that should never
                    // fail.  If it does,
                    // // there is a coding error somewhere.
                    // void run_consistency_check() {
                    //     // check that all elements in rows "matrix" equal
                    //     the corresponding
                    //     // elements in columns "matrix"
                    //     for ( size_t row = 0; row < rows(); row++ ) {
                    //         if ( _rows[ row ].size() != _ncols ) {
                    //             _rows[ row ].resize( _ncols );
                    //         }
                    //         std::map<size_t, ELEMENTTYPE> nonzero
                    //             = _rows[ row ].nonzero();
                    //         for ( auto mit = nonzero.cbegin();
                    //               mit != nonzero.cend(); ++mit ) {
                    //             if ( !NCPA::math::equals(
                    //                      mit->second,
                    //                      _cols[ mit->first ].get( row ) ) )
                    //                      {
                    //                 std::ostringstream oss;
                    //                 oss << "Matrix consistency failure:
                    //                 rows["
                    //                     << row << "][" << mit->first
                    //                     << "] does not equal columns["
                    //                     << mit->first << "][" << row << "]";
                    //                 throw std::logic_error( oss.str() );
                    //             }
                    //         }
                    //     }

                    //     // and vice versa
                    //     for (size_t col = 0; col < columns(); col++) {
                    //         if (_cols[ col ].size() != _nrows) {
                    //             _cols[ col ].resize( _nrows );
                    //         }
                    //         std::map<size_t, ELEMENTTYPE> nonzero
                    //             = _cols[ col ].nonzero();
                    //         for ( auto mit = nonzero.cbegin();
                    //               mit != nonzero.cend(); ++mit ) {
                    //             if ( !NCPA::math::equals(
                    //                      mit->second,
                    //                      _rows[ mit->first ].get( col ) ) )
                    //                      {
                    //                 std::ostringstream oss;
                    //                 oss << "Matrix consistency failure:
                    //                 columns["
                    //                     << col << "][" << mit->first
                    //                     << "] does not equal rows["
                    //                     << mit->first << "][" << col << "]";
                    //                 throw std::logic_error( oss.str() );
                    //             }
                    //         }
                    //     }
                    //     _finalized = true;
                    // }

                protected:
                    void _clear_column_cache() { _cols.clear(); }

                    void _build_column_cache() {
                        if ( _cols.size() == 0 ) {
                            _cols.resize( _ncols );
                            for ( size_t col = 0; col < _ncols; col++ ) {
                                _cols[ col ]
                                    = sparse_vector<ELEMENTTYPE>( _nrows );
                            }
                            for ( size_t row = 0; row < _nrows; row++ ) {
                                _rows[ row ].qc();
                                auto nzmap = _rows[ row ].nonzero();
                                for ( auto it = nzmap.cbegin();
                                      it != nzmap.cend(); ++it ) {
                                    _cols[ it->first ].set( row, it->second );
                                }
                            }
                        }
                    }

                    void _rebuild_column_cache() {
                        _clear_column_cache();
                        _build_column_cache();
                    }

                private:
                    size_t _nrows, _ncols;
                    // std::vector<std::unique_ptr<sparse_vector<ELEMENTTYPE>>>
                    std::vector<sparse_vector<ELEMENTTYPE>> _rows, _cols;
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::details::sparse_matrix<T>& a,
                  NCPA::linear::details::sparse_matrix<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::details::abstract_matrix<T>&>( a ),
            static_cast<NCPA::linear::details::abstract_matrix<T>&>( b ) );
    swap( a._rows, b._rows );
    swap( a._cols, b._cols );
    swap( a._nrows, b._nrows );
    swap( a._ncols, b._ncols );
}
