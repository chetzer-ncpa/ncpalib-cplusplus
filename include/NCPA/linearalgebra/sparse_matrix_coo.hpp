#pragma once

// Currently uses COO storage.  At some point write a CSR version

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
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

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::sparse_matrix, ELEMENTTYPE );

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

                    sparse_matrix( size_t nrows, size_t ncols ) :
                        _nrows { nrows }, _ncols { ncols } {
                    }

                    sparse_matrix() : sparse_matrix<ELEMENTTYPE>( 0, 0 ) {}

                    sparse_matrix( size_t nrows, size_t ncols,
                                   const std::vector<size_t>& row_inds,
                                   const std::vector<size_t>& col_inds,
                                   const std::vector<ELEMENTTYPE>& values ) :
                        sparse_matrix<ELEMENTTYPE>( nrows, ncols ) {
                        _row_inds = row_inds;
                        _col_inds = col_inds;
                        _values   = values;
                    }

                    sparse_matrix( const sparse_matrix<ELEMENTTYPE>& other ) :
                        sparse_matrix<ELEMENTTYPE>() {
                        _nrows    = other._nrows;
                        _ncols    = other._ncols;
                        _values   = other._values;
                        _row_inds = other._row_inds;
                        _col_inds = other._col_inds;
                    }

                    sparse_matrix(
                        const abstract_matrix<ELEMENTTYPE>& other ) :
                        sparse_matrix<ELEMENTTYPE>( other.rows(),
                                                    other.columns() ) {
                        for ( auto i = 0; i < other.rows(); i++ ) {
                            auto row  = other.get_row( i );
                            auto inds = row->nonzero_indices();
                            for ( auto it = inds.begin(); it != inds.end();
                                  ++it ) {
                                _row_inds.push_back( i );
                                _col_inds.push_back( *it );
                                _values.push_back( row->get( *it ) );
                            }
                        }
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
                        return *this;
                    }

                    sparse_matrix<ELEMENTTYPE>& operator=(
                        abstract_matrix<ELEMENTTYPE> other ) {
                        sparse_matrix<ELEMENTTYPE> copy( other );
                        swap( *this, copy );
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
                        if ( !_is_sparse( other ) ) {
                            return abstract_matrix<ELEMENTTYPE>::equals(
                                other );
                        }
                        const sparse_matrix<ELEMENTTYPE> *other_as_sparse
                            = _as_sparse( other );
                        for ( size_t i = 0; i < _values.size(); i++ ) {
                            if ( _row_inds[ i ]
                                     != other_as_sparse->_row_inds[ i ]
                                 || _col_inds[ i ]
                                        != other_as_sparse->_col_inds[ i ]
                                 || _values[ i ]
                                        != other_as_sparse->_values[ i ] ) {
                                return false;
                            }
                        }
                        return true;
                    }

                    virtual bool is_diagonal() const override {
                        if ( this->is_empty() ) {
                            return true;
                        }
                        for ( size_t i = 0; i < _values.size(); i++ ) {
                            if ( _row_inds[ i ] != _col_inds[ i ] ) {
                                return false;
                            }
                        }
                        return true;
                    }

                    virtual size_t rows() const override { return _nrows; }

                    virtual size_t columns() const override { return _ncols; }

                    virtual abstract_matrix<ELEMENTTYPE>& clear() override {
                        _values.clear();
                        _row_inds.clear();
                        _col_inds.clear();
                        _nrows = 0;
                        _ncols = 0;
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& zero(
                        size_t row, size_t col ) override {
                        assert_point_is_in_size( row, col );
                        size_t row_start = _col_inds[ col ];
                        size_t row_end   = _col_inds[ col + 1 ];
                        for ( size_t i = row_start; i < row_end; i++ ) {
                            if ( _row_inds[ i ] == row ) {
                                _row_inds.erase( _row_inds.begin() + i );
                                _values.erase( _values.begin() + i );
                                for ( size_t j = i; j < _col_inds.size();
                                      j++ ) {
                                    _col_inds[ j ]--;
                                }
                                return this->as_base_class();
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& zero() override {
                        _values.clear();
                        _row_inds.clear();
                        _col_inds.clear();
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& resize(
                        size_t newrows, size_t newcols ) override {
                        while ( _ncols > newcols ) {
                            _remove_column();
                        }
                        while ( _ncols < newcols ) {
                            _add_column();
                        }
                        while ( _nrows > newrows ) {
                            _remove_row();
                        }
                        while ( _nrows < newrows ) {
                            _add_row();
                        }
                        return this->as_base_class();
                    }

                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        multiply( const abstract_matrix<ELEMENTTYPE>& b )
                            const override {
                        if ( auto *derived = dynamic_cast<
                                 const sparse_matrix<ELEMENTTYPE> *>( &b ) ) {
                            return _sparse_multiply(
                                dynamic_cast<
                                    const sparse_matrix<ELEMENTTYPE>&>( b ) );
                        } else {
                            return abstract_matrix<ELEMENTTYPE>::multiply( b );
                        }
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& transpose()
                        override {
                        std::swap( _nrows, _ncols );
                        std::swap( _row_inds, _col_inds );
                        return this->reorder();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& reorder() {
                        auto perm = NCPA::arrays::sort_permutation_increasing(
                            _row_inds, _col_inds );
                        NCPA::arrays::apply_permutation_in_place( _row_inds,
                                                                  perm );
                        NCPA::arrays::apply_permutation_in_place( _col_inds,
                                                                  perm );
                        NCPA::arrays::apply_permutation_in_place( _values,
                                                                  perm );
                        return this->as_base_class();
                    }

                    virtual void assert_point_is_in_size( size_t row,
                                                          size_t col ) {
                        if ( row >= _nrows || col >= _ncols ) {
                            std::ostringstream oss;
                            oss << "Point [" << row << "," << col
                                << "] is out of range.";
                            throw std::range_error( oss.str() );
                        }
                    }

                    virtual void assert_same_size(
                        const abstract_matrix<ELEMENTTYPE>& b ) const {
                        if ( _nrows != b.rows() || _ncols != b.columns() ) {
                            throw std::range_error(
                                "Matrices are not the same size!" );
                        }
                    }

                    virtual size_t count() const {
                        if ( _row_inds.size() != _col_inds.size()
                             || _col_inds.size() != _values.size() ) {
                            throw std::logic_error(
                                "Matrix is not self-consistent!" );
                        }
                        return _values.size();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& add(
                        const abstract_matrix<ELEMENTTYPE>& b,
                        ELEMENTTYPE modifier = 1.0 ) override {
                        if ( _is_sparse( b ) ) {
                            return _sparse_add( b, modifier );
                        } else {
                            return abstract_matrix<ELEMENTTYPE>::add(
                                b, modifier );
                        }
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& _sparse_add(
                        const sparse_matrix<ELEMENTTYPE>& b,
                        ELEMENTTYPE modifier = 1.0 ) {
                        assert_same_size( b );
                        size_t a_ind = 0, b_ind = 0;
                        std::vector<size_t> sum_row_inds, sum_col_inds;
                        std::vector<ELEMENTTYPE> sum_values;
                        while ( a_ind < this->count() || b_ind < b.count() ) {
                            // if b's row and col is smaller
                            if ( _row_inds[ a_ind ] > b._row_inds[ b_ind ]
                                 || ( _row_inds[ a_ind ]
                                          == b._row_inds[ b_ind ]
                                      && _col_inds[ a_ind ]
                                             > b._col_inds[ b_ind ] ) ) {
                                // insert smaller value into result
                                sum_row_inds.push_back( b._row_inds[ b_ind ] );
                                sum_col_inds.push_back( b._col_inds[ b_ind ] );
                                sum_values.push_back( b._values[ b_ind ]
                                                      * modifier );
                                b_ind++;

                            } else if ( _row_inds[ a_ind ]
                                            < b._row_inds[ b_ind ]
                                        || ( _row_inds[ a_ind ]
                                                 == b._row_inds[ b_ind ]
                                             && _col_inds[ a_ind ]
                                                    < b._col_inds
                                                          [ b_ind ] ) ) {
                                // insert smaller value into result
                                sum_row_inds.push_back(
                                    this->_row_inds[ a_ind ] );
                                sum_col_inds.push_back(
                                    this->_col_inds[ a_ind ] );
                                sum_values.push_back( this->_values[ a_ind ] );
                                a_ind++;
                            } else {
                                // add the values as row and col is same
                                ELEMENTTYPE addedval
                                    = this->_values[ a_ind ]
                                    + b._values[ b_ind ] * modifier;

                                if ( !NCPA::math::is_zero( addedval ) ) {
                                    sum_row_inds.push_back(
                                        this->_row_inds[ a_ind ] );
                                    sum_col_inds.push_back(
                                        this->_col_inds[ a_ind ] );
                                    sum_values.push_back( addedval );
                                }
                                a_ind++;
                                b_ind++;
                            }
                        }

                        while ( a_ind < this->count() ) {
                            sum_row_inds.push_back( this->_row_inds[ a_ind ] );
                            sum_col_inds.push_back( this->_col_inds[ a_ind ] );
                            sum_values.push_back( this->_values[ a_ind++ ] );
                        }

                        while ( b_ind < b.count() ) {
                            sum_row_inds.push_back( b._row_inds[ b_ind ] );
                            sum_col_inds.push_back( b._col_inds[ b_ind ] );
                            sum_values.push_back( b._values[ b_ind++ ]
                                                  * modifier );
                        }
                        _row_inds = sum_row_inds;
                        _col_inds = sum_col_inds;
                        _values   = sum_values;
                        return this->as_base_class();
                    }

                    virtual bool is_tridiagonal() const override {
                        for ( size_t i = 0; i < this->count(); i++ ) {
                            if ( std::abs( (int)_row_inds[ i ]
                                           - (int)_col_inds[ i ] )
                                 > 1 ) {
                                return false;
                            }
                        }
                        return true;
                    }

                    virtual bool is_upper_triangular() const override {
                        for ( size_t i = 0; i < this->count(); i++ ) {
                            if ( _col_inds[ i ] < _row_inds[ i ] ) {
                                return false;
                            }
                        }
                        return true;
                    }

                    virtual bool is_lower_triangular() const override {
                        for ( size_t i = 0; i < this->count(); i++ ) {
                            if ( _col_inds[ i ] > _row_inds[ i ] ) {
                                return false;
                            }
                        }
                        return true;
                    }

                    std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        _sparse_multiply(
                            const sparse_matrix<ELEMENTTYPE>& b ) const {
                        this->check_size_for_mult( b );
                        std::cout << "Copying RHS matrix" << std::endl;
                        sparse_matrix<ELEMENTTYPE> bt = b;
                        // for this kind of matrix, transpose is cheap
                        std::cout << "Transposing RHS matrix" << std::endl;
                        bt.transpose();
                        // const sparse_matrix<ELEMENTTYPE> *mat1, *mat2;

                        std::vector<ELEMENTTYPE> new_values;
                        std::vector<size_t> new_row_inds, new_col_inds;
                        size_t a_ind = 0;
                        size_t b_ind = 0;
                        while ( a_ind < this->count() ) {
                            // rows in the LHS matrix will be used to calculate
                            // rows in the resultant matrix.  Since we store in
                            // row-first order, that's the order we want to go
                            // in
                            size_t a_row = this->_row_inds[ a_ind ];
                            std::cout << "a_row = " << a_row << " (a_ind = " << a_ind << ")" << std::endl;

                            // we want the b values where b._col_ind == row
                            while ( b_ind < b.count() ) {
                                std::cout << "b._row_inds[ b_ind ] = " << b._row_inds[ b_ind ] << " (b_ind = " << b_ind << ")" << std::endl;
                                while ( b._row_inds[ b_ind ] < a_row ) {
                                    b_ind++;
                                    std::cout << "b._row_inds[ b_ind ] = " << b._row_inds[ b_ind ] << " (b_ind = " << b_ind << ")" << std::endl;
                                }
                                while ( b._row_inds[ b_ind ] == a_row ) {
                                    if ( new_row_inds.back() == a_row
                                         && new_col_inds.back()
                                                == b._row_inds[ b_ind ] ) {
                                        // cell already populated, add to it
                                        new_values.back()
                                            += this->_values[ a_ind ]
                                             * b._values[ b_ind ];
                                    } else {
                                        new_row_inds.push_back( a_row );
                                        new_col_inds.push_back(
                                            b._row_inds[ b_ind ] );
                                        new_values.push_back(
                                            this->_values[ a_ind ]
                                            * b._values[ b_ind ] );
                                    }
                                    b_ind++;
                                }
                            }
                            a_ind++;
                            if ( this->_row_inds[ a_ind ]
                                 > this->_row_inds[ a_ind - 1 ] ) {
                                // new row, start b over
                                b_ind = 0;
                            }
                        }
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new sparse_matrix<ELEMENTTYPE>(
                                _nrows, b._ncols, new_row_inds, new_col_inds,
                                new_values ) );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& swap_rows(
                        size_t ind1, size_t ind2 ) override {
                        if ( ind1 == ind2 ) {
                            return this->as_base_class();
                        }
                        size_t i1 = std::min( ind1, ind2 );
                        size_t i2 = std::max( ind1, ind2 );
                        size_t i  = 0;
                        while ( _row_inds[ i ] <= i2 ) {
                            if ( _row_inds[ i ] == i1 ) {
                                _row_inds[ i ] = i2;
                            } else if ( _row_inds[ i ] == i2 ) {
                                _row_inds[ i ] = i1;
                            }
                            i++;
                        }
                        this->reorder();
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& swap_columns(
                        size_t ind1, size_t ind2 ) override {
                        if ( ind1 == ind2 ) {
                            return this->as_base_class();
                        }
                        for ( size_t i = 0; i < count(); i++ ) {
                            if ( _col_inds[ i ] == ind1 ) {
                                _col_inds[ i ] = ind2;
                            } else if ( _col_inds[ i ] == ind2 ) {
                                _col_inds[ i ] = ind1;
                            }
                        }
                        this->reorder();
                        return this->as_base_class();
                    }

                    // template<typename ANYTYPE,
                    //          ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                    //                        ELEMENTTYPE )>
                    abstract_matrix<ELEMENTTYPE>& add( ELEMENTTYPE b ) override {
                        std::vector<size_t> new_row_inds,
                                            new_col_inds,
                                            inds = NCPA::arrays::index_vector<size_t>( _ncols );
                        std::vector<ELEMENTTYPE> new_vals( _nrows * _ncols, b );
                        for (size_t i = 0; i < _nrows; i++) {
                            new_row_inds.insert( new_row_inds.end(), _ncols, i );
                            new_col_inds.insert( new_col_inds.end(), inds.begin(), inds.end() );
                        }
                        for (size_t i = 0; i < count(); i++) {
                            new_vals[ _row_inds[ i ] * _ncols + _col_inds[ i ] ] += _values[ i ];
                        }
                        _row_inds = new_row_inds;
                        _col_inds = new_col_inds;
                        _values = new_vals;
                        return this->as_base_class();
                    }

                    // template<typename ANYTYPE,
                    //          ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                    //                        ELEMENTTYPE )>
                    abstract_matrix<ELEMENTTYPE>& scale( ELEMENTTYPE val ) override {
                        if ( NCPA::math::is_zero( val ) ) {
                            this->zero();
                        } else {
                            for ( auto it = _values.begin();
                                  it != _values.end(); ++it ) {
                                *it *= val;
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& scale(
                        const abstract_matrix<ELEMENTTYPE>& b ) override {
                        this->assert_same_size( b );
                        for ( size_t i = 0; i < this->count(); i++ ) {
                            ELEMENTTYPE bval
                                = b.get( _row_inds[ i ], _col_inds[ i ] );
                            if ( NCPA::math::is_zero( bval ) ) {
                                this->zero( _row_inds[ i ], _col_inds[ i ] );
                            } else {
                                _values[ i ] *= bval;
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& identity(
                        size_t nrows, size_t ncols ) override {
                        ELEMENTTYPE one = NCPA::math::one<ELEMENTTYPE>();
                        this->clear().resize( nrows, ncols );
                        size_t ndiag = std::min( nrows, ncols );
                        for ( size_t i = 0; i < this->count(); i++ ) {
                            _row_inds.push_back( i );
                            _col_inds.push_back( i );
                            _values.push_back( one );
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        size_t row, size_t col, ELEMENTTYPE val ) override {
                        assert_point_is_in_size( row, col );
                        size_t i = 0;
                        // @todo speed up by searching from both ends?
                        while ( i < this->count() && _row_inds[ i ] < row
                                && _col_inds[ i ] < col ) {
                            i++;
                        }
                        if ( i == this->count() ) {
                            _row_inds.push_back( row );
                            _col_inds.push_back( col );
                            _values.push_back( val );
                        } else if ( row == _row_inds[ i ]
                                    && col == _col_inds[ i ] ) {
                            // have value already, replace
                            _values[ i ] = val;
                        } else {
                            // have to insert
                            _row_inds.insert( _row_inds.begin() + i, row );
                            _col_inds.insert( _col_inds.begin() + i, col );
                            _values.insert( _values.begin() + i, val );
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        ELEMENTTYPE val ) override {
                        this->zero();
                        if ( val != _zero ) {
                            for ( size_t i = 0; i < this->rows(); i++ ) {
                                for ( size_t j = 0; i < this->columns();
                                      j++ ) {
                                    _row_inds.push_back( i );
                                    _col_inds.push_back( j );
                                    _values.push_back( val );
                                }
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, size_t nvals, const size_t *column_inds,
                        const ELEMENTTYPE *vals ) override {
                        // std::cout
                        //     << "Called array-based set_row in sparse_matrix"
                        //     << std::endl;
                        assert_point_is_in_size( row, 0 );
                        size_t i = 0;
                        while ( i < this->count() && _row_inds[ i ] < row ) {
                            i++;
                        }
                        if ( i == this->count() ) {
                            _row_inds.insert( _row_inds.end(), nvals, row );
                            _col_inds.insert( _col_inds.end(), column_inds,
                                              column_inds + nvals );
                            _values.insert( _values.end(), vals,
                                            vals + nvals );
                        } else if ( _row_inds[ i ] > row ) {
                            // nothing in this row yet, insert
                            for ( size_t j = 0; j < nvals; j++ ) {
                                _row_inds.push_back( row );
                                _col_inds.push_back( column_inds[ j ] );
                                _values.push_back( vals[ j ] );
                            }
                        } else {
                            size_t j = 0, k = 0;
                            while ( _row_inds[ i ] == row && j < nvals ) {
                                if ( _col_inds[ i ] < column_inds[ j ] ) {
                                    i++;
                                } else if ( _col_inds[ i ]
                                            == column_inds[ j ] ) {
                                    _values[ i ] = vals[ j ];
                                    i++;
                                    j++;
                                } else {
                                    _row_inds.insert( _row_inds.begin() + i,
                                                      row );
                                    _col_inds.insert( _col_inds.begin() + i,
                                                      column_inds[ j ] );
                                    _values.insert( _values.begin() + i,
                                                    vals[ j ] );
                                    i++;
                                    j++;
                                }
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row,
                        const std::vector<ELEMENTTYPE>& vals ) override {
                        // std::cout
                        //     << "Called single-vector set_row in sparse_matrix"
                        //     << std::endl;
                        assert_point_is_in_size( row, 0 );
                        size_t i = 0;
                        while ( i < this->count() && _row_inds[ i ] < row ) {
                            i++;
                        }
                        if ( i == this->count() ) {
                            _row_inds.insert( _row_inds.end(), _ncols, row );
                            for ( size_t j = 0; j < this->columns(); j++ ) {
                                _col_inds.push_back( j );
                            }
                            _values.insert( _values.end(), vals.begin(),
                                            vals.end() );
                        } else {
                            while ( i << this->count()
                                    && _row_inds[ i ] == row ) {
                                _row_inds.erase( _row_inds.begin() + i );
                                _col_inds.erase( _col_inds.begin() + i );
                                _values.erase( _values.begin() + i );
                            }
                            _row_inds.insert( _row_inds.begin() + i, _ncols,
                                              row );
                            for ( size_t j = 0; j < this->columns(); j++ ) {
                                _col_inds.insert( _col_inds.begin() + i + j,
                                                  j );
                            }
                            _values.insert( _values.begin() + i, vals.cbegin(),
                                            vals.cend() );
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, size_t nvals, const size_t *row_inds,
                        const ELEMENTTYPE *vals ) override {
                        assert_point_is_in_size( 0, column );
                        size_t mat_ind = 0;
                        for ( size_t val_ind = 0; val_ind < nvals;
                              val_ind++ ) {
                            size_t row = row_inds[ val_ind ];
                            while ( mat_ind < count()
                                    && _row_inds[ mat_ind ] < row ) {
                                mat_ind++;
                            }

                            if ( mat_ind == count() ) {
                                // hit the end
                                while ( val_ind < nvals ) {
                                    _row_inds.push_back( row_inds[ val_ind ] );
                                    _col_inds.push_back( column );
                                    _values.push_back( vals[ val_ind++ ] );
                                }
                                return this->as_base_class();
                            } else if ( _row_inds[ mat_ind ] > row ) {
                                // overshot the row
                                _row_inds.insert( _row_inds.begin() + mat_ind,
                                                  row_inds[ val_ind ] );
                                _col_inds.insert( _col_inds.begin() + mat_ind,
                                                  column );
                                _values.insert( _values.begin() + mat_ind,
                                                vals[ val_ind ] );
                                mat_ind++;
                            } else {
                                while ( _row_inds[ mat_ind ] == row
                                        && _col_inds[ mat_ind ] < column ) {
                                    mat_ind++;
                                }
                                if ( _row_inds[ mat_ind ] == row
                                     && _col_inds[ mat_ind ] == column ) {
                                    _values[ mat_ind ] = vals[ val_ind ];
                                } else {
                                    _row_inds.insert( _row_inds.begin()
                                                          + mat_ind,
                                                      row_inds[ val_ind ] );
                                    _col_inds.insert(
                                        _col_inds.begin() + mat_ind, column );
                                    _values.insert( _values.begin() + mat_ind,
                                                    vals[ val_ind ] );
                                }
                                mat_ind++;
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column,
                        const std::vector<ELEMENTTYPE>& vals ) override {
                        size_t i = 0;
                        for ( size_t row = 0; row < vals.size(); row++ ) {
                            while ( i < count() && _row_inds[ i ] < row ) {
                                i++;
                            }
                            while ( i < count() && _row_inds[ i ] == row
                                    && _col_inds[ i ] < column ) {
                                i++;
                            }
                            if ( i == count() ) {
                                while ( row < vals.size() ) {
                                    _row_inds.push_back( row );
                                    _col_inds.push_back( column );
                                    _values.push_back( vals[ row++ ] );
                                }
                                return this->as_base_class();
                            } else if ( _row_inds[ i ] == row
                                        && _col_inds[ i ] == column ) {
                                _values[ i ] = vals[ row ];
                            } else {
                                _row_inds.insert( _row_inds.begin() + i, row );
                                _col_inds.insert( _col_inds.begin() + i,
                                                  column );
                                _values.insert( _values.begin() + i,
                                                vals[ row ] );
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual const ELEMENTTYPE& get(
                        size_t row, size_t col ) const override {
                        auto it = std::find( _row_inds.begin(),
                                             _row_inds.end(), row );
                        if ( it == _row_inds.end() ) {
                            return _zero;
                        }
                        size_t i = std::distance( _row_inds.begin(), it );
                        while ( i < count() && _row_inds[ i ] == row
                                && _col_inds[ i ] < col ) {
                            i++;
                        }
                        if ( i == count() ) {
                            return _zero;
                        }
                        if ( _row_inds[ i ] == row && _col_inds[ i ] == col ) {
                            return _values[ i ];
                        } else {
                            return _zero;
                        }
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_row( size_t row ) const override {
                        auto it = std::find( _row_inds.cbegin(),
                                             _row_inds.cend(), row );
                        if ( it == _row_inds.cend() ) {
                            return std::unique_ptr<
                                abstract_vector<ELEMENTTYPE>>(
                                new sparse_vector<ELEMENTTYPE>( _ncols ) );
                        }
                        std::vector<ELEMENTTYPE> vec( _ncols );
                        size_t i = std::distance( _row_inds.begin(), it );
                        while ( i < count() && _row_inds[ i ] == row ) {
                            vec[ _col_inds[ i ] ] = _values[ i ];
                            i++;
                        }
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new sparse_vector<ELEMENTTYPE>( vec ) );
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_column( size_t column ) const override {
                        sparse_vector<ELEMENTTYPE> vec( _nrows );
                        auto it = std::find( _col_inds.cbegin(),
                                             _col_inds.cend(), column );
                        while ( it != _col_inds.end() ) {
                            size_t i = std::distance( _col_inds.cbegin(), it );
                            vec.set( _row_inds[ i ], _values[ i ] );
                            it = std::find( it+1, _col_inds.cend(), column );
                        }
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new sparse_vector<ELEMENTTYPE>( vec ) );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& as_array(
                        size_t& nrows, size_t& ncols,
                        ELEMENTTYPE **& vals ) override {
                        NCPA::arrays::fill( vals, nrows, ncols, _zero );
                        for ( size_t i = 0; i < count(); i++ ) {
                            vals[ _row_inds[ i ] ][ _col_inds[ i ] ]
                                = _values[ i ];
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& as_base_class() override {
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>( this );
                    }

                private:
                    size_t _nrows, _ncols;

                    std::vector<ELEMENTTYPE> _values;
                    std::vector<size_t> _row_inds, _col_inds;
                    const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();

                    void _add_column() { _ncols++; }

                    void _add_row() { _nrows++; }

                    void _remove_column() {
                        if ( _ncols == 0 ) {
                            throw std::range_error(
                                "Can't remove columns from an empty matrix" );
                        }
                        _ncols--;
                        size_t i = 0;
                        while ( i < _values.size() ) {
                            if ( _col_inds[ i ] >= _ncols ) {
                                _values.erase( _values.begin() + i );
                                _row_inds.erase( _row_inds.begin() + i );
                                _col_inds.erase( _col_inds.begin() + i );
                            } else {
                                ++i;
                            }
                        }
                    }

                    void _remove_row() {
                        if ( _nrows == 0 ) {
                            throw std::range_error(
                                "Can't remove rows from an empty matrix" );
                        }
                        _nrows--;
                        auto it
                            = std::find_if( _row_inds.begin(), _row_inds.end(),
                                            [ this ]( size_t ind ) {
                                                return ( ind >= this->_nrows );
                                            } );
                        size_t offset = std::distance( _row_inds.begin(), it );
                        _row_inds.erase( _row_inds.begin() + offset,
                                         _row_inds.end() );
                        _col_inds.erase( _col_inds.begin() + offset,
                                         _col_inds.end() );
                        _values.erase( _values.begin() + offset,
                                       _values.end() );
                    }

                    bool _is_sparse(
                        const abstract_matrix<ELEMENTTYPE>& b ) const {
                        if ( auto *derived = dynamic_cast<
                                 const sparse_matrix<ELEMENTTYPE> *>( &b ) ) {
                            return true;
                        } else {
                            return false;
                        }
                    }

                    const sparse_matrix<ELEMENTTYPE> *_as_sparse(
                        const abstract_matrix<ELEMENTTYPE>& b ) const {
                        return dynamic_cast<
                            const sparse_matrix<ELEMENTTYPE> *>( &b );
                    }
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
    swap( a._nrows, b._nrows );
    swap( a._ncols, b._ncols );
    swap( a._row_inds, b._row_inds );
    swap( a._col_inds, b._col_inds );
    swap( a._values, b._values );
}
