#pragma once

#include "NCPA/exceptions.hpp"
#include "NCPA/linearalgebra/builders.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/functions.hpp"
#include "NCPA/linearalgebra/Matrix.hpp"
#include "NCPA/linearalgebra/Vector.hpp"
#include "NCPA/logging.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

template<typename T>
static void swap( NCPA::linear::BlockMatrix<T>& a,
                  NCPA::linear::BlockMatrix<T>& b ) noexcept;

namespace NCPA {
    namespace linear {

        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class BlockMatrix<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public Matrix<ELEMENTTYPE> {
                using Matrix<ELEMENTTYPE>::get;
                using Matrix<ELEMENTTYPE>::inverse;
                using Matrix<ELEMENTTYPE>::invert;
                using Matrix<ELEMENTTYPE>::is_column_matrix;
                using Matrix<ELEMENTTYPE>::is_row_matrix;
                using Matrix<ELEMENTTYPE>::lu;
                using Matrix<ELEMENTTYPE>::lu_decompose;
                using Matrix<ELEMENTTYPE>::set_row;
                using Matrix<ELEMENTTYPE>::set_column;
                using Matrix<ELEMENTTYPE>::set_diagonal;

            public:
                BlockMatrix() :
                    BlockMatrix<ELEMENTTYPE>( matrix_t::INVALID ) {}

                BlockMatrix( matrix_t mattype ) :
                    _rows_of_blocks { 0 },
                    _cols_of_blocks { 0 },
                    _rows_per_block { 0 },
                    _cols_per_block { 0 },
                    _blocktype { mattype } {}

                BlockMatrix( const BlockMatrix<ELEMENTTYPE>& other ) :
                    BlockMatrix<ELEMENTTYPE>( other._blocktype ) {
                    this->copy( other );
                    // for (auto it = other._elements.cbegin();
                    //      it != other._elements.cend(); ++it) {
                    //     _elements.push_back( *it );
                    // }
                    // _rows_of_blocks = other._rows_of_blocks;
                    // _cols_of_blocks = other._cols_of_blocks;
                    // _rows_per_block = other._rows_per_block;
                    // _cols_per_block = other._cols_per_block;
                }

                BlockMatrix( BlockMatrix<ELEMENTTYPE>&& source ) noexcept :
                    BlockMatrix<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                virtual ~BlockMatrix() {}

                friend void ::swap<ELEMENTTYPE>(
                    BlockMatrix<ELEMENTTYPE>& a,
                    BlockMatrix<ELEMENTTYPE>& b ) noexcept;

                /**
                 * Assignment operator.
                 * @param other The vector to assign to this.
                 */
                BlockMatrix<ELEMENTTYPE>& operator=(
                    BlockMatrix<ELEMENTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual std::unique_ptr<Matrix<ELEMENTTYPE>> clone() const override {
                    return std::unique_ptr<Matrix<ELEMENTTYPE>>( new BlockMatrix<ELEMENTTYPE>( *this ) );
                }

                virtual BlockMatrix<ELEMENTTYPE>& add(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    this->check();
                    other.check();
                    _checksizes( other );
                    for (size_t r = 0; r < this->block_rows(); ++r) {
                        for (size_t c = 0; c < this->block_columns(); ++c) {
                            // _elements[ _blockcoord2blockindex( r, c ) ]
                            this->get_block( r, c ) += other.get_block( r, c );
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& add(
                    const Matrix<ELEMENTTYPE>& other ) override {
                    if (auto derived
                        = dynamic_cast<const BlockMatrix<ELEMENTTYPE> *>(
                            &other )) {
                        return this->add( *derived );
                    }
                    this->check();
                    other.check();
                    check_same_size( other );
                    ELEMENTTYPE zero = NCPA::constants::zero<ELEMENTTYPE>();
                    for (size_t r = 0; r < this->rows(); ++r) {
                        for (size_t c = 0; c < this->columns(); ++c) {
                            ELEMENTTYPE otherval = other.get( r, c );
                            if (otherval != zero) {
                                this->set( r, c,
                                           this->get( r, c ) + otherval );
                            }
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& add(
                    ELEMENTTYPE other ) override {
                    if (other != NCPA::constants::zero<ELEMENTTYPE>()) {
                        for (auto it = _elements.begin();
                             it != _elements.end(); ++it) {
                            it->add( other );
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& as_array(
                    size_t& nrows, size_t& ncols,
                    ELEMENTTYPE **& vals ) override {
                    if (this->is_empty()) {
                        nrows = 0;
                        ncols = 0;
                        vals  = nullptr;
                    } else {
                        nrows = this->rows();
                        ncols = this->columns();
                        vals
                            = NCPA::arrays::zeros<ELEMENTTYPE>( nrows, ncols );
                        for (size_t r = 0; r < nrows; ++r) {
                            for (size_t c = 0; c < ncols; ++c) {
                                vals[ r ][ c ] = this->get( r, c );
                            }
                        }
                    }
                    return *this;
                }

                virtual size_t bandwidth() const override {
                    return lower_bandwidth() + upper_bandwidth() + 1;
                }

                virtual size_t block_columns() const {
                    return _cols_of_blocks;
                }

                virtual size_t block_rows() const { return _rows_of_blocks; }

                virtual BlockMatrix& block_type( matrix_t newtype ) {
                    if (newtype != _blocktype) {
                        _blocktype = newtype;
                        std::vector<Matrix<ELEMENTTYPE>> _newelements;
                        _newelements.reserve( _elements.size() );
                        for (auto it = _elements.cbegin();
                             it != _elements.cend(); ++it) {
                            _newelements.push_back(
                                MatrixFactory<ELEMENTTYPE>::build(
                                    _blocktype ) );
                            if (!it->is_zero()) {
                                _newelements.back().copy( *it );
                            }
                        }
                        std::swap( _elements, _newelements );
                    }
                    return *this;
                }

                virtual matrix_t block_type() const { return _blocktype; }

                virtual explicit operator bool() const override {
                    return ( _blocktype != matrix_t::INVALID );
                }

                virtual void check() const override { this->_check_this(); }

                virtual void check_same_size(
                    const Matrix<ELEMENTTYPE>& other ) const override {
                    if (auto *derived
                        = dynamic_cast<const BlockMatrix<ELEMENTTYPE> *>(
                            &other )) {
                        _checksizes( *derived );
                    } else {
                        Matrix<ELEMENTTYPE>::check_same_size( other );
                    }
                }

                virtual BlockMatrix<ELEMENTTYPE>& clear() override {
                    _elements.clear();
                    _rows_of_blocks = 0;
                    _cols_of_blocks = 0;
                    _rows_per_block = 0;
                    _cols_per_block = 0;
                    return *this;
                }

                virtual size_t columns() const override {
                    return _cols_of_blocks * _cols_per_block;
                }

                virtual size_t columns_per_block() const {
                    return _cols_per_block;
                }

                virtual BlockMatrix<ELEMENTTYPE>& copy(
                    const Matrix<ELEMENTTYPE>& other ) override {
                    if (auto *bptr
                        = dynamic_cast<const BlockMatrix<ELEMENTTYPE> *>(
                            &other )) {
                        _blocktype = bptr->_blocktype;
                        _elements.clear();
                        for (auto it = bptr->_elements.cbegin();
                             it != bptr->_elements.cend(); ++it) {
                            _elements.push_back( *it );
                        }
                        _rows_of_blocks = bptr->_rows_of_blocks;
                        _cols_of_blocks = bptr->_cols_of_blocks;
                        _rows_per_block = bptr->_rows_per_block;
                        _cols_per_block = bptr->_cols_per_block;
                    } else {
                        if (this->rows() != other.rows()
                            || this->columns() != other.columns()) {
                            throw std::logic_error(
                                "Cannot copy a non-block "
                                "matrix into a block matrix if the dimensions "
                                "do not agree" );
                        }
                        this->zero();
                        for (size_t i = 0; i < this->rows(); ++i) {
                            auto thisrow = this->get_row( i );
                            auto nzinds  = thisrow->nonzero_indices();
                            for (auto it = nzinds.begin(); it != nzinds.end();
                                 ++it) {
                                this->set( i, *it, thisrow->get( *it ) );
                            }
                        }
                    }
                    return *this;
                }

                virtual std::vector<int> diagonals() const override {
                    if (this->is_empty() || this->is_zero()) {
                        return std::vector<int>();
                    }
                    if (this->is_diagonal()) {
                        return std::vector<int> { 0 };
                    }
                    if (this->is_tridiagonal()) {
                        return std::vector<int> { -1, 0, 1 };
                    }

                    std::vector<int> diags;
                    bool isdiag = this->is_block_diagonal();
                    bool istri  = this->is_block_tridiagonal();
                    for (int br = 0; br < (int)this->block_rows(); ++br) {
                        int mincol, maxcol;
                        if (isdiag) {
                            mincol = br;
                            maxcol = br;
                        } else if (istri) {
                            mincol = std::max( 0, br - 1 );
                            maxcol = std::min( (int)this->block_rows() - 1,
                                               br + 1 );
                        } else {
                            mincol = 0;
                            maxcol = (int)this->block_columns() - 1;
                        }
                        for (int bc = mincol; bc <= maxcol; ++bc) {
                            const Matrix<ELEMENTTYPE> *element
                                = &this->get_block( br, bc );
                            int bdiag = bc - br;
                            if (!element->is_zero()) {
                                std::vector<int> blockdiags
                                    = element->diagonals();
                                for (auto bit = blockdiags.cbegin();
                                     bit != blockdiags.cend(); ++bit) {
                                    diags.push_back( bdiag * _cols_per_block + *bit );
                                    // if (bdiag < 0) {
                                    //     diags.push_back( bdiag * _cols_per_block + *bit );
                                    // } else if (bdiag > 0) {
                                    //     diags.push_back( br * _rows_per_block 
                                    //                      - *bit );
                                    // }
                                }
                            }
                        }
                    }
                    std::sort( diags.begin(), diags.end() );
                    auto last = std::unique( diags.begin(), diags.end() );
                    diags.erase( last, diags.end() );
                    return diags;
                }

                virtual bool equals(
                    const BlockMatrix<ELEMENTTYPE>& other ) const {
                    // compare structure
                    if (( (bool)*this ) != ( (bool)other )) {
                        return false;
                    }
                    if (!*this && !other) {
                        return true;
                    }
                    if (!( this->rows() == other.rows()
                           && this->columns() == other.columns()
                           && this->block_rows() == other.block_rows()
                           && this->block_columns()
                                  == other.block_columns() )) {
                        return false;
                    }

                    // compare contents
                    auto it1 = _elements.cbegin();
                    auto it2 = other._elements.cbegin();
                    for (; it1 != _elements.cend()
                           && it2 != other._elements.cend();
                         ++it1, ++it2) {
                        if (*it1 != *it2) {
                            return false;
                        }
                    }
                    return true;
                }

                virtual bool equals(
                    const Matrix<ELEMENTTYPE>& other ) const override {
                    if (auto derived
                        = dynamic_cast<const BlockMatrix<ELEMENTTYPE> *>(
                            &other )) {
                        return this->equals( *derived );
                    }
                    // compare structure
                    if (( (bool)*this ) != ( (bool)other )) {
                        return false;
                    }
                    if (!*this && !other) {
                        return true;
                    }
                    if (!( this->rows() == other.rows()
                           && this->columns() == other.columns() )) {
                        return false;
                    }

                    // compare contents
                    for (size_t r = 0; r < this->rows(); ++r) {
                        if (*this->get_row( r ) != *other.get_row( r )) {
                            return false;
                        }
                    }
                    return true;
                }

                // get the element by its overall coordinates
                virtual const ELEMENTTYPE& get( size_t row, size_t col ) {
                    return _elements[ _coord2blockindex( row, col ) ].get(
                        _r2rowinblock( row ), _c2colinblock( col ) );
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const override {
                    return _elements[ _coord2blockindex( row, col ) ].get(
                        _r2rowinblock( row ), _c2colinblock( col ) );
                }

                virtual Matrix<ELEMENTTYPE>& get_block( size_t row,
                                                        size_t col ) {
                    return _elements.at( _blockcoord2blockindex( row, col ) );
                }

                virtual const Matrix<ELEMENTTYPE>& get_block(
                    size_t row, size_t col ) const {
                    return _elements.at( _blockcoord2blockindex( row, col ) );
                }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_column(
                    size_t col ) const override {
                    this->check();
                    if (*this && this->columns() > col) {
                        size_t colstartblockind = _coord2blockindex( 0, col );
                        std::unique_ptr<Vector<ELEMENTTYPE>> c
                            = _elements.at( colstartblockind )
                                  .get_column( _c2colinblock( col ) );
                        c->resize( this->rows() );
                        for (size_t r = 1; r < this->block_rows(); ++r) {
                            std::unique_ptr<Vector<ELEMENTTYPE>> rc
                                = _elements
                                      .at( colstartblockind
                                           + r * _cols_of_blocks )
                                      .get_column( _c2colinblock( col ) );
                            std::vector<size_t> nz = rc->nonzero_indices();
                            for (auto nzit = nz.cbegin(); nzit != nz.cend();
                                 ++nzit) {
                                c->set( r * _rows_per_block + *nzit,
                                        rc->get( *nzit ) );
                            }
                        }
                        return c;
                    } else {
                        std::ostringstream oss;
                        oss << "Requested column index " << col
                            << " invalid for matrix with " << this->columns()
                            << " columns.";
                        throw std::invalid_argument( oss.str() );
                    }
                }

                virtual vector_ptr_t<ELEMENTTYPE> get_diagonal(
                    int offset = 0 ) const override {
                    this->check();
                    if (_elements.size() == 0) {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>();
                    } else {
                        // get the right kind of Vector
                        vector_ptr_t<ELEMENTTYPE> diag
                            = _elements[ 0 ].get_diagonal( 0 );
                        diag->zero();
                        int ndiag    = 0;
                        int startrow = 0;
                        // size properly.  account for non-square matrices
                        int sizediff
                            = (int)this->rows() - (int)this->columns();
                        int aoffset   = std::abs( offset );
                        int asizediff = std::abs( sizediff );

                        // case 1: square
                        if (sizediff == 0) {
                            ndiag    = (int)this->rows() - std::abs( offset );
                            startrow = offset < 0 ? -offset : 0;
                        } else if (sizediff > 0) {
                            if (offset >= 0) {
                                // case 2: row-dominant, main or upper diagonal
                                ndiag    = (int)this->columns() - aoffset;
                                startrow = 0;
                            } else if (aoffset <= asizediff) {
                                // case 3: row-dominant, lower diagonal within
                                // full-width
                                ndiag    = (int)this->columns();
                                startrow = (size_t)aoffset;
                            } else {
                                // case 4: row-dominant, lower diagonal outside
                                // full-width
                                ndiag    = (int)this->columns() - aoffset;
                                startrow = (size_t)aoffset;
                            }
                        } else {
                            if (offset <= 0) {
                                // case 5: column-dominant, main or lower
                                // diagonal
                                ndiag    = (int)this->rows() - aoffset;
                                startrow = aoffset;
                            } else if (aoffset <= asizediff) {
                                // case 6: row-dominant, upper diagonal within
                                // full-width
                                ndiag    = (int)this->rows();
                                startrow = 0;
                            } else {
                                // case 7: column-dominant, upper diagonal
                                // outside full-width
                                ndiag    = (int)this->rows() - aoffset;
                                startrow = 0;
                            }
                        }

                        if (ndiag <= 0) {
                            std::ostringstream oss;
                            oss << "Requested diagonal " << offset
                                << " invalid for " << this->rows() << " x "
                                << this->columns() << " matrix";
                            throw std::invalid_argument( oss.str() );
                        }
                        diag->resize( ndiag );

                        for (int i = 0; i < ndiag; ++i) {
                            int row = startrow + i;
                            int col = row + offset;
                            diag->set( i, this->get( row, col ) );
                        }

                        return diag;
                    }
                }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_row(
                    size_t row ) const override {
                    this->check();
                    if (*this && this->rows() > row) {
                        size_t rowstartblockind = _coord2blockindex( row, 0 );
                        std::unique_ptr<Vector<ELEMENTTYPE>> r
                            = _elements.at( rowstartblockind )
                                  .get_row( _r2rowinblock( row ) );
                        r->resize( this->columns() );
                        for (size_t c = 1; c < this->block_columns(); ++c) {
                            std::unique_ptr<Vector<ELEMENTTYPE>> cr
                                = _elements.at( rowstartblockind + c )
                                      .get_row( _r2rowinblock( row ) );
                            std::vector<size_t> nz = cr->nonzero_indices();
                            for (auto nzit = nz.cbegin(); nzit != nz.cend();
                                 ++nzit) {
                                r->set( c * _cols_per_block + *nzit,
                                        cr->get( *nzit ) );
                            }
                        }
                        return r;
                    } else {
                        std::ostringstream oss;
                        oss << "Requested row index " << row
                            << " invalid for matrix with " << this->rows()
                            << " rows.";
                        throw std::invalid_argument( oss.str() );
                    }
                }

                virtual BlockMatrix<ELEMENTTYPE>& identity() override {
                    if (!*this) {
                        throw std::logic_error(
                            "Block Matrix has not been initialized" );
                    }
                    if (!is_square()) {
                        throw std::invalid_argument(
                            "Cannot turn a non-square matrix into an "
                            "identity matrix" );
                    }
                    for (size_t i = 0; i < this->block_rows(); ++i) {
                        for (size_t j = 0; j < this->block_columns(); ++j) {
                            if (i == j) {
                                this->get_block( i, j ).identity();
                            } else {
                                this->get_block( i, j ).zero();
                            }
                        }
                    }
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& identity( size_t rows,
                                                       size_t cols ) override {
                    if (!*this) {
                        throw std::logic_error(
                            "Block Matrix has not been initialized" );
                    }
                    if (rows != cols) {
                        throw std::invalid_argument(
                            "Cannot create a non-square identity matrix!" );
                    }
                    return this->zero().resize( rows, cols ).identity();
                }

                virtual BlockMatrix<ELEMENTTYPE>& identity(
                    size_t rows, size_t cols, size_t blockrows,
                    size_t blockcols ) {
                    if (!*this) {
                        throw std::logic_error(
                            "Block Matrix has not been initialized" );
                    }
                    if (rows != cols || blockrows != blockcols) {
                        throw std::invalid_argument(
                            "Cannot create a non-square identity matrix!" );
                    }
                    this->clear();
                    this->resize( rows, cols, blockrows, blockcols );
                    this->identity();
                    return *this;
                }

                virtual bool is_band_diagonal() const override {
                    if (this->is_empty()) {
                        return true;
                    }
                    for (auto it = _elements.cbegin(); it != _elements.cend();
                         ++it) {
                        if (!it->is_band_diagonal()) {
                            return false;
                        }
                    }
                    return true;
                }

                virtual bool is_block_diagonal() const {
                    size_t ndiag = _maxblockdiag();
                    for (size_t di = 1; di < ndiag; ++di) {
                        // row index within main or upper diagonal.  Block
                        // coordinate is [i,i+di].  Matching lower diagonal
                        // coordinate is [i+di,i]
                        for (size_t i = 0; i < ndiag - di; ++i) {
                            if (!this->get_block( i, i + di ).is_zero()) {
                                return false;
                            }
                            if (!this->get_block( i + di, i ).is_zero()) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_block_matrix() const override { return true; }

                virtual bool is_block_tridiagonal() const {
                    for (int r = 0; r < this->block_rows(); ++r) {
                        for (int c = 0; c < ( r - 1 ); ++c) {
                            if (!this->get_block( r, c ).is_zero()) {
                                return false;
                            }
                        }
                        for (int c = r + 2; c < this->block_columns(); ++c) {
                            if (!this->get_block( r, c ).is_zero()) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_diagonal() const override {
                    for (size_t r = 0; r < _rows_of_blocks; ++r) {
                        for (size_t c = 0; c < _cols_of_blocks; ++c) {
                            if (r == c) {
                                if (!get_block( r, c ).is_diagonal()) {
                                    return false;
                                }
                            } else {
                                if (!get_block( r, c ).is_zero()) {
                                    return false;
                                }
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_empty() const override {
                    return ( rows() == 0 || columns() == 0 );
                }

                virtual bool is_identity() const {
                    if (this->is_empty()) {
                        return false;
                    }
                    for (size_t r = 0; r < _rows_of_blocks; ++r) {
                        for (size_t c = 0; c < _cols_of_blocks; ++c) {
                            if (r == c) {
                                if (!get_block( r, c ).is_identity()) {
                                    return false;
                                }
                            } else {
                                if (!get_block( r, c ).is_zero()) {
                                    return false;
                                }
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_lower_triangular() const override {
                    if (this->is_empty()) {
                        return false;
                    }
                    for (size_t r = 0; r < _rows_of_blocks; ++r) {
                        for (size_t c = r; c < _cols_of_blocks; ++c) {
                            if (r == c) {
                                if (!get_block( r, c ).is_lower_triangular()) {
                                    return false;
                                }
                            } else {
                                if (!get_block( r, c ).is_zero()) {
                                    return false;
                                }
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_square() const override {
                    return ( this->is_empty() || rows() == columns() );
                }

                virtual bool is_symmetric() const override {
                    if (this->block_rows() != this->block_columns()) {
                        return false;
                    }

                    // loop over diagonals
                    for (size_t di = 0; di < this->block_rows(); ++di) {
                        // row index within upper diagonal.  Block coordinate
                        // is [i,i+di].  Matching lower diagonal coordinate is
                        // [i+di,i]
                        if (!this->get_block( di, di ).is_symmetric()) {
                            return false;
                        }
                        for (size_t i = 1; i < this->block_columns() - di;
                             ++i) {
                            Matrix<ELEMENTTYPE> testmat
                                = this->get_block( i, i + di );
                            testmat.transpose();
                            if (this->get_block( i + di, i ) != testmat) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_tridiagonal() const override {
                    if (this->is_empty()) {
                        return false;
                    }
                    for (size_t r = 0; r < _rows_of_blocks; ++r) {
                        for (size_t c = 0; c < _cols_of_blocks; ++c) {
                            if (r == c) {
                                if (!get_block( r, c ).is_tridiagonal()) {
                                    return false;
                                }
                            } else {
                                if (!get_block( r, c ).is_zero()) {
                                    return false;
                                }
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_upper_triangular() const override {
                    if (this->is_empty()) {
                        return false;
                    }
                    for (size_t r = 0; r < _rows_of_blocks; ++r) {
                        for (size_t c = 0; c <= std::min( r, _cols_of_blocks );
                             ++c) {
                            if (r == c) {
                                if (!get_block( r, c ).is_upper_triangular()) {
                                    return false;
                                }
                            } else {
                                if (!get_block( r, c ).is_zero()) {
                                    return false;
                                }
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_zero() const override {
                    if (this->is_empty()) {
                        return true;
                    }
                    for (auto it = _elements.cbegin(); it != _elements.cend();
                         ++it) {
                        if (!it->is_zero()) {
                            return false;
                        }
                    }
                    return true;
                }

                virtual Vector<ELEMENTTYPE> left_multiply(
                    const Vector<ELEMENTTYPE>& v ) const override {
                    if (this->rows() != v.size()) {
                        std::ostringstream oss;
                        oss << "Size mismatch in vector-matrix "
                               "multiplication: "
                            << this->rows() << " rows in matrix vs "
                            << v.size() << " elements in vector";
                        throw std::invalid_argument( oss.str() );
                    }
                    Vector<ELEMENTTYPE> product
                        = VectorFactory<ELEMENTTYPE>::build(
                            ( (double)this->columns() / (double)v.size() < 0.5
                                  ? vector_t::SPARSE
                                  : vector_t::DENSE ),
                            this->columns() );
                    for (size_t i = 0; i < columns(); i++) {
                        ELEMENTTYPE sum = NCPA::math::zero<ELEMENTTYPE>();
                        for (size_t j = 0; j < rows(); j++) {
                            sum += get( j, i ) * v.get( j );
                        }
                        product.set( i, sum );
                    }
                    return product;
                }

                virtual size_t lower_bandwidth() const override {
                    if (this->is_empty() || this->is_zero()) {
                        return 0;
                    }
                    if (this->is_block_diagonal()) {
                        size_t lbw = 0;
                        for (size_t r = 0;
                             r < std::min( _rows_of_blocks, _cols_of_blocks );
                             ++r) {
                            lbw = std::max(
                                lbw, get_block( r, r ).lower_bandwidth() );
                        }
                        return lbw;
                    }
                    if (this->is_block_tridiagonal()) {
                        size_t lbw;
                        for (size_t r = 1;
                             r < std::min( _rows_of_blocks, _cols_of_blocks );
                             ++r) {
                            lbw = std::max(
                                lbw, get_block( r, r - 1 ).lower_bandwidth() );
                        }
                        return lbw + _rows_per_block;
                    }
                    size_t lbw      = 0;
                    int max_bdiag   = this->block_rows() - 1;
                    int n_max_bdiag = (int)std::min( this->block_rows(),
                                                     this->block_columns() );
                    int bcmin       = 0;
                    for (int bdiag = max_bdiag; bdiag > 1; --bdiag) {
                        int brmin  = bdiag;
                        int nbdiag = std::min(
                            n_max_bdiag, (int)this->block_rows() - bdiag );
                        int brmax = bdiag + nbdiag;
                        int bcmax = nbdiag;
                        for (int bc = bcmin; bc < bcmax; ++bc) {
                            int br = bdiag + bc;
                            if (!this->get_block( br, bc ).is_zero()) {
                                size_t blbw = this->get_block( br, bc )
                                                  .lower_bandwidth();
                                matrix_coordinate_t point = _blockcoord2coord(
                                    (size_t)br, (size_t)bc, blbw, 0 );
                                lbw = std::max( lbw, point.row );
                            }
                        }
                        if (lbw > 0) {
                            return lbw;
                        }
                    }
                    throw std::logic_error(
                        "Error calculating lower bandwidth, this shouldn't be "
                        "possible" );
                }

                virtual BlockMatrix<ELEMENTTYPE>& multiply(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    this->check();
                    if (this->block_columns() != other.block_rows()) {
                        std::ostringstream oss;
                        oss << "Block layout of first matrix ("
                            << this->block_rows() << " x "
                            << this->block_columns()
                            << " blocks) not multiplication-compatible with "
                               "second matrix ("
                            << other.block_rows() << " x "
                            << other.block_columns() << " blocks)";
                        throw std::invalid_argument( oss.str() );
                    }
                    if (this->columns_per_block() != other.rows_per_block()) {
                        std::ostringstream oss;
                        oss << "Block sizes of first matrix ("
                            << this->rows_per_block() << " x "
                            << this->columns_per_block()
                            << " blocks) not multiplication-compatible with "
                               "second matrix ("
                            << other.rows_per_block() << " x "
                            << other.columns_per_block() << " blocks)";
                        throw std::invalid_argument( oss.str() );
                    }
                    BlockMatrix<ELEMENTTYPE> product( _blocktype );
                    product.resize( this->block_rows(), other.block_columns(),
                                    this->rows_per_block(),
                                    other.columns_per_block() );

                    for (size_t r = 0; r < this->block_rows(); ++r) {
                        for (size_t c = 0; c < other.block_columns(); ++c) {
                            Matrix<ELEMENTTYPE> *P
                                = &product.get_block( r, c );
                            for (size_t k = 0; k < this->block_columns();
                                 ++k) {
                                const Matrix<ELEMENTTYPE> *A
                                    = &this->get_block( r, k );
                                const Matrix<ELEMENTTYPE> *B
                                    = &other.get_block( k, c );
                                if (!( A->is_zero() || B->is_zero() )) {
                                    // std::cout << "Computing A[ " << r << ", "
                                    //           << k << " ] * B[ " << k << ", "
                                    //           << c << " ] to add to C[ " << r
                                    //           << ", " << c << " ]"
                                    //           << std::endl;
                                    if (P->is_zero()) {
                                        *P = ( *A ) * ( *B );
                                    } else {
                                        *P += ( *A ) * ( *B );
                                    }
                                }
                            }
                        }
                    }

                    // for (size_t r = 0; r < product.block_rows(); ++r) {
                    //     for (size_t c = 0; c < product.block_columns(); ++c)
                    //     {
                    //         for (size_t k = 0; k < this->block_columns();
                    //              ++k) {
                    //             if (!( this->get_block( r, k ).is_zero()
                    //                    || other.get_block( k, c )
                    //                           .is_zero() )) {
                    //                 std::cout << "C[ " << r << ", " << c <<
                    //                 " ] += A[ " << r << ", " << k
                    //                             << " ] * B[ " << k << ", "
                    //                             << c << " ]" << std::endl;
                    //                 product.get_block( r, c )
                    //                     += this->get_block( r, k )
                    //                      * other.get_block( k, c );
                    //             }
                    //         }
                    //     }
                    // }
                    std::swap( *this, product );
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& multiply(
                    const Matrix<ELEMENTTYPE>& other ) override {
                    if (auto derived
                        = dynamic_cast<const BlockMatrix<ELEMENTTYPE> *>(
                            &other )) {
                        return this->multiply( *derived );
                    }
                    check();
                    other.check();
                    if (!( this->is_square() && other.is_square() )) {
                        throw NCPA::NotImplementedError(
                            "BlockMatrix.multiply(): Currently multiplication "
                            "of a block matrix by a non-block matrix requires "
                            "them both to be square." );
                    }
                    if (columns() != other.rows()) {
                        throw std::invalid_argument(
                            "Matrix size mismatch: cannot multiply" );
                    }

                    BlockMatrix<ELEMENTTYPE> product( _blocktype );
                    product.resize( _rows_of_blocks, _cols_of_blocks,
                                    _rows_per_block, _cols_per_block );
                    block_matrix_indexed_coordinate_t bcoord;
                    for (auto it = _elements.cbegin(); it != _elements.cend();
                         ++it) {
                        if (!it->is_zero()) {
                            bcoord.index
                                = std::distance( _elements.cbegin(), it );
                            bcoord.coordinates_in_block.row    = 0;
                            bcoord.coordinates_in_block.column = 0;
                            matrix_coordinate_t fullcoord
                                = _indexedblockcoord2coord( bcoord );
                            for (size_t br = 0; br < _rows_per_block; ++br) {
                                for (size_t bc = 0; bc < _cols_per_block;
                                     ++bc) {
                                    size_t i = fullcoord.row + br;
                                    size_t k = fullcoord.column + bc;
                                    // since Cij = Aik * Bkj, first see if Aik
                                    // is nonzero
                                    if (!it->is_zero( br, bc )) {
                                        auto krow = other.get_row( k );
                                        auto nonzero_j_inds
                                            = krow->internal()
                                                  ->nonzero_indices();
                                        for (auto jit = nonzero_j_inds.begin();
                                             jit != nonzero_j_inds.end();
                                             ++jit) {
                                            product.set(
                                                i, *jit,
                                                product.get( i, *jit )
                                                    + it->get( br, bc )
                                                          * krow->get(
                                                              *jit ) );
                                        }
                                    }
                                }
                            }
                        }
                    }
                    swap( *this, product );
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& resize(
                    size_t rows, size_t cols ) override {
                    if (rows % _rows_per_block != 0
                        || cols % _cols_per_block != 0) {
                        std::ostringstream oss;
                        oss << "Requested new size " << rows << " x " << cols
                            << " cannot be divided evenly into blocks of size "
                            << _rows_per_block << " x " << _cols_per_block
                            << "!";
                        throw std::range_error( oss.str() );
                    }
                    this->resize( rows / _rows_per_block,
                                  cols / _cols_per_block, _rows_per_block,
                                  _cols_per_block );
                    return *this;
                }

                virtual BlockMatrix& resize( size_t rows_of_blocks, size_t cols_of_blocks,
                                             size_t rows_per_block,
                                             size_t cols_per_block ) {
                    this->check();
                    if (rows_of_blocks == _rows_of_blocks && cols_of_blocks == _cols_of_blocks
                        && rows_per_block == _rows_per_block
                        && cols_per_block == _cols_per_block) {
                        return *this;
                    }
                    std::vector<Matrix<ELEMENTTYPE>> _newelements( rows_of_blocks
                                                                   * cols_of_blocks );
                    for (size_t rowind = 0; rowind < rows_of_blocks; ++rowind) {
                        for (size_t colind = 0; colind < cols_of_blocks; ++colind) {
                            size_t newind
                                = rc2index( rowind, colind, rows_of_blocks, cols_of_blocks );
                            if (rowind < _rows_of_blocks
                                && colind < _cols_of_blocks) {
                                std::swap( _newelements[ newind ],
                                           _elements[ _blockcoord2blockindex(
                                               rowind, colind ) ] );
                            } else {
                                _newelements[ newind ]
                                    = MatrixFactory<ELEMENTTYPE>::build(
                                        _blocktype );
                            }
                            _newelements[ newind ].resize( rows_per_block,
                                                           cols_per_block );
                        }
                    }
                    _rows_of_blocks = rows_of_blocks;
                    _cols_of_blocks = cols_of_blocks;
                    _rows_per_block = rows_per_block;
                    _cols_per_block = cols_per_block;
                    std::swap( _elements, _newelements );
                    return *this;
                }

                virtual Vector<ELEMENTTYPE> right_multiply(
                    const Vector<ELEMENTTYPE>& v ) const override {
                    if (this->columns() != v.size()) {
                        std::ostringstream oss;
                        oss << "Size mismatch in matrix-vector "
                               "multiplication: "
                            << this->columns() << " columns in matrix vs "
                            << v.size() << " elements in vector";
                        throw std::invalid_argument( oss.str() );
                    }
                    Vector<ELEMENTTYPE> product
                        = VectorFactory<ELEMENTTYPE>::build(vector_t::DENSE,
                            this->rows() );
                    for (size_t i = 0; i < rows(); i++) {
                        product.set( i, this->get_row( i )->dot( v ) );
                        // auto row = this->get_row( i );
                        // ELEMENTTYPE sum = NCPA::math::zero<ELEMENTTYPE>();
                        // for (size_t j = 0; j < columns(); j++) {
                        //     sum += get( i, j ) * v.get( j );
                        // }
                        // product.set( i, sum );
                    }
                    return product;
                }

                virtual size_t rows() const override {
                    return _rows_of_blocks * _rows_per_block;
                }

                virtual size_t rows_per_block() const {
                    return _rows_per_block;
                }

                virtual BlockMatrix<ELEMENTTYPE>& scale(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    this->check();
                    other.check();
                    _checksizes( other );
                    for (size_t r = 0; r < this->block_rows(); ++r) {
                        for (size_t c = 0; c < this->block_columns(); ++c) {
                            this->get_block( r, c ).scale(
                                other.get_block( r, c ) );
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& scale(
                    const Matrix<ELEMENTTYPE>& other ) override {
                    if (auto derived
                        = dynamic_cast<const BlockMatrix<ELEMENTTYPE> *>(
                            &other )) {
                        return this->scale( *derived );
                    }
                    this->check();
                    other.check();
                    this->check_same_size( other );
                    for (size_t i = 0; i < _elements.size(); ++i) {
                        if (!_elements[ i ].is_zero()) {
                            matrix_coordinate_t tl, br;
                            matrix_coordinate_span_t span
                                = _blockindex2span( i );
                            for (size_t r = 0; r < _rows_per_block; ++r) {
                                for (size_t c = 0; c < _cols_per_block; ++c) {
                                    ELEMENTTYPE product
                                        = _elements[ i ].get( r, c )
                                        * other.get( r + span.topleft.row,
                                                     c + span.topleft.column );
                                    if (!NCPA::math::is_zero( product )) {
                                        _elements[ i ].set( r, c, product );
                                    }
                                }
                            }
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& scale(
                    ELEMENTTYPE other ) override {
                    this->check();
                    for (auto it = _elements.begin(); it != _elements.end();
                         ++it) {
                        it->scale( other );
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& set(
                    size_t r, size_t c, ELEMENTTYPE val ) override {
                    this->check();
                    _elements[ _coord2blockindex( r, c ) ].set(
                        _r2rowinblock( r ), _c2colinblock( c ), val );
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& set(
                    ELEMENTTYPE val ) override {
                    this->check();
                    for (auto it = _elements.begin(); it != _elements.end();
                         ++it) {
                        it->set( val );
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_column(
                    size_t col, ELEMENTTYPE val ) override {
                    this->check();
                    for (size_t row  = 0; row < this->rows();
                         row        += _rows_per_block) {
                        block_matrix_indexed_coordinate_t coords
                            = _coord2indexedcoord( row, col );
                        _elements[ coords.index ].set_column(
                            coords.coordinates_in_block.column, val );
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_column(
                    size_t col, size_t nvals, const size_t *indices,
                    const ELEMENTTYPE *vals ) override {
                    this->check();
                    for (size_t i = 0; i < nvals; ++i) {
                        this->set( indices[ i ], col, vals[ i ] );
                    }
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::vector<size_t>& colinds,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    return this->set_column( col, colinds.size(),
                                             colinds.data(), vals.data() );
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    this->check();
                    if (vals.size() != this->rows()) {
                        std::ostringstream oss;
                        oss << "BlockMatrix.set_column(): Number of elements "
                               "in "
                               "vector ("
                            << vals.size()
                            << ") does not match number of rows ("
                            << this->rows() << ")!";
                        throw std::range_error( oss.str() );
                    }
                    std::vector<size_t> inds
                        = NCPA::arrays::index_vector<size_t>( vals.size() );
                    this->set_column( col, vals.size(), inds.data(),
                                      vals.data() );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::initializer_list<size_t> colinds,
                    const std::initializer_list<ELEMENTTYPE> vals ) override {
                    return this->set_column(
                        col, std::vector<size_t> { colinds },
                        std::vector<ELEMENTTYPE> { vals } );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const abstract_vector<ELEMENTTYPE>& vec ) override {
                    return this->set_column( col, vec.as_std() );
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_diagonal(
                    ELEMENTTYPE val, int offset = 0 ) override {
                    this->check();
                    matrix_coordinate_t coord
                        = matrix_diagonal_start( offset );
                    size_t ndiag = this->diagonal_size( offset );
                    for (size_t i = 0; i < ndiag; ++i) {
                        this->set( coord.row + i, coord.column + i, val );
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_diagonal(
                    size_t nvals, const ELEMENTTYPE *vals,
                    int offset = 0 ) override {
                    this->check();
                    if (nvals != this->diagonal_size( offset )) {
                        std::ostringstream oss;
                        oss << "set_diagonal(): Number of elements " << nvals
                            << " does not match diagonal size "
                            << this->diagonal_size( offset );
                        throw std::range_error( oss.str() );
                    }
                    matrix_coordinate_t coord
                        = matrix_diagonal_start( offset );
                    for (size_t i = 0; i < nvals; ++i) {
                        this->set( coord.row + i, coord.column + i,
                                   vals[ i ] );
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_row(
                    size_t row, ELEMENTTYPE val ) override {
                    this->check();
                    for (size_t col  = 0; col < this->columns();
                         col        += _cols_per_block) {
                        block_matrix_indexed_coordinate_t coords
                            = _coord2indexedcoord( row, col );
                        _elements[ coords.index ].set_row(
                            coords.coordinates_in_block.row, val );
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_row(
                    size_t row, size_t nvals, const size_t *indices,
                    const ELEMENTTYPE *vals ) override {
                    this->check();
                    for (size_t i = 0; i < nvals; ++i) {
                        this->set( row, indices[ i ], vals[ i ] );
                    }
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::vector<size_t>& rowinds,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    return this->set_row( row, rowinds.size(), rowinds.data(),
                                          vals.data() );
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    this->check();
                    if (vals.size() != this->columns()) {
                        std::ostringstream oss;
                        oss << "BlockMatrix.set_row(): Number of elements in "
                               "vector ("
                            << vals.size()
                            << ") does not match number of columns ("
                            << this->columns() << ")!";
                        throw std::range_error( oss.str() );
                    }
                    std::vector<size_t> inds
                        = NCPA::arrays::index_vector<size_t>( vals.size() );
                    this->set_row( row, vals.size(), inds.data(),
                                   vals.data() );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::initializer_list<size_t> rowinds,
                    const std::initializer_list<ELEMENTTYPE> vals ) override {
                    return this->set_row( row, std::vector<size_t> { rowinds },
                                          std::vector<ELEMENTTYPE> { vals } );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const abstract_vector<ELEMENTTYPE>& vec ) override {
                    return this->set_row( row, vec.as_std() );
                }

                virtual BlockMatrix<ELEMENTTYPE>& subtract(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    this->check();
                    other.check();
                    _checksizes( other );

                    for (size_t r = 0; r < this->block_rows(); ++r) {
                        for (size_t c = 0; c < this->block_columns(); ++c) {
                            if (!other.get_block( r, c ).is_zero()) {
                                this->get_block( r, c )
                                    -= other.get_block( r, c );
                            }
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& subtract(
                    const Matrix<ELEMENTTYPE>& other ) override {
                    if (auto derived
                        = dynamic_cast<const BlockMatrix<ELEMENTTYPE> *>(
                            &other )) {
                        return this->subtract( *derived );
                    }
                    this->check();
                    other.check();
                    check_same_size( other );

                    ELEMENTTYPE zero = NCPA::constants::zero<ELEMENTTYPE>();
                    for (size_t r = 0; r < this->rows(); ++r) {
                        for (size_t c = 0; c < this->columns(); ++c) {
                            ELEMENTTYPE otherval = other.get( r, c );
                            if (otherval != zero) {
                                this->set( r, c,
                                           this->get( r, c ) - otherval );
                            }
                        }
                    }
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& subtract( ELEMENTTYPE other ) {
                    return this->add( -other );
                }

                virtual Matrix<ELEMENTTYPE>& swap_rows(
                    size_t ind1, size_t ind2 ) override {
                    this->check();
                    auto row1 = this->get_row( ind1 );
                    auto row2 = this->get_row( ind2 );
                    return this->set_row( ind1, row2 ).set_row( ind2, row1 );
                }

                virtual Matrix<ELEMENTTYPE>& swap_columns(
                    size_t ind1, size_t ind2 ) override {
                    this->check();
                    auto col1 = this->get_row( ind1 );
                    auto col2 = this->get_row( ind2 );
                    return this->set_column( ind1, col2 )
                        .set_column( ind2, col1 );
                }

                BlockMatrix<ELEMENTTYPE> transpose() {
                    this->check();
                    BlockMatrix<ELEMENTTYPE> T( _blocktype );
                    T.resize( this->block_columns(), this->block_rows(),
                              this->columns_per_block(),
                              this->rows_per_block() );
                    for (size_t r = 0; r < this->block_rows(); ++r) {
                        for (size_t c = 0; c < this->block_columns(); ++c) {
                            auto Tt = this->get_block( r, c );
                            Tt.transpose();
                            T.set_block( c, r, Tt );
                        }
                    }
                    swap( *this, T );
                    return *this;
                }

                virtual size_t upper_bandwidth() const override {
                    if (this->is_zero()) {
                        return 0;
                    }
                    if (this->is_block_diagonal()) {
                        size_t ubw = 0;
                        for (size_t r = 0;
                             r < std::min( _rows_of_blocks, _cols_of_blocks );
                             ++r) {
                            ubw = std::max(
                                ubw, get_block( r, r ).upper_bandwidth() );
                        }
                        return ubw;
                    }
                    if (this->is_block_tridiagonal()) {
                        size_t ubw;
                        for (size_t r = 0; r > std::min( _rows_of_blocks - 1,
                                                         _cols_of_blocks - 1 );
                             ++r) {
                            ubw = std::max(
                                ubw, get_block( r, r + 1 ).upper_bandwidth() );
                        }
                        return ubw + _rows_per_block;
                    }
                    size_t ubw = 0;

                    int max_bdiag   = this->block_columns() - 1;
                    int n_max_bdiag = (int)std::min( this->block_rows(),
                                                     this->block_columns() );
                    int brmin       = 0;
                    for (int bdiag = max_bdiag; bdiag > 1; --bdiag) {
                        int bcmin  = bdiag;
                        int nbdiag = std::min(
                            n_max_bdiag, (int)this->block_columns() - bdiag );
                        int bcmax = bdiag + nbdiag;
                        int brmax = nbdiag;
                        for (int br = brmin; br < brmax; ++br) {
                            int bc = bdiag + br;
                            if (!this->get_block( br, bc ).is_zero()) {
                                size_t bubw = this->get_block( br, bc )
                                                  .upper_bandwidth();
                                matrix_coordinate_t point = _blockcoord2coord(
                                    (size_t)br, (size_t)bc, 0, bubw );
                                ubw = std::max( ubw, point.column );
                            }
                        }
                        if (ubw > 0) {
                            return ubw;
                        }
                    }
                    throw std::logic_error(
                        "Error calculating upper bandwidth. It shouldn't be "
                        "possible to reach this line." );
                }

                virtual BlockMatrix<ELEMENTTYPE>& zero( size_t r,
                                                        size_t c ) override {
                    this->check();
                    _elements[ _coord2blockindex( r, c ) ].zero(
                        _r2rowinblock( r ), _c2colinblock( c ) );
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& zero() override {
                    this->check();
                    for (auto it = _elements.begin(); it != _elements.end();
                         ++it) {
                        it->zero();
                    }
                    return *this;
                }

                ///////////////////////////////////////////////////
                // BELOW HERE NOT CHECKED


                virtual Matrix<ELEMENTTYPE> flatten() const {
                    this->check();
                    Matrix<ELEMENTTYPE> m
                        = MatrixFactory<ELEMENTTYPE>::build( _blocktype );
                    m.resize( this->rows(), this->columns() );
                    for (size_t r = 0; r < this->rows(); ++r) {
                        vector_ptr_t<ELEMENTTYPE> row = this->get_row( r );
                        std::vector<size_t> nz        = row->nonzero_indices();
                        for (auto nzit = nz.cbegin(); nzit != nz.cend();
                             ++nzit) {
                            m.set( r, *nzit, row->get( *nzit ) );
                        }
                    }
                    return m;
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_block(
                    size_t block_r, size_t block_c,
                    const Matrix<ELEMENTTYPE>& m ) {
                    this->check();
                    if (m.rows() != this->rows_per_block()
                        || m.columns() != this->columns_per_block()) {
                        std::ostringstream oss;
                        oss << "Matrix size " << m.rows() << " x "
                            << m.columns() << " does not match block size "
                            << this->rows_per_block() << " x "
                            << this->columns_per_block();
                        throw std::invalid_argument( oss.str() );
                    }

                    _elements.at( _blockcoord2blockindex( block_r, block_c ) )
                        = m;
                    return *this;
                }

                // virtual BlockMatrix<ELEMENTTYPE> transpose() {
                //     this->check();
                //     BlockMatrix<ELEMENTTYPE> T( _blocktype );
                //     T.resize( this->block_columns(), this->block_rows(),
                //               this->columns_per_block(),
                //               this->rows_per_block() );
                //     for (size_t r = 0; r < this->block_rows(); ++r) {
                //         for (size_t c = 0; c < this->block_columns(); ++c) {
                //             auto Tt = this->get_block( r, c );
                //             Tt.transpose();
                //             T.set_block( c, r, Tt );
                //         }
                //     }
                //     swap( *this, T );
                //     return *this;
                // }

                virtual BlockMatrix<ELEMENTTYPE>& like(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    this->resize( other.block_rows(), other.block_columns(),
                                  other.rows_per_block(),
                                  other.columns_per_block() );
                    return *this;
                }

                virtual bool same(
                    const BlockMatrix<ELEMENTTYPE>& other ) const {
                    return ( this->block_rows() == other.block_rows()
                             && this->block_columns() == other.block_columns()
                             && this->rows_per_block()
                                    == other.rows_per_block()
                             && this->columns_per_block()
                                    == other.columns_per_block() );
                }

                // virtual BlockMatrix<ELEMENTTYPE>& add(
                //     const BlockMatrix<ELEMENTTYPE>& other ) {
                //     this->check();
                //     other.check();
                //     _checksizes( other );
                //     for (size_t r = 0; r < this->block_rows(); ++r) {
                //         for (size_t c = 0; c < this->block_columns(); ++c) {
                //             _elements[ _blockcoord2blockindex( r, c ) ]
                //                 += other.get_block( r, c );
                //         }
                //     }
                //     return *this;
                // }

                // virtual BlockMatrix<ELEMENTTYPE>& subtract(
                //     const BlockMatrix<ELEMENTTYPE>& other ) {
                //     this->check();
                //     other.check();
                //     _checksizes( other );
                //     for (size_t r = 0; r < this->block_rows(); ++r) {
                //         for (size_t c = 0; c < this->block_columns(); ++c) {
                //             _elements[ _blockcoord2blockindex( r, c ) ]
                //                 -= other.get_block( r, c );
                //         }
                //     }
                //     return *this;
                // }

                // virtual BlockMatrix<ELEMENTTYPE>& scale(
                //     const BlockMatrix<ELEMENTTYPE>& other ) {
                //     this->check();
                //     other.check();
                //     _checksizes( other );
                //     for (size_t r = 0; r < this->block_rows(); ++r) {
                //         for (size_t c = 0; c < this->block_columns(); ++c) {
                //             _elements[ _blockcoord2blockindex( r, c )
                //             ].scale(
                //                 other.get_block( r, c ) );
                //         }
                //     }
                //     return *this;
                // }

                // virtual BlockMatrix<ELEMENTTYPE>& scale( ELEMENTTYPE other )
                // {
                //     this->check();
                //     for (size_t r = 0; r < this->block_rows(); ++r) {
                //         for (size_t c = 0; c < this->block_columns(); ++c) {
                //             _elements[ _blockcoord2blockindex( r, c )
                //             ].scale( other );
                //         }
                //     }
                //     return *this;
                // }

                // virtual BlockMatrix<ELEMENTTYPE>& multiply(
                //     const BlockMatrix<ELEMENTTYPE>& other ) {
                //     this->check();
                //     if (this->block_columns() != other.block_rows()) {
                //         std::ostringstream oss;
                //         oss << "Block layout of first matrix ("
                //             << this->block_rows() << " x "
                //             << this->block_columns()
                //             << " blocks) not multiplication-compatible with
                //             "
                //                "second matrix ("
                //             << other.block_rows() << " x "
                //             << other.block_columns() << " blocks)";
                //         throw std::invalid_argument( oss.str() );
                //     }
                //     if (this->columns_per_block() != other.rows_per_block())
                //     {
                //         std::ostringstream oss;
                //         oss << "Block sizes of first matrix ("
                //             << this->rows_per_block() << " x "
                //             << this->columns_per_block()
                //             << " blocks) not multiplication-compatible with
                //             "
                //                "second matrix ("
                //             << other.rows_per_block() << " x "
                //             << other.columns_per_block() << " blocks)";
                //         throw std::invalid_argument( oss.str() );
                //     }
                //     BlockMatrix<ELEMENTTYPE> product( _blocktype );
                //     product.resize( this->block_rows(),
                //     other.block_columns(),
                //                     this->rows_per_block(),
                //                     other.columns_per_block() );
                //     for (size_t r = 0; r < product.block_rows(); ++r) {
                //         for (size_t c = 0; c < product.block_columns(); ++c)
                //         {
                //             for (size_t k = 0; k < this->block_columns();
                //                  ++k) {
                //                 if (!( this->get_block( r, k ).is_zero()
                //                        || other.get_block( k, c )
                //                               .is_zero() )) {
                //                     product.get_block( r, c )
                //                         += this->get_block( r, k )
                //                          * other.get_block( k, c );
                //                 }
                //             }
                //         }
                //     }
                //     std::swap( *this, product );
                //     return *this;
                // }

                // virtual bool equals(
                //     const BlockMatrix<ELEMENTTYPE>& other ) const {
                //     // compare structure
                //     if (( (bool)*this ) != ( (bool)other )) {
                //         return false;
                //     }
                //     if (!*this && !other) {
                //         return true;
                //     }
                //     if (!( this->rows() == other.rows()
                //            && this->columns() == other.columns()
                //            && this->block_rows() == other.block_rows()
                //            && this->block_columns()
                //                   == other.block_columns() )) {
                //         return false;
                //     }

                //     // compare contents
                //     auto it1 = _elements.cbegin();
                //     auto it2 = other._elements.cbegin();
                //     for (; it1 != _elements.cend()
                //            && it2 != other._elements.cend();
                //          ++it1, ++it2) {
                //         if (*it1 != *it2) {
                //             return false;
                //         }
                //     }
                //     return true;
                // }

                // virtual bool equals( const Matrix<ELEMENTTYPE>& other )
                // const {
                //     // compare structure
                //     if (( (bool)*this ) != ( (bool)other )) {
                //         return false;
                //     }
                //     if (!*this && !other) {
                //         return true;
                //     }
                //     if (!( this->rows() == other.rows()
                //            && this->columns() == other.columns() )) {
                //         return false;
                //     }

                //     // compare contents
                //     for (size_t r = 0; r < this->rows(); ++r) {
                //         if (*this->get_row( r ) != *other.get_row( r )) {
                //             return false;
                //         }
                //     }
                //     return true;
                // }


                // friend operators
                friend bool operator==( const BlockMatrix<ELEMENTTYPE>& a,
                                        const BlockMatrix<ELEMENTTYPE>& b ) {
                    return a.equals( b );
                }

                friend bool operator==( const BlockMatrix<ELEMENTTYPE>& a,
                                        const Matrix<ELEMENTTYPE>& b ) {
                    return a.equals( b );
                }

                friend bool operator==( const Matrix<ELEMENTTYPE>& a,
                                        const BlockMatrix<ELEMENTTYPE>& b ) {
                    return b.equals( a );
                }

                friend bool operator!=( const BlockMatrix<ELEMENTTYPE>& a,
                                        const BlockMatrix<ELEMENTTYPE>& b ) {
                    return !( a.equals( b ) );
                }

                friend bool operator!=( const Matrix<ELEMENTTYPE>& a,
                                        const BlockMatrix<ELEMENTTYPE>& b ) {
                    return !( a.equals( b ) );
                }

                friend bool operator!=( const BlockMatrix<ELEMENTTYPE>& a,
                                        const Matrix<ELEMENTTYPE>& b ) {
                    return !( a.equals( b ) );
                }

                friend BlockMatrix<ELEMENTTYPE> operator+(
                    const BlockMatrix<ELEMENTTYPE>& c1,
                    const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend BlockMatrix<ELEMENTTYPE> operator+(
                    const BlockMatrix<ELEMENTTYPE>& c1,
                    const Matrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend BlockMatrix<ELEMENTTYPE> operator+(
                    const Matrix<ELEMENTTYPE>& c1,
                    const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c2 );
                    out += c1;
                    return out;
                }

                friend BlockMatrix<ELEMENTTYPE> operator+(
                    const BlockMatrix<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend BlockMatrix<ELEMENTTYPE> operator+(
                    ELEMENTTYPE c1, const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c2 );
                    out += c1;
                    return out;
                }

                friend BlockMatrix<ELEMENTTYPE> operator-(
                    const BlockMatrix<ELEMENTTYPE>& c1,
                    const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out -= c2;
                    return out;
                }

                BlockMatrix<ELEMENTTYPE> operator-() const {
                    BlockMatrix<ELEMENTTYPE> m = *this;
                    m.scale( -NCPA::math::one<ELEMENTTYPE>() );
                    return m;
                }

                friend BlockMatrix<ELEMENTTYPE> operator-(
                    const BlockMatrix<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out -= c2;
                    return out;
                }

                friend NCPA::linear::BlockMatrix<ELEMENTTYPE> operator*(
                    const BlockMatrix<ELEMENTTYPE>& c1,
                    const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out *= c2;
                    return out;
                }

                // friend NCPA::linear::Matrix<ELEMENTTYPE> operator*(
                //     const BlockMatrix<ELEMENTTYPE>& c1,
                //     const Matrix<ELEMENTTYPE>& c2 ) {
                //     // can we blockify c2?
                //     if (c1.columns() != c2.rows()) {
                //         std::ostringstream oss;
                //         oss << "Incompatible matrix sizes for
                //         multiplication: "
                //             << c1.rows() << " x " << c1.columns() << " vs "
                //             << c2.rows() << " x " << c2.columns();
                //         throw std::range_error( oss.str() );
                //     }

                //     Matrix<ELEMENTTYPE> out(c1);
                //     out                     *= c2;
                //     return out;
                // }

                // friend NCPA::linear::Matrix<ELEMENTTYPE> operator*(
                //     const Matrix<ELEMENTTYPE>& c1,
                //     const BlockMatrix<ELEMENTTYPE>& c2 ) {
                //     Matrix<ELEMENTTYPE> out( c1 );
                //     out *= c2.flatten();
                //     return out;
                // }

                friend BlockMatrix<ELEMENTTYPE> operator*(
                    const BlockMatrix<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out *= c2;
                    return out;
                }

                friend BlockMatrix<ELEMENTTYPE> operator*(
                    ELEMENTTYPE c1, const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c2 );
                    out *= c1;
                    return out;
                }

            private:
                size_t _rows_of_blocks, _cols_of_blocks, _rows_per_block,
                    _cols_per_block;
                matrix_t _blocktype;
                std::vector<Matrix<ELEMENTTYPE>> _elements;

                void _check_this() const {
                    if (!*this) {
                        throw std::logic_error(
                            "Matrix has not been initialized" );
                    }
                }

                void _checksizes( const BlockMatrix<ELEMENTTYPE>& m ) const {
                    if (this->block_rows() != m.block_rows()
                        || this->block_columns() != m.block_columns()) {
                        std::ostringstream oss;
                        oss << "Block layout of second matrix ("
                            << m.block_rows() << " x " << m.block_columns()
                            << " blocks) does not match first matrix ("
                            << this->block_rows() << " x "
                            << this->block_columns() << " blocks)";
                        throw std::invalid_argument( oss.str() );
                    }
                }

                size_t _maxblockdiag() const {
                    return std::max( _rows_of_blocks, _cols_of_blocks ) - 1;
                }

                // coordinate conversions

                // given a matrix-relative row, returns the block row
                // containing that row
                size_t _r2blockrow( size_t row ) const {
                    return row / _rows_per_block;
                }

                // given a matrix-relative column, returns the block column
                // containing that column
                size_t _c2blockcol( size_t col ) const {
                    return col / _cols_per_block;
                }

                // given a matrix-relative row, returns the row number within
                // its block
                size_t _r2rowinblock( size_t row ) const {
                    return row % _rows_per_block;
                }

                // given a matrix-relative column, returns the column number
                // within its block
                size_t _c2colinblock( size_t col ) const {
                    return col % _cols_per_block;
                }

                static matrix_coordinate_t _rc2coord( size_t row,
                                                      size_t col ) {
                    return matrix_coordinate_t { row, col };
                }

                // from block coordinates
                size_t _blockcoord2blockindex(
                    const matrix_coordinate_t& bcoord ) const {
                    return _blockcoord2blockindex( bcoord.row, bcoord.column );
                }

                size_t _blockcoord2blockindex( size_t brow,
                                               size_t bcol ) const {
                    return rc2index( brow, bcol, _rows_of_blocks,
                                     _cols_of_blocks );
                }

                matrix_coordinate_t _blockcoord2coord( size_t brow,
                                                       size_t bcol, size_t row,
                                                       size_t col ) const {
                    return matrix_coordinate_t { brow * _rows_per_block + row,
                                                 bcol * _cols_per_block
                                                     + col };
                }

                matrix_coordinate_t _blockcoord2coord(
                    matrix_coordinate_t bcoord,
                    matrix_coordinate_t coord_in_block ) const {
                    return _blockcoord2coord( bcoord.row, bcoord.column,
                                              coord_in_block.row,
                                              coord_in_block.column );
                }

                // from block index
                matrix_coordinate_t _blockindex2blockcoord(
                    size_t index ) const {
                    return matrix_coordinate_t { index / _cols_of_blocks,
                                                 index % _cols_of_blocks };
                }

                matrix_coordinate_span_t _blockindex2span(
                    size_t index ) const {
                    matrix_coordinate_t block
                        = _blockindex2blockcoord( index );
                    matrix_coordinate_span_t span;
                    span.topleft.row     = block.row * _rows_per_block;
                    span.topleft.column  = block.column * _cols_per_block;
                    span.bottomright.row = span.topleft.row + _rows_per_block;
                    span.bottomright.column
                        = span.topleft.column + _cols_per_block;
                    return span;
                }

                // from indexed block coordinates
                matrix_coordinate_t _indexedblockcoord2coord(
                    const block_matrix_indexed_coordinate_t& bmic ) const {
                    matrix_coordinate_t blockcoord
                        = _blockindex2blockcoord( bmic.index );
                    matrix_coordinate_t overall;
                    overall.row = blockcoord.row * _rows_per_block
                                + bmic.coordinates_in_block.row;
                    overall.column = blockcoord.column * _cols_per_block
                                   + bmic.coordinates_in_block.column;
                    return overall;
                }

                // from overall matrix coordinates
                matrix_coordinate_t _coord2blockcoord(
                    const matrix_coordinate_t& coord ) const {
                    return _coord2blockcoord( coord.row, coord.column );
                }

                matrix_coordinate_t _coord2blockcoord( size_t row,
                                                       size_t col ) const {
                    return matrix_coordinate_t { _r2blockrow( row ),
                                                 _c2blockcol( col ) };
                }

                size_t _coord2blockindex(
                    const matrix_coordinate_t& coord ) const {
                    return _coord2blockindex( coord.row, coord.column );
                }

                size_t _coord2blockindex( size_t row, size_t col ) const {
                    return _blockcoord2blockindex(
                        _coord2blockcoord( row, col ) );
                }

                matrix_coordinate_t _coord2coordinblock(
                    const matrix_coordinate_t& coord ) const {
                    return _coord2coordinblock( coord.row, coord.column );
                }

                matrix_coordinate_t _coord2coordinblock( size_t row,
                                                         size_t col ) const {
                    return matrix_coordinate_t { _r2rowinblock( row ),
                                                 _c2colinblock( col ) };
                }

                block_matrix_coordinate_t _coord2blockmatrixcoord(
                    const matrix_coordinate_t& coord ) const {
                    return _coord2blockmatrixcoord( coord.row, coord.column );
                }

                block_matrix_coordinate_t _coord2blockmatrixcoord(
                    size_t row, size_t col ) const {
                    block_matrix_coordinate_t bmc;
                    bmc.block_coordinates    = _coord2blockcoord( row, col );
                    bmc.coordinates_in_block = _coord2coordinblock( row, col );
                    return bmc;
                }

                block_matrix_indexed_coordinate_t _coord2indexedcoord(
                    size_t row, size_t col ) const {
                    block_matrix_indexed_coordinate_t icoord;
                    icoord.index = _coord2blockindex( row, col );
                    icoord.coordinates_in_block
                        = _coord2coordinblock( row, col );
                    return icoord;
                }
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::BlockMatrix<T>& a,
                  NCPA::linear::BlockMatrix<T>& b ) noexcept {
    using std::swap;
    swap( a._elements, b._elements );
    swap( a._rows_of_blocks, b._rows_of_blocks );
    swap( a._cols_of_blocks, b._cols_of_blocks );
    swap( a._rows_per_block, b._rows_per_block );
    swap( a._cols_per_block, b._cols_per_block );
    swap( a._blocktype, b._blocktype );
}
