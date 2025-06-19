#pragma once

#include "NCPA/linearalgebra/builders.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/vector.hpp"
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
            class BlockMatrix<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC> {
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
                    for (auto it = other._elements.cbegin();
                         it != other._elements.cend(); ++it) {
                        _elements.push_back( *it );
                    }
                    _rows_of_blocks = other._rows_of_blocks;
                    _cols_of_blocks = other._cols_of_blocks;
                    _rows_per_block = other._rows_per_block;
                    _cols_per_block = other._cols_per_block;
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

                // get the element by its overall coordinates
                virtual const ELEMENTTYPE& get( size_t row, size_t col ) {
                    return _elements[ _rc2blockindex( row, col ) ].get(
                        _r2blockrow( row ), _c2blockcol( col ) );
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const {
                    return _elements[ _rc2blockindex( row, col ) ].get(
                        _r2blockrow( row ), _c2blockcol( col ) );
                }

                virtual Matrix<ELEMENTTYPE>& get_block( size_t row,
                                                        size_t col ) {
                    return _elements.at( _blockindex( row, col ) );
                }

                virtual const Matrix<ELEMENTTYPE>& get_block(
                    size_t row, size_t col ) const {
                    return _elements.at( _blockindex( row, col ) );
                }

                explicit operator bool() const {
                    return ( _blocktype != matrix_t::INVALID );
                }

                virtual BlockMatrix& clear() {
                    _elements.clear();
                    _rows_of_blocks = 0;
                    _cols_of_blocks = 0;
                    _rows_per_block = 0;
                    _cols_per_block = 0;
                    return *this;
                }

                virtual size_t block_columns() const {
                    return _cols_of_blocks;
                }

                virtual size_t block_rows() const { return _rows_of_blocks; }

                virtual size_t rows_per_block() const {
                    return _rows_per_block;
                }

                virtual size_t columns_per_block() const {
                    return _cols_per_block;
                }

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
                            _newelements.back().copy( *it );
                        }
                        std::swap( _elements, _newelements );
                    }
                    return *this;
                }

                virtual matrix_t block_type() const { return _blocktype; }

                virtual BlockMatrix& resize( size_t rows, size_t cols,
                                             size_t blockrows,
                                             size_t blockcols ) {
                    _checkthis();
                    if (rows == _rows_of_blocks && cols == _cols_of_blocks
                        && blockrows == _rows_per_block
                        && blockcols == _cols_per_block) {
                        return *this;
                    }
                    std::vector<Matrix<ELEMENTTYPE>> _newelements( rows
                                                                   * cols );
                    for (size_t rowind = 0; rowind < rows; ++rowind) {
                        for (size_t colind = 0; colind < cols; ++colind) {
                            size_t newind
                                = rc2index( rowind, colind, rows, cols );
                            if (rowind < _rows_of_blocks
                                && colind < _cols_of_blocks) {
                                std::swap( _newelements[ newind ],
                                           _elements[ _blockindex(
                                               rowind, colind ) ] );
                            } else {
                                _newelements[ newind ]
                                    = MatrixFactory<ELEMENTTYPE>::build(
                                        _blocktype );
                            }
                            _newelements[ newind ].resize( blockrows,
                                                           blockcols );
                        }
                    }
                    _rows_of_blocks = rows;
                    _cols_of_blocks = cols;
                    _rows_per_block = blockrows;
                    _cols_per_block = blockcols;
                    std::swap( _elements, _newelements );
                    return *this;
                }

                virtual BlockMatrix& resize( size_t rows, size_t cols ) {
                    return this->resize( rows, cols, _rows_per_block,
                                         _cols_per_block );
                }

                virtual BlockMatrix<ELEMENTTYPE>& identity() {
                    if (!*this) {
                        throw std::logic_error(
                            "Block Matrix has not been initialized" );
                    }
                    if (!is_square()) {
                        throw std::invalid_argument(
                            "Cannot turn a non-square matrix into an "
                            "identity "
                            "matrix" );
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
                    return this->clear()
                        .resize( rows, cols, blockrows, blockcols )
                        .identity();
                }

                virtual bool is_square() const {
                    return ( this->rows() == this->columns()
                             && this->block_rows() == this->block_columns() );
                }

                virtual bool is_empty() const {
                    return ( this->rows() == 0 && this->columns() == 0 );
                }

                virtual bool is_zero() const {
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

                virtual bool is_identity() const {
                    if (!this->is_square()) {
                        return false;
                    }
                    for (size_t i = 0; i < this->block_rows(); ++i) {
                        if (!this->get_block( i, i ).is_identity()) {
                            return false;
                        }
                    }
                    return true;
                }

                virtual bool is_symmetric() const {
                    if (this->block_rows() != this->block_columns()) {
                        return false;
                    }

                    // loop over diagonals
                    for (size_t di = 1; di < this->block_rows(); ++di) {
                        // row index within upper diagonal.  Block coordinate
                        // is [i,i+di].  Matching lower diagonal coordinate is
                        // [i+di,i]
                        for (size_t i = 0; i < this->block_rows() - di; ++i) {
                            Matrix<ELEMENTTYPE> testmat
                                = this->get_block( i + di, i );
                            testmat.transpose();
                            if (this->get_block( i, i + di ) != testmat) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_diagonal() const {
                    size_t ndiag = _maxblockdiag();
                    for (size_t i = 0; i < ndiag; ++i) {
                        if (!this->get_block( i, i ).is_diagonal()) {
                            return false;
                        }
                    }
                    return this->is_block_diagonal();
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

                virtual bool is_tridiagonal() const {
                    if (_rows_per_block > 1 && _cols_per_block > 1) {
                        if (!this->is_block_diagonal()) {
                            return false;
                        }
                        size_t ndiag = _maxblockdiag();
                        for (size_t i = 0; i < ndiag; ++i) {
                            if (!this->get_block( i, i ).is_tridiagonal()) {
                                return false;
                            }
                        }
                        return true;
                    } else {
                        // special (stupid) case, block size is 1x1
                        return this->is_block_tridiagonal();
                    }
                }

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

                virtual size_t rows() const {
                    return _rows_of_blocks * _rows_per_block;
                }

                virtual size_t columns() const {
                    return _cols_of_blocks * _cols_per_block;
                }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_row(
                    size_t row ) const {
                    _checkthis();
                    if (*this && this->rows() > row) {
                        size_t rowstartblockind = _rc2blockindex( row, 0 );
                        std::unique_ptr<Vector<ELEMENTTYPE>> r
                            = _elements.at( rowstartblockind )
                                  .get_row( _r2blockrow( row ) );
                        r->resize( this->columns() );
                        for (size_t c = 1; c < this->block_columns(); ++c) {
                            std::unique_ptr<Vector<ELEMENTTYPE>> cr
                                = _elements.at( rowstartblockind + c )
                                      .get_row( _r2blockrow( row ) );
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

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_column(
                    size_t col ) const {
                    _checkthis();
                    if (*this && this->columns() > col) {
                        size_t colstartblockind = _rc2blockindex( 0, col );
                        std::unique_ptr<Vector<ELEMENTTYPE>> c
                            = _elements.at( colstartblockind )
                                  .get_column( _c2blockcol( col ) );
                        c->resize( this->rows() );
                        for (size_t r = 1; r < this->block_rows(); ++r) {
                            std::unique_ptr<Vector<ELEMENTTYPE>> rc
                                = _elements
                                      .at( colstartblockind
                                           + r * _cols_of_blocks )
                                      .get_column( _c2blockcol( col ) );
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

                virtual vector_ptr_t<ELEMENTTYPE> get_diagonal( int offset
                                                                = 0 ) const {
                    _checkthis();
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

                virtual Matrix<ELEMENTTYPE> as_matrix() const {
                    _checkthis();
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

                virtual BlockMatrix<ELEMENTTYPE>& set( size_t r, size_t c,
                                                       ELEMENTTYPE val ) {
                    _checkthis();
                    _elements[ _blockindex( r, c ) ].set(
                        _r2blockrow( r ), _c2blockcol( c ), val );
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& zero( size_t r, size_t c ) {
                    _checkthis();
                    _elements[ _rc2blockindex( r, c ) ].zero(
                        _r2blockrow( r ), _c2blockcol( c ) );
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& zero() {
                    _checkthis();
                    for (auto it = _elements.begin(); it != _elements.end();
                         ++it) {
                        it->zero();
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& set_block(
                    size_t r, size_t c, const Matrix<ELEMENTTYPE>& m ) {
                    _checkthis();
                    if (m.rows() != this->rows_per_block()
                        || m.columns() != this->columns_per_block()) {
                        std::ostringstream oss;
                        oss << "Matrix size " << m.rows() << " x "
                            << m.columns() << " does not match block size "
                            << this->rows_per_block() << " x "
                            << this->columns_per_block();
                        throw std::invalid_argument( oss.str() );
                    }

                    _elements.at( _blockindex( r, c ) ) = m;
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE> transpose() {
                    _checkthis();
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

                virtual BlockMatrix<ELEMENTTYPE>& like(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    this->resize( other.block_rows(), other.block_columns(),
                                  other.rows_per_block(),
                                  other.columns_per_block() );
                    return *this;
                }

                virtual bool same( const BlockMatrix<ELEMENTTYPE>& other ) const {
                    return (this->block_rows() == other.block_rows()
                            && this->block_columns() == other.block_columns()
                        && this->rows_per_block() == other.rows_per_block()
                        && this->columns_per_block() == other.columns_per_block() );
                }

                virtual BlockMatrix<ELEMENTTYPE>& add(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    _checkthis();
                    _checkother( other );
                    _checksizes( other );
                    for (size_t r = 0; r < this->block_rows(); ++r) {
                        for (size_t c = 0; c < this->block_columns(); ++c) {
                            _elements[ _blockindex( r, c ) ]
                                += other.get_block( r, c );
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& subtract(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    _checkthis();
                    _checkother( other );
                    _checksizes( other );
                    for (size_t r = 0; r < this->block_rows(); ++r) {
                        for (size_t c = 0; c < this->block_columns(); ++c) {
                            _elements[ _blockindex( r, c ) ]
                                -= other.get_block( r, c );
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& scale(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    _checkthis();
                    _checkother( other );
                    _checksizes( other );
                    for (size_t r = 0; r < this->block_rows(); ++r) {
                        for (size_t c = 0; c < this->block_columns(); ++c) {
                            _elements[ _blockindex( r, c ) ].scale(
                                other.get_block( r, c ) );
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& scale( ELEMENTTYPE other ) {
                    _checkthis();
                    for (size_t r = 0; r < this->block_rows(); ++r) {
                        for (size_t c = 0; c < this->block_columns(); ++c) {
                            _elements[ _blockindex( r, c ) ].scale( other );
                        }
                    }
                    return *this;
                }

                virtual BlockMatrix<ELEMENTTYPE>& multiply(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    _checkthis();
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
                    for (size_t r = 0; r < product.block_rows(); ++r) {
                        for (size_t c = 0; c < product.block_columns(); ++c) {
                            for (size_t k = 0; k < this->block_columns();
                                 ++k) {
                                product.get_block( r, c )
                                    += this->get_block( r, k )
                                     * other.get_block( k, c );
                            }
                        }
                    }
                    std::swap( *this, product );
                    return *this;
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

                virtual bool equals( const Matrix<ELEMENTTYPE>& other ) const {
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

                // assignment operators
                virtual BlockMatrix<ELEMENTTYPE>& operator+=(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    return this->add( other );
                }

                // virtual BlockMatrix<ELEMENTTYPE>& operator+=(
                //     const ELEMENTTYPE& other ) {
                //     return this->add( other );
                // }

                virtual BlockMatrix<ELEMENTTYPE>& operator-=(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    return this->subtract( other );
                }

                // virtual BlockMatrix<ELEMENTTYPE>& operator-=(
                //     const ELEMENTTYPE& other ) {
                //     return this->add( -other );
                // }

                virtual BlockMatrix<ELEMENTTYPE>& operator*=(
                    const BlockMatrix<ELEMENTTYPE>& other ) {
                    return this->multiply( other );
                }

                virtual BlockMatrix<ELEMENTTYPE>& operator*=(
                    const ELEMENTTYPE& other ) {
                    return this->scale( other );
                }

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

                friend BlockMatrix<ELEMENTTYPE> operator+(
                    const BlockMatrix<ELEMENTTYPE>& c1,
                    const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                // friend BlockMatrix<ELEMENTTYPE> operator+(
                //     const BlockMatrix<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                //     BlockMatrix<ELEMENTTYPE> out( c1 );
                //     out += c2;
                //     return out;
                // }

                // friend BlockMatrix<ELEMENTTYPE> operator+(
                //     ELEMENTTYPE c1, const BlockMatrix<ELEMENTTYPE>& c2 ) {
                //     BlockMatrix<ELEMENTTYPE> out( c2 );
                //     out += c1;
                //     return out;
                // }

                friend BlockMatrix<ELEMENTTYPE> operator-(
                    const BlockMatrix<ELEMENTTYPE>& c1,
                    const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out -= c2;
                    return out;
                }

                virtual BlockMatrix<ELEMENTTYPE> operator-() const {
                    BlockMatrix<ELEMENTTYPE> m = *this;
                    m.scale( -NCPA::math::one<ELEMENTTYPE>() );
                    return m;
                }

                // friend BlockMatrix<ELEMENTTYPE> operator-(
                //     const BlockMatrix<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                //     BlockMatrix<ELEMENTTYPE> out( c1 );
                //     out -= c2;
                //     return out;
                // }

                friend NCPA::linear::BlockMatrix<ELEMENTTYPE> operator*(
                    const BlockMatrix<ELEMENTTYPE>& c1,
                    const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    BlockMatrix<ELEMENTTYPE> out( c1 );
                    out *= c2;
                    return out;
                }

                friend NCPA::linear::Matrix<ELEMENTTYPE> operator*(
                    const BlockMatrix<ELEMENTTYPE>& c1,
                    const Matrix<ELEMENTTYPE>& c2 ) {
                    Matrix<ELEMENTTYPE> out = c1.as_matrix();
                    out *= c2;
                    return out;
                }

                friend NCPA::linear::Matrix<ELEMENTTYPE> operator*(
                    const Matrix<ELEMENTTYPE>& c1,
                    const BlockMatrix<ELEMENTTYPE>& c2 ) {
                    Matrix<ELEMENTTYPE> out( c1 );
                    out *= c2.as_matrix();
                    return out;
                }

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

            protected:
                size_t _rows_of_blocks, _cols_of_blocks, _rows_per_block,
                    _cols_per_block;
                matrix_t _blocktype;
                std::vector<Matrix<ELEMENTTYPE>> _elements;

                void _checkthis() const {
                    if (!*this) {
                        throw std::logic_error(
                            "Matrix has not been initialized" );
                    }
                }

                void _checkother( const BlockMatrix<ELEMENTTYPE>& m ) {
                    if (!m) {
                        throw std::logic_error(
                            "Other matrix has not been initialized" );
                    }
                }

                void _checksizes( const BlockMatrix<ELEMENTTYPE>& m ) {
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
                    return std::min( _rows_of_blocks, _cols_of_blocks );
                }

                // gets a block index based on its block coordinates
                size_t _blockindex( size_t row, size_t col ) const {
                    return rc2index( row, col, _rows_of_blocks,
                                     _cols_of_blocks );
                }

                // gets a block index based on the overall element coordinates
                size_t _rc2blockindex( size_t row, size_t col ) const {
                    return _blockindex( row / _rows_per_block,
                                        col / _cols_per_block );
                }

                size_t _r2blockrow( size_t row ) const {
                    return row % _rows_per_block;
                }

                size_t _c2blockcol( size_t col ) const {
                    return col % _cols_per_block;
                }

                void _rc2blockrc( size_t row, size_t col, size_t& index,
                                  size_t& blockrow, size_t& blockcol ) const {
                    index    = _blockindex( row, col );
                    blockrow = row % _rows_per_block;
                    blockcol = col % _cols_per_block;
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
