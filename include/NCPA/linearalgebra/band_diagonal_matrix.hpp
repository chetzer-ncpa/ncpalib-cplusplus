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

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::band_diagonal_matrix, ELEMENTTYPE );

namespace NCPA {
    namespace linear {

        enum class diagonal_index_status_t {
            INVALID,
            VALID_BUT_NOT_DEFINED,
            VALID
        };

        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class band_diagonal_matrix<ELEMENTTYPE,
                                       _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public abstract_matrix<ELEMENTTYPE> {
            public:
                using abstract_matrix<ELEMENTTYPE>::set_row;
                using abstract_matrix<ELEMENTTYPE>::set_column;
                using abstract_matrix<ELEMENTTYPE>::set_diagonal;
                friend class BandDiagonalLUDecomposition<ELEMENTTYPE>;
                friend class basic_band_diagonal_linear_system_solver<
                    ELEMENTTYPE>;

                band_diagonal_matrix( size_t nrows, size_t ncols,
                                      size_t n_lower_offdiags,
                                      size_t n_upper_offdiags ) :
                    _nrows { nrows },
                    _ncols { ncols },
                    _n_lower { n_lower_offdiags },
                    _n_upper { n_upper_offdiags } {
                    _prep_diagonals();
                }

                band_diagonal_matrix( size_t nrows, size_t ncols ) :
                    band_diagonal_matrix<ELEMENTTYPE>( nrows, ncols, 0, 0 ) {}

                band_diagonal_matrix() :
                    band_diagonal_matrix<ELEMENTTYPE>( 0, 0, 0, 0 ) {}

                band_diagonal_matrix(
                    const band_diagonal_matrix<ELEMENTTYPE>& other ) :
                    band_diagonal_matrix<ELEMENTTYPE>( other.rows(),
                                                       other.columns() ) {
                    _n_lower  = other._n_lower;
                    _n_upper  = other._n_upper;
                    _contents = other._contents;
                }

                band_diagonal_matrix(
                    const abstract_matrix<ELEMENTTYPE>& other ) :
                    band_diagonal_matrix<ELEMENTTYPE>() {
                    this->copy( other );
                }

                band_diagonal_matrix(
                    band_diagonal_matrix<ELEMENTTYPE>&& source ) noexcept :
                    band_diagonal_matrix<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                // constructors, destructors, copying, and assignment
                virtual ~band_diagonal_matrix() {}

                friend void ::swap<ELEMENTTYPE>(
                    band_diagonal_matrix<ELEMENTTYPE>& a,
                    band_diagonal_matrix<ELEMENTTYPE>& b ) noexcept;

                band_diagonal_matrix<ELEMENTTYPE>& operator=(
                    band_diagonal_matrix<ELEMENTTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& add(
                    const abstract_matrix<ELEMENTTYPE>& babs,
                    ELEMENTTYPE modifier = 1.0 ) override {
                    if (!is_this_subclass( babs )) {
                        return dynamic_cast<
                            band_diagonal_matrix<ELEMENTTYPE>&>(
                            abstract_matrix<ELEMENTTYPE>::add( babs,
                                                               modifier ) );
                    }
                    const band_diagonal_matrix<ELEMENTTYPE> *b
                        = downcast( &babs );

                    this->check_size( *b );
                    int sum_n_lower
                        = (int)std::max( this->_n_lower, b->_n_lower );
                    int sum_n_upper
                        = (int)std::max( this->_n_upper, b->_n_upper );
                    if (_n_lower < sum_n_lower) {
                        this->_add_subdiagonal( sum_n_lower - _n_lower );
                    }
                    if (_n_upper < sum_n_upper) {
                        this->_add_superdiagonal( sum_n_upper - _n_upper );
                    }
                    for (int i = -sum_n_lower; i <= sum_n_upper; ++i) {
                        if (b->has_diagonal( i )) {
                            int i1 = 0, i2 = 0;
                            if (!this->diag2internal( i, i1 )) {
                                throw std::logic_error(
                                    "Diagonal index mismatch in first "
                                    "matrix" );
                            }
                            if (!b->diag2internal( i, i2 )) {
                                throw std::logic_error(
                                    "Diagonal index mismatch in second "
                                    "matrix" );
                            }
                            if (this->internal( i1 ).empty()) {
                                this->internal( i1 )
                                    = NCPA::arrays::scale_vector(
                                        b->contents().at( i2 ), modifier );
                            } else {
                                this->internal( i1 )
                                    = NCPA::arrays::add_vectors(
                                        this->internal( i1 ),
                                        NCPA::arrays::scale_vector(
                                            b->contents().at( i2 ),
                                            modifier ) );
                            }
                        }
                    }


                    // // lower half
                    // size_t skip = 0, add_on = 0, both = 0;
                    // if (_n_lower > b->_n_lower) {
                    //     skip = _n_lower - b->_n_lower;
                    //     both = b->_n_lower;
                    // } else if (_n_lower < b->_n_lower) {
                    //     add_on = b->_n_lower - _n_lower;
                    //     both   = _n_lower;
                    // } else {
                    //     both = _n_lower;
                    // }
                    // for (size_t ind1 = 0; ind1 < both; ind1++) {
                    //     this->contents()[ ind1 + skip ]
                    //         = NCPA::arrays::add_vectors(
                    //             this->contents()[ ind1 + skip ],
                    //             NCPA::arrays::scale_vector(
                    //                 b->contents()[ ind1 + add_on ],
                    //                 modifier ) );
                    // }
                    // if (add_on > 0) {
                    //     this->contents().insert(
                    //         this->contents().begin(), b->contents().begin(),
                    //         b->contents().begin() + add_on );
                    //     for (size_t ind1 = 0; ind1 < add_on; ind1++) {
                    //         this->contents()[ ind1 ]
                    //             = NCPA::arrays::scale_vector(
                    //                 this->contents()[ ind1 ], modifier );
                    //     }
                    //     _n_lower += add_on;
                    // }

                    // // upper half
                    // skip   = 0;
                    // add_on = 0;
                    // both   = 0;
                    // if (_n_upper > b->_n_upper) {
                    //     both = b->_n_upper;
                    // } else if (_n_upper < b->_n_upper) {
                    //     both   = _n_upper;
                    //     add_on = b->_n_upper - _n_upper;
                    // } else {
                    //     both = _n_upper;
                    // }
                    // for (size_t ind1 = 0; ind1 < both; ind1++) {
                    //     this->contents()[ _n_lower + ind1 + 1 ]
                    //         = NCPA::arrays::add_vectors(
                    //             this->contents()[ _n_lower + ind1 + 1 ],
                    //             NCPA::arrays::scale_vector(
                    //                 b->contents()[ b->_n_lower + ind1 + 1 ],
                    //                 modifier ) );
                    // }

                    // // insert the last add_on diagonals from b, then negate
                    // // them
                    // if (add_on > 0) {
                    //     this->contents().insert( this->contents().end(),
                    //                              b->contents().end() -
                    //                              add_on, b->contents().end()
                    //                              );
                    //     _n_upper += add_on;
                    //     for (size_t ind1 = this->contents().size() - add_on;
                    //          ind1 < this->contents().size(); ind1++) {
                    //         this->contents()[ ind1 ]
                    //             = NCPA::arrays::scale_vector(
                    //                 this->contents()[ ind1 ], modifier );
                    //     }
                    // }

                    // // and the center diagonal
                    // this->contents()[ _n_lower ] =
                    // NCPA::arrays::add_vectors(
                    //     this->contents()[ _n_lower ],
                    //     NCPA::arrays::scale_vector(
                    //         b->contents()[ b->_n_lower ], modifier ) );
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& add(
                    ELEMENTTYPE b ) override {
                    for (size_t ind1 = 0; ind1 < this->contents().size();
                         ind1++) {
                        for (size_t ind2 = _min_ind2( ind1 );
                             ind2 < _max_ind2( ind1 ); ind2++) {
                            this->internal( ind1, ind2 ) += b;
                        }
                    }
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& as_array(
                    size_t& nrows, size_t& ncols,
                    ELEMENTTYPE **& vals ) override {
                    if (vals == nullptr) {
                        nrows = rows();
                        ncols = columns();
                        vals
                            = NCPA::arrays::zeros<ELEMENTTYPE>( nrows, ncols );
                    } else {
                        if (nrows == 0) {
                            nrows = rows();
                        } else if (nrows != rows()) {
                            throw std::invalid_argument(
                                "Wrong number of rows requested" );
                        }
                        if (ncols == 0) {
                            ncols = columns();
                        } else if (ncols != columns()) {
                            throw std::invalid_argument(
                                "Wrong number of columns requested" );
                        }
                        NCPA::arrays::fill( vals, nrows, ncols, _zero );
                    }
                    int row, col;
                    for (size_t ind1 = 0; ind1 < this->contents().size();
                         ind1++) {
                        for (size_t ind2 = _min_ind2( ind1 );
                             ind2 < _max_ind2( ind1 ); ind2++) {
                            // bool tf = internal2rowcol( ind1, ind2, row, col
                            // ); std::cout << "[ " << ind1 << "," << ind2
                            //           << " ] => [ " << row << "," << col
                            //           << " ] (" << tf << ")" << std::endl;
                            if (internal2rowcol( ind1, ind2, row, col )
                                == diagonal_index_status_t::VALID) {
                                vals[ row ][ col ]
                                    = this->internal( ind1, ind2 );
                            }
                        }
                    }
                    return *this;
                }

                virtual std::vector<size_t> band_column_indices(
                    size_t row ) const {
                    std::vector<size_t> inds;
                    for (int i = std::max( 0, (int)row - (int)_n_lower );
                         i < std::min( (int)columns(),
                                       (int)row + (int)_n_upper + 1 );
                         i++) {
                        inds.push_back( (size_t)i );
                    }
                    return inds;
                }

                virtual std::vector<size_t> band_row_indices( size_t col ) {
                    std::vector<size_t> inds;
                    for (int i = std::max( 0, (int)col - (int)_n_upper );
                         i < std::min( (int)rows(),
                                       (int)col + (int)_n_lower + 1 );
                         i++) {
                        inds.push_back( (size_t)i );
                    }
                    return inds;
                }

                virtual size_t bandwidth() const override {
                    std::vector<int> diags;
                    for (auto it = _contents.cbegin(); it != _contents.cend();
                         ++it) {
                        if (!it->empty()) {
                            diags.push_back(
                                (int)std::distance( _contents.cbegin(), it )
                                - (int)_n_lower );
                        }
                    }
                    if (diags.empty()) {
                        return 0;
                    } else {
                        return diags.back() - diags.front() + 1;
                    }
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    build_vector( size_t n = 0 ) const override {
                    return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                        new sparse_vector<ELEMENTTYPE>( n ) );
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& clear() override {
                    _nrows   = 0;
                    _ncols   = 0;
                    _n_lower = 0;
                    _n_upper = 0;
                    this->contents().clear();
                    this->contents().resize( 1 );
                    return *this;
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> clone()
                    const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new band_diagonal_matrix( *this ) );
                }

                virtual size_t columns() const override { return _ncols; }

                virtual std::vector<std::vector<ELEMENTTYPE>>& contents() {
                    return _contents;
                }

                virtual const std::vector<std::vector<ELEMENTTYPE>>& contents()
                    const {
                    return _contents;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& copy(
                    const abstract_matrix<ELEMENTTYPE>& other ) override {
                    this->clear().resize( other.rows(), other.columns() );
                    if (other.is_zero()) {
                        return *this;
                    }
                    std::vector<int> otherdiags = other.diagonals();
                    while (
                        other.get_diagonal( otherdiags.front() )->is_zero()) {
                        otherdiags.erase( otherdiags.begin() );
                    }
                    while (
                        other.get_diagonal( otherdiags.back() )->is_zero()) {
                        otherdiags.erase( otherdiags.begin()
                                          + otherdiags.size() - 1 );
                    }
                    for (auto it = otherdiags.begin(); it != otherdiags.end();
                         ++it) {
                        auto otherd = other.get_diagonal( *it );
                        if (!otherd->is_zero()) {
                            this->set_diagonal( otherd->as_std(), *it );
                        }
                    }
                    return *this;
                }

                virtual bool diag2internal( int diag, int& ind1 ) const {
                    ind1 = diag + (int)_n_lower;
                    return ( ind1 >= 0
                             && ind1 < (int)( this->contents().size() ) );
                }

                virtual std::vector<int> diagonals() const override {
                    std::vector<int> diag;
                    if (this->is_zero()) {
                        return diag;
                    }
                    for (size_t ii = 0; ii < _contents.size(); ++ii) {
                        if (!_contents.at( ii ).empty()) {
                            diag.push_back( (int)ii - _n_lower );
                        }
                    }
                    return diag;
                }

                const band_diagonal_matrix<ELEMENTTYPE> *downcast(
                    const abstract_matrix<ELEMENTTYPE> *in ) const {
                    return dynamic_cast<
                        const band_diagonal_matrix<ELEMENTTYPE> *>( in );
                }

                virtual bool equals( const abstract_matrix<ELEMENTTYPE>&
                                         other ) const override {
                    if (!is_this_subclass( other )) {
                        return abstract_matrix<ELEMENTTYPE>::equals( other );
                    }
                    const band_diagonal_matrix<ELEMENTTYPE> *b
                        = downcast( &other );

                    bool contents_equal
                        = ( this->contents() == b->contents() );
                    return ( rows() == b->rows() && columns() == b->columns()
                             && this->lower_bandwidth() == b->lower_bandwidth()
                             && this->upper_bandwidth() == b->upper_bandwidth()
                             && this->contents() == b->contents() );
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                    fresh_clone() const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new band_diagonal_matrix() );
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const override {
                    // this->check_size( row, col );
                    int ind1, ind2;
                    if (this->rowcol2internal( row, col, ind1, ind2 )
                        == diagonal_index_status_t::VALID) {
                        // && this->contents().at( ind1 ).size() > ind2) {
                        return this->internal( ind1, ind2 );
                    } else {
                        return _zero;
                    }
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    get_column( size_t column ) const override {
                    std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                        new sparse_vector<ELEMENTTYPE>( this->rows() ) );
                    if (this->is_zero()) {
                        return v;
                    }
                    int ind1, ind2;
                    for (size_t r = 0; r < this->rows(); r++) {
                        if (this->rowcol2internal( r, column, ind1, ind2 )
                            == diagonal_index_status_t::VALID) {
                            v->set( r, this->internal( ind1, ind2 ) );
                        }
                    }
                    return v;
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    get_diagonal( int offset = 0 ) const override {
                    std::vector<ELEMENTTYPE> diag(
                        this->diagonal_size( offset ), _zero );
                    if (this->has_diagonal( offset )) {
                        int ind1 = (int)_n_lower + offset;
                        diag.assign( this->internal( ind1 ).begin(),
                                     //  + _min_ind2( ind1 ),
                                     this->internal( ind1 ).begin()
                                         + _max_ind2( ind1 ) );
                    }
                    // int ind1 = (int)_n_lower + offset;
                    // if (ind1 >= 0 && ind1 < this->contents().size()) {
                    //     diag.assign( this->contents()[ ind1 ].begin()
                    //                      + _min_ind2( ind1 ),
                    //                  this->contents()[ ind1 ].begin()
                    //                      + _max_ind2( ind1 ) );
                    // }
                    return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                        new dense_vector<ELEMENTTYPE>( diag ) );
                }

                // @todo make read-write vector view for columns
                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>> get_row(
                    size_t row ) const override {
                    std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                        new sparse_vector<ELEMENTTYPE>( this->columns() ) );
                    for (int c = std::max( (int)row - (int)_n_lower, 0 );
                         c <= std::min( (int)row + (int)_n_upper,
                                        (int)this->columns() - 1 );
                         ++c) {
                        v->set( c, this->get( row, c ) );
                    }
                    // int r, c;
                    // for (size_t ind1 = 0; ind1 < this->contents().size();
                    //      ind1++) {
                    //     if (internal2rowcol( ind1, row, r, c )) {
                    //         v->set( c, this->contents()[ ind1 ][ row ] );
                    //     }
                    // }
                    return v;
                }

                virtual bool has_diagonal( int diag ) const {
                    int internal;
                    return ( this->diag2internal( diag, internal )
                             && !this->contents().at( internal ).empty() );
                }

                virtual std::string id() const override {
                    return "NCPA band-diagonal matrix";
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& identity(
                    size_t nrows, size_t ncols ) override {
                    this->clear().resize( nrows, ncols );
                    this->internal( 0 ) = std::vector<ELEMENTTYPE>(
                        this->diagonal_size( 0 ),
                        NCPA::math::one<ELEMENTTYPE>() );
                    return *this;
                }

                virtual std::vector<ELEMENTTYPE>& internal( size_t ind1 ) {
                    return this->contents().at( ind1 );
                }

                virtual const std::vector<ELEMENTTYPE>& internal(
                    size_t ind1 ) const {
                    return this->contents().at( ind1 );
                }

                virtual ELEMENTTYPE& internal( size_t ind1, size_t ind2 ) {
                    // return this->internal( ind1 ).at( ind2 );
                    return this->internal( ind1 )[ ind2 ];
                }

                virtual const ELEMENTTYPE& internal( size_t ind1,
                                                     size_t ind2 ) const {
                    return this->internal( ind1 ).at( ind2 );
                }

                virtual diagonal_index_status_t internal2rowcol(
                    const size_t& ind1, const size_t& ind2, int& row,
                    int& col ) const {
                    int diag = -(int)_n_lower + (int)ind1;
                    row      = (int)ind2 - std::min( diag, 0 );
                    col      = (int)ind2 + std::max( diag, 0 );
                    // col      = (int)ind1 + (int)ind2 - (int)_n_lower;
                    // row      = (int)ind2;
                    // return ( col >= 0 && col < this->columns()
                    //          && row < this->rows() );
                    if (ind1 < 0 || ind1 >= (int)( this->contents().size() )
                        || row >= this->rows() || col >= this->columns()) {
                        return diagonal_index_status_t::INVALID;
                    } else if (ind2 >= this->contents().at( ind1 ).size()) {
                        return diagonal_index_status_t::VALID_BUT_NOT_DEFINED;
                    } else {
                        return diagonal_index_status_t::VALID;
                    }
                }

                virtual bool is_band_diagonal() const override { return true; }

                virtual bool is_diagonal() const override {
                    return _n_lower == 0 && _n_upper == 0;
                }

                virtual bool is_identity() const override {
                    if (this->is_zero() || !this->is_square()
                        || !this->is_diagonal()) {
                        return false;
                    }
                    for (auto it = this->internal( 0 ).cbegin();
                         it != this->internal( 0 ).cend(); ++it) {
                        if (!NCPA::math::within(
                                *it, NCPA::math::one<ELEMENTTYPE>(),
                                1.0e-12 )) {
                            return false;
                        }
                    }
                    return true;
                }

                virtual bool is_lower_triangular() const override {
                    if (this->is_empty()) {
                        return true;
                    }
                    if (_n_upper == 0) {
                        return true;
                    }
                    for (size_t i = 0; i < _n_upper; i++) {
                        for (auto it = this->contents()
                                           .at( _n_lower + i + 1 )
                                           .cbegin();
                             it != this->internal( _n_lower + i + 1 ).cend();
                             ++it) {
                            if (*it != _zero) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_this_subclass(
                    const abstract_matrix<ELEMENTTYPE>& b ) const override {
                    if (auto *derived = dynamic_cast<
                            const band_diagonal_matrix<ELEMENTTYPE> *>( &b )) {
                        return true;
                    } else {
                        return false;
                    }
                }

                virtual bool is_tridiagonal() const override {
                    return _n_lower <= 1 && _n_upper <= 1;
                }

                virtual bool is_upper_triangular() const override {
                    if (this->is_empty()) {
                        return true;
                    }
                    if (_n_lower == 0) {
                        return true;
                    }
                    for (size_t i = 0; i < _n_lower; i++) {
                        for (auto it = this->internal( i ).cbegin();
                             it != this->internal( i ).cend(); ++it) {
                            if (*it != _zero) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_zero( double tol = 1.0e-12 ) const override {
                    if (this->is_empty()) {
                        return true;
                    }
                    for (auto it1 = _contents.cbegin();
                         it1 != _contents.cend(); ++it1) {
                        for (auto it2 = it1->cbegin(); it2 != it1->cend();
                             ++it2) {
                            if (std::abs( *it2 ) > tol) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_zero( size_t r, size_t c,
                                      double tol = 1.0e-12 ) const override {
                    if (this->is_empty()) {
                        return true;
                    }
                    int diag       = 0;
                    size_t element = 0;
                    rc2diag( r, c, diag, element );
                    return !( this->has_diagonal( diag )
                              && std::abs( this->get( r, c ) ) > tol );
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    left_multiply( const abstract_vector<ELEMENTTYPE>& x )
                        const override {
                    if (this->rows() != x.size()) {
                        std::ostringstream oss;
                        oss << "Size mismatch in vector-matrix "
                               "multiplication: "
                            << x.size() << " elements in vector vs "
                            << columns() << " rows in matrix";
                        throw std::invalid_argument( oss.str() );
                    }
                    std::unique_ptr<abstract_vector<ELEMENTTYPE>> b(
                        new dense_vector<ELEMENTTYPE>( this->columns() ) );
                    int n = (int)columns();
                    for (int i = 0; i < n; i++) {
                        int minloop
                            = std::max( i - (int)this->upper_bandwidth(), 0 ),
                            maxloop
                            = std::min( i + this->lower_bandwidth() + 1,
                                        this->rows() );
                        ELEMENTTYPE bval = _zero;
                        for (int k = minloop; k < maxloop; k++) {
                            bval += x.get( k ) * this->get( k, i );
                        }
                        b->set( i, bval );
                        // int k            = i - (int)_n_upper;
                        // int tmploop      = std::min( bw, n - k );
                        // ELEMENTTYPE bval = _zero;
                        // for ( int j = std::max( 0, -k ); j < tmploop; j++ )
                        // {
                        //     bval += this->contents()[ j ][ i ] * x.get( j +
                        //     k );
                        // }
                        // b->set( i, bval );
                    }
                    return b;
                }

                virtual size_t lower_bandwidth() const override {
                    int lband = (int)_n_lower;
                    int i     = 0;
                    while (i < _n_lower && _contents.at( i++ ).empty()) {
                        lband--;
                    }
                    return lband;
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> multiply(
                    const abstract_matrix<ELEMENTTYPE>& other )
                    const override {
                    if (!is_this_subclass( other )) {
                        return abstract_matrix<ELEMENTTYPE>::multiply( other );
                    }
                    const band_diagonal_matrix<ELEMENTTYPE> *b
                        = downcast( &other );

                    // size_t new_n_lower
                    //     = this->lower_bandwidth() + b->lower_bandwidth(),
                    //     new_n_upper
                    //     = this->upper_bandwidth() + b->upper_bandwidth();
                    std::unique_ptr<abstract_matrix<ELEMENTTYPE>> product(
                        new band_diagonal_matrix<ELEMENTTYPE>(
                            this->rows(), b->columns(), 0, 0 ) );
                    band_diagonal_matrix<ELEMENTTYPE> *bproduct
                        = dynamic_cast<band_diagonal_matrix<ELEMENTTYPE> *>(
                            product.get() );
                    band_diagonal_matrix<
                        ELEMENTTYPE>::_do_band_diagonal_multiply( *this, *b,
                                                                  *bproduct );
                    // int prows = (int)( product->rows() );
                    // int pcols = (int)( product->columns() );
                    // for ( int r = 0; r < prows; r++ ) {
                    //     for ( int c = std::max( 0, r - (int)new_n_lower );
                    //           c < std::min( (int)( product->columns() ),
                    //                         r + (int)new_n_upper + 1 );
                    //           c++ ) {
                    //         ELEMENTTYPE val = _zero;
                    //         int kmin        = std::max(
                    //             std::max( 0, r -
                    //             (int)this->lower_bandwidth() ), std::max( 0,
                    //             c - (int)( b->upper_bandwidth() ) ) );
                    //         int kmax = std::min(
                    //             std::min( (int)columns(),
                    //                       r + (int)this->upper_bandwidth() +
                    //                       1 ),
                    //             std::min( (int)( b->rows() ),
                    //                       c + (int)( b->lower_bandwidth() )
                    //                       + 1 ) );
                    //         for ( int k = kmin; k < kmax; k++ ) {
                    //             val += this->get( r, k ) * b->get( k, c );
                    //         }
                    //         product->set( (size_t)r, (size_t)c, val );
                    //     }
                    // }
                    return product;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& resize(
                    size_t r, size_t c ) override {
                    if (this->contents().size() == 0) {
                        // starting from scratch
                        _nrows   = r;
                        _ncols   = c;
                        _n_lower = 0;
                        _n_upper = 0;
                        _prep_diagonals();
                    } else {
                        // did the diagonal size increase or decrease?
                        int ddiag = (int)std::min( r, c )
                                  - (int)std::min( rows(), columns() );
                        if (ddiag > 0) {
                            // increase each diagonal by ddiag elements
                            for (auto it = this->contents().begin();
                                 it != this->contents().end(); ++it) {
                                if (it->size() > 0) {
                                    it->insert( it->cend(), (size_t)ddiag,
                                                _zero );
                                }
                            }
                        } else if (ddiag < 0) {
                            for (auto it = this->contents().begin();
                                 it != this->contents().end(); ++it) {
                                if (it->size() > 0) {
                                    it->erase( it->cend() - size_t( -ddiag ),
                                               it->cend() );
                                }
                            }
                        }
                        _nrows = r;
                        _ncols = c;
                    }

                    return *this;
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    right_multiply( const abstract_vector<ELEMENTTYPE>& x )
                        const override {
                    if (columns() != x.size()) {
                        std::ostringstream oss;
                        oss << "Size mismatch in matrix-vector "
                               "multiplication: "
                            << columns() << " columns in matrix vs "
                            << x.size() << " elements in vector";
                        throw std::invalid_argument( oss.str() );
                    }
                    std::unique_ptr<abstract_vector<ELEMENTTYPE>> b(
                        new dense_vector<ELEMENTTYPE>( rows() ) );
                    int n = (int)rows();  // , bw = (int)bandwidth();
                    for (int i = 0; i < n; i++) {
                        int minloop
                            = std::max( 0, i - (int)this->lower_bandwidth() ),
                            maxloop
                            = std::min( i + (int)this->upper_bandwidth() + 1,
                                        (int)this->columns() );
                        ELEMENTTYPE bval = _zero;
                        for (int k = minloop; k < maxloop; k++) {
                            bval += this->get( i, k ) * x.get( k );
                        }
                        b->set( i, bval );
                    }
                    // for ( int i = 0; i < n; i++ ) {
                    //     int k            = i - (int)_n_lower;
                    //     int tmploop      = std::min( bw, n - k );
                    //     ELEMENTTYPE bval = _zero;
                    //     for ( int j = std::max( 0, -k ); j < tmploop; j++ )
                    //     {
                    //         bval += this->contents()[ j ][ i ] * x.get( j +
                    //         k );
                    //     }
                    //     b->set( i, bval );
                    // }
                    return b;
                }

                virtual bool row_in_range( size_t row ) const {
                    return ( row < this->rows() );
                }

                virtual diagonal_index_status_t rowcol2internal(
                    const size_t& row, const size_t& col, int& ind1,
                    int& ind2 ) const {
                    int diag = (int)col - (int)row;
                    ind1     = diag + (int)_n_lower;
                    ind2     = (int)row + std::min( 0, diag );
                    if (ind1 < 0 || ind1 >= (int)( this->contents().size() )
                        || row >= this->rows() || col >= this->columns()) {
                        return diagonal_index_status_t::INVALID;
                    } else if (ind2 >= this->internal( ind1 ).size()) {
                        return diagonal_index_status_t::VALID_BUT_NOT_DEFINED;
                    } else {
                        return diagonal_index_status_t::VALID;
                    }
                    // return ( ind1 >= 0
                    //          && ind1 < (int)( this->contents().size() )
                    //          && row < this->rows() && col < this->columns()
                    //          );
                }

                virtual size_t rows() const override { return _nrows; }

                virtual band_diagonal_matrix<ELEMENTTYPE>& scale(
                    const abstract_matrix<ELEMENTTYPE>& b ) override {
                    this->check_size( b );
                    int r, c;
                    for (size_t ind1 = 0; ind1 < this->contents().size();
                         ind1++) {
                        for (size_t ind2
                             = 0;  // size_t ind2 = _min_ind2( ind1 );
                             ind2 < this->internal( ind1 ).size(); ind2++) {
                            if (internal2rowcol( ind1, ind2, r, c )
                                == diagonal_index_status_t::VALID) {
                                this->internal( ind1, ind2 ) *= b.get( r, c );
                            }
                        }
                    }
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& scale(
                    ELEMENTTYPE val ) override {
                    for (auto it1 = this->contents().begin();
                         it1 != this->contents().end(); ++it1) {
                        for (auto it2 = it1->begin(); it2 != it1->end();
                             ++it2) {
                            *it2 *= val;
                        }
                    }
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set_internal(
                    size_t ind1, size_t ind2, ELEMENTTYPE val ) {
                    if (ind1 >= this->contents().size()) {
                        std::ostringstream oss;
                        oss << "Requested internal index " << ind1
                            << " out of range for bandwidth "
                            << this->bandwidth();
                        throw std::range_error( oss.str() );
                    }
                    if (ind2 >= this->internal( ind1 ).size()) {
                        this->internal( ind1 ).resize( ind2 + 1 );
                    }
                    this->internal( ind1, ind2 ) = val;
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set(
                    size_t row, size_t col, ELEMENTTYPE val ) override {
                    this->check_size( row, col );
                    int ind1, ind2;
                    while (rowcol2internal( row, col, ind1, ind2 )
                           == diagonal_index_status_t::INVALID) {
                        if (ind1 < 0) {
                            this->_add_subdiagonal( (size_t)( -ind1 ) );
                        } else {
                            this->_add_superdiagonal(
                                (size_t)( ind1 - _contents.size() + 1 ) );
                        }
                    }
                    if (this->internal( ind1 ).empty()) {
                        this->internal( ind1 ).resize(
                            this->diagonal_size( (int)col - (int)row ) );
                    }
                    this->internal( ind1, ind2 ) = val;
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set(
                    ELEMENTTYPE val ) override {
                    for (auto it1 = this->contents().begin();
                         it1 != this->contents().end(); ++it1) {
                        it1->assign( this->diagonal_size( 0 ), val );
                    }
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set_column(
                    size_t column, size_t nvals, const size_t *row_inds,
                    const ELEMENTTYPE *vals ) override {
                    for (size_t i = 0; i < nvals; i++) {
                        this->set_safe( row_inds[ i ], column, vals[ i ] );
                    }
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    if (vals.size() == rows()) {
                        return dynamic_cast<
                            band_diagonal_matrix<ELEMENTTYPE>&>(
                            this->set_column(
                                col,
                                NCPA::arrays::index_vector<size_t>(
                                    vals.size() ),
                                vals ) );
                    }
                    if (vals.size() == bandwidth()) {
                        return dynamic_cast<
                            band_diagonal_matrix<ELEMENTTYPE>&>(
                            this->set_column(
                                col, this->band_row_indices( col ), vals ) );
                    }
                    throw std::invalid_argument(
                        "Value vector size matches neither number of rows nor "
                        "bandwidth" );
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set_column(
                    size_t col, ELEMENTTYPE val ) override {
                    return this->set_column(
                        col, std::vector<ELEMENTTYPE>( bandwidth(), val ) );
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set_diagonal(
                    size_t nvals, const ELEMENTTYPE *vals,
                    int offset = 0 ) override {
                    if (nvals > this->diagonal_size( offset )) {
                        throw std::out_of_range(
                            "Too many values for requested diagonal" );
                    }
                    while (offset > (int)_n_upper) {
                        this->_add_superdiagonal( offset - (int)_n_upper );
                    }
                    while (-offset > (int)_n_lower) {
                        this->_add_subdiagonal( -offset - (int)_n_lower );
                    }
                    int ind1 = (int)_n_lower + offset;
                    this->internal( ind1 ).resize(
                        this->diagonal_size( offset ) );
                    for (size_t i = 0; i < nvals; i++) {
                        this->internal( ind1, i ) = vals[ i ];
                    }

                    // if (offset > 0) {
                    //     while (offset > _n_upper) {
                    //         this->_add_superdiagonal();
                    //     }
                    //     int ind1 = _n_lower + offset;
                    //     this->contents().at( ind1 ).resize(
                    //         this->diagonal_size( offset ) );
                    //     for (size_t i = 0; i < nvals; i++) {
                    //         this->internal( ind1,  i ) = vals[ i ];
                    //     }
                    // } else if (offset < 0) {
                    //     while (-offset > (int)_n_lower) {
                    //         this->_add_subdiagonal();
                    //     }
                    //     int ind1 = _n_lower + offset;
                    //     this->contents().at( ind1 ).resize(
                    //         this->diagonal_size( offset ) );
                    //     for (size_t i = 0; i < nvals; i++) {
                    //         int ind2 = _min_ind2( ind1 ) + i;
                    //         this->internal( ind1,  ind2 ) = vals[ i
                    //         ];
                    //     }
                    // } else {
                    //     this->contents()
                    //         .at( _n_lower )
                    //         .resize( this->diagonal_size( 0 ) );
                    //     for (size_t i = 0; i < nvals; i++) {
                    //         this->contents().at( _n_lower ).at( i )
                    //             = vals[ i ];
                    //     }
                    // }
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set_row(
                    size_t row, size_t nvals, const size_t *column_inds,
                    const ELEMENTTYPE *vals ) override {
                    for (size_t i = 0; i < nvals; i++) {
                        this->set_safe( row, column_inds[ i ], vals[ i ] );
                    }
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    if (vals.size() == columns()) {
                        return dynamic_cast<
                            band_diagonal_matrix<ELEMENTTYPE>&>( this->set_row(
                            row,
                            NCPA::arrays::index_vector<size_t>( vals.size() ),
                            vals ) );
                    }
                    if (vals.size() == bandwidth()) {
                        return dynamic_cast<
                            band_diagonal_matrix<ELEMENTTYPE>&>( this->set_row(
                            row, band_column_indices( row ), vals ) );
                    }
                    throw std::invalid_argument(
                        "Value vector size matches neither number of columns "
                        "nor bandwidth" );
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set_row(
                    size_t row, ELEMENTTYPE val ) override {
                    return set_row(
                        row, std::vector<ELEMENTTYPE>( bandwidth(), val ) );
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& set_safe(
                    size_t row, size_t col, ELEMENTTYPE val ) {
                    this->check_size( row, col );
                    int ind1, ind2;
                    if (rowcol2internal( row, col, ind1, ind2 )
                        != diagonal_index_status_t::INVALID) {
                        return this->set( row, col, val );
                        // this->internal( ind1,  row ) = val;
                    } else {
                        std::ostringstream oss;
                        oss << "Element [" << row << ", " << col
                            << "] is out of band.";
                        throw std::out_of_range( oss.str() );
                    }
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& swap_columns(
                    size_t ind1, size_t ind2 ) override {
                    throw std::logic_error(
                        "Cannot swap columns of a band-diagonal matrix!" );
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& swap_rows(
                    size_t ind1, size_t ind2 ) override {
                    throw std::logic_error(
                        "Cannot swap rows of a band-diagonal matrix!" );
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& transpose()
                    override {
                    std::vector<std::vector<ELEMENTTYPE>> newcontents
                        = this->contents();
                    std::reverse( newcontents.begin(), newcontents.end() );
                    std::swap( _n_lower, _n_upper );
                    std::swap( _nrows, _ncols );
                    // if (_n_lower > 0) {
                    //     // take the zeros from the end and put them up front
                    //     for (size_t ind1 = 0; ind1 < _n_lower; ind1++) {
                    //         if (!this->contents().at( ind1 ).empty()) {
                    //             size_t invalid = _n_lower - ind1;

                    //             for (size_t i = 0; i < invalid; i++) {
                    //                 newcontents[ ind1 ].pop_back();
                    //                 newcontents[ ind1 ].insert(
                    //                     newcontents[ ind1 ].begin(), 1,
                    //                     _zero );
                    //             }
                    //         }
                    //     }
                    // }
                    // if (_n_upper > 0) {
                    //     // take the zeros from the front and put them at the
                    //     // end
                    //     for (size_t ind1 = _n_lower + 1;
                    //          ind1 < newcontents.size(); ind1++) {
                    //         if (!this->contents().at( ind1 ).empty()) {
                    //             size_t invalid = ind1 - _n_lower;
                    //             for (size_t i = 0; i < invalid; i++) {
                    //                 newcontents[ ind1 ].erase(
                    //                     newcontents[ ind1 ].cbegin() );
                    //                 newcontents[ ind1 ].push_back( _zero );
                    //             }
                    //         }
                    //     }
                    // }
                    _contents = newcontents;
                    return *this;
                }

                virtual size_t upper_bandwidth() const override {
                    int uband = (int)_n_upper;
                    int i     = 0;
                    auto rit  = _contents.crbegin();
                    while (i < _n_upper && ( rit + ( i++ ) )->empty()) {
                        uband--;
                    }
                    return uband;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& zero() override {
                    _n_lower = 0;
                    _n_upper = 0;
                    this->contents().clear();
                    this->contents().resize( 1 );
                    return *this;
                }

                virtual band_diagonal_matrix<ELEMENTTYPE>& zero(
                    size_t row, size_t col ) override {
                    return this->set( row, col, _zero );
                }


            protected:
                size_t _nrows, _ncols, _n_lower, _n_upper;

                std::vector<std::vector<ELEMENTTYPE>> _contents;

                const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();

                void _prep_diagonals() {
                    this->contents().clear();
                    this->contents().resize( _n_lower + _n_upper + 1 );
                }

                void _add_subdiagonal( size_t n = 1 ) {
                    if (n == 0) return;
                    bool was_empty = this->contents().empty();
                    this->contents().insert( this->contents().cbegin(), n,
                                             std::vector<ELEMENTTYPE>() );
                    if (was_empty) {
                        // one of the things set was the main diagonal so don't
                        // count that one
                        --n;
                    }
                    _n_lower += n;
                }

                void _add_superdiagonal( size_t n = 1 ) {
                    if (n == 0) return;
                    bool was_empty = this->contents().empty();
                    this->contents().insert( this->contents().end(), n,
                                             std::vector<ELEMENTTYPE>() );
                    if (was_empty) {
                        // one of the things set was the main diagonal so don't
                        // count that one
                        --n;
                    }
                    _n_upper += n;
                }

                size_t _min_ind2( size_t ind1 ) const {
                    return 0;
                    // return (size_t)std::max( (int)_n_lower - (int)ind1, 0 );
                }

                size_t _max_ind2( size_t ind1 ) const {
                    return this->diagonal_size( (int)ind1 - (int)_n_lower );
                    // return _min_ind2( ind1 ) + this->diagonal_size( ind1 ) -
                    // 1; return (size_t)std::min( this->rows(),
                    //                          this->rows() + _n_lower - ind1
                    //                          );
                }

                // virtual band_diagonal_matrix<ELEMENTTYPE>& _add(
                //     const band_diagonal_matrix<ELEMENTTYPE>& b,
                //     ELEMENTTYPE modifier = 1.0 ) {
                //     this->check_size( b );
                //     size_t new_n_lower = std::max( _n_lower, b._n_lower );
                //     size_t new_n_upper = std::max( _n_upper, b._n_upper );
                //     std::vector<std::vector<ELEMENTTYPE>> newcontents(
                //         new_n_lower + new_n_upper + 1 );

                //     // diagonals first
                //     newcontents[ new_n_lower ] = NCPA::arrays::add_vectors(
                //         this->contents()[ _n_lower ],
                //         NCPA::arrays::scale_vector( b.contents()[ b._n_lower
                //         ],
                //                                     modifier ) );

                //     // subdiagonals
                //     for (size_t n = 1; n < new_n_lower; n++) {
                //         if (n <= _n_lower && n <= b._n_lower) {
                //             newcontents[ new_n_lower - n ]
                //                 = NCPA::arrays::add_vectors(
                //                     this->contents()[ _n_lower - n ],
                //                     NCPA::arrays::scale_vector(
                //                         b.contents()[ b._n_lower - n ],
                //                         modifier ) );
                //         } else if (n <= _n_lower) {
                //             newcontents[ new_n_lower - n ]
                //                 = this->contents()[ _n_lower - n ];
                //         } else if (n <= b._n_lower) {
                //             newcontents[ new_n_lower - n ]
                //                 = NCPA::arrays::scale_vector(
                //                     b.contents()[ b._n_lower - n ], modifier
                //                     );
                //         }
                //     }

                //     // superdiagonals
                //     for (size_t n = 1; n < new_n_upper; n++) {
                //         if (n <= _n_upper && n <= b._n_upper) {
                //             newcontents[ new_n_lower + n ]
                //                 = NCPA::arrays::add_vectors(
                //                     this->contents()[ _n_lower + n ],
                //                     NCPA::arrays::scale_vector(
                //                         b.contents()[ b._n_lower + n ],
                //                         modifier ) );
                //         } else if (n <= _n_upper) {
                //             newcontents[ new_n_lower + n ]
                //                 = this->contents()[ _n_lower + n ];
                //         } else if (n <= b._n_upper) {
                //             newcontents[ new_n_lower + n ]
                //                 = NCPA::arrays::scale_vector(
                //                     b.contents()[ b._n_lower + n ], modifier
                //                     );
                //         }
                //     }
                //     _contents = newcontents;
                //     _n_lower  = new_n_lower;
                //     _n_upper  = new_n_upper;
                //     return *this;
                //     ;
                // }

                static void _do_band_diagonal_multiply(
                    const band_diagonal_matrix<ELEMENTTYPE>& a,
                    const band_diagonal_matrix<ELEMENTTYPE>& b,
                    band_diagonal_matrix<ELEMENTTYPE>& product ) {
                    int prows      = (int)( product.rows() );
                    int pcols      = (int)( product.columns() );
                    // int new_n_lower = product.lower_bandwidth();
                    // int new_n_upper = product.upper_bandwidth();
                    size_t counter = 0;
                    for (int r = 0; r < (int)product.rows(); ++r) {
                        for (int c = 0; c < (int)product.columns(); ++c) {
                            int min_k_a
                                = std::max( r - (int)a.lower_bandwidth(), 0 );
                            int max_k_a
                                = std::min( r + (int)a.upper_bandwidth(),
                                            (int)a.columns() - 1 );
                            int min_k_b
                                = std::min( c - (int)b.upper_bandwidth(), 0 );
                            int max_k_b
                                = std::max( c + (int)b.lower_bandwidth(),
                                            (int)b.rows() - 1 );
                            int min_k = std::max( min_k_a, min_k_b );
                            int max_k = std::min( max_k_a, max_k_b );

                            // for (int k = 0; k < a.columns(); ++k) {
                            for (int k = min_k; k <= max_k; ++k) {
                                if (a.has_diagonal( k - r )
                                    && b.has_diagonal( c - k )) {
                                    product.set( r, c,
                                                 product.get( r, c )
                                                     + a.get( r, k )
                                                           * b.get( k, c ) );
                                }
                            }
                        }
                    }
                    return;


                    // for (int r = 0; r < prows; r++) {
                    //     for (int c = std::max( 0, r - (int)new_n_lower );
                    //          c < std::min( pcols, r + (int)new_n_upper + 1
                    //          ); c++) {
                    //         ELEMENTTYPE val = a._zero;
                    //         int kmin        = std::max(
                    //             std::max( 0,
                    //                              r - (int)(
                    //                              a.lower_bandwidth() ) ),
                    //             std::max( 0,
                    //                              c - (int)(
                    //                              b.upper_bandwidth() ) ) );
                    //         int kmax = std::min(
                    //             std::min( (int)( a.columns() ),
                    //                       r + (int)a.upper_bandwidth() + 1
                    //                       ),
                    //             std::min( (int)( b.rows() ),
                    //                       c + (int)( b.lower_bandwidth() )
                    //                           + 1 ) );
                    //         for (int k = kmin; k < kmax; k++) {
                    //             val += a.get( r, k ) * b.get( k, c );
                    //             counter++;
                    //         }
                    //         product.set( (size_t)r, (size_t)c, val );
                    //     }
                    // }
                    // std::cout << counter << " multiplications" << std::endl;
                }
        };

    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::band_diagonal_matrix<T>& a,
                  NCPA::linear::band_diagonal_matrix<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::abstract_matrix<T>&>( a ),
            static_cast<NCPA::linear::abstract_matrix<T>&>( b ) );
    swap( a._nrows, b._nrows );
    swap( a._ncols, b._ncols );
    swap( a._n_lower, b._n_lower );
    swap( a._n_upper, b._n_upper );
    swap( a._contents, b._contents );
}
