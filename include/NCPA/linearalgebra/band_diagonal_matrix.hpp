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

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> clone()
                    const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new band_diagonal_matrix( *this ) );
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                    fresh_clone() const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new band_diagonal_matrix() );
                }

                virtual abstract_matrix<ELEMENTTYPE>& copy(
                    const abstract_matrix<ELEMENTTYPE>& other ) override {
                    resize( other.rows(), other.columns() );
                    set_diagonal( *other.get_diagonal() );
                    int countdown = other.max_off_diagonal();
                    int ndiag     = 1;
                    // std::cout << "Getting off-diagonal " << -ndiag <<
                    // std::endl;
                    auto diag     = other.get_diagonal( -ndiag );
                    while (ndiag <= countdown && !diag->is_zero()) {
                        // std::cout << "Off-diagonal " << -ndiag << " is
                        // nonzero" << std::endl;
                        set_diagonal( *diag, -ndiag );
                        ndiag++;
                        // std::cout << "Getting off-diagonal " << -ndiag <<
                        // std::endl;
                        diag = other.get_diagonal( -ndiag );
                    }
                    ndiag = 1;
                    //  std::cout << "Getting off-diagonal " << ndiag <<
                    //  std::endl;
                    diag  = other.get_diagonal( ndiag );
                    while (ndiag <= countdown && !diag->is_zero()) {
                        // std::cout << "Off-diagonal " << ndiag << " is
                        // nonzero" << std::endl;
                        set_diagonal( *diag, ndiag );
                        ndiag++;
                        // std::cout << "Getting off-diagonal " << ndiag <<
                        // std::endl;
                        diag = other.get_diagonal( ndiag );
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                const band_diagonal_matrix<ELEMENTTYPE> *downcast(
                    const abstract_matrix<ELEMENTTYPE> *in ) const {
                    return dynamic_cast<
                        const band_diagonal_matrix<ELEMENTTYPE> *>( in );
                }

                // virtual abstract_matrix<ELEMENTTYPE>& upcast()
                //     override {
                //     return *static_cast<abstract_matrix<ELEMENTTYPE> *>(
                //         this );
                // }

                virtual std::string id() const override {
                    return "NCPA band-diagonal matrix";
                }

                virtual size_t rows() const override { return _nrows; }

                virtual size_t columns() const override { return _ncols; }

                virtual abstract_matrix<ELEMENTTYPE>& clear() override {
                    _nrows   = 0;
                    _ncols   = 0;
                    _n_lower = 0;
                    _n_upper = 0;
                    this->contents().clear();
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual bool rowcol2internal( const size_t& row,
                                              const size_t& col, int& ind1,
                                              int& ind2 ) const {
                    ind1 = (int)col - (int)row + (int)_n_lower;
                    ind2 = (int)row;
                    return ( ind1 >= 0
                             && ind1 < (int)( this->contents().size() )
                             && row < this->rows() && col < this->columns() );
                }

                virtual bool internal2rowcol( const size_t& ind1,
                                              const size_t& ind2, int& row,
                                              int& col ) const {
                    col = (int)ind1 + (int)ind2 - (int)_n_lower;
                    row = (int)ind2;
                    return ( col >= 0 && col < this->columns()
                             && row < this->rows() );
                    // if ( icol < 0 || icol <= columns() ) {
                    //     std::ostringstream oss;
                    //     oss << "Attempted to retrieve invalid internal value
                    //     "
                    //            "at ind1="
                    //         << ind1 << ", ind2=" << ind2
                    //         << " (_n_lower = " << _n_lower << ")";
                    //     throw std::out_of_range( oss.str() );
                    // }
                    // return (size_t)icol;
                }

                virtual bool row_in_range( size_t row ) const {
                    return ( row >= 0 && row < this->rows() );
                }

                virtual abstract_matrix<ELEMENTTYPE>& set(
                    size_t row, size_t col, ELEMENTTYPE val ) override {
                    this->check_size( row, col );
                    int ind1, ind2;
                    while (!rowcol2internal( row, col, ind1, ind2 )) {
                        if (ind1 < 0) {
                            this->_add_subdiagonal();
                        } else {
                            this->_add_superdiagonal();
                        }
                    }
                    this->contents().at( ind1 ).at( row ) = val;
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_safe(
                    size_t row, size_t col, ELEMENTTYPE val ) {
                    this->check_size( row, col );
                    int ind1, ind2;
                    if (rowcol2internal( row, col, ind1, ind2 )) {
                        this->contents().at( ind1 ).at( row ) = val;
                    } else {
                        std::ostringstream oss;
                        oss << "Element [" << row << ", " << col
                            << "] is out of band.";
                        throw std::out_of_range( oss.str() );
                    }

                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set(
                    ELEMENTTYPE val ) override {
                    for (auto it1 = this->contents().begin();
                         it1 != this->contents().end(); ++it1) {
                        for (auto it2 = it1->begin(); it2 != it1->end();
                             ++it2) {
                            *it2 = val;
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_row(
                    size_t row, size_t nvals, const size_t *column_inds,
                    const ELEMENTTYPE *vals ) override {
                    for (size_t i = 0; i < nvals; i++) {
                        this->set_safe( row, column_inds[ i ], vals[ i ] );
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    if (vals.size() == columns()) {
                        return set_row(
                            row,
                            NCPA::arrays::index_vector<size_t>( vals.size() ),
                            vals );
                    }
                    if (vals.size() == bandwidth()) {
                        return set_row( row, band_column_indices( row ),
                                        vals );
                    }
                    throw std::invalid_argument(
                        "Value vector size matches neither number of columns "
                        "nor bandwidth" );
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_row(
                    size_t row, ELEMENTTYPE val ) override {
                    return set_row(
                        row, std::vector<ELEMENTTYPE>( bandwidth(), val ) );
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_column(
                    size_t column, size_t nvals, const size_t *row_inds,
                    const ELEMENTTYPE *vals ) override {
                    // NCPA_DEBUG
                    //     << "Made it to band_diagonal_matrix.set_column()"
                    //     << std::endl;
                    for (size_t i = 0; i < nvals; i++) {
                        this->set_safe( row_inds[ i ], column, vals[ i ] );
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const std::vector<ELEMENTTYPE>& vals ) override {
                    if (vals.size() == rows()) {
                        return set_column(
                            col,
                            NCPA::arrays::index_vector<size_t>( vals.size() ),
                            vals );
                    }
                    if (vals.size() == bandwidth()) {
                        return set_column( col, band_row_indices( col ),
                                           vals );
                    }
                    throw std::invalid_argument(
                        "Value vector size matches neither number of rows nor "
                        "bandwidth" );
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_column(
                    size_t col, ELEMENTTYPE val ) override {
                    return set_column(
                        col, std::vector<ELEMENTTYPE>( bandwidth(), val ) );
                }

                virtual abstract_matrix<ELEMENTTYPE>& as_array(
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
                            if (internal2rowcol( ind1, ind2, row, col )) {
                                vals[ row ][ col ]
                                    = this->contents()[ ind1 ][ ind2 ];
                            }
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                // @todo make read-write vector view for columns
                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>> get_row(
                    size_t row ) const override {
                    std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                        new sparse_vector<ELEMENTTYPE>( this->columns() ) );
                    int r, c;
                    for (size_t ind1 = 0; ind1 < this->contents().size();
                         ind1++) {
                        if (internal2rowcol( ind1, row, r, c )) {
                            v->set( c, this->contents()[ ind1 ][ row ] );
                        }
                    }
                    return v;
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    get_column( size_t column ) const override {
                    std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                        new sparse_vector<ELEMENTTYPE>( this->rows() ) );
                    int ind1, ind2;
                    for (size_t r = 0; r < this->rows(); r++) {
                        if (rowcol2internal( r, column, ind1, ind2 )) {
                            v->set( r, this->contents()[ ind1 ][ ind2 ] );
                        }
                    }
                    // for ( size_t ind1 = 0; ind1 < this->contents().size();
                    // ind1++ )
                    // {
                    //     size_t ind2 = (size_t)std::max(
                    //         (int)( column + _n_lower ) - (int)ind1, 0 );
                    //     ind2 = std::min( ind2, this->contents()[ ind1
                    //     ].size() - 1
                    //     ); v->set( ind2, this->contents()[ ind1 ][ ind2 ] );
                    // }
                    return v;
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const override {
                    // this->check_size( row, col );
                    int ind1, ind2;
                    if (rowcol2internal( row, col, ind1, ind2 )) {
                        return this->contents()[ ind1 ][ ind2 ];
                    } else {
                        return _zero;
                    }
                }

                virtual abstract_matrix<ELEMENTTYPE>& resize(
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
                                it->insert( it->cend(), (size_t)ddiag, _zero );
                            }
                        } else if (ddiag < 0) {
                            for (auto it = this->contents().begin();
                                 it != this->contents().end(); ++it) {
                                it->erase( it->cend() - size_t( -ddiag ),
                                           it->cend() );
                            }
                        }
                        _nrows = r;
                        _ncols = c;
                    }

                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& transpose() override {
                    std::vector<std::vector<ELEMENTTYPE>> newcontents
                        = this->contents();
                    std::reverse( newcontents.begin(), newcontents.end() );
                    // NCPA_DEBUG << "Reversed contents" << std::endl;
                    // NCPA_DEBUG << "Before, had " << _n_lower
                    //            << " subdiagonals and " << _n_upper
                    //            << " superdiagonals." << std::endl;
                    std::swap( _n_lower, _n_upper );
                    // NCPA_DEBUG << "Now, has " << _n_lower
                    //            << " subdiagonals and " << _n_upper
                    //            << " superdiagonals." << std::endl;
                    std::swap( _nrows, _ncols );
                    if (_n_lower > 0) {
                        // take the zeros from the end and put them up front
                        for (size_t ind1 = 0; ind1 < _n_lower; ind1++) {
                            size_t invalid = _n_lower - ind1;
                            // NCPA_DEBUG << "For ind1 = " << ind1 << ": moving
                            // "
                            //            << invalid
                            //            << " zeros from back to front"
                            //            << std::endl;
                            for (size_t i = 0; i < invalid; i++) {
                                newcontents[ ind1 ].pop_back();
                                newcontents[ ind1 ].insert(
                                    newcontents[ ind1 ].begin(), 1, _zero );
                                // NCPA_DEBUG << "Cycled one." << std::endl;
                            }
                        }
                    }
                    if (_n_upper > 0) {
                        // take the zeros from the front and put them at the
                        // end
                        for (size_t ind1 = _n_lower + 1;
                             ind1 < newcontents.size(); ind1++) {
                            size_t invalid = ind1 - _n_lower;
                            // NCPA_DEBUG << "For ind1 = " << ind1 << ": moving
                            // "
                            //            << invalid
                            //            << " zeros from back to front"
                            //            << std::endl;
                            for (size_t i = 0; i < invalid; i++) {
                                newcontents[ ind1 ].erase(
                                    newcontents[ ind1 ].cbegin() );
                                newcontents[ ind1 ].push_back( _zero );
                                // NCPA_DEBUG << "Cycled one." << std::endl;
                            }
                        }
                    }
                    _contents = newcontents;
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& swap_rows(
                    size_t ind1, size_t ind2 ) override {
                    throw std::logic_error(
                        "Cannot swap rows of a band-diagonal matrix!" );
                }

                virtual abstract_matrix<ELEMENTTYPE>& swap_columns(
                    size_t ind1, size_t ind2 ) override {
                    throw std::logic_error(
                        "Cannot swap columns of a band-diagonal matrix!" );
                }

                virtual abstract_matrix<ELEMENTTYPE>& zero(
                    size_t row, size_t col ) override {
                    return set( row, col, _zero );
                }

                virtual abstract_matrix<ELEMENTTYPE>& zero() override {
                    _n_lower = 0;
                    _n_upper = 0;
                    this->contents().clear();
                    this->contents().push_back(
                        std::vector<ELEMENTTYPE>( this->diagonal_size() ) );
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    build_vector( size_t n = 0 ) const override {
                    return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                        new sparse_vector<ELEMENTTYPE>( n ) );
                }

                // overrides/specializations of non-pure virtual methods
                virtual bool equals( const abstract_matrix<ELEMENTTYPE>&
                                         other ) const override {
                    if (!is_this_subclass( other )) {
                        return abstract_matrix<ELEMENTTYPE>::equals( other );
                    }
                    const band_diagonal_matrix<ELEMENTTYPE> *b
                        = downcast( &other );

                    return ( rows() == b->rows() && columns() == b->columns()
                             && this->contents() == b->contents() );
                }

                virtual bool is_identity() const override {
                    if (!this->is_square() || !this->is_diagonal()) {
                        return false;
                    }
                    for (size_t i = 0; i < this->contents()[ 0 ].size(); i++) {
                        if (this->contents()[ 0 ][ i ]
                            != NCPA::math::one<ELEMENTTYPE>()) {
                            return false;
                        }
                    }
                    return true;
                }

                virtual bool is_zero( double tol = 1.0e-12 ) const override {
                    if ( this->is_empty() ) {
                        return true;
                    }
                    for (auto it1 = _contents.cbegin(); it1 != _contents.cend(); ++it1) {
                        for (auto it2 = it1->cbegin(); it2 != it1->cend(); ++it2) {
                            if (std::abs( *it2 - _zero) > tol) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual bool is_diagonal() const override {
                    if (this->is_empty()) {
                        return true;
                    }
                    return _n_lower == 0 && _n_upper == 0;
                }

                virtual bool is_tridiagonal() const override {
                    if (this->is_empty()) {
                        return true;
                    }
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
                        for (auto it = this->contents()[ i ].cbegin();
                             it != this->contents()[ i ].cend(); ++it) {
                            if (*it != _zero) {
                                return false;
                            }
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
                        for (auto it
                             = this->contents()[ _n_lower + i + 1 ].cbegin();
                             it != this->contents()[ _n_lower + i + 1 ].cend();
                             ++it) {
                            if (*it != _zero) {
                                return false;
                            }
                        }
                    }
                    return true;
                }

                virtual abstract_matrix<ELEMENTTYPE>& add(
                    const abstract_matrix<ELEMENTTYPE>& babs,
                    ELEMENTTYPE modifier = 1.0 ) override {
                    if (!is_this_subclass( babs )) {
                        return abstract_matrix<ELEMENTTYPE>::add( babs,
                                                                  modifier );
                    }
                    const band_diagonal_matrix<ELEMENTTYPE> *b
                        = downcast( &babs );

                    this->check_size( *b );
                    // lower half
                    size_t skip = 0, add_on = 0, both = 0;
                    if (_n_lower > b->_n_lower) {
                        skip = _n_lower - b->_n_lower;
                        both = b->_n_lower;
                    } else if (_n_lower < b->_n_lower) {
                        add_on = b->_n_lower - _n_lower;
                        both   = _n_lower;
                    } else {
                        both = _n_lower;
                    }
                    for (size_t ind1 = 0; ind1 < both; ind1++) {
                        this->contents()[ ind1 + skip ]
                            = NCPA::arrays::add_vectors(
                                this->contents()[ ind1 + skip ],
                                NCPA::arrays::scale_vector(
                                    b->contents()[ ind1 + add_on ],
                                    modifier ) );
                    }
                    if (add_on > 0) {
                        this->contents().insert(
                            this->contents().begin(), b->contents().begin(),
                            b->contents().begin() + add_on );
                        for (size_t ind1 = 0; ind1 < add_on; ind1++) {
                            this->contents()[ ind1 ]
                                = NCPA::arrays::scale_vector(
                                    this->contents()[ ind1 ], modifier );
                        }
                        _n_lower += add_on;
                    }

                    // upper half
                    skip   = 0;
                    add_on = 0;
                    both   = 0;
                    if (_n_upper > b->_n_upper) {
                        both = b->_n_upper;
                    } else if (_n_upper < b->_n_upper) {
                        both   = _n_upper;
                        add_on = b->_n_upper - _n_upper;
                    } else {
                        both = _n_upper;
                    }
                    for (size_t ind1 = 0; ind1 < both; ind1++) {
                        this->contents()[ _n_lower + ind1 + 1 ]
                            = NCPA::arrays::add_vectors(
                                this->contents()[ _n_lower + ind1 + 1 ],
                                NCPA::arrays::scale_vector(
                                    b->contents()[ b->_n_lower + ind1 + 1 ],
                                    modifier ) );
                    }

                    // insert the last add_on diagonals from b, then negate
                    // them
                    if (add_on > 0) {
                        this->contents().insert( this->contents().end(),
                                                 b->contents().end() - add_on,
                                                 b->contents().end() );
                        _n_upper += add_on;
                        for (size_t ind1 = this->contents().size() - add_on;
                             ind1 < this->contents().size(); ind1++) {
                            this->contents()[ ind1 ]
                                = NCPA::arrays::scale_vector(
                                    this->contents()[ ind1 ], modifier );
                        }
                    }

                    // and the center diagonal
                    this->contents()[ _n_lower ] = NCPA::arrays::add_vectors(
                        this->contents()[ _n_lower ],
                        NCPA::arrays::scale_vector(
                            b->contents()[ b->_n_lower ], modifier ) );
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& add(
                    ELEMENTTYPE b ) override {
                    for (size_t ind1 = 0; ind1 < this->contents().size();
                         ind1++) {
                        for (size_t ind2 = _min_ind2( ind1 );
                             ind2 < _max_ind2( ind1 ); ind2++) {
                            this->contents()[ ind1 ][ ind2 ] += b;
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual bool is_band_diagonal() const override { return true; }

                virtual bool is_this_subclass(
                    const abstract_matrix<ELEMENTTYPE>& b ) const override {
                    if (auto *derived = dynamic_cast<
                            const band_diagonal_matrix<ELEMENTTYPE> *>( &b )) {
                        return true;
                    } else {
                        return false;
                    }
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

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> multiply(
                    const abstract_matrix<ELEMENTTYPE>& other )
                    const override {
                    if (!is_this_subclass( other )) {
                        return abstract_matrix<ELEMENTTYPE>::multiply( other );
                    }
                    const band_diagonal_matrix<ELEMENTTYPE> *b
                        = downcast( &other );

                    size_t new_n_lower
                        = this->lower_bandwidth() + b->lower_bandwidth(),
                        new_n_upper
                        = this->upper_bandwidth() + b->upper_bandwidth();
                    std::unique_ptr<abstract_matrix<ELEMENTTYPE>> product(
                        new band_diagonal_matrix<ELEMENTTYPE>(
                            this->rows(), b->columns(), new_n_lower,
                            new_n_upper ) );
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

                virtual abstract_matrix<ELEMENTTYPE>& scale(
                    const abstract_matrix<ELEMENTTYPE>& b ) override {
                    this->check_size( b );
                    int r, c;
                    for (size_t ind1 = 0; ind1 < this->contents().size();
                         ind1++) {
                        for (size_t ind2 = _min_ind2( ind1 );
                             ind2 < _max_ind2( ind1 ); ind2++) {
                            if (internal2rowcol( ind1, ind2, r, c )) {
                                this->contents()[ ind1 ][ ind2 ]
                                    *= b.get( r, c );
                            }
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& scale(
                    ELEMENTTYPE val ) override {
                    for (auto it1 = this->contents().begin();
                         it1 != this->contents().end(); ++it1) {
                        for (auto it2 = it1->begin(); it2 != it1->end();
                             ++it2) {
                            *it2 *= val;
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& identity(
                    size_t nrows, size_t ncols ) override {
                    this->clear();
                    _nrows   = nrows;
                    _ncols   = ncols;
                    _n_lower = 0;
                    _n_upper = 0;
                    std::vector<ELEMENTTYPE> diag(
                        this->diagonal_size( 0 ),
                        NCPA::math::one<ELEMENTTYPE>() );
                    this->contents().push_back( diag );
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_diagonal(
                    size_t nvals, const ELEMENTTYPE *vals,
                    int offset = 0 ) override {
                    if (nvals > this->diagonal_size( offset )) {
                        throw std::out_of_range(
                            "Too many values for requested diagonal" );
                    }
                    if (offset > 0) {
                        while (offset > _n_upper) {
                            this->_add_superdiagonal();
                        }
                        for (size_t i = 0; i < nvals; i++) {
                            this->contents()[ _n_lower + offset ][ i ]
                                = vals[ i ];
                        }
                    } else if (offset < 0) {
                        while (-offset > (int)_n_lower) {
                            this->_add_subdiagonal();
                        }
                        int ind1 = _n_lower + offset;
                        for (size_t i = 0; i < nvals; i++) {
                            int ind2 = _min_ind2( ind1 ) + i;
                            this->contents()[ ind1 ][ ind2 ] = vals[ i ];
                        }
                    } else {
                        for (size_t i = 0; i < nvals; i++) {
                            this->contents()[ _n_lower ][ i ] = vals[ i ];
                        }
                    }
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                    ;
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    get_diagonal( int offset = 0 ) const override {
                    // std::cout << "Calling get_diagonal( " << offset << " )"
                    //           << std::endl;
                    std::vector<ELEMENTTYPE> diag(
                        this->diagonal_size( offset ), _zero );
                    // if ( offset > _n_upper || -offset > _n_lower ) {
                    //     diag.resize( this->diagonal_size( offset ), _zero );
                    // } else {
                    int ind1 = (int)_n_lower + offset;
                    if (ind1 >= 0 && ind1 < this->contents().size()) {
                        // std::cout
                        //     << "ind1 = " << ind1
                        //     << " is in range for _n_upper = " << _n_upper
                        //     << " and _n_lower = " << _n_lower
                        //     << ", this->contents().size() = " <<
                        //     this->contents().size()
                        //     << std::endl;
                        // std::cout << "Fetching internal index " << ind1
                        //           << std::endl;
                        // std::cout << "Fetching ind2 = " << _min_ind2( ind1 )
                        //           << " to " << _max_ind2( ind1 ) <<
                        //           std::endl;
                        diag.assign( this->contents()[ ind1 ].begin()
                                         + _min_ind2( ind1 ),
                                     this->contents()[ ind1 ].begin()
                                         + _max_ind2( ind1 ) );
                        // } else {
                        //     std::cout << "ind1 = " << ind1
                        //               << " out of rangefor _n_upper = " <<
                        //               _n_upper
                        //               << " and _n_lower = " << _n_lower
                        //               << ", returning zeros" << std::endl;
                    }
                    return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                        new dense_vector<ELEMENTTYPE>( diag ) );
                    // std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    //     vdiag = build_vector();
                    // vdiag->set( diag );
                    // return vdiag;
                }

                virtual size_t bandwidth() const override {
                    return this->contents().size();
                }

                virtual size_t lower_bandwidth() const override {
                    return _n_lower;
                }

                virtual size_t upper_bandwidth() const override {
                    return _n_upper;
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
                    // NCPA_DEBUG << "Getting nonzero indices for column " <<
                    // col << std::endl;
                    for (int i = std::max( 0, (int)col - (int)_n_upper );
                         i < std::min( (int)rows(),
                                       (int)col + (int)_n_lower + 1 );
                         i++) {
                        inds.push_back( (size_t)i );
                        // NCPA_DEBUG << "Adding index " << i << std::endl;
                    }
                    return inds;
                }

                // virtual abstract_matrix<ELEMENTTYPE>& operator-() const override {
                //     band_diagonal_matrix<ELEMENTTYPE> neg(*this);
                //     ELEMENTTYPE negone = -NCPA::math::one<ELEMENTTYPE>();
                //     for (auto it = neg._contents.begin(); it != neg._contents.end(); ++it) {
                //         for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
                //             NCPA::arrays::scale_vector( *it2, negone );
                //         }
                //     }
                //     return neg;
                // }

                virtual std::vector<std::vector<ELEMENTTYPE>>& contents() {
                    return _contents;
                }

                virtual const std::vector<std::vector<ELEMENTTYPE>>& contents()
                    const {
                    return _contents;
                }

            protected:
                size_t _nrows, _ncols, _n_lower, _n_upper;

                std::vector<std::vector<ELEMENTTYPE>> _contents;

                const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();

                void _prep_diagonals() {
                    this->contents().clear();
                    size_t n_diag = _n_lower + _n_upper + 1;
                    // std::cout << "Adding " << n_diag << " diagonals" <<
                    // std::endl;
                    std::vector<ELEMENTTYPE> diag( this->diagonal_size( 0 ) );
                    for (size_t i = 0; i < n_diag; i++) {
                        this->contents().push_back( diag );
                    }
                    // std::cout << "Done" << std::endl;
                }

                void _add_subdiagonal() {
                    std::vector<ELEMENTTYPE> diag( this->diagonal_size( 0 ) );
                    this->contents().insert( this->contents().cbegin(), diag );
                    _n_lower++;
                }

                void _add_superdiagonal() {
                    std::vector<ELEMENTTYPE> diag( this->diagonal_size( 0 ) );
                    this->contents().push_back( diag );
                    _n_upper++;
                }

                size_t _min_ind2( size_t ind1 ) const {
                    return (size_t)std::max( (int)_n_lower - (int)ind1, 0 );
                }

                size_t _max_ind2( size_t ind1 ) const {
                    return (size_t)std::min( this->rows(),
                                             this->rows() + _n_lower - ind1 );
                }

                virtual abstract_matrix<ELEMENTTYPE>& _add(
                    const band_diagonal_matrix<ELEMENTTYPE>& b,
                    ELEMENTTYPE modifier = 1.0 ) {
                    this->check_size( b );
                    size_t new_n_lower = std::max( _n_lower, b._n_lower );
                    size_t new_n_upper = std::max( _n_upper, b._n_upper );
                    std::vector<std::vector<ELEMENTTYPE>> newcontents(
                        new_n_lower + new_n_upper + 1 );

                    // diagonals first
                    newcontents[ new_n_lower ] = NCPA::arrays::add_vectors(
                        this->contents()[ _n_lower ],
                        NCPA::arrays::scale_vector( b.contents()[ b._n_lower ],
                                                    modifier ) );

                    // subdiagonals
                    for (size_t n = 1; n < new_n_lower; n++) {
                        if (n <= _n_lower && n <= b._n_lower) {
                            newcontents[ new_n_lower - n ]
                                = NCPA::arrays::add_vectors(
                                    this->contents()[ _n_lower - n ],
                                    NCPA::arrays::scale_vector(
                                        b.contents()[ b._n_lower - n ],
                                        modifier ) );
                        } else if (n <= _n_lower) {
                            newcontents[ new_n_lower - n ]
                                = this->contents()[ _n_lower - n ];
                        } else if (n <= b._n_lower) {
                            newcontents[ new_n_lower - n ]
                                = NCPA::arrays::scale_vector(
                                    b.contents()[ b._n_lower - n ], modifier );
                        }
                    }

                    // superdiagonals
                    for (size_t n = 1; n < new_n_upper; n++) {
                        if (n <= _n_upper && n <= b._n_upper) {
                            newcontents[ new_n_lower + n ]
                                = NCPA::arrays::add_vectors(
                                    this->contents()[ _n_lower + n ],
                                    NCPA::arrays::scale_vector(
                                        b.contents()[ b._n_lower + n ],
                                        modifier ) );
                        } else if (n <= _n_upper) {
                            newcontents[ new_n_lower + n ]
                                = this->contents()[ _n_lower + n ];
                        } else if (n <= b._n_upper) {
                            newcontents[ new_n_lower + n ]
                                = NCPA::arrays::scale_vector(
                                    b.contents()[ b._n_lower + n ], modifier );
                        }
                    }
                    _contents = newcontents;
                    _n_lower  = new_n_lower;
                    _n_upper  = new_n_upper;
                    RETURN_THIS_AS_ABSTRACT_MATRIX;
                }

                static void _do_band_diagonal_multiply(
                    const band_diagonal_matrix<ELEMENTTYPE>& a,
                    const band_diagonal_matrix<ELEMENTTYPE>& b,
                    band_diagonal_matrix<ELEMENTTYPE>& product ) {
                    int prows       = (int)( product.rows() );
                    int pcols       = (int)( product.columns() );
                    int new_n_lower = product.lower_bandwidth();
                    int new_n_upper = product.upper_bandwidth();
                    size_t counter  = 0;
                    for (int r = 0; r < prows; r++) {
                        for (int c = std::max( 0, r - (int)new_n_lower );
                             c < std::min( pcols, r + (int)new_n_upper + 1 );
                             c++) {
                            ELEMENTTYPE val = a._zero;
                            int kmin        = std::max(
                                std::max( 0,
                                                 r - (int)( a.lower_bandwidth() ) ),
                                std::max( 0,
                                                 c - (int)( b.upper_bandwidth() ) ) );
                            int kmax = std::min(
                                std::min( (int)( a.columns() ),
                                          r + (int)a.upper_bandwidth() + 1 ),
                                std::min( (int)( b.rows() ),
                                          c + (int)( b.lower_bandwidth() )
                                              + 1 ) );
                            for (int k = kmin; k < kmax; k++) {
                                val += a.get( r, k ) * b.get( k, c );
                                counter++;
                            }
                            product.set( (size_t)r, (size_t)c, val );
                        }
                    }
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
