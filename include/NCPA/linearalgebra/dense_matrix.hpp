#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/defines.hpp"
#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/dense_vector.hpp"
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

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::dense_matrix,
                                             ELEMENTTYPE );

namespace NCPA {
    namespace linear {


        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class dense_matrix<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public abstract_matrix<ELEMENTTYPE> {
            public:
                dense_matrix() : abstract_matrix<ELEMENTTYPE>() {}

                dense_matrix( size_t nrows, size_t ncols ) : dense_matrix() {
                    if (nrows > 0 && ncols > 0) {
                        resize( nrows, ncols );
                    }
                }

                dense_matrix( const dense_matrix<ELEMENTTYPE>& other ) :
                    dense_matrix<ELEMENTTYPE>( other.rows(),
                                               other.columns() ) {
                    _contents = other._contents;
                    // std::memcpy( _contents, other._contents,
                    //              other.rows() * other.columns()
                    //                  * sizeof( ELEMENTTYPE ) );
                }

                dense_matrix( const abstract_matrix<ELEMENTTYPE>& other ) :
                    dense_matrix<ELEMENTTYPE>() {
                    this->copy( other );
                }

                virtual ~dense_matrix() {
                    // if (_contents != nullptr) {
                    //     delete[] _contents;
                    //     _contents = nullptr;
                    // }
                }

                dense_matrix<ELEMENTTYPE>& operator=(
                    dense_matrix<ELEMENTTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                dense_matrix<ELEMENTTYPE>& operator=(
                    abstract_matrix<ELEMENTTYPE> other ) {
                    dense_matrix<ELEMENTTYPE> othercopy( other );
                    swap( *this, othercopy );
                    return *this;
                }

                friend void ::swap<ELEMENTTYPE>(
                    dense_matrix<ELEMENTTYPE>& a,
                    dense_matrix<ELEMENTTYPE>& b ) noexcept;


                // bring some methods in unchanged
                using abstract_matrix<ELEMENTTYPE>::add;
                using abstract_matrix<ELEMENTTYPE>::scale;
                using abstract_matrix<ELEMENTTYPE>::set_row;
                using abstract_matrix<ELEMENTTYPE>::set_column;

                template<typename ANYTYPE,
                         ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                                       ELEMENTTYPE )>
                abstract_matrix<ELEMENTTYPE>& add( ANYTYPE b ) {
                    for (auto it = _contents.begin(); it != _contents.end();
                         ++it) {
                        *it += b;
                        // for (size_t i = 0; i < _rows * _cols; i++) {
                        //     _contents[ i ] += b;
                    }
                    return *this;
                }

                virtual dense_matrix<ELEMENTTYPE>& as_array(
                    size_t& nrows, size_t& ncols,
                    ELEMENTTYPE **& vals ) override {
                    if (nrows == 0 && ncols == 0) {
                        nrows = rows();
                        ncols = columns();
                        vals
                            = NCPA::arrays::zeros<ELEMENTTYPE>( nrows, ncols );
                    } else if (nrows != rows() || ncols != columns()) {
                        throw std::invalid_argument(
                            "Size mismatch between vector and target "
                            "array" );
                    }
                    for (size_t i = 0; i < rows(); i++) {
                        for (size_t j = 0; j < columns(); j++) {
                            std::memcpy( vals[ i ],
                                         _contents.data() + _rc2ind( i, 0 ),
                                         _cols * sizeof( ELEMENTTYPE ) );
                        }
                    }
                    return *this;
                }

                virtual size_t bandwidth() const override {
                    return this->upper_bandwidth() + this->lower_bandwidth()
                         + 1;
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    build_vector( size_t n = 0 ) const override {
                    return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                        new dense_vector<ELEMENTTYPE>( n ) );
                }

                virtual dense_matrix<ELEMENTTYPE>& clean(
                    ELEMENTTYPE tol ) override {
                    // auto atol = std::abs( tol );
                    for (auto it = _contents.begin(); it != _contents.end();
                         ++it) {
                        NCPA::constants::zero_out( *it, tol );
                        // if (std::abs( *it ) < atol) {
                        //     *it = _zero;
                        // }
                    }
                    return *this;
                }

                virtual dense_matrix<ELEMENTTYPE>& clear() override {
                    _contents.clear();
                    _rows = 0;
                    _cols = 0;
                    return *this;
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> clone()
                    const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new dense_matrix( *this ) );
                }

                virtual size_t columns() const override { return _cols; }

                virtual dense_matrix<ELEMENTTYPE>& copy(
                    const abstract_matrix<ELEMENTTYPE>& other ) override {
                    resize( other.rows(), other.columns() );
                    for (auto i = 0; i < other.rows(); i++) {
                        for (auto j = 0; j < other.columns(); j++) {
                            set( i, j, other.get( i, j ) );
                        }
                    }
                    return *this;
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                    fresh_clone() const override {
                    return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        new dense_matrix() );
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const override {
                    // return _contents.[ _rc2ind( row, col ) ];
                    return _contents.at( _rc2ind( row, col ) );
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    get_column( size_t column ) const override {
                    this->check_size( 0, column );
                    std::vector<ELEMENTTYPE> vcol( rows() );
                    for (size_t row = 0; row < rows(); row++) {
                        vcol[ row ] = this->get( row, column );
                        // vcol[ row ] = _elements[ row ][ column ];
                    }
                    return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                        new dense_vector<ELEMENTTYPE>( vcol ) );
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>> get_row(
                    size_t row ) const override {
                    this->check_size( row, 0 );
                    return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                        new dense_vector<ELEMENTTYPE>(
                            this->columns(),
                            _contents.data() + _rc2ind( row, 0 ) ) );
                }

                virtual std::string id() const override {
                    return "NCPA Basic Dense Matrix";
                }

                virtual bool is_this_subclass(
                    const abstract_matrix<ELEMENTTYPE>& b ) const override {
                    if (auto *derived
                        = dynamic_cast<const dense_matrix<ELEMENTTYPE> *>(
                            &b )) {
                        return true;
                    } else {
                        return false;
                    }
                }

                virtual bool is_zero( double tol = 1.0e-12 ) const override {
                    if (this->is_empty()) {
                        return true;
                    }
                    double atol = std::abs( tol );
                    for (auto it = _contents.begin(); it != _contents.end();
                         ++it) {
                        if (std::abs( *it ) > atol) {
                            return false;
                        }
                    }
                    return true;
                }

                virtual bool is_zero( size_t r, size_t c,
                                      double tol = 1.0e-12 ) const override {
                    return ( std::abs( this->get( r, c ) ) <= tol );
                }

                virtual size_t lower_bandwidth() const override {
                    int off = this->max_off_diagonal();
                    while (off > 0) {
                        if (this->get_diagonal( -off )->count_nonzero_indices()
                            > 0) {
                            return (size_t)off;
                        }
                        off--;
                    }
                    return 0;
                }

                virtual dense_matrix<ELEMENTTYPE>& resize(
                    size_t rows, size_t cols ) override {
                    if (_contents.empty() || _rows == 0 || _cols == 0) {
                        this->clear();  // make sure we're clean
                        _contents = std::vector<ELEMENTTYPE>( rows * cols );
                        _rows     = rows;
                        _cols     = cols;
                    } else {
                        // ELEMENTTYPE *newcontents
                        //     = NCPA::arrays::zeros<ELEMENTTYPE>( rows * cols
                        //     );
                        std::vector<ELEMENTTYPE> newcontents( rows * cols );
                        for (size_t r = 0; r < std::min( rows, _rows ); ++r) {
                            std::memcpy( newcontents.data() + r * cols,
                                         _contents.data() + _rc2ind( r, 0 ),
                                         std::min( _cols, cols )
                                             * sizeof( ELEMENTTYPE ) );
                        }
                        this->clear();
                        _contents = newcontents;
                        _rows     = rows;
                        _cols     = cols;
                    }
                    return *this;
                }

                virtual size_t rows() const override { return _rows; }

                template<typename ANYTYPE,
                         ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                                       ELEMENTTYPE )>
                dense_matrix<ELEMENTTYPE>& scale( ANYTYPE val ) {
                    for (auto it = _contents.begin(); it != _contents.end();
                         ++it) {
                        *it *= val;
                    }
                    return *this;
                }

                // set individual elements
                virtual dense_matrix<ELEMENTTYPE>& set(
                    size_t row, size_t col, ELEMENTTYPE val ) override {
                    _contents.at( _rc2ind( row, col ) ) = val;
                    return *this;
                }

                virtual dense_matrix<ELEMENTTYPE>& set(
                    ELEMENTTYPE val ) override {
                    _contents = std::vector<ELEMENTTYPE>( _rows * _cols, val );
                    return *this;
                }

                // set multiple elements in one column
                virtual dense_matrix<ELEMENTTYPE>& set_column(
                    size_t column, size_t nvals, const size_t *row_inds,
                    const ELEMENTTYPE *vals ) override {
                    for (size_t i = 0; i < nvals; i++) {
                        this->set( row_inds[ i ], column, vals[ i ] );
                    }
                    return *this;
                }

                // set multiple elements in one row
                virtual dense_matrix<ELEMENTTYPE>& set_row(
                    size_t row, size_t nvals, const size_t *column_inds,
                    const ELEMENTTYPE *vals ) override {
                    size_t start = _rc2ind( row, 0 );
                    for (size_t i = 0; i < nvals; i++) {
                        _contents.at( start + column_inds[ i ] ) = vals[ i ];
                    }
                    return *this;
                }

                virtual dense_matrix<ELEMENTTYPE>& swap_columns(
                    size_t ind1, size_t ind2 ) override {
                    auto col1 = this->get_column( ind1 );
                    auto col2 = this->get_column( ind2 );
                    this->set_column( ind1, *col2 );
                    this->set_column( ind2, *col1 );
                    return *this;
                }

                virtual dense_matrix<ELEMENTTYPE>& swap_rows(
                    size_t ind1, size_t ind2 ) override {
                    auto row1 = this->get_row( ind1 );
                    auto row2 = this->get_row( ind2 );
                    this->set_row( ind2, *row1 );
                    this->set_row( ind1, *row2 );
                    return *this;
                }

                virtual dense_matrix<ELEMENTTYPE>& transpose() override {
                    std::unique_ptr<abstract_matrix<ELEMENTTYPE>> orig
                        = this->clone();
                    this->clear().resize( orig->columns(), orig->rows() );
                    for (size_t i = 0; i < orig->rows(); ++i) {
                        for (size_t j = 0; j < orig->columns(); ++j) {
                            this->set( j, i, orig->get( i, j ) );
                        }
                    }
                    return *this;
                }

                virtual size_t upper_bandwidth() const override {
                    int off = this->max_off_diagonal();
                    while (off > 0) {
                        if (this->get_diagonal( off )->count_nonzero_indices()
                            > 0) {
                            return (size_t)off;
                        }
                        off--;
                    }
                    return 0;
                }

                // virtual dense_matrix<ELEMENTTYPE>& operator-() const
                // override {
                //     dense_matrix<ELEMENTTYPE> neg( this->rows(),
                //     this->columns() ); ELEMENTTYPE one =
                //     NCPA::math::one<ELEMENTTYPE>(); for (size_t r = 0; r <
                //     this->rows(); ++r) {
                //         for (size_t c = 0; c < this->columns(); ++c) {
                //             neg.set( r, c, this->get( r, c ) * -one );
                //         }
                //     }
                //     return neg;
                // }

            protected:
                // ELEMENTTYPE *_contents = nullptr;
                std::vector<ELEMENTTYPE> _contents;
                size_t _rows = 0, _cols = 0;
                bool _finalized = false;

                const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();

                size_t _rc2ind( size_t row, size_t col ) const {
                    return row * _cols + col;
                }
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::dense_matrix<T>& a,
                  NCPA::linear::dense_matrix<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::abstract_matrix<T>&>( a ),
            static_cast<NCPA::linear::abstract_matrix<T>&>( b ) );
    swap( a._contents, b._contents );
    swap( a._rows, b._rows );
    swap( a._cols, b._cols );
    swap( a._finalized, b._finalized );
}
