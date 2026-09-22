#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/matrix_view.hpp"

#include <memory>

namespace NCPA {
    namespace linear {
        template<typename ELEMENTTYPE>
        class transpose_view : public matrix_view<ELEMENTTYPE> {
            public:
                transpose_view( const abstract_matrix<ELEMENTTYPE>& m ) :
                    matrix_view<ELEMENTTYPE>( m ) {}

                transpose_view(
                    const std::unique_ptr<abstract_matrix<ELEMENTTYPE>>& m ) :
                    matrix_view<ELEMENTTYPE>( m ) {}

                virtual ~transpose_view() {}

                virtual abstract_matrix<ELEMENTTYPE>& as_array(
                    size_t& nrows, size_t& ncols, ELEMENTTYPE **& vals ) {
                    if (_is_transposed) {
                        if (nrows == 0 && ncols == 0) {
                            nrows = this->rows();
                            ncols = this->columns();
                            vals  = zeros<ELEMENTTYPE>( nrows, ncols );
                        } else if (nrows != rows() || ncols != columns()) {
                            throw std::invalid_argument(
                                "Size mismatch between vector and target "
                                "array" );
                        }
                        for (size_t i = 0; i < rows(); i++) {
                            for (size_t j = 0; j < columns(); j++) {
                                vals[ i ][ j ] = this->get( i, j );
                            }
                        }
                    } else {
                        this->internal()->as_array( nrows, ncols, vals );
                    }
                    return *this;
                }

                virtual size_t bandwidth() const {
                    return this->internal()->bandwidth();
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    build_vector( size_t n = 0 ) const {
                    return this->internal()->build_vector( n );
                }

                virtual abstract_matrix<ELEMENTTYPE>& clean(
                    ELEMENTTYPE tol ) {
                    return this->internal()->clean( tol );
                }

                virtual abstract_matrix<ELEMENTTYPE>& clear() {
                    return this->internal()->clear();
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>> clone()
                    const {
                    if (_is_transposed) {
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new transpose_view( this->internal()->clone() ) );
                    } else {
                        return this->internal()->clone();
                    }
                }

                virtual size_t columns() const {
                    if (_is_transposed) {
                        return this->internal()->rows();
                    } else {
                        return this->internal()->columns();
                    }
                }

                virtual abstract_matrix<ELEMENTTYPE>& copy(
                    const abstract_matrix<ELEMENTTYPE>& other ) {
                    throw std::logic_error( "Not implemented" );
                }

                virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                    fresh_clone() const {
                    if (_is_transposed) {
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new transpose_view(
                                this->internal()->fresh_clone() ) );
                    } else {
                        return this->internal()->fresh_clone();
                    }
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const {
                    if (_is_transposed) {
                        return this->internal()->get( col, row );
                    } else {
                        return this->internal()->get( row, col );
                    }
                }

                // @todo make read-write vector view for columns
                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    get_column( size_t column ) const {
                    if (_is_transposed) {
                        return this->internal()->get_row( column );
                    } else {
                        return this->internal()->get_column( column );
                    }
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>> get_row(
                    size_t row ) const {
                    if (_is_transposed) {
                        return this->internal()->get_column( row );
                    } else {
                        return this->internal()->get_row( row );
                    }
                }

                virtual std::string id() const {
                    return this->internal()->id();
                }

                virtual bool is_this_subclass(
                    const abstract_matrix<ELEMENTTYPE>& b ) const {
                    return this->internal()->is_this_subclass( b );
                }

                virtual bool is_zero( double tol = 1.0e-12 ) const {
                    return this->internal()->is_zero( tol );
                }

                virtual bool is_zero( size_t r, size_t c,
                                      double tol = 1.0e-12 ) const {
                    if (_is_transposed) {
                        return this->internal()->is_zero( c, r, tol );
                    } else {
                        return this->internal()->is_zero( r, c, tol );
                    }
                }

                virtual size_t lower_bandwidth() const {
                    if (_is_transposed) {
                        return this->internal()->upper_bandwidth();
                    } else {
                        return this->internal()->lower_bandwidth();
                    }
                }

                virtual abstract_matrix<ELEMENTTYPE>& resize( size_t rows,
                                                              size_t cols ) {
                    if (_is_transposed) {
                        this->internal()->resize( cols, rows );
                    } else {
                        this->internal()->resize( rows, cols );
                    }
                    return *this;
                }

                virtual size_t rows() const {
                    if (_is_transposed) {
                        return this->internal()->columns();
                    } else {
                        return this->internal()->rows();
                    }
                }

                virtual abstract_matrix<ELEMENTTYPE>& set( size_t row,
                                                           size_t col,
                                                           ELEMENTTYPE val ) {
                    if (_is_transposed) {
                        this->internal()->set( col, row, val );
                    } else {
                        this->internal()->set( row, col, val );
                    }
                    return *this;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set( ELEMENTTYPE val ) {
                    this->internal()->set( val );
                    return *this;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_column(
                    size_t column, size_t nvals, const size_t *row_inds,
                    const ELEMENTTYPE *vals ) {
                    if (_is_transposed) {
                        this->internal()->set_row( column, nvals, row_inds,
                                                   vals );
                    } else {
                        this->internal()->set_column( column, nvals, row_inds,
                                                      vals );
                    }
                    return *this;
                }

                virtual abstract_matrix<ELEMENTTYPE>& set_row(
                    size_t row, size_t nvals, const size_t *column_inds,
                    const ELEMENTTYPE *vals ) {
                    if (_is_transposed) {
                        this->internal()->set_column( row, nvals, column_inds,
                                                      vals );
                    } else {
                        this->internal()->set_row( row, nvals, column_inds,
                                                   vals );
                    }
                    return *this;
                }

                virtual abstract_matrix<ELEMENTTYPE>& swap_columns(
                    size_t ind1, size_t ind2 ) {
                    if (_is_transposed) {
                        this->internal()->swap_rows( ind1, ind2 );
                    } else {
                        this->internal()->swap_columns( ind1, ind2 );
                    }
                    return *this;
                }

                virtual abstract_matrix<ELEMENTTYPE>& swap_rows(
                    size_t ind1, size_t ind2 ) {
                    if (_is_transposed) {
                        this->internal()->swap_columns( ind1, ind2 );
                    } else {
                        this->internal()->swap_rows( ind1, ind2 );
                    }
                    return *this;
                }

                virtual abstract_matrix<ELEMENTTYPE>& transpose() {
                    _is_transposed = !_is_transposed;
                    return *this;
                }

                virtual size_t upper_bandwidth() const {
                    if (_is_transposed) {
                        return this->internal()->lower_bandwidth();
                    } else {
                        return this->internal()->upper_bandwidth();
                    }
                }

            private:
                bool _is_transposed = true;
        };
    }  // namespace linear
}  // namespace NCPA
