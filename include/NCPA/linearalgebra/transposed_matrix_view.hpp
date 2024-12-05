#pragma once

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

#define _MATRIX_VIEW_OVERRIDE_THROW                                     \
    override {                                                          \
        throw std::logic_error( "Matrix transpose view is read-only" ); \
    }

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::transposed_matrix_view, ELEMENTTYPE );


#define RETURN_SUPERCLASS_REFERENCE \
    return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>( this );

namespace NCPA {
    namespace linear {

        namespace details {
            NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
                class transposed_matrix_view<ELEMENTTYPE,
                                             _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
                : public abstract_matrix<ELEMENTTYPE> {
                public:
                    transposed_matrix_view(
                        const abstract_matrix<ELEMENTTYPE>& mat ) {
                        _mat = &mat;
                    }

                    transposed_matrix_view(
                        const transposed_matrix_view<ELEMENTTYPE>& other ) :
                        transposed_matrix_view<ELEMENTTYPE>() {
                        _mat = other._mat;
                    }

                    virtual ~transposed_matrix_view() {}

                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        clone() const override {
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new transposed_matrix_view( *this ) );
                    }

                    sparse_matrix<ELEMENTTYPE>& operator=(
                        sparse_matrix<ELEMENTTYPE> other ) {
                        swap( *this, other );
                        return *this;
                    }

                    virtual std::string id() const override {
                        return "NCPA Transposed Matrix View";
                    }

                    friend void ::swap<ELEMENTTYPE>(
                        sparse_matrix<ELEMENTTYPE>& a,
                        sparse_matrix<ELEMENTTYPE>& b ) noexcept;

                    virtual bool equals( const abstract_matrix<ELEMENTTYPE>&
                                             other ) const override {
                        if ( rows() != other.columns() ) {
                            return false;
                        }
                        if ( columns() != other.rows() ) {
                            return false;
                        }

                        for ( auto r = 0; r < rows(); r++ ) {
                            for ( auto c = 0; c < columns(); c++ ) {
                                if ( !NCPA::math::equals(
                                         get( r, c ), other.get( c, r ) ) ) {
                                    return false;
                                }
                            }
                        }
                        return true;
                    }

                    virtual bool is_upper_triangular() const override {
                        return _mat->is_lower_triangular();
                    }

                    virtual bool is_lower_triangular() const override {
                        return _mat->is_upper_triangular();
                    }

                    virtual size_t rows() const override {
                        return _mat->columns();
                    }

                    virtual size_t columns() const override {
                        return _mat->rows();
                    }

                    virtual const ELEMENTTYPE& get(
                        size_t row, size_t col ) const override {
                        return _mat->get( col, row );
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_row( size_t row ) const override {
                        return _mat->get_column( row );
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_column( size_t column ) const override {
                        return _mat->get_row( column );
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_diagonal( int offset = 0 ) const override {
                        return _mat->get_diagonal( -offset );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& as_array(
                        size_t& nrows, size_t& ncols,
                        ELEMENTTYPE **& vals ) override {
                        _mat->as_array( ncols, nrows, vals );
                        RETURN_SUPERCLASS_REFERENCE;
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        build_vector( size_t n = 0 ) const override {
                        return _mat->build_vector( n );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& operator[](
                        size_t i ) _MATRIX_VIEW_OVERRIDE_THROW;


                    virtual abstract_matrix<ELEMENTTYPE>& clear()
                        _MATRIX_VIEW_OVERRIDE_THROW;

                    virtual abstract_matrix<ELEMENTTYPE>& resize(
                        size_t rows, size_t cols ) _MATRIX_VIEW_OVERRIDE_THROW;

                    virtual abstract_matrix<ELEMENTTYPE>& transpose()
                        _MATRIX_VIEW_OVERRIDE_THROW;

                    // @todo make swapped-row and swapped-column views?
                    virtual abstract_matrix<ELEMENTTYPE>& swap_rows(
                        size_t ind1, size_t ind2 ) _MATRIX_VIEW_OVERRIDE_THROW;
                    virtual abstract_matrix<ELEMENTTYPE>& swap_columns(
                        size_t ind1, size_t ind2 ) _MATRIX_VIEW_OVERRIDE_THROW;
                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        size_t row, size_t col,
                        ELEMENTTYPE val ) _MATRIX_VIEW_OVERRIDE_THROW;
                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        ELEMENTTYPE val ) _MATRIX_VIEW_OVERRIDE_THROW;
                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, size_t nvals, const size_t *column_inds,
                        const ELEMENTTYPE *vals ) _MATRIX_VIEW_OVERRIDE_THROW;

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, size_t nvals, const size_t *row_inds,
                        const ELEMENTTYPE *vals ) _MATRIX_VIEW_OVERRIDE_THROW;

                    virtual abstract_matrix<ELEMENTTYPE>& set_diagonal(
                        size_t nvals, const ELEMENTTYPE *vals,
                        int offset = 0 ) _MATRIX_VIEW_OVERRIDE_THROW;



                private:
                    const abstract_matrix<ELEMENTTYPE> *_mat;
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::details::transposed_matrix_view<T>& a,
                  NCPA::linear::details::transposed_matrix_view<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::details::abstract_matrix<T>&>( a ),
            static_cast<NCPA::linear::details::abstract_matrix<T>&>( b ) );
    swap( a._mat, b._mat );
}