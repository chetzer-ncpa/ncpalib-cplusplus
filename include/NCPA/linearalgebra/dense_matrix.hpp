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

// namespace NCPA {
//     namespace linear {
//         namespace details {
//             NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( dense_matrix,
//                                                          abstract_matrix );
//         }
//     }  // namespace linear
// }  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::dense_matrix, ELEMENTTYPE );

namespace NCPA {
    namespace linear {

        namespace details {
            NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
                class dense_matrix<ELEMENTTYPE,
                                   _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
                : public abstract_matrix<ELEMENTTYPE> {
                public:
                    dense_matrix() : abstract_matrix<ELEMENTTYPE>() {}

                    dense_matrix( size_t nrows, size_t ncols ) :
                        dense_matrix() {
                        resize( nrows, ncols );
                    }

                    dense_matrix( const dense_matrix<ELEMENTTYPE>& other ) :
                        dense_matrix<ELEMENTTYPE>() {
                        _elements.clear();
                        for ( auto it = other._elements.cbegin();
                              it != other._elements.cend(); ++it ) {
                            _elements.emplace_back(
                                dense_vector<ELEMENTTYPE>( *it ) );
                        }
                    }

                    dense_matrix( const abstract_matrix<ELEMENTTYPE>& other ) :
                        dense_matrix<ELEMENTTYPE>() {
                            this->copy( other );
                        // resize( other.rows(), other.columns() );
                        // for ( auto i = 0; i < other.rows(); i++ ) {
                        //     for ( auto j = 0; j < other.columns(); j++ ) {
                        //         set( i, j, other.get( i, j ) );
                        //     }
                        // }
                    }

                    virtual ~dense_matrix() {}

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

                    // bring some methods in unchanged
                    using abstract_matrix<ELEMENTTYPE>::set;
                    using abstract_matrix<ELEMENTTYPE>::set_row;
                    using abstract_matrix<ELEMENTTYPE>::set_column;
                    using abstract_matrix<ELEMENTTYPE>::add;
                    using abstract_matrix<ELEMENTTYPE>::scale;

                    virtual abstract_matrix<ELEMENTTYPE>& copy(
                        const abstract_matrix<ELEMENTTYPE>& other ) override {
                        resize( other.rows(), other.columns() );
                        for ( auto i = 0; i < other.rows(); i++ ) {
                            for ( auto j = 0; j < other.columns(); j++ ) {
                                set( i, j, other.get( i, j ) );
                            }
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual std::string id() const override {
                        return "NCPA Basic Dense Matrix";
                    }

                    friend void ::swap<ELEMENTTYPE>(
                        dense_matrix<ELEMENTTYPE>& a,
                        dense_matrix<ELEMENTTYPE>& b ) noexcept;

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        build_vector( size_t n = 0 ) const override {
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new dense_vector<ELEMENTTYPE>( n ) );
                    }

                    virtual size_t rows() const override {
                        return _elements.size();
                    }

                    virtual size_t columns() const override {
                        if ( _elements.size() > 0 ) {
                            return _elements.at( 0 ).size();
                        } else {
                            return 0;
                        }
                    }

                    // virtual ELEMENTTYPE& get( size_t row,
                    //                           size_t col ) override {
                    //     this->check_size( row, col );
                    //     return _elements[ row ][ col ];
                    // }

                    virtual const ELEMENTTYPE& get(
                        size_t row, size_t col ) const override {
                        this->check_size( row, col );
                        return _elements[ row ][ col ];
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_row( size_t row ) const override {
                        this->check_size( row, 0 );
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new dense_vector<ELEMENTTYPE>(
                                _elements[ row ] ) );
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_column( size_t column ) const override {
                        this->check_size( 0, column );
                        std::vector<ELEMENTTYPE> vcol( rows() );
                        for ( size_t row = 0; row < rows(); row++ ) {
                            vcol[ row ] = _elements[ row ][ column ];
                        }
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new dense_vector<ELEMENTTYPE>( vcol ) );
                    }

                    // virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                    //     get_diagonal( int offset = 0 ) const override {
                    //     size_t absoffset = (size_t)( std::abs( offset ) );
                    //     this->check_size( absoffset, absoffset );
                    //     // size_t ndiag = (size_t)( std::min( rows(),
                    //     columns()
                    //     // )
                    //     //                          - absoffset );
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
                    //         new dense_vector<ELEMENTTYPE>( diag ) );
                    // }

                    virtual abstract_matrix<ELEMENTTYPE>& clear() override {
                        _elements.clear();
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        clone() const override {
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new dense_matrix( *this ) );
                    }

                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        fresh_clone() const override {
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new dense_matrix() );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& resize(
                        size_t rows, size_t cols ) override {
                        _elements.resize( rows );
                        for ( size_t i = 0; i < rows; i++ ) {
                            _elements[ i ].resize( cols );
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

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
                        for ( size_t i = 0; i < rows(); i++ ) {
                            for ( size_t j = 0; j < columns(); j++ ) {
                                vals[ i ][ j ] = _elements[ i ][ j ];
                            }
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    // set individual elements
                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        size_t row, size_t col, ELEMENTTYPE val ) override {
                        this->check_size( row, col );
                        _elements[ row ][ col ] = val;
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        ELEMENTTYPE val ) override {
                        for ( auto it = _elements.begin();
                              it != _elements.end(); ++it ) {
                            it->set( val );
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    // set multiple elements in one row
                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, size_t nvals, const size_t *column_inds,
                        const ELEMENTTYPE *vals ) override {
                        this->check_size(
                            row, NCPA::math::max( column_inds, nvals ) );
                        for ( size_t i = 0; i < nvals; i++ ) {
                            _elements[ row ][ column_inds[ i ] ] = vals[ i ];
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    // set multiple elements in one column
                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, size_t nvals, const size_t *row_inds,
                        const ELEMENTTYPE *vals ) override {
                        this->check_size( NCPA::math::max( row_inds, nvals ),
                                          column );
                        for ( size_t i = 0; i < nvals; i++ ) {
                            _elements[ row_inds[ i ] ][ column ] = vals[ i ];
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    template<typename ANYTYPE,
                             ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                                           ELEMENTTYPE )>
                    abstract_matrix<ELEMENTTYPE>& scale( ANYTYPE val ) {
                        for ( auto it = _elements.begin();
                              it != _elements.end(); ++it ) {
                            it->scale( val );
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    // virtual abstract_matrix<ELEMENTTYPE>& scale(
                    //     const abstract_matrix<ELEMENTTYPE>& b ) override {
                    //     this->check_size( b );
                    //     for ( size_t row = 0; row < rows(); row++ ) {
                    //         for ( size_t col = 0; col < columns(); col++ ) {
                    //             _elements[ row ][ col ] *= b.get( row, col
                    //             );
                    //         }
                    //     }
                    //     return *dynamic_cast<abstract_matrix<ELEMENTTYPE>
                    //     *>(
                    //         this );
                    // }

                    // virtual abstract_matrix<ELEMENTTYPE>& add(
                    //     const abstract_matrix<ELEMENTTYPE>& b,
                    //     ELEMENTTYPE modifier = 1.0 ) override {
                    //     this->check_size( b );
                    //     for ( size_t row = 0; row < rows(); row++ ) {
                    //         for ( size_t col = 0; col < columns(); col++ ) {
                    //             _elements[ row ][ col ] += b.get( row, col
                    //             );
                    //         }
                    //     }
                    //     return *dynamic_cast<abstract_matrix<ELEMENTTYPE>
                    //     *>(
                    //         this );
                    // }

                    template<typename ANYTYPE,
                             ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                                           ELEMENTTYPE )>
                    abstract_matrix<ELEMENTTYPE>& add( ELEMENTTYPE b ) {
                        for ( auto it = _elements.begin();
                              it != _elements.end(); ++it ) {
                            it->add( b );
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& transpose()
                        override {
                        std::unique_ptr<abstract_matrix<ELEMENTTYPE>> orig
                            = clone();
                        resize( orig->columns(), orig->rows() );
                        for ( size_t i = 0; i < orig->rows(); i++ ) {
                            set_column( i, orig->get_row( i )->as_std() );
                        }
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& swap_rows(
                        size_t ind1, size_t ind2 ) override {
                        std::swap( _elements[ ind1 ], _elements[ ind2 ] );
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& swap_columns(
                        size_t ind1, size_t ind2 ) override {
                        auto col1 = get_column( ind1 );
                        auto col2 = get_column( ind2 );
                        set_column( ind1, *col2 );
                        set_column( ind2, *col1 );
                        return *dynamic_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual bool is_this_subclass(
                        const abstract_matrix<ELEMENTTYPE>& b )
                        const override {
                        if ( auto *derived
                             = dynamic_cast<const dense_matrix<ELEMENTTYPE> *>(
                                 &b ) ) {
                            return true;
                        } else {
                            return false;
                        }
                    }

                    // virtual abstract_vector<ELEMENTTYPE>& operator[](
                    //     size_t i ) override {
                    //     return _elements[ i ];
                    // }

                    // virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                    // transpose()
                    //     const override {
                    //     std::unique_ptr<abstract_matrix<ELEMENTTYPE>> trans
                    //         = fresh_clone();
                    //     trans->resize( columns(), rows() );
                    //     for ( size_t i = 0; i < rows(); i++ ) {
                    //         trans->set_column( i, get_row( i )->as_std() );
                    //     }
                    //     return trans;
                    // }

                    typename std::vector<dense_vector<ELEMENTTYPE>>::iterator
                        begin() noexcept {
                        return _elements.begin();
                    }

                    typename std::vector<dense_vector<ELEMENTTYPE>>::iterator
                        end() noexcept {
                        return _elements.end();
                    }

                    typename std::vector<
                        dense_vector<ELEMENTTYPE>>::const_iterator
                        cbegin() const noexcept {
                        return _elements.cbegin();
                    }

                    typename std::vector<
                        dense_vector<ELEMENTTYPE>>::const_iterator
                        cend() const noexcept {
                        return _elements.cend();
                    }

                private:
                    bool _finalized = false;
                    std::vector<dense_vector<ELEMENTTYPE>> _elements;
                    // std::vector<WrapperVector<ELEMENTTYPE>> _wrappers;
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::details::dense_matrix<T>& a,
                  NCPA::linear::details::dense_matrix<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::details::abstract_matrix<T>&>( a ),
            static_cast<NCPA::linear::details::abstract_matrix<T>&>( b ) );
    swap( a._elements, b._elements );
    swap( a._finalized, b._finalized );
}
