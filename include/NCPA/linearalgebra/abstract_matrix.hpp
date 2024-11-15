#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/vector.hpp"
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

// forward declarations for operators and friend functions
namespace NCPA {
    namespace linear {
        namespace details {
            template<typename ELEMENTTYPE>
            class abstract_matrix;
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::abstract_matrix, ELEMENTTYPE );

namespace NCPA {
    namespace linear {

        namespace details {
            template<typename ELEMENTTYPE>
            class abstract_matrix {
                public:
                    virtual ~abstract_matrix()     = default;
                    virtual std::string id() const = 0;
                    friend void ::swap<ELEMENTTYPE>(
                        abstract_matrix<ELEMENTTYPE>& a,
                        abstract_matrix<ELEMENTTYPE>& b ) noexcept;

                    virtual size_t rows() const                        = 0;
                    virtual size_t columns() const                     = 0;
                    virtual ELEMENTTYPE& get( size_t row, size_t col ) = 0;
                    virtual const ELEMENTTYPE& get( size_t row,
                                                    size_t col ) const
                        = 0;

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_row( size_t row ) const = 0;
                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_column( size_t column ) const = 0;
                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_diagonal( int offset ) const = 0;

                    virtual abstract_matrix<ELEMENTTYPE>& clear() = 0;
                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        clone() const = 0;
                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        fresh_clone() const = 0;
                    virtual abstract_matrix<ELEMENTTYPE>& resize( size_t rows,
                                                                  size_t cols )
                        = 0;
                    virtual abstract_matrix<ELEMENTTYPE>& as_array(
                        size_t& nrows, size_t& ncols, ELEMENTTYPE **& vals )
                        = 0;

                    // set individual elements
                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        size_t row, size_t col, ELEMENTTYPE val )
                        = 0;

                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        ELEMENTTYPE val )
                        = 0;


                    // set multiple elements in one row
                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, size_t nvals, const size_t *column_inds,
                        const ELEMENTTYPE *vals )
                        = 0;

                    // set multiple elements in one column
                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, size_t nvals, const size_t *row_inds,
                        const ELEMENTTYPE *vals )
                        = 0;

                    virtual abstract_matrix<ELEMENTTYPE>& transpose() = 0;

                    // implementations, not abstract
                    // metaconstructors
                    virtual abstract_matrix<ELEMENTTYPE>& identity( size_t nrows, size_t ncols ) {
                        clear();
                        resize( nrows, ncols );
                        ELEMENTTYPE one = NCPA::math::one<ELEMENTTYPE>();
                        for (auto i = 0; i < diagonal_size(0); i++) {
                            set(i,i,one);
                        }
                        return *this;
                    }

                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        multiply(
                            const abstract_matrix<ELEMENTTYPE>& b ) const {
                        check_size_for_mult( b );
                        std::unique_ptr<abstract_matrix<ELEMENTTYPE>> product
                            = fresh_clone();
                        product->resize( rows(), b.columns() );
                        for ( size_t row = 0; row < product->rows(); row++ ) {
                            for ( size_t col = 0; col < product->columns();
                                  col++ ) {
                                product->set( row, col,
                                              get_row( row )->dot(
                                                  *( b.get_column( col ) ) ) );
                            }
                        }
                        return product;
                    }

                    template<typename ANYTYPE,
                             ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                                           ELEMENTTYPE )>
                     abstract_matrix<ELEMENTTYPE>& scale(
                        ANYTYPE val ) {
                        for ( size_t row = 0; row < rows(); row++ ) {
                            for ( size_t col = 0; col < columns(); col++ ) {
                                get( row, col ) *= val;
                            }
                        }
                        return *this;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& scale(
                        const abstract_matrix<ELEMENTTYPE>& b ) {
                        check_size( b );
                        for ( size_t row = 0; row < rows(); row++ ) {
                            for ( size_t col = 0; col < columns(); col++ ) {
                                get( row, col ) *= b.get( row, col );
                            }
                        }
                        return *this;
                    }

                    
                    virtual abstract_matrix<ELEMENTTYPE>& add(
                        const abstract_matrix<ELEMENTTYPE>& b,
                        ELEMENTTYPE modifier = 1.0 ) {
                        check_size( b );
                        for ( size_t row = 0; row < rows(); row++ ) {
                            for ( size_t col = 0; col < columns(); col++ ) {
                                get( row, col )
                                    += b.get( row, col ) * modifier;
                            }
                        }
                        return *this;
                    }

                    template<typename ANYTYPE,
                             ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                                           ELEMENTTYPE )>
                     abstract_matrix<ELEMENTTYPE>& add( ANYTYPE b ) {
                        for ( size_t row = 0; row < rows(); row++ ) {
                            for ( size_t col = 0; col < columns(); col++ ) {
                                get( row, col ) += b;
                            }
                        }
                        return *this;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& zero( size_t row,
                                                                size_t col ) {
                        return set( row, col,
                                    NCPA::math::zero<ELEMENTTYPE>() );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& zero() {
                        return scale( NCPA::math::zero<ELEMENTTYPE>() );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, const std::vector<size_t>& column_inds,
                        const std::vector<ELEMENTTYPE>& vals ) {
                        return set_row( row, column_inds.size(),
                                        &column_inds[ 0 ], &vals[ 0 ] );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, std::initializer_list<size_t> column_inds,
                        std::initializer_list<ELEMENTTYPE>& vals ) {
                        return set_row( row,
                                        std::vector<size_t>( column_inds ),
                                        std::vector<ELEMENTTYPE>( vals ) );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, const std::vector<ELEMENTTYPE>& vals ) {
                        return set_row(
                            row,
                            NCPA::arrays::index_vector<size_t>( vals.size() ),
                            vals );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, ELEMENTTYPE val ) {
                        return set_row(
                            row, std::vector<ELEMENTTYPE>( columns(), val ) );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, const std::vector<size_t>& row_inds,
                        const std::vector<ELEMENTTYPE>& vals ) {
                        return set_column( column, row_inds.size(),
                                           &row_inds[ 0 ], &vals[ 0 ] );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, std::initializer_list<size_t> row_inds,
                        std::initializer_list<ELEMENTTYPE>& vals ) {
                        return set_column( column,
                                           std::vector<size_t>( row_inds ),
                                           std::vector<ELEMENTTYPE>( vals ) );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, const std::vector<ELEMENTTYPE>& vals ) {
                        return set_column(
                            column,
                            NCPA::arrays::index_vector<size_t>( vals.size() ),
                            vals );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, ELEMENTTYPE val ) {
                        return set_column(
                            column, std::vector<ELEMENTTYPE>( rows(), val ) );
                    }

                    // set the diagonal or off-diagonal
                    virtual abstract_matrix<ELEMENTTYPE>& set_diagonal(
                        size_t nvals, const ELEMENTTYPE *vals,
                        int offset = 0 ) {
                        size_t absoffset = (size_t)std::abs( offset );
                        if ( absoffset > max_off_diagonal() ) {
                            std::ostringstream oss;
                            oss << absoffset
                                << "'th off-diagonal requested, but"
                                   " matrix of size "
                                << rows() << "x" << columns() << " only has "
                                << max_off_diagonal() << " off-diagonals";
                            throw std::range_error( oss.str() );
                        }
                        if ( nvals > diagonal_size( offset ) ) {
                            std::ostringstream oss;
                            oss << nvals << " samples is too many for the "
                                << absoffset << "'th diagonal of a " << rows()
                                << "x" << columns() << " matrix; max size is "
                                << diagonal_size( offset ) << ".";
                            throw std::invalid_argument( oss.str() );
                        }
                        for ( size_t i = 0; i < nvals; i++ ) {
                            if ( offset >= 0 ) {
                                // upper triangle
                                set( i, i + absoffset, vals[ i ] );
                            } else {
                                // lower triangle
                                set( i + absoffset, i, vals[ i ] );
                            }
                        }
                        return *this;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_diagonal(
                        const std::vector<ELEMENTTYPE>& vals,
                        int offset = 0 ) {
                        return set_diagonal(
                            std::min( vals.size(), diagonal_size( offset ) ),
                            NCPA::arrays::as_array( vals ), offset );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_diagonal(
                        std::initializer_list<ELEMENTTYPE>& vals,
                        int offset = 0 ) {
                        return set_diagonal( std::vector<ELEMENTTYPE>( vals ),
                                             offset );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_diagonal(
                        ELEMENTTYPE val, int offset = 0 ) {
                        return set_diagonal(
                            std::vector<ELEMENTTYPE>( diagonal_size( offset ),
                                                      val ),
                            offset );
                    }

                    // template<typename OTHERTYPE>
                    // bool equals(
                    //     const abstract_matrix<OTHERTYPE>& other ) const {
                    //     return false;
                    // }
                    template<typename ANYTYPE,
                             ENABLE_IF_TU( std::is_convertible, ANYTYPE,
                                           ELEMENTTYPE )>
                    bool equals(
                        const abstract_matrix<ANYTYPE>& other ) const {
                        if ( rows() != other.rows() ) {
                            return false;
                        }
                        if ( columns() != other.columns() ) {
                            return false;
                        }

                        for ( auto r = 0; r < rows(); r++ ) {
                            for ( auto c = 0; c < columns(); c++ ) {
                                if ( get( r, c ) != other.get( r, c ) ) {
                                    return false;
                                }
                            }
                        }
                        return true;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& operator+=(
                        const abstract_matrix<ELEMENTTYPE>& other ) {
                        this->add( other );
                        return *this;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& operator+=(
                        const ELEMENTTYPE& other ) {
                        this->add( other );
                        return *this;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& operator-=(
                        const abstract_matrix<ELEMENTTYPE>& other ) {
                        this->add( other, -1.0 );
                        return *this;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& operator-=(
                        const ELEMENTTYPE& other ) {
                        this->add( -other );
                        return *this;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& operator*=(
                        const abstract_matrix<ELEMENTTYPE>& other ) {
                        this->scale( other );
                        return *this;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& operator*=(
                        const ELEMENTTYPE& other ) {
                        this->scale( other );
                        return *this;
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& operator/=(
                        const ELEMENTTYPE& other ) {
                        this->scale( 1.0 / other );
                        return *this;
                    }

                    // virtual abstract_vector<ELEMENTTYPE>& operator[](
                    //     size_t i ) {
                    //     return *get_row( i );
                    // }

                    // virtual const abstract_vector<ELEMENTTYPE>&
                    //     operator[]( size_t i ) const {
                    //     return *get_row( i );
                    // }

                    friend bool operator==(
                        const abstract_matrix<ELEMENTTYPE>& a,
                        const abstract_matrix<ELEMENTTYPE>& b ) {
                        return a.equals( b );
                    }

                    friend bool operator!=(
                        const abstract_matrix<ELEMENTTYPE>& a,
                        const abstract_matrix<ELEMENTTYPE>& b ) {
                        return !( a.equals( b ) );
                    }

                    virtual void check_size(
                        const abstract_matrix<ELEMENTTYPE>& b ) const {
                        if ( rows() != b.rows() || columns() != b.columns() ) {
                            std::ostringstream oss;
                            oss << "Size mismatch between Matrices: ["
                                << rows() << "," << columns() << "] vs ["
                                << b.rows() << "," << b.columns() << "]";
                            throw std::invalid_argument( oss.str() );
                        }
                    }

                    virtual void check_size_for_mult(
                        const abstract_matrix<ELEMENTTYPE>& b ) const {
                        if ( columns() != b.rows() ) {
                            std::ostringstream oss;
                            oss << "Multiplication size mismatch between "
                                   "Matrices: ["
                                << rows() << "," << columns() << "] vs ["
                                << b.rows() << "," << b.columns() << "]";
                            throw std::invalid_argument( oss.str() );
                        }
                    }

                    virtual void check_size( size_t nrows,
                                             size_t ncols ) const {
                        std::ostringstream oss;
                        if ( rows() == 0 || columns() == 0 ) {
                            throw std::logic_error(
                                "Matrix has not been initialized!" );
                        }

                        if ( rows() <= nrows ) {
                            oss << "Index " << nrows
                                << " too large for Matrix with " << rows()
                                << " rows";
                            throw std::range_error( oss.str() );
                        }

                        if ( columns() <= ncols ) {
                            oss << "Index " << ncols
                                << " too large for Matrix with " << columns()
                                << " columns";
                            throw std::range_error( oss.str() );
                        }
                    }

                    virtual size_t diagonal_size( int offset ) const {
                        return max_off_diagonal() + 1
                             - (size_t)std::abs( offset );
                    }

                    virtual size_t max_off_diagonal() const {
                        check_size( 0, 0 );
                        return std::max( rows(), columns() ) - 1;
                    }

                    virtual bool is_empty() const {
                        return ( rows() == 0 && columns() == 0 );
                    }

                    virtual bool is_diagonal() const {
                        if ( is_empty() ) {
                            return true;
                        }
                        ELEMENTTYPE zero = NCPA::math::zero<ELEMENTTYPE>();
                        for ( auto i = 0; i < rows(); i++ ) {
                            for ( auto j = 0; j < columns(); j++ ) {
                                if ( i != j && get( i, j ) != zero ) {
                                    return false;
                                }
                            }
                        }
                        return true;
                    }

                    virtual bool is_tridiagonal() const {
                        if ( is_empty() ) {
                            return true;
                        }
                        ELEMENTTYPE zero = NCPA::math::zero<ELEMENTTYPE>();
                        int nr           = (int)rows();
                        int nc           = (int)columns();
                        for ( int i = 0; i < rows(); i++ ) {
                            for ( int j = 0; j < columns(); j++ ) {
                                if ( std::abs( i - j ) > 1
                                     && get( i, j ) != zero ) {
                                    return false;
                                }
                            }
                        }
                        return true;
                    }

                    virtual bool is_upper_triangular() const {
                        if ( is_empty() ) {
                            return true;
                        }
                        ELEMENTTYPE zero = NCPA::math::zero<ELEMENTTYPE>();
                        for ( size_t i = 1; i < rows(); i++ ) {
                            for ( size_t j = 0; j < i; j++ ) {
                                if ( get( i, j ) != zero ) {
                                    return false;
                                }
                            }
                        }
                        return true;
                    }

                    virtual bool is_lower_triangular() const {
                        if ( is_empty() ) {
                            return true;
                        }
                        ELEMENTTYPE zero = NCPA::math::zero<ELEMENTTYPE>();
                        for ( size_t i = 0; i < rows() - 1; i++ ) {
                            for ( size_t j = i + 1; j < columns(); j++ ) {
                                if ( get( i, j ) != zero ) {
                                    return false;
                                }
                            }
                        }
                        return true;
                    }
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::details::abstract_matrix<T>& a,
                  NCPA::linear::details::abstract_matrix<T>& b ) noexcept {}
