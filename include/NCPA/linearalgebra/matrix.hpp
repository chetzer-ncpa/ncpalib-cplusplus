#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_matrix.hpp"
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

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE(
            Matrix, details::abstract_matrix );
    }
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::Matrix,
                                             ELEMENTTYPE );

template<typename ELEMENTTYPE>
NCPA::linear::Matrix<ELEMENTTYPE> operator+(
    const NCPA::linear::Matrix<ELEMENTTYPE>& c1,
    const NCPA::linear::Matrix<ELEMENTTYPE>& c2 );
template<typename ELEMENTTYPE>
NCPA::linear::Matrix<ELEMENTTYPE> operator-(
    const NCPA::linear::Matrix<ELEMENTTYPE>& c1,
    const NCPA::linear::Matrix<ELEMENTTYPE>& c2 );
template<typename ELEMENTTYPE>
NCPA::linear::Matrix<ELEMENTTYPE> operator*(
    const NCPA::linear::Matrix<ELEMENTTYPE>& c1,
    const NCPA::linear::Matrix<ELEMENTTYPE>& c2 );

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class Matrix<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC> {
            public:
                Matrix() {}

                Matrix( std::unique_ptr<details::abstract_matrix<ELEMENTTYPE>>
                            ptr ) :
                    Matrix<ELEMENTTYPE>() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                Matrix( const Matrix<ELEMENTTYPE>& other ) :
                    Matrix<ELEMENTTYPE>() {
                    _ptr = std::move( other._ptr->clone() );
                }

                /**
                 * Move constructor.
                 * @param source The vector to assimilate.
                 */
                Matrix( Matrix<ELEMENTTYPE>&& source ) noexcept :
                    Matrix<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                virtual ~Matrix() {}

                friend void ::swap<ELEMENTTYPE>(
                    Matrix<ELEMENTTYPE>& a, Matrix<ELEMENTTYPE>& b ) noexcept;

                /**
                 * Assignment operator.
                 * @param other The vector to assign to this.
                 */
                Matrix<ELEMENTTYPE>& operator=( Matrix<ELEMENTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual bool is_row_matrix() const {
                    return ( columns() == 1 );
                }

                virtual bool is_column_matrix() const {
                    return ( rows() == 1 );
                }

                virtual bool is_empty() const {
                    return ( _ptr ? _ptr->is_empty() : true );
                }

                virtual bool is_diagonal() const {
                    return ( _ptr ? _ptr->is_diagonal() : true );
                }

                virtual bool is_tridiagonal() const {
                    return ( _ptr ? _ptr->is_tridiagonal() : true );
                }

                virtual bool is_lower_triangular() const {
                    return ( _ptr ? _ptr->is_lower_triangular() : true );
                }

                virtual bool is_upper_triangular() const {
                    return ( _ptr ? _ptr->is_upper_triangular() : true );
                }

                virtual size_t rows() const {
                    return ( _ptr ? _ptr->rows() : 0 );
                }

                virtual size_t columns() const {
                    return ( _ptr ? _ptr->columns() : 0 );
                }

                virtual ELEMENTTYPE& get( size_t row, size_t col ) {
                    check_pointer();
                    return _ptr->get( row, col );
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const {
                    return ( _ptr ? _ptr->get( row, col ) : _zero );
                }

                virtual std::unique_ptr<Matrix<ELEMENTTYPE>> get_row(
                    size_t row ) const {
                    if ( _ptr ) {
                        std::unique_ptr<Matrix<ELEMENTTYPE>> rowmat(
                            new Matrix<ELEMENTTYPE>( *this ) );
                        rowmat->resize( 1, columns() )
                            .zero()
                            .set_row( 0, _ptr->get_row( row ) );
                        return rowmat;
                    } else {
                        return std::unique_ptr<Matrix<ELEMENTTYPE>>();
                    }
                }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_row_vector(
                    size_t row ) const {
                    if ( _ptr ) {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>(
                            new Vector<ELEMENTTYPE>( _ptr->get_row( row ) ) );
                    } else {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>();
                    }
                }

                virtual std::unique_ptr<Matrix<ELEMENTTYPE>> get_column(
                    size_t col ) const {
                    if ( _ptr ) {
                        std::unique_ptr<Matrix<ELEMENTTYPE>> colmat(
                            new Matrix<ELEMENTTYPE>( *this ) );
                        colmat->resize( rows(), 1 )
                            .zero()
                            .set_column( 0, _ptr->get_column( col ) );
                        return colmat;
                    } else {
                        return std::unique_ptr<Matrix<ELEMENTTYPE>>();
                    }
                }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_column_vector(
                    size_t col ) const {
                    if ( _ptr ) {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>(
                            new Vector<ELEMENTTYPE>(
                                _ptr->get_column( col ) ) );
                    } else {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>();
                    }
                }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_diagonal(
                    int offset = 0 ) const {
                    if ( _ptr ) {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>(
                            new Vector<ELEMENTTYPE>( _ptr->get_diagonal( offset ) ) );
                    } else {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>();
                    }
                }

                virtual Matrix<ELEMENTTYPE>& clear() {
                    if ( _ptr ) {
                        _ptr->clear();
                    }
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& resize( size_t rows,
                                                     size_t cols ) {
                    check_pointer();
                    _ptr->resize( rows, cols );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& as_array( size_t& nrows,
                                                       size_t& ncols,
                                                       ELEMENTTYPE **& vals ) {
                    check_pointer();
                    _ptr->as_array( nrows, ncols, vals );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set( size_t row, size_t col,
                                                  ELEMENTTYPE val ) {
                    check_pointer();
                    _ptr->set( row, col, val );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set( ELEMENTTYPE val ) {
                    check_pointer();
                    _ptr->set( val );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& zero( size_t row, size_t col ) {
                    check_pointer();
                    _ptr->zero( row, col );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& zero() {
                    check_pointer();
                    _ptr->zero();
                    return *this;
                }

                /*
                set_row options:
                    0. constant
                    1. array
                    2. std::vector with indices
                    3. std::vector without indices
                    4. initializer list with indices
                    5. initializer list without indices
                    6. abstract_vector
                    7. Vector
                    8. Matrix, specify row
                    9. Row matrix
                */
                virtual Matrix<ELEMENTTYPE>& set_row( size_t row,
                                                      ELEMENTTYPE val ) {
                    check_pointer();
                    _ptr->set_row( row, val );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, size_t nvals, const size_t *indices,
                    const ELEMENTTYPE *vals ) {
                    check_pointer();
                    _ptr->set_row( row, nvals, indices, vals );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::vector<size_t>& rowinds,
                    const std::vector<ELEMENTTYPE>& vals ) {
                    check_pointer();
                    _ptr->set_row( row, rowinds, vals );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::vector<ELEMENTTYPE>& vals ) {
                    check_pointer();
                    _ptr->set_row( row, vals );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::initializer_list<size_t> rowinds,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    check_pointer();
                    _ptr->set_row( row, rowinds, vals );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    return set_row( row, std::vector<ELEMENTTYPE>( vals ) );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const details::abstract_vector<ELEMENTTYPE>& vec ) {
                    return set_row( row, vec.as_std() );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const Vector<ELEMENTTYPE>& vec ) {
                    return set_row( row, vec.as_std() );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const Matrix<ELEMENTTYPE>& mat,
                    size_t matrow ) {
                    return set_row( row,
                                    mat.get_row_vector( matrow )->as_std() );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const Matrix<ELEMENTTYPE>& mat ) {
                    if ( !mat.is_row_matrix() ) {
                        throw std::logic_error( "Matrix is not a row matrix" );
                    }
                    return set_row( row, mat, 0 );
                }

                /*
                set_column options:
                    0. constant
                    1. array
                    2. std::vector with indices
                    3. std::vector without indices
                    4. initializer list with indices
                    5. initializer list without indices
                    6. abstract_vector
                    7. Vector
                    8. Matrix, specify row
                    9. Row matrix
                */
                virtual Matrix<ELEMENTTYPE>& set_column( size_t col,
                                                         ELEMENTTYPE val ) {
                    check_pointer();
                    _ptr->set_column( col, val );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, size_t nvals, const size_t *indices,
                    const ELEMENTTYPE *vals ) {
                    check_pointer();
                    _ptr->set_column( col, nvals, indices, vals );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::vector<size_t>& colinds,
                    const std::vector<ELEMENTTYPE>& vals ) {
                    check_pointer();
                    _ptr->set_column( col, colinds, vals );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::vector<ELEMENTTYPE>& vals ) {
                    check_pointer();
                    _ptr->set_column( col, vals );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::initializer_list<size_t> colinds,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    check_pointer();
                    _ptr->set_column( col, colinds, vals );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    return set_column( col, std::vector<ELEMENTTYPE>( vals ) );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const details::abstract_vector<ELEMENTTYPE>& vec ) {
                    return set_column( col, vec.as_std() );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const Vector<ELEMENTTYPE>& vec ) {
                    return set_column( col, vec.as_std() );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const Matrix<ELEMENTTYPE>& mat,
                    size_t matcol ) {
                    return set_column(
                        col, mat.get_column_vector( matcol )->as_std() );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const Matrix<ELEMENTTYPE>& mat ) {
                    if ( !mat.is_column_matrix() ) {
                        throw std::logic_error(
                            "Matrix is not a column matrix" );
                    }
                    return set_column( col, mat, 0 );
                }

                /*
                    set_diagonal options:
                        0. constant
                        1. array
                        2. std::vector with indices
                        3. std::vector without indices
                        4. initializer list with indices
                        5. initializer list without indices
                        6. abstract_vector
                        7. Vector
                */
                virtual Matrix<ELEMENTTYPE>& set_diagonal( ELEMENTTYPE val,
                                                           int offset = 0 ) {
                    check_pointer();
                    _ptr->set_diagonal( val, offset );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    size_t nvals, const ELEMENTTYPE *vals, int offset = 0 ) {
                    check_pointer();
                    _ptr->set_diagonal( nvals, vals, offset );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const std::vector<ELEMENTTYPE>& vals, int offset = 0 ) {
                    check_pointer();
                    _ptr->set_diagonal( vals, offset );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const std::initializer_list<ELEMENTTYPE> vals,
                    int offset = 0 ) {
                    check_pointer();
                    _ptr->set_diagonal( vals, offset );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const details::abstract_vector<ELEMENTTYPE>& vec,
                    int offset = 0 ) {
                    return set_diagonal( vec.as_std(), offset );
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const Vector<ELEMENTTYPE>& vec, int offset = 0 ) {
                    return set_diagonal( vec.as_std(), offset );
                }

                virtual Matrix<ELEMENTTYPE>& transpose() {
                    check_pointer();
                    _ptr->transpose();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& add(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check_pointer();
                    other.check_pointer();
                    check_same_size( other );
                    *_ptr += *( other._ptr );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& add( ELEMENTTYPE other ) {
                    check_pointer();
                    *_ptr += other;
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& subtract(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check_pointer();
                    other.check_pointer();
                    check_same_size( other );
                    *_ptr -= *( other._ptr );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& subtract( ELEMENTTYPE other ) {
                    check_pointer();
                    *_ptr -= other;
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& scale(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check_pointer();
                    other.check_pointer();
                    check_same_size( other );
                    *_ptr *= *( other._ptr );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& scale( ELEMENTTYPE other ) {
                    check_pointer();
                    *_ptr *= other;
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& multiply(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check_pointer();
                    other.check_pointer();
                    if ( columns() != other.rows() ) {
                        throw std::invalid_argument(
                            "Matrix size mismatch: cannot multiply" );
                    }
                    _ptr->multiply( *( other._ptr ) );
                    return *this;
                }

                virtual bool equals( const Matrix<ELEMENTTYPE>& other ) const {
                    return ( is_empty() || other.is_empty()
                                 ? false
                                 : _ptr->equals( *( other._ptr ) ) );
                }

                // assignment operators
                virtual Matrix<ELEMENTTYPE>& operator+=(
                    const Matrix<ELEMENTTYPE>& other ) {
                    this->add( other );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& operator+=(
                    const ELEMENTTYPE& other ) {
                    this->add( other );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& operator-=(
                    const Matrix<ELEMENTTYPE>& other ) {
                    this->subtract( other );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& operator-=(
                    const ELEMENTTYPE& other ) {
                    this->add( -other );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& operator*=(
                    const Matrix<ELEMENTTYPE>& other ) {
                    this->multiply( other );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& operator*=(
                    const ELEMENTTYPE& other ) {
                    this->scale( other );
                    return *this;
                }

                // friend binary operators
                friend bool operator==( const Matrix<ELEMENTTYPE>& a,
                                        const Matrix<ELEMENTTYPE>& b ) {
                    return a.equals( b );
                }

                friend bool operator!=( const Matrix<ELEMENTTYPE>& a,
                                        const Matrix<ELEMENTTYPE>& b ) {
                    return !( a.equals( b ) );
                }

                friend Matrix<ELEMENTTYPE> operator+(
                    const Matrix<ELEMENTTYPE>& c1,
                    const Matrix<ELEMENTTYPE>& c2 ) {
                    Matrix<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend Matrix<ELEMENTTYPE> operator-(
                    const Matrix<ELEMENTTYPE>& c1,
                    const Matrix<ELEMENTTYPE>& c2 ) {
                    Matrix<ELEMENTTYPE> out( c1 );
                    out -= c2;
                    return out;
                }

                friend NCPA::linear::Matrix<ELEMENTTYPE> operator*(
                    const Matrix<ELEMENTTYPE>& c1,
                    const Matrix<ELEMENTTYPE>& c2 ) {
                    Matrix<ELEMENTTYPE> out( c1 );
                    out *= c2;
                    return out;
                }

                virtual void check_pointer() const {
                    if ( !_ptr ) {
                        throw std::logic_error(
                            "Matrix: Internal pointer has not been set!" );
                    }
                }

                virtual void check_same_size(
                    const Matrix<ELEMENTTYPE>& other ) const {
                    if ( rows() != other.rows()
                         || columns() != other.columns() ) {
                        throw std::invalid_argument(
                            "Matrices are not the same size!" );
                    }
                }

            private:
                std::unique_ptr<details::abstract_matrix<ELEMENTTYPE>> _ptr;
                const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::Matrix<T>& a,
                  NCPA::linear::Matrix<T>& b ) noexcept {
    // using std::swap;
    a._ptr.swap( b._ptr );
}
