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
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Matrix );
    }
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::Matrix,
                                             ELEMENTTYPE );

template<typename ELEMENTTYPE>
std::ostream& operator<<( std::ostream& os,
                          const NCPA::linear::Matrix<ELEMENTTYPE>& obj );

NCPA_LINEARALGEBRA_DECLARE_FRIEND_BINARY_OPERATORS( NCPA::linear::Matrix,
                                                    ELEMENTTYPE )
// template<typename ELEMENTTYPE>
// NCPA::linear::Vector<ELEMENTTYPE> operator*(
//     const NCPA::linear::Matrix<ELEMENTTYPE>& c1,
//     const NCPA::linear::Vector<ELEMENTTYPE>& c2 );

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

                // std::unique_ptr<Matrix<ELEMENTTYPE>> clone() const {
                //     return std::unique_ptr<Matrix<ELEMENTTYPE>>(
                //         new Matrix<ELEMENTTYPE>( *this ) );
                // }

                // std::unique_ptr<Matrix<ELEMENTTYPE>> fresh_clone() const {
                //     std::unique_ptr<Matrix<ELEMENTTYPE>> fresh(
                //         new Matrix<ELEMENTTYPE>(
                //             std::unique_ptr<
                //                 details::abstract_matrix<ELEMENTTYPE>>() )
                //                 );
                //     // fresh->clear();
                //     return fresh;
                // }

                Matrix<ELEMENTTYPE>& identity() {
                    check_pointer();
                    if ( !is_square() ) {
                        throw std::invalid_argument(
                            "Cannot turn a non-square matrix into an identity "
                            "matrix" );
                    }
                    size_t ndiag = _ptr->diagonal_size( 0 );
                    _ptr->zero();
                    for ( size_t i = 0; i < ndiag; i++ ) {
                        _ptr->set( i, i, NCPA::math::one<ELEMENTTYPE>() );
                    }
                    return *this;
                }

                Matrix<ELEMENTTYPE>& identity( size_t rows, size_t cols = 0 ) {
                    check_pointer();
                    if ( cols == 0 ) {
                        cols = rows;
                    }
                    resize( rows, cols );
                    return identity();
                }

                virtual bool is_square() const {
                    return ( _ptr ? _ptr->is_square() : true );
                }

                virtual bool is_row_matrix() const { return ( rows() == 1 ); }

                virtual bool is_column_matrix() const {
                    return ( columns() == 1 );
                }

                virtual bool is_empty() const {
                    return ( _ptr ? _ptr->is_empty() : true );
                }

                virtual bool is_identity() const {
                    return ( _ptr ? _ptr->is_identity() : true );
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

                explicit operator bool() const {
                    return ( _ptr ? true : false );
                }

                virtual size_t rows() const {
                    return ( _ptr ? _ptr->rows() : 0 );
                }

                virtual size_t columns() const {
                    return ( _ptr ? _ptr->columns() : 0 );
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const {
                    return ( _ptr ? _ptr->get( row, col ) : _zero );
                }

                virtual const ELEMENTTYPE& get( size_t ind ) const {
                    if ( _ptr ) {
                        if ( is_row_matrix() ) {
                            return _ptr->get( 0, ind );
                        } else if ( is_column_matrix() ) {
                            return _ptr->get( ind, 0 );
                        } else {
                            throw std::range_error(
                                "Both dimensions must be specified for "
                                "non-row and non-column matrices" );
                        }
                    } else {
                        return _zero;
                    }
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
                            new Vector<ELEMENTTYPE>(
                                _ptr->get_diagonal( offset ) ) );
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

                virtual Matrix<ELEMENTTYPE>& copy(
                    const Matrix<ELEMENTTYPE>& M ) {
                    resize( M.rows(), M.columns() );
                    for ( size_t r = 0; r < rows(); r++ ) {
                        set_row( r, *( M.get_row( r ) ) );
                    }
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
                    // return set_row( row, vec.as_std() );
                    check_pointer();
                    _ptr->set_row( row, vec );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const Vector<ELEMENTTYPE>& vec ) {
                    return set_row( row, vec.as_std() );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::unique_ptr<Vector<ELEMENTTYPE>>& vec ) {
                    return set_row( row, *vec );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const Matrix<ELEMENTTYPE>& mat,
                    size_t matrow ) {
                    return set_row( row,
                                    mat.get_row_vector( matrow )->as_std() );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::unique_ptr<Matrix<ELEMENTTYPE>>& mat,
                    size_t matrow ) {
                    return set_row( row, *mat );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const Matrix<ELEMENTTYPE>& mat ) {
                    if ( !mat.is_row_matrix() ) {
                        throw std::logic_error( "Matrix is not a row matrix" );
                    }
                    return set_row( row, mat, 0 );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::unique_ptr<Matrix<ELEMENTTYPE>>& mat ) {
                    return set_row( row, *mat, 0 );
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
                    size_t col,
                    const std::unique_ptr<Vector<ELEMENTTYPE>>& vec ) {
                    return set_column( col, *vec );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const Matrix<ELEMENTTYPE>& mat,
                    size_t matcol ) {
                    return set_column(
                        col, mat.get_column_vector( matcol )->as_std() );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const std::unique_ptr<Matrix<ELEMENTTYPE>>& mat,
                    size_t matcol ) {
                    return set_column( col, *mat, matcol );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const Matrix<ELEMENTTYPE>& mat ) {
                    if ( !mat.is_column_matrix() ) {
                        throw std::logic_error(
                            "Matrix is not a column matrix" );
                    }
                    return set_column( col, mat, 0 );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const std::unique_ptr<Matrix<ELEMENTTYPE>>& mat ) {
                    return set_column( col, *mat );
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

                virtual Matrix<ELEMENTTYPE>& swap_rows( size_t ind1,
                                                        size_t ind2 ) {
                    check_pointer();
                    _ptr->swap_rows( ind1, ind2 );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& swap_columns( size_t ind1,
                                                           size_t ind2 ) {
                    check_pointer();
                    _ptr->swap_columns( ind1, ind2 );
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
                    auto newmat = Matrix<ELEMENTTYPE>(
                        _ptr->multiply( *( other._ptr ) ) );
                    swap( *this, newmat );
                    return *this;
                }

                // virtual Vector<ELEMENTTYPE>& multiply(
                //     const Vector<ELEMENTTYPE>& other ) {
                //     check_pointer();
                //     other.check_pointer();
                //     if ( columns() != size.rows() ) {
                //         throw std::invalid_argument(
                //             "Matrix-vector size mismatch: cannot multiply" );
                //     }
                //     auto newvec = Vector<ELEMENTTYPE>(
                //         _ptr->multiply( *( other.internal() ) ) );
                //     return newvec;
                // }

                virtual bool equals( const Matrix<ELEMENTTYPE>& other ) const {
                    return ( is_empty() || other.is_empty()
                                 ? false
                                 : _ptr->equals( *( other._ptr ) ) );
                }

                virtual Vector<ELEMENTTYPE>& operator[]( size_t ind ) {
                    if ( ind >= _ptr->rows() ) {
                        _ptr->resize( ind + 1, _ptr->columns() );
                    }
                    _wrappers[ ind ]
                        = WrapperVector<ELEMENTTYPE>( ( *_ptr )[ ind ] );
                    return _wrappers[ ind ];
                }

                virtual const details::abstract_vector<ELEMENTTYPE>&
                    operator[]( size_t ind ) const {
                    if ( ind >= _ptr->rows() ) {
                        std::ostringstream oss;
                        oss << "Index " << ind << " too large for matrix of "
                            << rows() << " rows.";
                        throw std::out_of_range( oss.str() );
                    }
                    return ( *_ptr )[ ind ];
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

                // friend operators
                friend std::ostream& operator<<(
                    std::ostream& os, const Matrix<ELEMENTTYPE>& mat ) {
                    os << "[";
                    for ( size_t r = 0; r < mat.rows(); r++ ) {
                        if ( r > 0 ) {
                            os << " ";
                        }
                        os << " [ ";
                        for ( size_t c = 0; c < mat.columns(); c++ ) {
                            if ( c > 0 ) {
                                os << ", ";
                            }
                            os << mat[ r ][ c ];
                        }
                        if ( r != mat.rows() - 1 ) {
                            os << " ]," << std::endl;
                        } else {
                            os << " ] ]";
                        }
                    }
                    return os;
                }

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

                friend Matrix<ELEMENTTYPE> operator+(
                    const Matrix<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    Matrix<ELEMENTTYPE> out( c1 );
                    out += c2;
                    return out;
                }

                friend Matrix<ELEMENTTYPE> operator+(
                    ELEMENTTYPE c1, const Matrix<ELEMENTTYPE>& c2 ) {
                    Matrix<ELEMENTTYPE> out( c2 );
                    out += c1;
                    return out;
                }

                friend Matrix<ELEMENTTYPE> operator-(
                    const Matrix<ELEMENTTYPE>& c1,
                    const Matrix<ELEMENTTYPE>& c2 ) {
                    Matrix<ELEMENTTYPE> out( c1 );
                    out -= c2;
                    return out;
                }

                friend Matrix<ELEMENTTYPE> operator-(
                    const Matrix<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
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

                // friend NCPA::linear::Matrix<ELEMENTTYPE> operator*(
                //     const Matrix<ELEMENTTYPE>& c1,
                //     const Vector<ELEMENTTYPE>& c2 ) {
                //     return c1.multiply( c2 );
                // }

                friend Matrix<ELEMENTTYPE> operator*(
                    const Matrix<ELEMENTTYPE>& c1, ELEMENTTYPE c2 ) {
                    Matrix<ELEMENTTYPE> out( c1 );
                    out *= c2;
                    return out;
                }

                friend Matrix<ELEMENTTYPE> operator*(
                    ELEMENTTYPE c1, const Matrix<ELEMENTTYPE>& c2 ) {
                    Matrix<ELEMENTTYPE> out( c2 );
                    out *= c1;
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
                std::map<size_t, WrapperVector<ELEMENTTYPE>> _wrappers;
                const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::Matrix<T>& a,
                  NCPA::linear::Matrix<T>& b ) noexcept {
    // using std::swap;
    a._ptr.swap( b._ptr );
    swap( a._wrappers, b._wrappers );
}
