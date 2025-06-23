#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/Vector.hpp"
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

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class Matrix<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC> {
            public:
                friend class Vector<ELEMENTTYPE>;
                friend class LUDecomposition<ELEMENTTYPE>;
                friend class BlockMatrix<ELEMENTTYPE>;

                Matrix() {}

                Matrix( std::unique_ptr<abstract_matrix<ELEMENTTYPE>> ptr ) :
                    Matrix<ELEMENTTYPE>() {
                    _ptr = std::move( ptr );
                }

                Matrix( const abstract_matrix<ELEMENTTYPE>& mat ) :
                    Matrix<ELEMENTTYPE>() {
                    _ptr = std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                        mat.clone() );
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

                Matrix<ELEMENTTYPE>& copy( const Matrix<ELEMENTTYPE>& other ) {
                    if ( _ptr ) {
                        _ptr->copy( *other._ptr );
                    } else {
                        _ptr = std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            other._ptr->clone() );
                    }
                    if ( other._lu ) {
                        _lu = std::unique_ptr<LUDecomposition<ELEMENTTYPE>>(
                            new LUDecomposition<ELEMENTTYPE>( *other._lu ) );
                    }
                    return *this;
                }

                Matrix<ELEMENTTYPE>& identity() {
                    check_pointer();
                    if ( !is_square() ) {
                        throw std::invalid_argument(
                            "Cannot turn a non-square matrix into an "
                            "identity "
                            "matrix" );
                    }
                    this->_clear_cache();
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

                virtual bool is_zero() const {
                    return ( _ptr ? _ptr->is_zero() : true );
                }

                virtual bool is_identity() const {
                    return ( _ptr ? _ptr->is_identity() : true );
                }

                virtual bool is_symmetric() const {
                    return ( _ptr ? _ptr->is_symmetric() : true );
                }

                virtual bool is_diagonal() const {
                    return ( _ptr ? _ptr->is_diagonal() : true );
                }

                virtual bool is_band_diagonal() const {
                    return ( _ptr ? _ptr->is_band_diagonal() : true );
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

                virtual size_t lower_bandwidth() const {
                    return ( _ptr ? _ptr->lower_bandwidth() : 0 );
                }

                virtual size_t upper_bandwidth() const {
                    return ( _ptr ? _ptr->upper_bandwidth() : 0 );
                }

                virtual size_t bandwidth() const {
                    return ( _ptr ? _ptr->bandwidth() : 1 );
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

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_row(
                    size_t row ) const {
                    if ( _ptr ) {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>(
                            new Vector<ELEMENTTYPE>( _ptr->get_row( row ) ) );
                    } else {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>();
                    }
                }

                // virtual std::unique_ptr<Matrix<ELEMENTTYPE>> get_column(
                //     size_t col ) const {
                //     if ( _ptr ) {
                //         std::unique_ptr<Matrix<ELEMENTTYPE>> colmat(
                //             new Matrix<ELEMENTTYPE>( *this ) );
                //         NCPA_DEBUG << "Created column matrix wrapper" <<
                //         std::endl; colmat->resize( rows(), 1 );
                //         NCPA_DEBUG
                //         << "Resize complete" << std::endl;
                //         colmat->zero(); NCPA_DEBUG << "Zero complete" <<
                //         std::endl; auto colvec = _ptr->get_column( col
                //         ); NCPA_DEBUG << "Got column vector" <<
                //         std::endl; for (auto i = 0; i < colvec->size();
                //         i++) {
                //             NCPA_DEBUG << "colvec[" << i << "] = " <<
                //             colvec->get( i ) << std::endl;
                //         }
                //         NCPA_DEBUG << "Calling set_column()" <<
                //         std::endl; colmat->set_column( 0, *colvec );
                //         NCPA_DEBUG << "Set column done" << std::endl;
                //         return colmat;
                //     } else {
                //         return std::unique_ptr<Matrix<ELEMENTTYPE>>();
                //     }
                // }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_column(
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
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& resize( size_t rows,
                                                     size_t cols ) {
                    check_pointer();
                    if ( rows != this->rows() || cols != this->columns() ) {
                        this->_clear_cache();
                        _ptr->resize( rows, cols );
                    }
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
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set( ELEMENTTYPE val ) {
                    check_pointer();
                    _ptr->set( val );
                    this->_clear_cache();
                    return *this;
                }


                virtual Matrix<ELEMENTTYPE>& zero( size_t row, size_t col ) {
                    check_pointer();
                    _ptr->zero( row, col );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& zero() {
                    check_pointer();
                    _ptr->zero();
                    this->_clear_cache();
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
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, size_t nvals, const size_t *indices,
                    const ELEMENTTYPE *vals ) {
                    check_pointer();
                    _ptr->set_row( row, nvals, indices, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::vector<size_t>& rowinds,
                    const std::vector<ELEMENTTYPE>& vals ) {
                    check_pointer();
                    _ptr->set_row( row, rowinds, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::vector<ELEMENTTYPE>& vals ) {
                    check_pointer();
                    _ptr->set_row( row, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::initializer_list<size_t> rowinds,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    check_pointer();
                    _ptr->set_row( row, rowinds, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    return set_row( row, std::vector<ELEMENTTYPE>( vals ) );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const abstract_vector<ELEMENTTYPE>& vec ) {
                    // return set_row( row, vec.as_std() );
                    check_pointer();
                    _ptr->set_row( row, vec );
                    this->_clear_cache();
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
                    return set_row( row, mat.get_row( matrow )->as_std() );
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
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, size_t nvals, const size_t *indices,
                    const ELEMENTTYPE *vals ) {
                    check_pointer();
                    _ptr->set_column( col, nvals, indices, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::vector<size_t>& colinds,
                    const std::vector<ELEMENTTYPE>& vals ) {
                    check_pointer();
                    _ptr->set_column( col, colinds, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::vector<ELEMENTTYPE>& vals ) {
                    check_pointer();
                    _ptr->set_column( col, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::initializer_list<size_t> colinds,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    check_pointer();
                    _ptr->set_column( col, colinds, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    return set_column( col, std::vector<ELEMENTTYPE>( vals ) );
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const abstract_vector<ELEMENTTYPE>& vec ) {
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
                    return set_column( col,
                                       mat.get_column( matcol )->as_std() );
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
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    size_t nvals, const ELEMENTTYPE *vals, int offset = 0 ) {
                    check_pointer();
                    _ptr->set_diagonal( nvals, vals, offset );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const std::vector<ELEMENTTYPE>& vals, int offset = 0 ) {
                    check_pointer();
                    _ptr->set_diagonal( vals, offset );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const std::initializer_list<ELEMENTTYPE> vals,
                    int offset = 0 ) {
                    check_pointer();
                    _ptr->set_diagonal( vals, offset );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const abstract_vector<ELEMENTTYPE>& vec, int offset = 0 ) {
                    return set_diagonal( vec.as_std(), offset );
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const Vector<ELEMENTTYPE>& vec, int offset = 0 ) {
                    return set_diagonal( vec.as_std(), offset );
                }

                virtual Matrix<ELEMENTTYPE>& transpose() {
                    check_pointer();
                    _ptr->transpose();
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& swap_rows( size_t ind1,
                                                        size_t ind2 ) {
                    check_pointer();
                    _ptr->swap_rows( ind1, ind2 );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& swap_columns( size_t ind1,
                                                           size_t ind2 ) {
                    check_pointer();
                    _ptr->swap_columns( ind1, ind2 );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& add(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check_pointer();
                    other.check_pointer();
                    check_same_size( other );
                    *_ptr += *( other._ptr );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& add( ELEMENTTYPE other ) {
                    check_pointer();
                    *_ptr += other;
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& subtract(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check_pointer();
                    other.check_pointer();
                    check_same_size( other );
                    *_ptr -= *( other._ptr );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& subtract( ELEMENTTYPE other ) {
                    check_pointer();
                    *_ptr -= other;
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& scale(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check_pointer();
                    other.check_pointer();
                    check_same_size( other );
                    _ptr->scale( *other._ptr );
                    // *_ptr *= *( other._ptr );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& scale( ELEMENTTYPE other ) {
                    check_pointer();
                    *_ptr *= other;
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& lu_decompose() {
                    check_pointer();
                    _lu.reset();
                    _lu = std::unique_ptr<LUDecomposition<ELEMENTTYPE>>(
                        new LUDecomposition<ELEMENTTYPE>() );
                    _lu->decompose( *this );
                    return *this;
                }

                virtual LUDecomposition<ELEMENTTYPE>& lu() {
                    check_pointer();
                    if ( !_lu ) {
                        this->lu_decompose();
                    }
                    return *_lu;
                }

                virtual Matrix<ELEMENTTYPE>& invert() {
                    Solver<ELEMENTTYPE> solver
                        = SolverFactory<ELEMENTTYPE>::build( solver_t::BASIC );
                    solver.set_system_matrix( *this );
                    Matrix<ELEMENTTYPE> inv
                        = MatrixFactory<ELEMENTTYPE>::build( matrix_t::DENSE );
                    inv.resize( this->rows(), this->columns() );
                    Vector<ELEMENTTYPE> vec
                        = VectorFactory<ELEMENTTYPE>::build(
                            vector_t::SPARSE );
                    vec.resize( this->columns() );

                    for ( size_t i = 0; i < this->columns(); i++ ) {
                        vec.zero().set( i, NCPA::math::one<ELEMENTTYPE>() );
                        inv.set_column( i, solver.solve( vec ) );
                    }
                    swap( *this, inv );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE> inverse() const {
                    Matrix<ELEMENTTYPE> inv = *this;
                    return inv.invert();
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
                    this->_clear_cache();
                    return *this;
                }

                virtual Vector<ELEMENTTYPE> right_multiply(
                    const Vector<ELEMENTTYPE>& other ) const {
                    check_pointer();
                    other.check_pointer();
                    if ( columns() != other.size() ) {
                        throw std::invalid_argument(
                            "Matrix-vector size mismatch: cannot "
                            "multiply" );
                    }
                    // _ptr->right_multiply( *(other._ptr) );
                    return Vector<ELEMENTTYPE>(
                        _ptr->right_multiply( *( other._ptr ) ) );
                }

                virtual Vector<ELEMENTTYPE> left_multiply(
                    const Vector<ELEMENTTYPE>& other ) const {
                    check_pointer();
                    other.check_pointer();
                    if ( columns() != other.size() ) {
                        throw std::invalid_argument(
                            "Matrix-vector size mismatch: cannot "
                            "multiply" );
                    }
                    return Vector<ELEMENTTYPE>(
                        _ptr->left_multiply( *( other._ptr ) ) );
                }

                virtual bool equals( const Matrix<ELEMENTTYPE>& other ) const {
                    return ( is_empty() || other.is_empty()
                                 ? false
                                 : _ptr->equals( *( other._ptr ) ) );
                }

                // assignment operators
                virtual Matrix<ELEMENTTYPE>& operator+=(
                    const Matrix<ELEMENTTYPE>& other ) {
                    return this->add( other );
                }

                virtual Matrix<ELEMENTTYPE>& operator+=(
                    const ELEMENTTYPE& other ) {
                    return this->add( other );
                }

                virtual Matrix<ELEMENTTYPE>& operator-=(
                    const Matrix<ELEMENTTYPE>& other ) {
                    return this->subtract( other );
                }

                virtual Matrix<ELEMENTTYPE>& operator-=(
                    const ELEMENTTYPE& other ) {
                    return this->add( -other );
                }

                virtual Matrix<ELEMENTTYPE>& operator*=(
                    const Matrix<ELEMENTTYPE>& other ) {
                    return this->multiply( other );
                }

                virtual Matrix<ELEMENTTYPE>& operator*=(
                    const ELEMENTTYPE& other ) {
                    return this->scale( other );
                }

                virtual Matrix<ELEMENTTYPE> operator-() const {
                    Matrix<ELEMENTTYPE> m = *this;
                    m.scale( -NCPA::math::one<ELEMENTTYPE>() );
                    return m;
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
                            os << mat.get( r, c );
                        }
                        if ( r != mat.rows() - 1 ) {
                            os << " ] ;" << std::endl;
                        } else {
                            os << " ] ]";
                        }
                    }
                    return os;
                }

                template<typename U = ELEMENTTYPE,
                         ENABLE_FUNCTION_IF_COMPLEX( U )>
                void print_nonzero( std::ostream& os,
                                    const std::string& sep = " " ) {
                    ELEMENTTYPE element;
                    os << this->rows() << std::endl;
                    for ( size_t r = 0; r < this->rows(); r++ ) {
                        auto nzinds = this->get_row( r )
                                          ->internal()
                                          ->nonzero_indices();
                        for ( auto cit = nzinds.cbegin(); cit != nzinds.cend();
                              ++cit ) {
                            element = this->get( r, *cit );
                            os << r << sep << *cit << sep << element.real()
                               << sep << element.imag() << std::endl;
                        }
                    }
                }

                template<typename U = ELEMENTTYPE,
                         ENABLE_FUNCTION_IF_REAL( U )>
                void print_nonzero( std::ostream& os,
                                    const std::string& sep = " " ) {
                    os << this->rows() << std::endl;
                    for ( size_t r = 0; r < this->rows(); r++ ) {
                        auto nzinds = this->get_row( r )
                                          ->internal()
                                          ->nonzero_indices();
                        for ( auto cit = nzinds.cbegin(); cit != nzinds.cend();
                              ++cit ) {
                            os << r << sep << *cit << sep
                               << this->get( r, *cit ) << std::endl;
                        }
                    }
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

                friend Vector<ELEMENTTYPE> operator*(
                    const Vector<ELEMENTTYPE>& vec,
                    const Matrix<ELEMENTTYPE>& mat ) {
                    return mat.left_multiply( vec );
                }

                friend Vector<ELEMENTTYPE> operator*(
                    const Matrix<ELEMENTTYPE>& mat,
                    const Vector<ELEMENTTYPE>& vec ) {
                    return mat.right_multiply( vec );
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

                const abstract_matrix<ELEMENTTYPE>& internal() const {
                    return *_ptr;
                }

            private:
                std::unique_ptr<abstract_matrix<ELEMENTTYPE>> _ptr;
                // std::map<size_t, WrapperVector<ELEMENTTYPE>> _wrappers;
                const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();

                std::unique_ptr<LUDecomposition<ELEMENTTYPE>> _lu;

                void _clear_cache() { _lu.reset( nullptr ); }
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::Matrix<T>& a,
                  NCPA::linear::Matrix<T>& b ) noexcept {
    // using std::swap;
    a._ptr.swap( b._ptr );
    // swap( a._wrappers, b._wrappers );
    a._lu.swap( b._lu );
}

// template<typename ELEMENTTYPE>
// NCPA::linear::Vector<ELEMENTTYPE> operator*(
//     const NCPA::linear::Vector<ELEMENTTYPE>& vec,
//     const NCPA::linear::Matrix<ELEMENTTYPE>& mat ) {
//     return mat.left_multiply( vec );
// }

// template<typename ELEMENTTYPE>
// NCPA::linear::Vector<ELEMENTTYPE> operator*(
//     const NCPA::linear::Matrix<ELEMENTTYPE>& mat,
//     const NCPA::linear::Vector<ELEMENTTYPE>& vec ) {
//     return mat.right_multiply( vec );
// }
