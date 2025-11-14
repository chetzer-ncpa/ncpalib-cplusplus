#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/functions.hpp"
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
                    if (other._ptr) {
                        _ptr = std::move( other._ptr->clone() );
                    }
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

                virtual std::unique_ptr<Matrix<ELEMENTTYPE>> clone() const {
                    return std::unique_ptr<Matrix<ELEMENTTYPE>>(
                        new Matrix<ELEMENTTYPE>( *this ) );
                }

                virtual Matrix<ELEMENTTYPE>& add(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check();
                    other.check();
                    check_same_size( other );
                    *_ptr += *( other._ptr );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& add( ELEMENTTYPE other ) {
                    check();
                    *_ptr += other;
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& as_array( size_t& nrows,
                                                       size_t& ncols,
                                                       ELEMENTTYPE **& vals ) {
                    check();
                    _ptr->as_array( nrows, ncols, vals );
                    return *this;
                }

                virtual explicit operator bool() const {
                    return ( _ptr ? true : false );
                }

                virtual size_t bandwidth() const {
                    return ( _ptr ? _ptr->bandwidth() : 1 );
                }

                virtual Matrix<ELEMENTTYPE>& clean( ELEMENTTYPE tol ) {
                    if (_ptr) {
                        _ptr->clean( tol );
                    }
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& clear() {
                    if (_ptr) {
                        _ptr->clear();
                    }
                    this->_clear_cache();
                    return *this;
                }

                virtual size_t columns() const {
                    return ( _ptr ? _ptr->columns() : 0 );
                }

                virtual Matrix<ELEMENTTYPE>& copy(
                    const Matrix<ELEMENTTYPE>& other ) {
                    if (_ptr) {
                        _ptr->copy( *other._ptr );
                    } else {
                        _ptr = std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            other._ptr->clone() );
                    }
                    if (other._lu) {
                        _lu = std::unique_ptr<LUDecomposition<ELEMENTTYPE>>(
                            new LUDecomposition<ELEMENTTYPE>( *other._lu ) );
                    }
                    return *this;
                }

                virtual ELEMENTTYPE determinant() const {
                    if (this->is_empty()) {
                        return NCPA::math::one<ELEMENTTYPE>();
                    }
                    if (this->is_zero()) {
                        return 0;
                    }
                    if (this->is_identity()) {
                        return NCPA::math::one<ELEMENTTYPE>();
                    }
                    if (this->is_diagonal() || this->is_upper_triangular() || this->is_lower_triangular()) {
                        ELEMENTTYPE prod = NCPA::math::one<ELEMENTTYPE>();
                        for (size_t i = 0; i < this->diagonal_size(0); ++i) {
                            prod *= this->get(i,i);
                        }
                        return prod;
                    }
                    std::unique_ptr<LUDecomposition<ELEMENTTYPE>> lu( 
                        (this->is_band_diagonal() 
                        ? new BandDiagonalLUDecomposition<ELEMENTTYPE>()
                        : new LUDecomposition<ELEMENTTYPE>())
                    );
                    lu->decompose( *this );
                    return lu->lower().determinant() * lu->upper().determinant();
                }

                virtual size_t diagonal_size( int offset = 0 ) const {
                    return matrix_diagonal_size( this->rows(), this->columns(),
                                                 offset );
                }

                virtual std::vector<int> diagonals() const {
                    return ( _ptr ? _ptr->diagonals() : std::vector<int>() );
                }

                virtual bool equals( const Matrix<ELEMENTTYPE>& other ) const {
                    return ( is_empty() || other.is_empty()
                                 ? false
                                 : _ptr->equals( *( other._ptr ) ) );
                }

                virtual const ELEMENTTYPE& get( size_t row,
                                                size_t col ) const {
                    return ( _ptr ? _ptr->get( row, col ) : _zero );
                }

                virtual const ELEMENTTYPE& get( size_t ind ) const {
                    if (*this) {
                        if (is_row_matrix()) {
                            return get( 0, ind );
                        } else if (is_column_matrix()) {
                            return get( ind, 0 );
                        } else {
                            throw std::range_error(
                                "Both dimensions must be specified for "
                                "non-row and non-column matrices" );
                        }
                    } else {
                        return _zero;
                    }
                }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_column(
                    size_t col ) const {
                    if (_ptr) {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>(
                            new Vector<ELEMENTTYPE>(
                                _ptr->get_column( col ) ) );
                    } else {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>();
                    }
                }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_diagonal(
                    int offset = 0 ) const {
                    if (_ptr) {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>(
                            new Vector<ELEMENTTYPE>(
                                _ptr->get_diagonal( offset ) ) );
                    } else {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>();
                    }
                }

                virtual std::unique_ptr<Vector<ELEMENTTYPE>> get_row(
                    size_t row ) const {
                    if (_ptr) {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>(
                            new Vector<ELEMENTTYPE>( _ptr->get_row( row ) ) );
                    } else {
                        return std::unique_ptr<Vector<ELEMENTTYPE>>();
                    }
                }

                virtual Matrix<ELEMENTTYPE>& identity() {
                    check();
                    this->_clear_cache();
                    _ptr->identity();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& identity( size_t rows,
                                                       size_t cols = 0 ) {
                    check();
                    if (cols == 0) {
                        cols = rows;
                    }
                    _ptr->identity( rows, cols );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE> inverse() const {
                    Matrix<ELEMENTTYPE> inv = *this;
                    return inv.invert();
                }

                virtual Matrix<ELEMENTTYPE>& invert() {
                    if (this->is_empty()) {
                        return *this;
                    }
                    if (this->is_zero()) {
                        throw std::runtime_error( "Zero matrices are singular "
                                                  "and cannot be inverted" );
                    }
                    if (this->is_diagonal()) {
                        return this->invert_diagonal();
                    }
                    if (this->is_tridiagonal()) {
                        return this->invert_tridiagonal();
                    }
                    return this->invert_general();
                }

                virtual Matrix<ELEMENTTYPE>& invert_diagonal() {
                    if (!this->is_square()) {
                        throw std::range_error(
                            "Cannot invert a non-square matrix" );
                    }
                    size_t ds = this->diagonal_size( 0 );
                    for (size_t i = 0; i < ds; ++i) {
                        if (this->is_zero( i, i )) {
                            std::ostringstream oss;
                            oss << "Element " << i
                                << " of input matrix is zero, matrix is "
                                   "singular";
                            throw std::range_error( oss.str() );
                        }
                        this->set( i, i,
                                   NCPA::math::inverse( this->get( i, i ) ) );
                    }
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& invert_general() {
                    if (!this->is_square()) {
                        throw std::logic_error(
                            "Cannot invert a non-square matrix" );
                    }
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

                    for (size_t i = 0; i < this->columns(); i++) {
                        vec.zero().set( i, NCPA::math::one<ELEMENTTYPE>() );
                        inv.set_column( i, solver.solve( vec ) );
                    }
                    swap( *this, inv );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& invert_tridiagonal() {
                    if (!this->is_square()) {
                        throw std::logic_error(
                            "Cannot invert a non-square matrix" );
                    }
                    Matrix T
                        = MatrixFactory<ELEMENTTYPE>::build( matrix_t::DENSE );
                    T.resize( this->rows(), this->columns() );
                    size_t n = T.diagonal_size( 0 );
                    std::vector<ELEMENTTYPE> theta( n + 1 ), phi( n + 2 );
                    const ELEMENTTYPE one = NCPA::math::one<ELEMENTTYPE>();

                    // ICs
                    theta[ 0 ]   = one;
                    theta[ 1 ]   = this->get( 0, 0 );
                    phi[ n ]     = this->get( n - 1, n - 1 );
                    phi[ n + 1 ] = one;

                    for (size_t i = 2; i <= n; ++i) {
                        theta[ i ] = this->get( i - 1, i - 1 ) * theta[ i - 1 ]
                                   - this->get( i - 2, i - 1 )
                                         * this->get( i - 1, i - 2 )
                                         * theta[ i - 2 ];
                        // NCPA_DEBUG << "theta[ " << i << " ] = " << theta[ i
                        // ] << std::endl;
                        size_t iprime = n - i + 1;
                        phi[ iprime ] = this->get( iprime - 1, iprime - 1 )
                                          * phi[ iprime + 1 ]
                                      - this->get( iprime - 1, iprime )
                                            * this->get( iprime, iprime - 1 )
                                            * phi[ iprime + 2 ];
                        // NCPA_DEBUG << "phi[ " << iprime << " ] = " << phi[
                        // iprime ] << std::endl;
                    }
                    int rows = (int)this->rows();
                    int cols = (int)this->columns();

                    // use i=1..nr and j=1..nc to match source notation
                    for (size_t i = 1; i <= this->rows(); ++i) {
                        size_t r = i - 1;
                        T.set( r, r,
                               theta[ i - 1 ] * phi[ i + 1 ] / theta[ n ] );
                        // NCPA_DEBUG << "T[ " << r << ", " << r << " ] = " <<
                        // T.get(r,r) << std::endl;

                        if (i < rows) {
                            // i < j case
                            ELEMENTTYPE bprod = one;
                            for (int j = i + 1; j <= cols; ++j) {
                                size_t c          = j - 1;
                                ELEMENTTYPE sign  = std::pow( -one, i + j );
                                bprod            *= this->get( j - 2, j - 1 );
                                // NCPA_DEBUG << "bprod *= " << this->get( j-2,
                                // j-1 ) << " = " << bprod << std::endl;
                                T.set( r, c,
                                       sign * bprod * theta[ i - 1 ]
                                           * phi[ j + 1 ] / theta[ n ] );
                                // NCPA_DEBUG << "T[ " << r << ", " << c << " ]
                                // = " << T.get(r,c) << std::endl;
                            }
                        }

                        if (i > 1) {
                            // i > j case
                            ELEMENTTYPE cprod = one;
                            for (int j = i - 1; j >= 1; --j) {
                                size_t c          = j - 1;
                                ELEMENTTYPE sign  = std::pow( -one, i + j );
                                cprod            *= this->get( j, j - 1 );
                                // NCPA_DEBUG << "cprod *= " << this->get( j,
                                // j-1 ) << " = " << cprod << std::endl;
                                T.set( r, c,
                                       sign * cprod * theta[ j - 1 ]
                                           * phi[ i + 1 ] / theta[ n ] );
                                // NCPA_DEBUG << "T[ " << r << ", " << c << " ]
                                // = " << T.get(r,c) << std::endl;
                            }
                        }
                    }
                    swap( *this, T );
                    return *this;
                }

                virtual bool is_band_diagonal() const {
                    return ( _ptr ? _ptr->is_band_diagonal() : true );
                }

                virtual bool is_block_matrix() const { return false; }

                virtual bool is_column_matrix() const {
                    return ( columns() == 1 );
                }

                virtual bool is_diagonal() const {
                    return ( _ptr ? _ptr->is_diagonal() : true );
                }

                virtual bool is_empty() const {
                    return ( _ptr ? _ptr->is_empty() : true );
                }

                virtual bool is_identity() const {
                    return ( _ptr ? _ptr->is_identity() : false );
                }

                virtual bool is_lower_triangular() const {
                    return ( _ptr ? _ptr->is_lower_triangular() : true );
                }

                virtual bool is_null() const {
                    return this->is_zero();
                }

                virtual bool is_row_matrix() const { return ( rows() == 1 ); }

                virtual bool is_square() const {
                    return ( _ptr ? _ptr->is_square() : true );
                }

                virtual bool is_symmetric() const {
                    return ( _ptr ? _ptr->is_symmetric() : true );
                }

                virtual bool is_tridiagonal() const {
                    return ( _ptr ? _ptr->is_tridiagonal() : true );
                }

                virtual bool is_upper_triangular() const {
                    return ( _ptr ? _ptr->is_upper_triangular() : true );
                }

                virtual bool is_zero() const {
                    return ( _ptr ? _ptr->is_zero() : true );
                }

                virtual bool is_zero( size_t r, size_t c ) const {
                    return ( _ptr ? _ptr->is_zero( r, c ) : true );
                }

                virtual Vector<ELEMENTTYPE> left_multiply(
                    const Vector<ELEMENTTYPE>& other ) const {
                    check();
                    // other.check();
                    if (columns() != other.size()) {
                        throw std::invalid_argument(
                            "Matrix-vector size mismatch: cannot "
                            "multiply" );
                    }
                    return Vector<ELEMENTTYPE>(
                        _ptr->left_multiply( *( other._ptr ) ) );
                }

                virtual size_t lower_bandwidth() const {
                    return ( _ptr ? _ptr->lower_bandwidth() : 0 );
                }

                virtual LUDecomposition<ELEMENTTYPE>& lu() {
                    check();
                    if (!_lu) {
                        this->lu_decompose();
                    }
                    return *_lu;
                }

                virtual Matrix<ELEMENTTYPE>& lu_decompose() {
                    check();
                    _lu.reset();
                    _lu = std::unique_ptr<LUDecomposition<ELEMENTTYPE>>(
                        new LUDecomposition<ELEMENTTYPE>() );
                    _lu->decompose( *this );
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& multiply(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check();
                    other.check();
                    if (columns() != other.rows()) {
                        throw std::invalid_argument(
                            "Matrix size mismatch: cannot multiply" );
                    }
                    auto newmat = Matrix<ELEMENTTYPE>(
                        _ptr->multiply( *( other._ptr ) ) );
                    swap( *this, newmat );
                    this->_clear_cache();
                    return *this;
                }

                virtual void print( std::ostream& os = std::cout ) const {
                    os << *this << std::endl;
                }

                virtual Matrix<ELEMENTTYPE>& resize( size_t rows,
                                                     size_t cols ) {
                    check();
                    if (rows != this->rows() || cols != this->columns()) {
                        this->_clear_cache();
                        _ptr->resize( rows, cols );
                    }
                    return *this;
                }

                virtual Vector<ELEMENTTYPE> right_multiply(
                    const Vector<ELEMENTTYPE>& other ) const {
                    check();
                    // other.check();
                    if (columns() != other.size()) {
                        throw std::invalid_argument(
                            "Matrix-vector size mismatch: cannot "
                            "multiply" );
                    }
                    return Vector<ELEMENTTYPE>(
                        _ptr->right_multiply( *( other._ptr ) ) );
                }

                virtual size_t rows() const {
                    return ( _ptr ? _ptr->rows() : 0 );
                }

                virtual Matrix<ELEMENTTYPE>& scale(
                    const Matrix<ELEMENTTYPE>& other ) {
                    check();
                    other.check();
                    check_same_size( other );
                    _ptr->scale( *other._ptr );
                    // *_ptr *= *( other._ptr );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& scale( ELEMENTTYPE other ) {
                    check();
                    *_ptr *= other;
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set( size_t row, size_t col,
                                                  ELEMENTTYPE val ) {
                    check();
                    _ptr->set( row, col, val );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set( ELEMENTTYPE val ) {
                    check();
                    _ptr->set( val );
                    this->_clear_cache();
                    return *this;
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
                    check();
                    _ptr->set_column( col, val );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, size_t nvals, const size_t *indices,
                    const ELEMENTTYPE *vals ) {
                    check();
                    _ptr->set_column( col, nvals, indices, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::vector<size_t>& colinds,
                    const std::vector<ELEMENTTYPE>& vals ) {
                    check();
                    _ptr->set_column( col, colinds, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::vector<ELEMENTTYPE>& vals ) {
                    check();
                    _ptr->set_column( col, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_column(
                    size_t col, const std::initializer_list<size_t> colinds,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    check();
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
                    if (!mat.is_column_matrix()) {
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
                    check();
                    _ptr->set_diagonal( val, offset );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    size_t nvals, const ELEMENTTYPE *vals, int offset = 0 ) {
                    check();
                    _ptr->set_diagonal( nvals, vals, offset );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const std::vector<ELEMENTTYPE>& vals, int offset = 0 ) {
                    return this->set_diagonal( vals.size(), vals.data(),
                                               offset );
                    // check();
                    // _ptr->set_diagonal( vals, offset );
                    // this->_clear_cache();
                    // return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const std::initializer_list<ELEMENTTYPE> vals,
                    int offset = 0 ) {
                    return this->set_diagonal(
                        std::vector<ELEMENTTYPE>( vals ), offset );
                    // check();
                    // _ptr->set_diagonal( vals, offset );
                    // this->_clear_cache();
                    // return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const abstract_vector<ELEMENTTYPE>& vec, int offset = 0 ) {
                    return set_diagonal( vec.as_std(), offset );
                }

                virtual Matrix<ELEMENTTYPE>& set_diagonal(
                    const Vector<ELEMENTTYPE>& vec, int offset = 0 ) {
                    return set_diagonal( vec.as_std(), offset );
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
                    check();
                    _ptr->set_row( row, val );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, size_t nvals, const size_t *indices,
                    const ELEMENTTYPE *vals ) {
                    check();
                    _ptr->set_row( row, nvals, indices, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::vector<size_t>& rowinds,
                    const std::vector<ELEMENTTYPE>& vals ) {
                    check();
                    _ptr->set_row( row, rowinds, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::vector<ELEMENTTYPE>& vals ) {
                    check();
                    _ptr->set_row( row, vals );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row, const std::initializer_list<size_t> rowinds,
                    const std::initializer_list<ELEMENTTYPE> vals ) {
                    check();
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
                    check();
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
                    if (!mat.is_row_matrix()) {
                        throw std::logic_error( "Matrix is not a row matrix" );
                    }
                    return set_row( row, mat, 0 );
                }

                virtual Matrix<ELEMENTTYPE>& set_row(
                    size_t row,
                    const std::unique_ptr<Matrix<ELEMENTTYPE>>& mat ) {
                    return set_row( row, *mat, 0 );
                }

                virtual Matrix<ELEMENTTYPE>& subtract(
                    const Matrix<ELEMENTTYPE>& other ) {
                    this->check();
                    other.check();
                    check_same_size( other );
                    *_ptr -= *( other._ptr );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& subtract( ELEMENTTYPE other ) {
                    check();
                    *_ptr -= other;
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& swap_columns( size_t ind1,
                                                           size_t ind2 ) {
                    check();
                    _ptr->swap_columns( ind1, ind2 );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& swap_rows( size_t ind1,
                                                        size_t ind2 ) {
                    check();
                    _ptr->swap_rows( ind1, ind2 );
                    this->_clear_cache();
                    return *this;
                }

                Matrix<ELEMENTTYPE>& transpose() {
                    check();
                    _ptr->transpose();
                    this->_clear_cache();
                    return *this;
                }

                virtual size_t upper_bandwidth() const {
                    return ( _ptr ? _ptr->upper_bandwidth() : 0 );
                }

                virtual Matrix<ELEMENTTYPE>& zero( size_t row, size_t col ) {
                    check();
                    _ptr->zero( row, col );
                    this->_clear_cache();
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE>& zero() {
                    check();
                    _ptr->zero();
                    this->_clear_cache();
                    return *this;
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

                // NOT virtual because of invalid covariant return type; would
                // have to be a smart pointer
                Matrix<ELEMENTTYPE> operator-() const {
                    Matrix<ELEMENTTYPE> m = *this;
                    m.scale( -NCPA::math::one<ELEMENTTYPE>() );
                    return m;
                }

                // friend operators
                friend std::ostream& operator<<(
                    std::ostream& os, const Matrix<ELEMENTTYPE>& mat ) {
                    os << "[";
                    for (size_t r = 0; r < mat.rows(); r++) {
                        if (r > 0) {
                            os << " ";
                        }
                        os << " [ ";
                        for (size_t c = 0; c < mat.columns(); c++) {
                            if (c > 0) {
                                os << ", ";
                            }
                            os << mat.get( r, c );
                        }
                        if (r != mat.rows() - 1) {
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
                                    const std::string& sep = " " ) const {
                    ELEMENTTYPE element;
                    os << this->rows() << std::endl;
                    for (size_t r = 0; r < this->rows(); r++) {
                        auto nzinds = this->get_row( r )
                                          ->internal()
                                          ->nonzero_indices();
                        for (auto cit = nzinds.cbegin(); cit != nzinds.cend();
                             ++cit) {
                            element = this->get( r, *cit );
                            os << r << sep << *cit << sep << element.real()
                               << sep << element.imag() << std::endl;
                        }
                    }
                }

                template<typename U = ELEMENTTYPE,
                         ENABLE_FUNCTION_IF_REAL( U )>
                void print_nonzero( std::ostream& os,
                                    const std::string& sep = " " ) const {
                    os << this->rows() << std::endl;
                    for (size_t r = 0; r < this->rows(); r++) {
                        auto nzinds = this->get_row( r )
                                          ->internal()
                                          ->nonzero_indices();
                        for (auto cit = nzinds.cbegin(); cit != nzinds.cend();
                             ++cit) {
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

                virtual void check() const { this->_check_pointer(); }

                virtual void check_same_size(
                    const Matrix<ELEMENTTYPE>& other ) const {
                    if (rows() != other.rows()
                        || columns() != other.columns()) {
                        throw std::invalid_argument(
                            "Matrices are not the same size!" );
                    }
                }

                const abstract_matrix<ELEMENTTYPE> *internal() const {
                    return ( _ptr ? _ptr.get() : nullptr );
                }


            private:
                std::unique_ptr<abstract_matrix<ELEMENTTYPE>> _ptr;
                // std::map<size_t, WrapperVector<ELEMENTTYPE>> _wrappers;
                const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();

                std::unique_ptr<LUDecomposition<ELEMENTTYPE>> _lu;

                matrix_t _type;

                void _clear_cache() { _lu.reset( nullptr ); }

                virtual void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "Matrix: Internal pointer has not been set!" );
                    }
                }
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
