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
#include <cstring>
#include <initializer_list>
#include <map>
#include <memory>
#include <sstream>
#include <vector>


NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::sparse_matrix_compressed, ELEMENTTYPE );

namespace NCPA {
    namespace linear {
        namespace details {

            enum class compression_axis_t { ROWS, COLUMNS };

            NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
                class sparse_matrix_compressed<
                    ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
                : public abstract_matrix<ELEMENTTYPE> {
                public:
                    using abstract_matrix<ELEMENTTYPE>::set_row;
                    using abstract_matrix<ELEMENTTYPE>::set_column;

                    sparse_matrix_compressed( size_t nrows, size_t ncols,
                                              compression_axis_t axis
                                              = compression_axis_t::ROWS ) :
                         _axis { axis } {
                            this->zero( nrows, ncols );
                        }

                    sparse_matrix_compressed( compression_axis_t axis
                                              = compression_axis_t::ROWS ) :
                        sparse_matrix_compressed<ELEMENTTYPE>( 0, 0, axis ) {}

                    sparse_matrix_compressed(
                        const sparse_matrix_compressed<ELEMENTTYPE>& other ) :
                        sparse_matrix_compressed<ELEMENTTYPE>(
                            other.rows(), other.columns(), other.axis() ) {}

                    sparse_matrix_compressed(
                        sparse_matrix_compressed<ELEMENTTYPE>&&
                            source ) noexcept :
                        sparse_matrix_compressed<ELEMENTTYPE>() {
                        ::swap( *this, source );
                    }

                    // constructors, destructors, copying, and assignment
                    virtual ~sparse_matrix_compressed() {}

                    friend void ::swap<ELEMENTTYPE>(
                        sparse_matrix_compressed<ELEMENTTYPE>& a,
                        sparse_matrix_compressed<ELEMENTTYPE>& b ) noexcept;

                    sparse_matrix_compressed<ELEMENTTYPE>& operator=(
                        sparse_matrix_compressed<ELEMENTTYPE> other ) {
                        swap( *this, other );
                        return *this;
                    }

                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        clone() const override {
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new sparse_matrix_compressed( *this ) );
                    }

                    virtual std::unique_ptr<abstract_matrix<ELEMENTTYPE>>
                        fresh_clone() const override {
                        return std::unique_ptr<abstract_matrix<ELEMENTTYPE>>(
                            new sparse_matrix_compressed() );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& as_base_class() override {
                        return *static_cast<abstract_matrix<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual std::string id() const override {
                        return "NCPA compressed-storage sparse matrix";
                    }

                    virtual size_t rows() const override {
                        return ( this->axis() == compression_axis_t::ROWS
                                     ? _n_pointer_axis
                                     : _n_indexed_axis );
                    }

                    virtual size_t columns() const override {
                        return ( this->axis() == compression_axis_t::COLUMNS
                                     ? _n_pointer_axis
                                     : _n_indexed_axis );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& clear() override {
                        _n_pointer_axis = 0;
                        _n_indexed_axis = 0;
                        _indexed_array.clear();
                        _pointer_array.clear();
                        _value_array.clear();
                        return this->as_base_class();
                    }

                    virtual compression_axis_t axis() const { return _axis; }

                    virtual abstract_matrix<ELEMENTTYPE>& _set_indexed_value(
                        size_t pointer_axis_index, size_t indexed_axis_index,
                        ELEMENTTYPE val ) {
                        // for row_compressed, this is the row_index.
                        // _pointer_array[ n ] contains first index of the nth
                        // row values
                        size_t i_start = _pointer_array[ pointer_axis_index ];
                        size_t i_end
                            = _pointer_array[ pointer_axis_index + 1 ];
                        for ( size_t i = i_start; i < i_end; i++ ) {
                            if ( _indexed_array[ i ] > indexed_axis_index ) {
                                return _insert_indexed_value(
                                    i, indexed_axis_index, val );
                            } else if ( _indexed_array[ i ]
                                        == indexed_axis_index ) {
                                _value_array[ i ] = val;
                                return this->as_base_class();
                            }
                        }
                        return _insert_indexed_value(
                            i_end, indexed_axis_index, val );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& _insert_indexed_value(
                        size_t offset, size_t new_index, ELEMENTTYPE val ) {
                        _indexed_array.insert( _indexed_array.begin() + offset,
                                               new_index );
                        _value_array.insert( _value_array.begin() + offset,
                                             val );
                        for ( size_t j = offset + 1; j < _pointer_array.size();
                              j++ ) {
                            _pointer_array[ j ]++;
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        size_t row, size_t col, ELEMENTTYPE val ) override {
                        if ( this->axis() == compression_axis_t::ROWS ) {
                            return this->_set_indexed_value( row, col, val );
                        } else {
                            return this->_set_indexed_value( col, row, val );
                        }
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set(
                        ELEMENTTYPE val ) override {
                        this->zero();
                        for ( size_t i = 0; i < rows(); i++ ) {
                            for ( size_t j = 0; j < columns(); j++ ) {
                                this->set( i, j, val );
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_row(
                        size_t row, size_t nvals, const size_t *column_inds,
                        const ELEMENTTYPE *vals ) override {
                        for ( size_t i = 0; i < nvals; i++ ) {
                            this->set( row, column_inds[ i ], vals[ i ] );
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& set_column(
                        size_t column, size_t nvals, const size_t *row_inds,
                        const ELEMENTTYPE *vals ) override {
                        for ( size_t i = 0; i < nvals; i++ ) {
                            this->set( row_inds[ i ], column, vals[ i ] );
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& as_array(
                        size_t& nrows, size_t& ncols, ELEMENTTYPE **& vals ) override {
                        nrows = this->rows();
                        ncols = this->columns();
                        NCPA::arrays::fill( vals, nrows, ncols, _zero );
                        for ( size_t i = 0; i < _pointer_array.size(); i++ ) {
                            for ( size_t j = _pointer_array[ i ];
                                  j < _pointer_array[ j + 1 ]; j++ ) {
                                vals[ i ][ _indexed_array[ j ] ]
                                    = _value_array[ j ];
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        _get_pointer_axis_values( size_t index ) const {
                        std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                            new sparse_vector<ELEMENTTYPE>(
                                _n_indexed_axis ) );
                        for ( size_t j = _pointer_array[ index ];
                              j < _pointer_array[ index + 1 ]; j++ ) {
                            v->set( _indexed_array[ j ], _value_array[ j ] );
                        }
                        return v;
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        _get_indexed_axis_values( size_t index ) const {
                        std::unique_ptr<abstract_vector<ELEMENTTYPE>> v(
                            new sparse_vector<ELEMENTTYPE>(
                                _n_pointer_axis ) );
                        for ( size_t i = 0; i < _n_pointer_axis; i++ ) {
                            for ( size_t j = _pointer_array[ i ];
                                  j < _pointer_array[ i + 1 ]; j++ ) {
                                if ( _indexed_array[ j ] == index ) {
                                    v->set( i, _value_array[ j ] );
                                }
                            }
                        }
                        return v;
                    }

                    // @todo make read-write vector view for columns
                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_row( size_t row ) const override {
                        return ( this->axis() == compression_axis_t::ROWS
                                     ? this->_get_pointer_axis_values( row )
                                     : this->_get_indexed_axis_values( row ) );
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        get_column( size_t column ) const override {
                        return (
                            this->axis() == compression_axis_t::COLUMNS
                                ? this->_get_pointer_axis_values( column )
                                : this->_get_indexed_axis_values( column ) );
                    }

                    virtual const ELEMENTTYPE& _get_value(
                        size_t pointer_index, size_t indexed_index ) const {
                        for ( size_t j = _pointer_array[ pointer_index ];
                              j < _pointer_array[ pointer_index + 1 ]; j++ ) {
                            if ( _indexed_array[ j ] == indexed_index ) {
                                return _value_array[ j ];
                            }
                        }
                        return _zero;
                    }

                    virtual const ELEMENTTYPE& get(
                        size_t row, size_t col ) const override {
                        return ( this->axis() == compression_axis_t::ROWS
                                     ? this->_get_value( row, col )
                                     : this->_get_value( col, row ) );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& resize(
                        size_t rows, size_t cols ) override {
                        sparse_matrix_compressed<ELEMENTTYPE> newmat(
                            rows, cols, this->axis() );
                        size_t counter = 0;
                        for ( size_t i = 0;
                              i < std::min( rows, _n_pointer_axis ); i++ ) {
                            for ( size_t j = _pointer_array[ i ];
                                  j
                                  < std::min( cols, _pointer_array[ i + 1 ] );
                                  j++ ) {
                                newmat._indexed_array[ counter ]
                                    = this->_indexed_array[ j ];
                                newmat._value_array[ counter ]
                                    = this->_value_array[ j ];
                                newmat._pointer_array[ i + 1 ] = ++counter;
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& transpose()
                        override {
                        _axis = ( this->axis() == compression_axis_t::ROWS
                                      ? compression_axis_t::COLUMNS
                                      : compression_axis_t::ROWS );
                        return this->as_base_class();
                    }

                    abstract_matrix<ELEMENTTYPE>& _swap_pointer_indices(
                        size_t ind1, size_t ind2 ) {
                        if ( ind1 == ind2 ) {
                            return this->as_base_class();
                        }
                        size_t lower = std::min( ind1, ind2 );
                        size_t upper = std::max( ind1, ind2 );
                        size_t j;
                        std::vector<ELEMENTTYPE> valbuffer_lower,
                            valbuffer_upper;
                        std::vector<size_t> indbuffer_lower,
                            indbuffer_upper;
                        int size_diff  = (int)upper - (int)lower;
                        size_t n_upper = _pointer_array[ upper + 1 ]
                                       - _pointer_array[ upper ];
                        for ( j = _pointer_array[ upper ];
                              j < _pointer_array[ upper + 1 ]; j++ ) {
                            valbuffer_upper.push_back( _value_array[ j ] );
                            indbuffer_upper.push_back( _indexed_array[ j ] );
                        }
                        for ( j = _pointer_array[ lower ];
                              j < _pointer_array[ lower + 1 ]; j++ ) {
                            valbuffer_lower.push_back( _value_array[ j ] );
                            indbuffer_lower.push_back( _indexed_array[ j ] );
                        }
                        // replace higher-indexed values
                        _indexed_array.erase(
                            _indexed_array.begin() + _pointer_array[ upper ],
                            _indexed_array.begin()
                                + _pointer_array[ upper + 1 ] );
                        _indexed_array.insert(
                            _indexed_array.cbegin() + _pointer_array[ upper ],
                            indbuffer_lower.begin(), indbuffer_lower.end() );
                        _value_array.erase(
                            _value_array.begin() + _pointer_array[ upper ],
                            _value_array.begin()
                                + _pointer_array[ upper + 1 ] );
                        _value_array.insert(
                            _value_array.cbegin() + _pointer_array[ upper ],
                            valbuffer_lower.begin(), valbuffer_lower.end() );

                        // replace lower-indexed values
                        _indexed_array.erase(
                            _indexed_array.begin() + _pointer_array[ lower ],
                            _indexed_array.begin()
                                + _pointer_array[ lower + 1 ] );
                        _indexed_array.insert(
                            _indexed_array.cbegin() + _pointer_array[ lower ],
                            indbuffer_upper.begin(), indbuffer_upper.end() );
                        _value_array.erase(
                            _value_array.begin() + _pointer_array[ lower ],
                            _value_array.begin()
                                + _pointer_array[ lower + 1 ] );
                        _value_array.insert(
                            _value_array.cbegin() + _pointer_array[ lower ],
                            valbuffer_upper.begin(), valbuffer_upper.end() );

                        // finally adjust intermediate pointers
                        for ( j = lower + 1; j <= upper; j++ ) {
                            _pointer_array[ j ]
                                = (size_t)( (int)_pointer_array[ j ]
                                            + size_diff );
                        }
                        return this->as_base_class();
                    }

                    abstract_matrix<ELEMENTTYPE>& _swap_indexed_indices(
                        size_t ind1, size_t ind2 ) {
                        if ( ind1 == ind2 ) {
                            return this->as_base_class();
                        }
                        size_t lower = std::min( ind1, ind2 );
                        size_t upper = std::max( ind1, ind2 );

                        for ( size_t i = 0; i < _n_pointer_axis; i++ ) {
                            bool havelower = false, haveupper = false;
                            size_t lower_ind = 0, upper_ind = 0;
                            ELEMENTTYPE lower_val, upper_val;
                            for ( size_t j = _pointer_array[ i ];
                                  j < _pointer_array[ i + 1 ]; j++ ) {
                                if ( _indexed_array[ j ] < lower ) {
                                    lower_ind = j;
                                    upper_ind = j;
                                } else if ( _indexed_array[ j ] == lower ) {
                                    lower_ind = j;
                                    upper_ind = j;
                                    havelower = true;
                                    lower_val = _value_array[ j ];
                                } else if ( _indexed_array[ j ] < upper ) {
                                    upper_ind = j;
                                } else if ( _indexed_array[ j ] == upper ) {
                                    upper_ind = j;
                                    haveupper = true;
                                    upper_val = _value_array[ j ];
                                }
                            }
                            if ( havelower && haveupper ) {
                                std::swap( _value_array[ lower_ind ],
                                           _value_array[ upper_ind ] );
                            } else if ( havelower ) {
                                _indexed_array.insert( _indexed_array.begin()
                                                           + upper_ind + 1,
                                                       upper );
                                _indexed_array.erase( _indexed_array.begin()
                                                      + lower_ind );
                                _value_array.insert( _value_array.begin()
                                                         + upper_ind + 1,
                                                     upper_val );
                                _value_array.erase( _value_array.begin()
                                                    + lower_ind );
                            } else if ( haveupper ) {
                                _indexed_array.erase( _indexed_array.begin()
                                                      + upper_ind );
                                _indexed_array.insert( _indexed_array.begin()
                                                           + lower_ind + 1,
                                                       lower );
                                _value_array.erase( _value_array.begin()
                                                    + upper_ind );
                                _value_array.insert( _value_array.begin()
                                                         + lower_ind + 1,
                                                     lower );
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& swap_rows(
                        size_t ind1, size_t ind2 ) override {
                        return (
                            this->axis() == compression_axis_t::ROWS
                                ? this->_swap_pointer_indices( ind1, ind2 )
                                : this->_swap_indexed_indices( ind1, ind2 ) );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& swap_columns(
                        size_t ind1, size_t ind2 ) override {
                        return (
                            this->axis() == compression_axis_t::COLUMNS
                                ? this->_swap_pointer_indices( ind1, ind2 )
                                : this->_swap_indexed_indices( ind1, ind2 ) );
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& change_axis(
                        compression_axis_t new_axis ) {
                        if ( this->axis() == new_axis ) {
                            return this->as_base_class();
                        }
                        sparse_matrix_compressed<ELEMENTTYPE> newmat(
                            this->rows(), this->columns(), new_axis );
                        for ( size_t i = 0; i < _pointer_array.size() - 1;
                              i++ ) {
                            for ( size_t j = _pointer_array[ i ];
                                  j < _pointer_array[ i + 1 ]; j++ ) {
                                if ( this->axis()
                                     == compression_axis_t::ROWS ) {
                                    newmat.set( i, _indexed_array[ j ],
                                                _value_array[ j ] );
                                } else {
                                    newmat.set( _indexed_array[ j ], i,
                                                _value_array[ j ] );
                                }
                            }
                        }
                        std::swap( *this, newmat );
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& zero(
                        size_t row, size_t col ) override {
                        return ( this->axis() == compression_axis_t::ROWS
                                     ? this->_remove_element( row, col )
                                     : this->_remove_element( col, row ) );
                    }

                    abstract_matrix<ELEMENTTYPE>& _remove_element(
                        size_t pointer_index, size_t indexed_index ) {
                        for ( size_t j = _pointer_array[ pointer_index ];
                              j < _pointer_array[ pointer_index + 1 ]; j++ ) {
                            if ( _indexed_array[ j ] == indexed_index ) {
                                _indexed_array.erase( _indexed_array.begin()
                                                      + j );
                                _value_array.erase( _value_array.begin() + j );
                                for ( size_t i = pointer_index + 1;
                                      i < _pointer_array.size(); i++ ) {
                                    _pointer_array[ i ]--;
                                }
                            } else if ( _indexed_array[ j ] > indexed_index ) {
                                return this->as_base_class();
                            }
                        }
                        return this->as_base_class();
                    }

                    virtual abstract_matrix<ELEMENTTYPE>& zero() override {
                        if ( this->axis() == compression_axis_t::ROWS ) {
                            this->_set_dimensions( this->rows(),
                                                   this->columns() );
                        } else {
                            this->_set_dimensions( this->columns(),
                                                   this->rows() );
                        }
                        return this->as_base_class();
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        build_vector( size_t n = 0 ) const override {
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new sparse_vector<ELEMENTTYPE>( n ) );
                    }

                protected:
                    compression_axis_t _axis;
                    size_t _n_pointer_axis,  // == _nrows if row-compressed
                        _n_indexed_axis;     // == _ncols if row-compressed
                    

                    std::vector<size_t> _pointer_array, _indexed_array;
                    std::vector<ELEMENTTYPE> _value_array;
                    const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();

                    void _set_dimensions( size_t _n_ptr, size_t _n_ind ) {
                        _n_pointer_axis = _n_ptr;
                        _n_indexed_axis = _n_ind;
                        _pointer_array.clear();
                        _pointer_array.resize( _n_ptr + 1 );
                        _indexed_array.clear();
                        _value_array.clear();
                    }
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap(
    NCPA::linear::details::sparse_matrix_compressed<T>& a,
    NCPA::linear::details::sparse_matrix_compressed<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::details::abstract_matrix<T>&>( a ),
            static_cast<NCPA::linear::details::abstract_matrix<T>&>( b ) );
    swap( a._n_pointer_axis, b._n_pointer_axis );
    swap( a._n_indexed_axis, b._n_indexed_axis );
    swap( a._axis, b._axis );
    swap( a._indexed_array, b._indexed_array );
    swap( a._pointer_array, b._pointer_array );
    swap( a._value_array, b._value_array );
}
