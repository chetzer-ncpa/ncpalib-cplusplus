#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
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

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::sparse_vector,
                                             ELEMENTTYPE );

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class sparse_vector<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
            : public abstract_vector<ELEMENTTYPE> {
            public:
                sparse_vector() {}

                sparse_vector( const std::vector<ELEMENTTYPE>& v ) :
                    sparse_vector() {
                    set( v );
                }

                sparse_vector( size_t n ) : sparse_vector<ELEMENTTYPE>() {
                    resize( n );
                }

                sparse_vector( const std::initializer_list<ELEMENTTYPE>& v ) :
                    sparse_vector( std::vector<ELEMENTTYPE>( v ) ) {
                    // set( std::vector<ELEMENTTYPE>( v ) );
                }

                // copy constructor
                sparse_vector( const sparse_vector<ELEMENTTYPE>& other ) :
                    sparse_vector<ELEMENTTYPE>() {
                    _elements = other._elements;
                    _capacity = other._capacity;
                }

                /**
                 * Move constructor.
                 * @param source The vector to assimilate.
                 */
                sparse_vector( sparse_vector<ELEMENTTYPE>&& source ) noexcept :
                    sparse_vector<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                virtual ~sparse_vector() {}

                friend void ::swap<ELEMENTTYPE>(
                    sparse_vector<ELEMENTTYPE>& a,
                    sparse_vector<ELEMENTTYPE>& b ) noexcept;

                /**
                 * Assignment operator.
                 * @param other The vector to assign to this.
                 */
                sparse_vector<ELEMENTTYPE>& operator=(
                    sparse_vector<ELEMENTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual abstract_vector<ELEMENTTYPE>& set(
                    size_t n, const ELEMENTTYPE *val ) override {
                    clear();
                    for ( size_t i = 0; i < n; i++ ) {
                        if ( !NCPA::math::is_zero( val[ i ] ) ) {
                            // std::cout << "Setting _elements[" << i << "]
                            // = " << val[i] << std::endl;
                            _elements[ i ] = val[ i ];
                        }
                    }

                    _capacity = n;
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& set(
                    const std::vector<ELEMENTTYPE>& v ) override {
                    return set( v.size(), &v[ 0 ] );
                }

                virtual abstract_vector<ELEMENTTYPE>& set(
                    ELEMENTTYPE val ) override {
                    if ( NCPA::math::is_zero( val ) ) {
                        _elements.clear();
                    } else {
                        for ( size_t i = 0; i < _capacity; i++ ) {
                            set( i, val );
                        }
                    }
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual size_t size() const override { return _capacity; }

                virtual const ELEMENTTYPE& get( size_t n ) const override {
                    _check_out_of_range( n );
                    auto it = _elements.find( n );
                    if ( it != _elements.cend() ) {
                        return it->second;
                    } else {
                        return _zero;
                    }
                }

                virtual ELEMENTTYPE& get( size_t n ) override {
                    _check_out_of_range( n );
                    auto it = _elements.find( n );
                    if ( it == _elements.end() ) {
                        _elements[ n ] = _zero;
                    }
                    return _elements[ n ];
                }

                virtual std::vector<ELEMENTTYPE> as_std() const override {
                    std::vector<ELEMENTTYPE> v( _capacity );
                    for ( auto it = _elements.cbegin(); it != _elements.cend();
                          ++it ) {
                        v[ it->first ] = it->second;
                    }
                    return v;
                }

                virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>> clone()
                    override {
                    // return std::make_unique<Implementation>(*this);
                    return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                        new sparse_vector( *this ) );
                }

                virtual abstract_vector<ELEMENTTYPE>& clear() override {
                    _elements.clear();
                    _capacity = 0;
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& resize(
                    size_t n ) override {
                    if ( n < _capacity ) {
                        auto last = _elements.upper_bound( n - 1 );
                        while ( last != _elements.end() ) {
                            _elements.erase( last );
                            last = _elements.upper_bound( n - 1 );
                        }
                    }
                    _capacity = n;
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& as_array(
                    size_t& n, ELEMENTTYPE *& vals ) override {
                    if ( n == 0 ) {
                        n = size();
                        if ( vals != nullptr ) {
                            delete[] vals;
                            vals = nullptr;
                        }
                        vals = NCPA::arrays::zeros<ELEMENTTYPE>( n );
                    } else if ( n != size() ) {
                        throw std::invalid_argument(
                            "Size mismatch between vector and target "
                            "array" );
                    }
                    std::fill( vals, vals + n, _zero );
                    for ( auto it = _elements.cbegin(); it != _elements.cend();
                          ++it ) {
                        vals[ it->first ] = it->second;
                    }
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& set(
                    size_t n, ELEMENTTYPE val ) override {
                    _check_out_of_range( n );
                    if ( NCPA::math::is_zero( val ) ) {
                        _elements.erase( n );
                    } else {
                        _elements[ n ] = val;
                    }
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& scale(
                    ELEMENTTYPE val ) override {
                    for ( auto it = _elements.begin(); it != _elements.end();
                          ++it ) {
                        it->second *= val;
                    }
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& scale(
                    const abstract_vector<ELEMENTTYPE>& b ) override {
                    _check_same_size( b );
                    std::vector<size_t> keep = this->nonzero_index_union( b );

                    sparse_vector<ELEMENTTYPE> newv;
                    newv.resize( size() );
                    for ( auto uit = keep.cbegin(); uit != keep.cend();
                          ++uit ) {
                        newv.set( *uit, get( *uit ) * b.get( *uit ) );
                    }
                    swap( *this, newv );
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual ELEMENTTYPE dot(
                    const abstract_vector<ELEMENTTYPE>& b ) const override {
                    if ( auto *derived
                         = dynamic_cast<const sparse_vector<ELEMENTTYPE> *>(
                             &b ) ) {
                        return _sparse_dot(
                            dynamic_cast<const sparse_vector<ELEMENTTYPE>&>(
                                b ) );
                    } else {
                        _check_same_size( b );
                        ELEMENTTYPE product = _zero;
                        for ( auto it = _elements.cbegin();
                              it != _elements.cend(); ++it ) {
                            product += it->second * b.get( it->first );
                        }
                        return product;
                    }
                }

                ELEMENTTYPE _sparse_dot(
                    const sparse_vector<ELEMENTTYPE>& b ) const {
                    _check_same_size( b );
                    ELEMENTTYPE product = _zero;
                    auto a_it           = this->_elements.cbegin();
                    auto b_it           = b._elements.cbegin();
                    while ( a_it != this->_elements.cend()
                            && b_it != b._elements.cend() ) {
                        if ( a_it->first == b_it->first ) {
                            product += a_it->second * b_it->second;
                            ++a_it;
                            ++b_it;
                        } else if ( a_it->first < b_it->first ) {
                            ++a_it;
                        } else {
                            ++b_it;
                        }
                    }
                    return product;
                }

                // ELEMENTTYPE _sparse_dot( const
                // sparse_vector<ELEMENTTYPE>& b ) const {
                //     _check_same_size( b );
                //     ELEMENTTYPE product = _zero;
                //     auto inds = this->nonzero_index_intersection( b );
                //     for (auto it = inds.cbegin(); it != inds.cend();
                //     ++it) {
                //         product += this->get( *it ) * b.get( *it );
                //     }
                //     return product;
                // }

                virtual std::map<size_t, ELEMENTTYPE> nonzero()
                    const override {
                    return std::map<size_t, ELEMENTTYPE>( _elements );
                }

                virtual std::vector<size_t> nonzero_indices() const override {
                    std::vector<size_t> nz;
                    for ( auto it = _elements.cbegin(); it != _elements.cend();
                          ++it ) {
                        if ( !NCPA::math::is_zero( it->second ) ) {
                            nz.push_back( it->first );
                        }
                    }
                    return nz;
                }

                virtual size_t count_nonzero_indices() const override {
                    return _elements.size();
                }

                virtual void qc() override { _verify_nonzero(); }

                virtual bool is_zero() const override {
                    return _elements.size() == 0;
                }

                virtual abstract_vector<ELEMENTTYPE>& zero() override {
                    _elements.clear();
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& zero(
                    size_t n ) override {
                    _erase( n );
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& zero(
                    const std::vector<size_t>& n ) override {
                    for ( auto it = n.cbegin(); it != n.cend(); ++it ) {
                        _erase( *it );
                    }
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& zero(
                    std::initializer_list<size_t> n ) override {
                    return zero( std::vector<size_t>( n ) );
                }

                virtual abstract_vector<ELEMENTTYPE>& add(
                    const abstract_vector<ELEMENTTYPE>& b,
                    ELEMENTTYPE modifier = 1.0 ) override {
                    _check_same_size( b );
                    for ( size_t i = 0; i < b.size(); i++ ) {
                        set( i, get( i ) + b.get( i ) * modifier );
                    }
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                virtual abstract_vector<ELEMENTTYPE>& add(
                    ELEMENTTYPE b ) override {
                    for ( auto i = 0; i < size(); i++ ) {
                        set( i, get( i ) + b );
                        // _elements[ i ] += b;
                    }
                    return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                        this );
                }

                typename std::map<size_t, ELEMENTTYPE>::iterator
                    begin() noexcept {
                    return _elements.begin();
                }

                typename std::map<size_t, ELEMENTTYPE>::iterator
                    end() noexcept {
                    return _elements.end();
                }

                typename std::map<size_t, ELEMENTTYPE>::const_iterator cbegin()
                    const noexcept {
                    return _elements.cbegin();
                }

                typename std::map<size_t, ELEMENTTYPE>::const_iterator cend()
                    const noexcept {
                    return _elements.cend();
                }

                friend bool operator==( const sparse_vector<ELEMENTTYPE>& a,
                                        const sparse_vector<ELEMENTTYPE>& b ) {
                    return a.equals( b );
                }

                friend bool operator!=( const sparse_vector<ELEMENTTYPE>& a,
                                        const sparse_vector<ELEMENTTYPE>& b ) {
                    return !( a.equals( b ) );
                }

                friend class band_diagonal_matrix<ELEMENTTYPE>;

            protected:
                void _verify_nonzero() {
                    std::vector<size_t> bad_inds;
                    for ( auto it = _elements.begin(); it != _elements.end();
                          ++it ) {
                        if ( NCPA::math::is_zero( it->second ) ) {
                            bad_inds.push_back( it->first );
                        }
                    }
                    for ( auto it = bad_inds.begin(); it != bad_inds.end();
                          ++it ) {
                        _elements.erase( *it );
                    }
                }

                void _erase( size_t n ) {
                    auto it = _elements.find( n );
                    if ( it != _elements.cend() ) {
                        _elements.erase( it );
                    }
                }

                void _check_same_size(
                    const abstract_vector<ELEMENTTYPE>& b ) const {
                    if ( b.size() != size() ) {
                        throw std::invalid_argument(
                            "Operation requires vectors to be the same "
                            "size" );
                    }
                }

                void _check_out_of_range( size_t n ) const {
                    if ( n >= _capacity ) {
                        std::ostringstream oss;
                        oss << "Index " << n
                            << " exceeds range for vector of size "
                            << _capacity;
                        throw std::range_error( oss.str() );
                    }
                }

            private:
                size_t _capacity = 0;
                std::map<size_t, ELEMENTTYPE> _elements;
                const ELEMENTTYPE _zero = 0;
                // std::vector<ELEMENTTYPE> _elements;
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::sparse_vector<T>& a,
                  NCPA::linear::sparse_vector<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::abstract_vector<T>&>( a ),
            static_cast<NCPA::linear::abstract_vector<T>&>( b ) );
    swap( a._elements, b._elements );
    swap( a._capacity, b._capacity );
}
