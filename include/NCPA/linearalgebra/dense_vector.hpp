#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
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

namespace NCPA {
    namespace linear {
        namespace details {
            NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE( dense_vector,
                                                         abstract_vector );
        }
    }  // namespace linear
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::dense_vector, ELEMENTTYPE );

namespace NCPA {
    namespace linear {

        namespace details {
            NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
                class dense_vector<ELEMENTTYPE,
                                   _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
                : public abstract_vector<ELEMENTTYPE> {
                public:
                    dense_vector() {}

                    dense_vector( const std::vector<ELEMENTTYPE>& v ) :
                        dense_vector<ELEMENTTYPE>() {
                        set( v );
                    }

                    dense_vector( size_t n ) : dense_vector<ELEMENTTYPE>() {
                        resize( n );
                    }

                    dense_vector(
                        const std::initializer_list<ELEMENTTYPE>& v ) :
                        dense_vector<ELEMENTTYPE>() {
                        set( std::vector<ELEMENTTYPE>( v ) );
                    }

                    // copy constructor
                    dense_vector( const dense_vector<ELEMENTTYPE>& other ) :
                        dense_vector<ELEMENTTYPE>() {
                        set( other._elements );
                    }

                    /**
                     * Move constructor.
                     * @param source The vector to assimilate.
                     */
                    dense_vector( dense_vector<ELEMENTTYPE>&& source ) noexcept
                        :
                        dense_vector<ELEMENTTYPE>() {
                        ::swap( *this, source );
                    }

                    virtual ~dense_vector() {}

                    friend void ::swap<ELEMENTTYPE>(
                        dense_vector<ELEMENTTYPE>& a,
                        dense_vector<ELEMENTTYPE>& b ) noexcept;

                    /**
                     * Assignment operator.
                     * @param other The vector to assign to this.
                     */
                    dense_vector<ELEMENTTYPE>& operator=(
                        dense_vector<ELEMENTTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual abstract_vector<ELEMENTTYPE>& set(
                        size_t n, const ELEMENTTYPE *val ) override {
                        resize( n );
                        std::copy( val, val + n, begin() );
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& set(
                        const std::vector<ELEMENTTYPE>& v ) override {
                        _elements = v;
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& set(
                        ELEMENTTYPE val ) override {
                        for ( size_t i = 0; i < size(); i++ ) {
                            set( i, val );
                        }
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual size_t size() const override {
                        return _elements.size();
                    }

                    virtual abstract_vector<ELEMENTTYPE>& zero(
                        size_t n ) override {
                        return set( n, _zero );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& zero(
                        const std::vector<size_t>& n ) override {
                        for ( auto it = n.cbegin(); it != n.cend(); ++it ) {
                            _elements.at( *it ) = _zero;
                        }
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& zero(
                        std::initializer_list<size_t> n ) override {
                        return zero( std::vector<size_t>( n ) );
                    }

                    virtual std::vector<size_t> nonzero_indices()
                        const override {
                        std::vector<size_t> nz;
                        for ( auto i = 0; i < size(); i++ ) {
                            if ( !NCPA::math::is_zero( _elements[ i ] ) ) {
                                nz.push_back( i );
                            }
                        }
                        return nz;
                    }

                    virtual size_t count_nonzero_indices() const override {
                        return nonzero_indices().size();
                    }

                    virtual const ELEMENTTYPE& get( size_t n ) const override {
                        return _elements.at( n );
                    }

                    virtual ELEMENTTYPE& get( size_t n ) override {
                        return _elements.at( n );
                    }

                    virtual std::vector<ELEMENTTYPE> as_std() const override {
                        std::vector<ELEMENTTYPE> v( _elements );
                        return v;
                    }

                    virtual std::unique_ptr<abstract_vector<ELEMENTTYPE>>
                        clone() override {
                        return std::unique_ptr<abstract_vector<ELEMENTTYPE>>(
                            new dense_vector( *this ) );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& clear() override {
                        _elements.clear();
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& resize(
                        size_t n ) override {
                        _elements.resize( n );
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
                        std::copy( cbegin(), cend(), vals );
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& set(
                        size_t n, ELEMENTTYPE val ) override {
                        _elements[ n ] = val;
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& scale(
                        ELEMENTTYPE val ) override {
                        for ( auto it = begin(); it != end(); ++it ) {
                            *it *= val;
                        }
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& scale(
                        const abstract_vector<ELEMENTTYPE>& b ) override {
                        _check_same_size( b );
                        for ( auto i = 0; i < size(); i++ ) {
                            _elements[ i ] *= b.get( i );
                        }
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual ELEMENTTYPE dot(
                        const abstract_vector<ELEMENTTYPE>& b )
                        const override {
                        _check_same_size( b );
                        ELEMENTTYPE product = 0.0;
                        for ( auto i = 0; i < size(); i++ ) {
                            product += _elements[ i ] * b.get( i );
                        }
                        return product;
                    }

                    virtual abstract_vector<ELEMENTTYPE>& add(
                        const abstract_vector<ELEMENTTYPE>& b,
                        ELEMENTTYPE modifier = 1.0 ) override {
                        _check_same_size( b );
                        for ( auto i = 0; i < size(); i++ ) {
                            _elements[ i ] += b.get( i ) * modifier;
                        }
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_vector<ELEMENTTYPE>& add(
                        ELEMENTTYPE b ) override {
                        for ( auto it = begin(); it != end(); ++it ) {
                            *it += b;
                        }
                        return *dynamic_cast<abstract_vector<ELEMENTTYPE> *>(
                            this );
                    }

                    typename std::vector<ELEMENTTYPE>::iterator
                        begin() noexcept {
                        return _elements.begin();
                    }

                    typename std::vector<ELEMENTTYPE>::iterator
                        end() noexcept {
                        return _elements.end();
                    }

                    typename std::vector<ELEMENTTYPE>::const_iterator cbegin()
                        const noexcept {
                        return _elements.cbegin();
                    }

                    typename std::vector<ELEMENTTYPE>::const_iterator cend()
                        const noexcept {
                        return _elements.cend();
                    }

                    friend bool operator==(
                        const dense_vector<ELEMENTTYPE>& a,
                        const dense_vector<ELEMENTTYPE>& b ) {
                        return a.equals( b );
                    }

                    friend bool operator!=(
                        const dense_vector<ELEMENTTYPE>& a,
                        const dense_vector<ELEMENTTYPE>& b ) {
                        return !( a.equals( b ) );
                    }

                protected:
                    void _check_same_size(
                        const abstract_vector<ELEMENTTYPE>& b ) const {
                        if ( b.size() != size() ) {
                            throw std::invalid_argument(
                                "Dot product requires vectors to be the same "
                                "size" );
                        }
                    }

                private:
                    std::vector<ELEMENTTYPE> _elements;
                    const ELEMENTTYPE _zero = 0;
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::details::dense_vector<T>& a,
                  NCPA::linear::details::dense_vector<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::linear::details::abstract_vector<T>&>( a ),
            static_cast<NCPA::linear::details::abstract_vector<T>&>( b ) );
    swap( a._elements, b._elements );
}