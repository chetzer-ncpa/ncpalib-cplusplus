#pragma once

#include "NCPA/arrays.hpp"
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

// forward declarations for operators and friend functions
namespace NCPA {
    namespace linear {
        namespace details {
            template<typename ELEMENTTYPE>
            class abstract_vector;
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::abstract_vector, ELEMENTTYPE );

namespace NCPA {
    namespace linear {

        namespace details {
            template<typename ELEMENTTYPE>
            class abstract_vector {
                public:
                    virtual ~abstract_vector() = default;
                    friend void ::swap<ELEMENTTYPE>(
                        abstract_vector<ELEMENTTYPE>& a,
                        abstract_vector<ELEMENTTYPE>& b ) noexcept;

                    virtual size_t size() const                      = 0;
                    virtual ELEMENTTYPE& get( size_t n )             = 0;
                    virtual const ELEMENTTYPE& get( size_t n ) const = 0;
                    virtual std::vector<ELEMENTTYPE> as_std() const  = 0;
                    virtual abstract_vector<ELEMENTTYPE>& clear()    = 0;
                    virtual std::unique_ptr<abstract_vector> clone() = 0;
                    virtual abstract_vector<ELEMENTTYPE>& resize( size_t n )
                        = 0;
                    virtual abstract_vector<ELEMENTTYPE>& as_array(
                        size_t& n, ELEMENTTYPE *& vals )
                        = 0;
                    virtual abstract_vector<ELEMENTTYPE>& set(
                        size_t n, ELEMENTTYPE val )
                        = 0;
                    virtual abstract_vector<ELEMENTTYPE>& set(
                        size_t n, const ELEMENTTYPE *val )
                        = 0;
                    virtual abstract_vector<ELEMENTTYPE>& set(
                        const std::vector<ELEMENTTYPE>& v )
                        = 0;
                    virtual abstract_vector<ELEMENTTYPE>& set(
                        ELEMENTTYPE val )
                        = 0;

                    virtual abstract_vector<ELEMENTTYPE>& zero( size_t n ) = 0;
                    virtual abstract_vector<ELEMENTTYPE>& zero(
                        const std::vector<size_t>& n )
                        = 0;
                    virtual abstract_vector<ELEMENTTYPE>& zero(
                        std::initializer_list<size_t> n )
                        = 0;

                    virtual abstract_vector<ELEMENTTYPE>& scale(
                        ELEMENTTYPE val )
                        = 0;

                    virtual abstract_vector<ELEMENTTYPE>& scale(
                        const abstract_vector<ELEMENTTYPE>& b )
                        = 0;

                    virtual abstract_vector<ELEMENTTYPE>& add(
                        const abstract_vector<ELEMENTTYPE>& b,
                        ELEMENTTYPE modifier = 1.0 )
                        = 0;

                    virtual abstract_vector<ELEMENTTYPE>& add( ELEMENTTYPE b )
                        = 0;

                    virtual ELEMENTTYPE dot(
                        const abstract_vector<ELEMENTTYPE>& b ) const
                        = 0;

                    virtual size_t count_nonzero_indices() const        = 0;
                    virtual std::vector<size_t> nonzero_indices() const = 0;

                    // implementations, not abstract
                    template<typename OTHERTYPE>
                    bool equals(
                        const abstract_vector<OTHERTYPE>& other ) const {
                        return false;
                    }

                    virtual bool equals(
                        const abstract_vector<ELEMENTTYPE>& other ) const {
                        if ( size() != other.size() ) {
                            return false;
                        }
                        for ( auto i = 0; i < size(); i++ ) {
                            if ( get( i ) != other.get( i ) ) {
                                return false;
                            }
                        }
                        return true;
                    }

                    virtual abstract_vector<ELEMENTTYPE>& operator+=(
                        const abstract_vector<ELEMENTTYPE>& other ) {
                        this->add( other );
                        return *this;
                    }

                    virtual abstract_vector<ELEMENTTYPE>& operator+=(
                        const ELEMENTTYPE& other ) {
                        this->add( other );
                        return *this;
                    }

                    virtual abstract_vector<ELEMENTTYPE>& operator-=(
                        const abstract_vector<ELEMENTTYPE>& other ) {
                        this->add( other, -1.0 );
                        return *this;
                    }

                    virtual abstract_vector<ELEMENTTYPE>& operator-=(
                        const ELEMENTTYPE& other ) {
                        this->add( -other );
                        return *this;
                    }

                    virtual abstract_vector<ELEMENTTYPE>& operator*=(
                        const abstract_vector<ELEMENTTYPE>& other ) {
                        this->scale( other );
                        return *this;
                    }

                    virtual abstract_vector<ELEMENTTYPE>& operator*=(
                        const ELEMENTTYPE& other ) {
                        this->scale( other );
                        return *this;
                    }

                    virtual abstract_vector<ELEMENTTYPE>& operator/=(
                        const ELEMENTTYPE& other ) {
                        this->scale( 1.0 / other );
                        return *this;
                    }

                    virtual ELEMENTTYPE& operator[]( size_t i ) {
                        return get( i );
                    }

                    virtual const ELEMENTTYPE& operator[]( size_t i ) const {
                        return get( i );
                    }

                    friend bool operator==(
                        const abstract_vector<ELEMENTTYPE>& a,
                        const abstract_vector<ELEMENTTYPE>& b ) {
                        return a.equals( b );
                    }

                    friend bool operator!=(
                        const abstract_vector<ELEMENTTYPE>& a,
                        const abstract_vector<ELEMENTTYPE>& b ) {
                        return !( a.equals( b ) );
                    }
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::details::abstract_vector<T>& a,
                  NCPA::linear::details::abstract_vector<T>& b ) noexcept {}
