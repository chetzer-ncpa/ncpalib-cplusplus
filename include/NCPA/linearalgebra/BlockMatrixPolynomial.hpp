#pragma once

#include "NCPA/linearalgebra/BlockMatrix.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/Matrix.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

template<typename T>
static void swap( NCPA::linear::BlockMatrixPolynomial<T>& a,
                  NCPA::linear::BlockMatrixPolynomial<T>& b ) noexcept;

namespace NCPA {
    namespace linear {

        template<typename T>
        class BlockMatrixPolynomial : public std::vector<BlockMatrix<T>> {
            public:
                BlockMatrixPolynomial() : std::vector<BlockMatrix<T>>() {}

                BlockMatrixPolynomial( const BlockMatrix<T>& base,
                                       size_t order = 0 ) :
                    BlockMatrixPolynomial<T>() {
                    this->set_base( base );
                    this->compute( order );
                }

                BlockMatrixPolynomial(
                    const BlockMatrixPolynomial<T>& other ) :
                    BlockMatrixPolynomial<T>() {
                    this->assign( other.cbegin(), other.cend() );
                }

                BlockMatrixPolynomial(
                    BlockMatrixPolynomial<T>&& source ) noexcept :
                    BlockMatrixPolynomial<T>() {
                    ::swap( *this, source );
                }

                virtual ~BlockMatrixPolynomial() {}

                BlockMatrixPolynomial<T>& operator=(
                    BlockMatrixPolynomial<T> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void ::swap<T>( BlockMatrixPolynomial<T>& a,
                                       BlockMatrixPolynomial<T>& b ) noexcept;

                virtual int order() const { return this->size() - 1; }

                virtual BlockMatrixPolynomial<T>& compute( size_t order ) {
                    if (this->empty()) {
                        throw std::logic_error(
                            "BlockMatrixPolynomial.compute(): "
                            "No base matrix set!" );
                    }
                    size_t n_mats = order + 1;
                    while (this->size() < n_mats) {
                        this->_compute_next_power();
                    }
                    return *this;
                }

                virtual BlockMatrixPolynomial<T>& set_base(
                    const BlockMatrix<T>& base ) {
                    this->assign( 1, base );
                    return *this;
                }

                virtual BlockMatrix<T>& get( size_t power ) {
                    this->compute( power );
                    return this->at( power - 1 );
                }

            protected:
                void _compute_next_power() {
                    this->push_back( this->front() * this->back() );
                }
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::BlockMatrixPolynomial<T>& a,
                  NCPA::linear::BlockMatrixPolynomial<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<std::vector<NCPA::linear::BlockMatrix<T>>&>( a ),
            static_cast<std::vector<NCPA::linear::BlockMatrix<T>>&>( b ) );
}
