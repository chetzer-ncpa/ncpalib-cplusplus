#pragma once

#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/Matrix.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

template<typename T>
static void swap( NCPA::linear::MatrixPolynomial<T>& a,
                  NCPA::linear::MatrixPolynomial<T>& b ) noexcept;

namespace NCPA {
    namespace linear {

        template<typename T>
        class MatrixPolynomial : public std::vector<Matrix<T>> {
            public:
                MatrixPolynomial() : std::vector<Matrix<T>>() {}

                MatrixPolynomial( const Matrix<T>& base, size_t order = 0 ) :
                    MatrixPolynomial<T>() {
                    this->set_base( base );
                    this->compute( order );
                }

                MatrixPolynomial( const MatrixPolynomial<T>& other ) :
                    MatrixPolynomial<T>() {
                    this->assign( other.cbegin(), other.cend() );
                }

                MatrixPolynomial( MatrixPolynomial<T>&& source ) noexcept :
                    MatrixPolynomial<T>() {
                    ::swap( *this, source );
                }

                virtual ~MatrixPolynomial() {}

                MatrixPolynomial<T>& operator=( MatrixPolynomial<T> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void ::swap<T>( MatrixPolynomial<T>& a,
                                       MatrixPolynomial<T>& b ) noexcept;

                virtual int order() const { return this->size() - 1; }

                virtual MatrixPolynomial<T>& compute( size_t order ) {
                    if (this->empty()) {
                        throw std::logic_error( "MatrixPolynomial.compute(): "
                                                "No base matrix set!" );
                    }
                    size_t n_mats = order + 1;
                    while (this->size() < n_mats) {
                        this->_compute_next_power();
                    }
                    return *this;
                }

                virtual MatrixPolynomial<T>& set_base(
                    const Matrix<T>& base ) {
                    this->assign( 1, base );
                    return *this;
                }

                virtual Matrix<T>& get( size_t power ) {
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
static void swap( NCPA::linear::MatrixPolynomial<T>& a,
                  NCPA::linear::MatrixPolynomial<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<std::vector<NCPA::linear::Matrix<T>>&>( a ),
            static_cast<std::vector<NCPA::linear::Matrix<T>>&>( b ) );
}
