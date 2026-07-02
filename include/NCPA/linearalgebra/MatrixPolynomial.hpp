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

                // contents: this->at(0) = Q, this->at(1) = Q^2, etc.
                // So order() = this->size(), 
                virtual int order() const { return this->size(); }

                virtual MatrixPolynomial<T>& compute( size_t order ) {
                    if (this->empty()) {
                        throw std::logic_error( "MatrixPolynomial.compute(): "
                                                "No base matrix set!" );
                    }
                    // size_t n_mats = order + 1;
                    while (this->size() < order) {
                        this->_compute_next_power();
                    }
                    return *this;
                }

                virtual MatrixPolynomial<T>& set_base(
                    const Matrix<T>& base ) {
                    this->assign( 1, base );
                    return *this;
                }

                virtual Matrix<T>& power( size_t pow ) {
                    this->compute( pow );
                    return this->at( pow - 1 );
                }

                virtual Matrix<T> scale_and_sum(
                    const std::vector<T>& coeffs ) {
                    // if (coeffs.size() != this->size()) {
                    //     std::ostringstream oss;
                    //     oss << "Matrix polynomial order " << this->size()
                    //         << " and coefficient vector size " << coeffs.size()
                    //         << " do not match!";
                    //     throw std::out_of_range( oss.str() );
                    // }
                    Matrix<T> sum = this->front();
                    sum.identity().scale( coeffs[ 0 ] );

                    for (size_t i = 1; i < coeffs.size(); ++i) {
                        sum += this->power( i ) * coeffs[ i ];
                    }
                    return sum;
                }

                virtual Matrix<T> scale_and_sum(
                    const Vector<T>& coeffs ) {
                    return this->scale_and_sum( coeffs.as_std() );
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
