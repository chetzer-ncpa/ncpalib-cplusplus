#pragma once

#include "NCPA/linearalgebra.hpp"
// #include "NCPA/types.hpp"

#include <complex>
// #include <type_traits>
#include <vector>

namespace NCPA {
    namespace approx {

        template<typename T = std::complex<double>>
        class PadeApproximator {
            public:
                PadeApproximator() {}

                PadeApproximator& calculate( const std::vector<T>& taylor_coefficients,
                                size_t n_numerator, size_t n_denominator ) {
                    PadeApproximator<T>::calculate(
                        taylor_coefficients, n_numerator, n_denominator,
                        _numerator_coefficients, _denominator_coefficients );
                    return *this;
                }

                PadeApproximator& apply( const NCPA::linear::MatrixPolynomial<T>& Q,
                            NCPA::linear::Matrix<T>& B,
                            NCPA::linear::Matrix<T>& C ) {
                    size_t NQ = Q.at( 0 ).rows();

                    B.clear().identity( NQ, NQ );
                    C.clear().identity( NQ, NQ );
                    for (auto i = 1; i < _denominator_coefficients.size();
                         i++) {
                        C += Q[ i - 1 ] * _denominator_coefficients[ i ];
                    }
                    for (auto i = 1; i < _numerator_coefficients.size(); i++) {
                        B += Q[ i - 1 ] * _numerator_coefficients[ i ];
                    }
                    return *this;
                }

                static void calculate(
                    const std::vector<T>& taylor_coefficients,
                    size_t n_numerator, size_t n_denominator,
                    std::vector<T>& numerator_coefficients,
                    std::vector<T>& denominator_coefficients ) {
                    // sanity checks
                    if (n_denominator < n_numerator) {
                        throw std::invalid_argument(
                            "Denominator count must be >= numerator "
                            "count for Pade calculation" );
                    }
                    size_t n        = n_numerator - 1;    // numerator order
                    size_t m        = n_denominator - 1;  // denominator order
                    size_t N        = n + m;
                    size_t n_taylor = taylor_coefficients.size();
                    if (n_taylor < ( N + 1 )) {
                        std::ostringstream oss;
                        oss << "Count of Taylor series must be at least "
                            << ( N + 1 ) << " for numerator count "
                            << n_numerator << " and denominator count "
                            << n_denominator;
                        throw std::invalid_argument( oss.str() );
                    }

                    NCPA::linear::Solver<T> solver
                        = NCPA::linear::SolverFactory<T>::build(
                            NCPA::linear::solver_t::BASIC );
                    NCPA::linear::Matrix<T> A
                        = NCPA::linear::MatrixFactory<T>::build(
                            NCPA::linear::matrix_t::DENSE );
                    NCPA::linear::Vector<T> y
                        = NCPA::linear::VectorFactory<T>::build(
                            NCPA::linear::vector_t::DENSE );
                    A.resize( N, N );
                    y.resize( N );
                    T tempsc( -1.0, 0.0 );
                    for (size_t ii = 0; ii < n; ii++) {
                        A.set( ii, ii, tempsc );
                    }
                    for (size_t ii = 0; ii < N; ii++) {
                        for (size_t jj = n;
                             jj <= std::min( A.columns() - 1, ii + n ); jj++) {
                            A.set( ii, jj,
                                   taylor_coefficients.at( ii - jj + n ) );
                            y[ ii ] = -taylor_coefficients[ ii + 1 ];
                        }
                    }

                    solver.set_system_matrix( A );
                    NCPA::linear::Vector<T> x = solver.solve( y );

                    numerator_coefficients.clear();
                    denominator_coefficients.clear();
                    numerator_coefficients.push_back(
                        taylor_coefficients[ 0 ] );
                    for (auto ii = 0; ii < n; ii++) {
                        numerator_coefficients.push_back( x[ ii ] );
                    }
                    denominator_coefficients.push_back( T( 1.0, 0.0 ) );
                    for (auto ii = n; ii < N; ii++) {
                        denominator_coefficients.push_back( x[ ii ] );
                    }
                }

                std::vector<T> numerator() const {
                    return _numerator_coefficients;
                }

                std::vector<T> denominator() const {
                    return _denominator_coefficients;
                }


            protected:
                std::vector<T> _numerator_coefficients,
                    _denominator_coefficients;
        };
    }  // namespace approx
}  // namespace NCPA
