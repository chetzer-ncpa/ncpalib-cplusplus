#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra.hpp"
// #include "NCPA/types.hpp"

#include <complex>
// #include <type_traits>
#include <unordered_map>
#include <vector>

namespace NCPA {
    namespace approx {


        template<typename T = std::complex<double>>
        class PadeApproximator {
            public:
                PadeApproximator() {}

                PadeApproximator& calculate(
                    const std::vector<T>& taylor_coefficients,
                    size_t n_numerator, size_t n_denominator ) {
                    PadeApproximator<T>::calculate_from_taylor(
                        taylor_coefficients, n_numerator, n_denominator,
                        _numerator_coefficients, _denominator_coefficients );
                    return *this;
                }

                PadeApproximator& calculate_e( size_t n_numerator,
                                               size_t n_denominator ) {
                    PadeApproximator<T>::calculate_exponential(
                        n_numerator, n_denominator, _numerator_coefficients,
                        _denominator_coefficients );
                    return *this;
                }

                PadeApproximator& apply(
                    const NCPA::linear::MatrixPolynomial<T>& Q,
                    NCPA::linear::Matrix<T>& B, NCPA::linear::Matrix<T>& C ) {
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

                static void calculate_from_taylor(
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

                static void calculate_exponential(
                    size_t n_numerator, size_t n_denominator, 
                    std::vector<T>& numerator_coefficients,
                    std::vector<T>& denominator_coefficients ) {
                    size_t numerator_order   = n_numerator - 1;
                    size_t denominator_order = n_denominator - 1;
                    PadeApproximator::_exponential_coefficients(
                        numerator_order, denominator_order,
                        numerator_coefficients, denominator_coefficients );
                }


            protected:
                std::vector<T> _numerator_coefficients,
                    _denominator_coefficients;

                static void _exponential_coefficients(
                    size_t n_order, size_t d_order, std::vector<T>& num,
                    std::vector<T>& den ) {
                    switch (n_order) {
                        case 0:
                            switch (d_order) {
                                case 0:
                                    num.assign( { 1. } );
                                    den.assign( { 1. } );
                                    break;
                                case 1:
                                    num.assign( { 1. } );
                                    den.assign( { 1., -1. } );
                                    break;
                                case 2:
                                    num.assign( { 2. } );
                                    den.assign( { 2., -2., 1. } );
                                    break;
                                case 3:
                                    num.assign( { 6. } );
                                    den.assign( { 6., -6., 3., -1. } );
                                    break;
                                default:
                                    throw std::range_error(
                                        "Tabulations for exponential "
                                        "coefficients only "
                                        "calculated for n <= 4; calculate "
                                        "manually or "
                                        "reduce order" );
                            }
                            break;
                        case 1:
                            switch (d_order) {
                                case 0:
                                    num.assign( { 1., 1. } );
                                    den.assign( { 1. } );
                                    break;
                                case 1:
                                    num.assign( { 2., 1. } );
                                    den.assign( { 2., -1. } );
                                    break;
                                case 2:
                                    num.assign( { 6., 2. } );
                                    den.assign( { 6., -4., 1. } );
                                    break;
                                case 3:
                                    num.assign( { 24., 6. } );
                                    den.assign( { 24., -18., 6., -1. } );
                                    break;
                                default:
                                    throw std::range_error(
                                        "Tabulations for exponential "
                                        "coefficients only "
                                        "calculated for n <= 4; calculate "
                                        "manually or "
                                        "reduce order" );
                            }
                            break;
                        case 2:
                            switch (d_order) {
                                case 0:
                                    num.assign( { 2., 2., 1. } );
                                    den.assign( { 1. } );
                                    break;
                                case 1:
                                    num.assign( { 6., 4., 1. } );
                                    den.assign( { 6., -2. } );
                                    break;
                                case 2:
                                    num.assign( { 12., 6., 1. } );
                                    den.assign( { 12., -6., 1. } );
                                    break;
                                case 3:
                                    num.assign( { 60., 24., 3. } );
                                    den.assign( { 60., -36., 9., -1. } );
                                    break;
                                default:
                                    throw std::range_error(
                                        "Tabulations for exponential "
                                        "coefficients only "
                                        "calculated for n <= 4; calculate "
                                        "manually or "
                                        "reduce order" );
                            }
                            break;
                        case 3:
                            switch (d_order) {
                                case 0:
                                    num.assign( { 6., 6., 3., 1. } );
                                    den.assign( { 6. } );
                                    break;
                                case 1:
                                    num.assign( { 24., 18., 6., 1. } );
                                    den.assign( { 24., -6. } );
                                    break;
                                case 2:
                                    num.assign( { 60., 36., 9., 1. } );
                                    den.assign( { 60., -24., 3. } );
                                    break;
                                case 3:
                                    num.assign( { 120., 60., 12., 1. } );
                                    den.assign( { 120., -60., 12., -1. } );
                                    break;
                                default:
                                    throw std::range_error(
                                        "Tabulations for exponential "
                                        "coefficients only "
                                        "calculated for n <= 4; calculate "
                                        "manually or "
                                        "reduce order" );
                            }
                            break;
                        default:
                            throw std::range_error(
                                "Tabulations for exponential "
                                "coefficients only "
                                "calculated for n <= 4; calculate "
                                "manually or "
                                "reduce order" );
                            break;
                    }
                    T norm = NCPA::math::one<T>() / num.front();
                    for (auto it = num.begin(); it != num.end(); ++it) {
                        *it *= norm;
                    }
                    for (auto it = den.begin(); it != den.end(); ++it) {
                        *it *= norm;
                    }
                }
        };
    }  // namespace approx
}  // namespace NCPA
