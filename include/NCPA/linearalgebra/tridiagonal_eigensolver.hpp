#pragma once

#ifndef NCPA_LINEAR_TRIDIAGONAL_EIGENSOLVER_TOLERANCE
#  define NCPA_LINEAR_TRIDIAGONAL_EIGENSOLVER_TOLERANCE 1.0e-8
#endif

#include "NCPA/linearalgebra/Matrix.hpp"
#include "NCPA/linearalgebra/Vector.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

static void swap( NCPA::linear::TridiagonalEigensolver& a,
                  NCPA::linear::TridiagonalEigensolver& b ) noexcept;

namespace NCPA {
    namespace linear {
        class TridiagonalEigensolver {
            public:
                TridiagonalEigensolver() {
                    _Q = MatrixFactory<std::complex<double>>::build(
                        matrix_t::BAND_DIAGONAL );
                }

                TridiagonalEigensolver(
                    const NCPA::linear::Matrix<std::complex<double>>& m ) :
                    TridiagonalEigensolver() {
                    _Q.copy( m );
                }

                // copy constructor
                TridiagonalEigensolver( const TridiagonalEigensolver& other ) :
                    TridiagonalEigensolver() {
                    _lmax  = other._lmax;
                    _lmin  = other._lmin;
                    _eigs  = other._eigs;
                    _eigvs = other._eigvs;
                    _Q     = other._Q;
                }

                // move constructor
                TridiagonalEigensolver(
                    TridiagonalEigensolver&& source ) noexcept :
                    TridiagonalEigensolver() {
                    ::swap( *this, source );
                }

                virtual ~TridiagonalEigensolver() {}

                friend void ::swap( TridiagonalEigensolver& a,
                                    TridiagonalEigensolver& b ) noexcept;

                TridiagonalEigensolver& operator=(
                    TridiagonalEigensolver other ) {
                    ::swap( *this, other );
                    return *this;
                }

                size_t count() {
                    size_t cmax = _sturm_count( _lmax );
                    size_t cmin = _sturm_count( _lmin );
                    return ( cmax >= cmin ) ? cmax - cmin : 0;
                }

                size_t count( double lambda ) {
                    return _sturm_count( lambda );
                }

                size_t count( double lmin, double lmax ) {
                    size_t maxcount = _sturm_count( lmax );
                    size_t mincount = _sturm_count( lmin );
                    if (maxcount < mincount) {
                        throw std::logic_error(
                            "Eigenvalue counter returns negative count!" );
                    }
                    return maxcount - mincount;
                }

                virtual const std::vector<double>& eigenvalues() const {
                    return _eigs;
                }

                virtual const std::vector<std::vector<std::complex<double>>>&
                    eigenvectors() const {
                    return _eigvs;
                }

                virtual TridiagonalEigensolver& set(
                    const NCPA::linear::Matrix<std::complex<double>>& m ) {
                    if (!m.is_square() || !m.is_tridiagonal()
                        || !m.is_symmetric()) {
                        throw std::invalid_argument(
                            "Matrix must be square, tridiagonal and "
                            "symmetric" );
                    }

                    _Q = m;
                    return *this;
                }

                virtual TridiagonalEigensolver& set_range( double kkmin,
                                                           double kkmax ) {
                    _lmin = kkmin;
                    _lmax = kkmax;
                    return *this;
                }

                virtual TridiagonalEigensolver& solve( bool find_eigenvecs
                                                       = true ) {
                    _calculate_eigenvalues( _lmin, _lmax );
                    if (find_eigenvecs) {
                        _calculate_eigenvectors();
                    }
                    return *this;
                }

                virtual TridiagonalEigensolver& solve( double kmin,
                                                       double kmax,
                                                       bool find_eigenvecs
                                                       = true ) {
                    set_range( kmin, kmax );
                    return solve( find_eigenvecs );
                }

                virtual TridiagonalEigensolver& solve(
                    const Matrix<std::complex<double>>& M, double kmin,
                    double kmax, bool find_eigenvecs = true ) {
                    this->set( M );
                    return this->solve( kmin, kmax, find_eigenvecs );
                }

            protected:
                // size_t eig_count  = 0;
                double _lmax = -1.0, _lmin = -1.0;
                std::vector<double> _eigs;
                std::vector<std::vector<std::complex<double>>> _eigvs;

                Matrix<std::complex<double>> _Q;

                static constexpr double ERR = 1.0e-8;

                void _calculate_eigenvalues( double lambda1, double lambda2 ) {
                    double k1 = std::min( lambda1, lambda2 );
                    double k2 = std::max( lambda1, lambda2 );
                    int nphi
                        = (int)_sturm_count( k2 ) - (int)_sturm_count( k1 );
                    std::vector<double> phi( nphi, 0.0 );
                    double kstep = 1.0e-8 * (double)( k2 - k1 );
                    for (size_t i = 0; i < nphi; ++i) {
                        if (!_get_first_eigenvalue( k1, k2, phi[ i ] )) {
                            throw std::logic_error(
                                "Unexpected lack of eigenvalues in range!" );
                        }
                        k1 = phi[ i ] + kstep;
                    }
                    _eigs = phi;
                }

                std::vector<std::complex<double>> _calculate_eigenvector(
                    double lambda ) {
                    double tol      = 1e-10;
                    size_t max_iter = 10;
                    size_t rows     = _Q.rows();
                    int NN          = rows - 1;
                    double dNNp1    = (double)( NN + 1 );
                    Vector<std::complex<double>> yy
                        = VectorFactory<std::complex<double>>::build(
                            vector_t::DENSE );
                    yy.resize( NN + 1 ).set( 1.0 );
                    Matrix<std::complex<double>> Qpert = _Q;
                    Qpert.set_diagonal( *( _Q.get_diagonal( 0 ) ) - lambda );
                    Solver<std::complex<double>> solver
                        = SolverFactory<std::complex<double>>::build(
                            solver_t::TRIDIAGONAL );
                    solver.set_system_matrix( Qpert );
                    double dy = 0.0;
                    Vector<std::complex<double>> xx;
                    Vector<std::complex<double>> lasty = yy;

                    size_t counter = 0;
                    do {
                        xx = solver.solve( yy );
                        std::complex<double> norm( 0.0, 0.0 );
                        dy = 0.0;
                        for (size_t i = 0; i <= NN; ++i) {
                            norm += xx.get( i ) * xx.get( i );
                        }
                        norm = std::pow( norm, -0.5 );
                        yy   = xx * norm;
                        for (size_t i = 0; i <= NN; ++i) {
                            dy += ( std::abs( yy.get( i ) )
                                    - std::abs( lasty.get( i ) ) )
                                / dNNp1;
                        }
                        lasty = yy;
                    } while (++counter <= max_iter && std::abs( dy ) > tol);
                    return yy.as_std();
                }

                void _calculate_eigenvectors() {
                    double cup, pot;
                    _eigvs.clear();
                    _eigvs.resize( _eigs.size() );
                    for (size_t i = 0; i < _eigs.size(); ++i) {
                        std::vector<std::complex<double>> evec
                            = _calculate_eigenvector( _eigs[ i ] );
                        _eigvs[ i ] = evec;
                    }
                }

                bool _get_first_eigenvalue( double k1, double k2,
                                            double& phi_n,
                                            double tol = 1e-8 ) {
                    phi_n = 0.0;
                    if (k2 < k1) {
                        std::swap( k1, k2 );
                    }
                    size_t lower_count = _sturm_count( k1 );
                    size_t upper_count = _sturm_count( k2 );
                    if (upper_count <= lower_count) {
                        return false;
                    }
                    double dk = 0.2 * ( k2 - k1 );
                    phi_n     = k2;
                    while (dk >= tol) {
                        dk          *= 0.5;
                        phi_n       -= dk;
                        upper_count  = _sturm_count( phi_n );
                        while (upper_count > lower_count) {
                            phi_n       -= dk;
                            upper_count  = _sturm_count( phi_n );
                        }
                        phi_n += dk;
                    }
                    return true;
                }

                size_t _sturm_count( double lambda ) {
                    if (!_Q.is_square() || !_Q.is_tridiagonal()) {
                        throw std::invalid_argument(
                            "Matrix must be square and tridiagonal" );
                    }
                    size_t n            = 0;
                    size_t rows         = _Q.rows();
                    double last_nonzero = 1.0;
                    std::vector<double> P( rows + 1 );
                    P[ 0 ] = 1.0;
                    P[ 1 ] = _Q.get( 0, 0 ).real() - lambda;
                    if (P[ 1 ] != 0.0) {
                        last_nonzero = P[ 1 ];
                        if (P[ 1 ] < 0.0) {
                            n++;
                        }
                    }
                    for (int i = 2; i <= rows; ++i) {
                        P[ i ] = ( _Q.get( i - 1, i - 1 ).real() - lambda )
                                   * P[ i - 1 ]
                               - _Q.get( i - 1, i - 2 ).real()
                                     * _Q.get( i - 1, i - 2 ).real()
                                     * P[ i - 2 ];
                        if (P[ i ] != 0.0) {
                            if (P[ i ] * last_nonzero < 0.0) {
                                n++;
                            }
                            last_nonzero  = P[ i ];
                            P[ i - 1 ]   /= std::abs( P[ i ] );
                            P[ i ]       /= std::abs( P[ i ] );
                        }
                    }
                    return n;
                }
        };
    }  // namespace linear
}  // namespace NCPA

static void swap( NCPA::linear::TridiagonalEigensolver& a,
                  NCPA::linear::TridiagonalEigensolver& b ) noexcept {
    using std::swap;
    swap( a._lmin, b._lmin );
    swap( a._lmax, b._lmax );
    swap( a._eigs, b._eigs );
    swap( a._eigvs, b._eigvs );
    swap( a._Q, b._Q );
}
