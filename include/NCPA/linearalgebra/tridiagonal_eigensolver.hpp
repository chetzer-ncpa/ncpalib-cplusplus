#pragma once

#ifndef NCPA_LINEAR_TRIDIAGONAL_EIGENSOLVER_TOLERANCE
#  define NCPA_LINEAR_TRIDIAGONAL_EIGENSOLVER_TOLERANCE 1.0e-12
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

namespace NCPA {
    namespace linear {

        template<typename ELEMENTTYPE>
        class TridiagonalEigensolver {
            public:
                TridiagonalEigensolver() {}

                TridiagonalEigensolver(
                    const NCPA::linear::Matrix<ELEMENTTYPE>& m ) {}

                virtual ~TridiagonalEigensolver() {}

                virtual const std::vector<std::complex<ELEMENTTYPE>>&
                    eigenvalues() const {
                    return eigs;
                }

                virtual void set( const NCPA::linear::Matrix<ELEMENTTYPE>& m,
                                  bool check = true ) {
                    if (!m.is_square() || !m.is_tridiagonal()
                        || !m.is_symmetric()) {
                        throw std::invalid_argument(
                            "Matrix must be square, tridiagonal and "
                            "symmetric" );
                    }

                    if (check) {
                        check_offdiags( m );
                    }

                    set_diagonal( m.get_diagonal( 0 )->as_std() );
                    set_offdiagonal( m.get( 0, 1 ) );
                }

                virtual void set_diagonal(
                    const std::vector<ELEMENTTYPE>& d ) {
                    diag = d;
                }

                virtual void set_wavenumber_limits( ELEMENTTYPE kmin,
                                                    ELEMENTTYPE kmax ) {
                    kk_min = kmin * kmin;
                    kk_max = kmax * kmax;
                }

                virtual void set_offdiagonal( ELEMENTTYPE od ) { o_d_2 = od; }

                virtual void set_dz( ELEMENTTYPE dz ) { del_z = dz; }

                virtual void solve() { find_eigs(); }

                virtual void solve( ELEMENTTYPE kmin, ELEMENTTYPE kmax ) {
                    set_wavenumber_limits( kmin, kmax );
                    solve();
                }

            protected:
                size_t eig_count = 0;
                ELEMENTTYPE o_d_2 = 0.0, kk_max=-1.0, kk_min = -1.0, del_z = -1.0;

                std::complex<ELEMENTTYPE> admittence = 0.0;
                std::vector<std::complex<ELEMENTTYPE>> eigs;
                std::vector<ELEMENTTYPE> diag, Pn;

                const ELEMENTTYPE ERR = 1.0e-8;

                void charpoly( ELEMENTTYPE lambda ) {
                    size_t N = diag.size() + 1;
                    Pn.resize( N) ;
                    Pn[0] = 1.0;
                    Pn[1] = diag[0] - lambda;
                    for (size_t i = 2; i <= N; ++i) {
                        Pn[i] = (diag[i-1] - lambda) * Pn[i-1] - o_d_2 * o_d_2 * Pn[i-2];
                    }
                }

                

                void check_offdiags(
                    const NCPA::linear::Matrix<ELEMENTTYPE>& m,
                    ELEMENTTYPE tol
                    = NCPA_LINEAR_TRIDIAGONAL_EIGENSOLVER_TOLERANCE ) {
                    size_t noffdiag = m.rows() - 1;
                    auto udiagptr   = m.get_diagonal( 1 );
                    auto ldiagptr   = m.get_diagonal( -1 );
                    ELEMENTTYPE val = udiagptr->get( 0 );
                    for (size_t i = 0; i < noffdiag; ++i) {
                        if (std::abs( udiagptr->get( i ) - val ) > tol) {
                            throw std::invalid_argument(
                                "Matrix off-diagonals not constant within "
                                "tolerance" );
                        }
                        if (std::abs( ldiagptr->get( i ) - val ) > tol) {
                            throw std::invalid_argument(
                                "Matrix off-diagonals not constant within "
                                "tolerance" );
                        }
                    }
                }

                int eig_counter() {
                    return sturm_count( kk_max ) - sturm_count( kk_min );
                }

                void find_eigs() {
                    ELEMENTTYPE diff, eig_cup;
                    int i;
                    ELEMENTTYPE R_bc = admittence.real();

                    // size_t eig_count = eig_counter();
                    size_t eig_count = diag.size();
                    eigs.resize( eig_count );
                    diff = 1.0;
                    kk_min = 1e6;
                    kk_max = -1e6;
                    // diff    = kk_max - kk_min;
                    // eig_cup = kk_min + 5.0 * R_bc;
                    eig_cup = kk_min;
                    for (i = 0; i < eig_count; i++) {
                        eigs[ i ] = find_first_eig( eig_cup, kk_max - R_bc );
                        eig_cup   = eigs[ i ].real() + ERR * diff;
                    }
                }

                ELEMENTTYPE find_first_eig( ELEMENTTYPE lower,
                                            ELEMENTTYPE upper ) {
                    ELEMENTTYPE lambda, d_lambda;
                    int lower_count, upper_count, int_cup;

                    lower_count = sturm_count( lower );
                    upper_count = sturm_count( upper );
                    if (upper_count <= lower_count) {
                        std::ostringstream oss;
                        oss << "No eigenvalues in specified interval " << upper
                            << " to " << lower << "   " << eig_count;
                        throw std::range_error( oss.str() );
                    }
                    lambda = upper;

                    d_lambda = 0.2 * ( lower - upper );
                    while (d_lambda >= ERR * ( lower - upper )) {
                        d_lambda = 0.5 * d_lambda;
                        lambda   = lambda + d_lambda;
                        while (( int_cup = sturm_count( lambda ) )
                               > lower_count) {
                            lambda = lambda + d_lambda;
                        }
                        lambda = lambda - d_lambda;
                    }

                    return ( lambda );
                }

                int sturm_count( ELEMENTTYPE lambda ) {
                    ELEMENTTYPE pot, cup0, cup1, cup2;
                    int i, pm;
                    size_t NN = diag.size() - 1;
                    
                    pm   = 0;
                    cup0 = diag[ NN ] + 0.5 * del_z * del_z * lambda
                         - 0.5 * del_z * std::sqrt( lambda );
                    pot  = diag[ NN - 1 ] + 0.5 * del_z * del_z * lambda;
                    cup1 = cup0 * pot;
                    if (cup0 * cup1 < 0.0) {
                        pm++;
                    }
                    cup0 = cup0 / std::fabs( cup1 );
                    cup1 = cup1 / std::fabs( cup1 );

                    for (i = NN - 2; i >= 0; i--) {
                        pot  = diag[ i ] + 0.5 * del_z * del_z * lambda;
                        cup2 = pot * cup1 - o_d_2 * cup0;
                        if (cup1 * cup2 < 0.0) {
                            pm++;
                        }
                        cup0 = cup1 / std::fabs( cup2 );
                        cup1 = cup2 / std::fabs( cup2 );
                    }

                    return pm;
                }
        };
    }  // namespace linear
}  // namespace NCPA
