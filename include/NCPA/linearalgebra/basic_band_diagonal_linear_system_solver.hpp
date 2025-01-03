#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_linear_system_solver.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/lu.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/vector.hpp"
#include "NCPA/logging.hpp"
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

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS(
    NCPA::linear::details::basic_band_diagonal_linear_system_solver,
    ELEMENTTYPE );

namespace NCPA {
    namespace linear {
        namespace details {

            NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
                class basic_band_diagonal_linear_system_solver<
                    ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC>
                : public abstract_linear_system_solver<ELEMENTTYPE> {
                public:
                    basic_band_diagonal_linear_system_solver() :
                        abstract_linear_system_solver<ELEMENTTYPE>() {}

                    // copy constructor
                    basic_band_diagonal_linear_system_solver(
                        const basic_band_diagonal_linear_system_solver<
                            ELEMENTTYPE>& other ) :
                        abstract_linear_system_solver<ELEMENTTYPE>() {
                        _mat  = other._mat;
                        _au   = other._au;
                        _al   = other._al;
                        n     = other.n;
                        m1    = other.m1;
                        m2    = other.m2;
                        d     = other.d;
                        _indx = other._indx;
                    }

                    /**
                     * Move constructor.
                     * @param source The vector to assimilate.
                     */
                    basic_band_diagonal_linear_system_solver(
                        basic_band_diagonal_linear_system_solver<ELEMENTTYPE>&&
                            source ) noexcept :
                        abstract_linear_system_solver<ELEMENTTYPE>() {
                        ::swap( *this, source );
                    }

                    virtual ~basic_band_diagonal_linear_system_solver() {
                        this->clear();
                    }

                    friend void ::swap<ELEMENTTYPE>(
                        basic_band_diagonal_linear_system_solver<ELEMENTTYPE>&
                            a,
                        basic_band_diagonal_linear_system_solver<ELEMENTTYPE>&
                            b ) noexcept;

                    /**
                     * Assignment operator.
                     * @param other The vector to assign to this.
                     */
                    basic_band_diagonal_linear_system_solver<ELEMENTTYPE>&
                        operator=( basic_band_diagonal_linear_system_solver<
                                   ELEMENTTYPE>
                                       other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual abstract_linear_system_solver<ELEMENTTYPE>& clear()
                        override {
                        _mat.clear();
                        if ( _al != nullptr ) {
                            NCPA::arrays::free_array( _al, n, m1 );
                            _al = nullptr;
                        }
                        if ( _au != nullptr ) {
                            NCPA::arrays::free_array( _au, n, m1 + m2 + 1 );
                            _au = nullptr;
                        }
                        if ( _indx != nullptr ) {
                            NCPA::arrays::free_array( _indx, n );
                            _indx = nullptr;
                        }
                        n  = 0;
                        m1 = 0;
                        m2 = 0;
                        d  = 0.0;
                        return *static_cast<
                            abstract_linear_system_solver<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_linear_system_solver<ELEMENTTYPE>&
                        set_system_matrix(
                            const Matrix<ELEMENTTYPE>& M ) override {
                        if ( !M.is_band_diagonal() ) {
                            throw std::logic_error(
                                "System matrix must be band-diagonal" );
                        }
                        if ( !M.is_square() ) {
                            throw std::logic_error(
                                "System matrix must be square!" );
                        }
                        this->clear();
                        _mat.copy( M.internal() );
                        _decompose();
                        return *static_cast<
                            abstract_linear_system_solver<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual abstract_linear_system_solver<ELEMENTTYPE>&
                        _decompose( bool pivot = true ) {
                        int i, j, k, l, mm;
                        ELEMENTTYPE dum;
                        _pivot = pivot;

                        n  = (int)_mat.rows();
                        m1 = (int)_mat._n_lower;
                        m2 = (int)_mat._n_upper;
                        mm = m1 + m2 + 1;

                        _indx = NCPA::arrays::zeros<int>( n );

                        std::vector<std::vector<ELEMENTTYPE>> a
                            = _mat._contents;

                        // NCPA_DEBUG << "Before decomposition a = (" <<
                        // a.size()
                        //            << "," << a[ 0 ].size()
                        //            << "):" << std::endl;
                        // for ( i = 0; i < mm; i++ ) {
                        //     NCPA_DEBUG << "[ ";
                        //     for ( j = 0; j < n; j++ ) {
                        //         if ( j != 0 ) {
                        //             NCPA_DEBUG << ", ";
                        //         }
                        //         NCPA_DEBUG << a[ i ][ j ];
                        //     }
                        //     NCPA_DEBUG << "]" << std::endl;
                        // }

                        _au = NCPA::arrays::zeros<ELEMENTTYPE>( n, mm );
                        if ( m1 > 0 ) {
                            _al = NCPA::arrays::zeros<ELEMENTTYPE>( n, m1 );
                        }
                        for ( i = 0; i < mm; i++ ) {
                            for ( j = 0; j < n; j++ ) {
                                // reverse storage indices to
                                // match Numerical Recipes
                                _au[ j ][ i ] = a[ i ][ j ];
                            }
                        }

                        // NCPA_DEBUG << "Before decomposition _au ="
                        //            << std::endl;
                        // for ( i = 0; i < n; i++ ) {
                        //     NCPA_DEBUG << "[ ";
                        //     for ( j = 0; j < mm; j++ ) {
                        //         if ( j != 0 ) {
                        //             NCPA_DEBUG << ", ";
                        //         }
                        //         NCPA_DEBUG << _au[ i ][ j ];
                        //     }
                        //     NCPA_DEBUG << "]" << std::endl;
                        // }

                        // _al.clear();
                        // for ( i = 0; i < m1; i++ ) {
                        //     _al.push_back( std::vector<ELEMENTTYPE>( n ) );
                        // }

                        // if ( m1 > 0 ) {
                        //     NCPA_DEBUG << "Before decomposition _al ="
                        //                << std::endl;
                        //     for ( i = 0; i < n; i++ ) {
                        //         NCPA_DEBUG << "[ ";
                        //         for ( j = 0; j < m1; j++ ) {
                        //             if ( j != 0 ) {
                        //                 NCPA_DEBUG << ", ";
                        //             }
                        //             NCPA_DEBUG << _al[ i ][ j ];
                        //         }
                        //         NCPA_DEBUG << "]" << std::endl;
                        //     }
                        // } else {
                        //     NCPA_DEBUG << "_al is empty" << std::endl;
                        // }

                        l = m1;
                        for ( i = 0; i < m1; i++ ) {
                            for ( j = m1 - i; j < mm; j++ ) {
                                _au[ i ][ j - 1 ] = _au[ i ][ j ];
                            }
                            l--;
                            for ( j = mm - l - 1; j < mm; j++ ) {
                                _au[ i ][ j ] = 0.0;
                            }
                        }

                        // NCPA_DEBUG << "After rearranging _au =" <<
                        // std::endl; for ( i = 0; i < n; i++ ) {
                        //     NCPA_DEBUG << "[ ";
                        //     for ( j = 0; j < mm; j++ ) {
                        //         if ( j != 0 ) {
                        //             NCPA_DEBUG << ", ";
                        //         }
                        //         NCPA_DEBUG << _au[ i ][ j ];
                        //     }
                        //     NCPA_DEBUG << "]" << std::endl;
                        // }


                        d = 1.0;
                        l = m1;
                        for ( k = 0; k < n; k++ ) {
                            dum = _au[ k ][ 0 ];
                            i   = k;
                            if ( l < n ) {
                                l++;
                            }
                            // if ( _pivot ) {
                                for ( j = k + 1; j < l; j++ ) {
                                    if ( std::abs( _au[ j ][ 0 ] )
                                         > std::abs( dum ) ) {
                                        dum = _au[ j ][ 0 ];
                                        i   = j;
                                    }
                                }
                                // _indx[ k ] = i + 1;
                                if ( dum == 0.0 ) {
                                    // matrix is algorithmically singular but
                                    // keep going with tiny pivot
                                    _au[ k ][ 0 ] = 1e-40;
                                    std::cout << "Singular!" << std::endl;
                                }
                            if ( _pivot ) {
                                _indx[ k ] = i + 1;
                                if ( i != k ) {
                                    // interchange rows
                                    d = -d;
                                    for ( j = 0; j < mm; j++ ) {
                                        auto tmp      = _au[ k ][ j ];
                                        _au[ k ][ j ] = _au[ i ][ j ];
                                        _au[ i ][ j ] = tmp;
                                    }
                                }
                            } else {
                                _indx[k] = k+1;
                            }

                            for ( i = k + 1; i < l; i++ ) {
                                dum = _au[ i ][ 0 ] / _au[ k ][ 0 ];
                                _al[ k ][ i - k - 1 ] = dum;
                                for ( j = 1; j < mm; j++ ) {
                                    _au[ i ][ j - 1 ]
                                        = _au[ i ][ j ] - dum * _au[ k ][ j ];
                                }
                                _au[ i ][ mm - 1 ] = _zero;
                            }
                        }

                        NCPA_DEBUG << "After decomposition _au = :"
                                   << std::endl;
                        for ( i = 0; i < n; i++ ) {
                            NCPA_DEBUG << "[ ";
                            for ( j = 0; j < mm; j++ ) {
                                if ( j != 0 ) {
                                    NCPA_DEBUG << ", ";
                                }
                                NCPA_DEBUG << _au[ i ][ j ];
                            }
                            NCPA_DEBUG << "]" << std::endl;
                        }

                        NCPA_DEBUG << "After decomposition _al = :"
                                   << std::endl;
                        for ( i = 0; i < n; i++ ) {
                            NCPA_DEBUG << "[ ";
                            for ( j = 0; j < m1; j++ ) {
                                if ( j != 0 ) {
                                    NCPA_DEBUG << ", ";
                                }
                                NCPA_DEBUG << _al[ i ][ j ];
                            }
                            NCPA_DEBUG << "]" << std::endl;
                        }

                        NCPA_DEBUG
                            << "After decomposition, _indx = " << std::endl
                            << "[ ";
                        for ( i = 0; i < n; i++ ) {
                            if ( i != 0 ) {
                                NCPA_DEBUG << ", ";
                            }
                            NCPA_DEBUG << _indx[ i ];
                        }
                        NCPA_DEBUG << " ]" << std::endl;

                        return *static_cast<
                            abstract_linear_system_solver<ELEMENTTYPE> *>(
                            this );
                    }

                    virtual NCPA::linear::Vector<ELEMENTTYPE> _solve_using_lu(
                        const NCPA::linear::Vector<ELEMENTTYPE>& b ) {
                        // use nomenclature from Numerical Recipes except
                        // indexing of _al and _au is reversed.
                        int i, j, k, l, mm;
                        ELEMENTTYPE dum;
                        mm = m1 + m2 + 1;
                        l  = m1;

                        size_t N = b.size();
                        std::vector<ELEMENTTYPE> x( N, _zero );

                        for ( k = 0; k < n; k++ ) {
                            x[ k ] = b[ k ];
                        }
                        for ( k = 0; k < n; k++ ) {
                            // if ( _pivot ) {
                                j = _indx[ k ] - 1;
                                if ( j != k ) {
                                    std::swap( x[ k ], x[ j ] );
                                }
                            // }
                            if ( l < n ) {
                                l++;
                            }
                            for ( j = k + 1; j < l; j++ ) {
                                x[ j ] -= _al[ k ][ j - k - 1 ] * x[ k ];
                            }
                        }

                        l = 1;
                        for ( i = n - 1; i >= 0; i-- ) {
                            dum = x[ i ];
                            for ( k = 1; k < l; k++ ) {
                                dum -= _au[ i ][ k ] * x[ k + i ];
                            }
                            x[ i ] = dum / _au[ i ][ 0 ];
                            if ( l < mm ) {
                                l++;
                            }
                        }

                        Vector<ELEMENTTYPE> xvec
                            = VectorFactory<ELEMENTTYPE>::build(
                                family_t::NCPA_DENSE );
                        xvec.resize( x.size() );
                        xvec.set( x );
                        return xvec;
                    }

                    virtual NCPA::linear::Vector<ELEMENTTYPE> solve(
                        const NCPA::linear::Vector<ELEMENTTYPE>& b ) override {
                        size_t N = b.size();
                        if ( N != _mat.rows() ) {
                            std::ostringstream oss;
                            oss << "solver: size mismatch between system "
                                   "matrix "
                                   "size ["
                                << _mat.rows() << "x" << _mat.columns()
                                << "] and input vector size " << b.size();
                            throw std::logic_error( oss.str() );
                        }
                        if ( !_lu_ready() ) {
                            _decompose();
                        }

                        return _solve_using_lu( b );
                    }

                    virtual NCPA::linear::Vector<ELEMENTTYPE> solve(
                        const NCPA::linear::Matrix<ELEMENTTYPE>& b ) override {
                        if ( b.is_column_matrix() ) {
                            return solve( *b.get_column( 0 ) );
                        } else if ( b.is_row_matrix() ) {
                            return solve( *b.get_row( 0 ) );
                        } else {
                            throw std::logic_error(
                                "solve(): input Matrix is neither a column "
                                "matrix nor a row matrix!" );
                        }
                    }

                    // protected:
                    //     void _build_lu() {
                    //         _lu = std::unique_ptr<
                    //             NCPA::linear::LUDecomposition<ELEMENTTYPE>>(
                    //             new
                    //             NCPA::linear::BandDiagonalLUDecomposition<
                    //                 ELEMENTTYPE>() );
                    //     }

                private:
                    bool _lu_ready() const {
                        return ( _au != nullptr && _al != nullptr
                                 && _indx != nullptr );
                    }

                    NCPA::linear::band_diagonal_matrix<ELEMENTTYPE> _mat;
                    // std::unique_ptr<
                    //     NCPA::linear::BandDiagonalLUDecomposition<ELEMENTTYPE>>
                    //     _lu;
                    int n, m1, m2;
                    ELEMENTTYPE d;
                    int *_indx        = nullptr;
                    ELEMENTTYPE **_au = nullptr, **_al = nullptr;

                    // std::vector<int> _indx;
                    // std::vector<std::vector<ELEMENTTYPE>> _au, _al;
                    const ELEMENTTYPE _zero = NCPA::math::zero<ELEMENTTYPE>();
                    bool _pivot             = false;
            };
        }  // namespace details
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap(
    NCPA::linear::details::basic_band_diagonal_linear_system_solver<T>& a,
    NCPA::linear::details::basic_band_diagonal_linear_system_solver<T>&
        b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::linear::details::abstract_linear_system_solver<T>&>(
            a ),
        static_cast<NCPA::linear::details::abstract_linear_system_solver<T>&>(
            b ) );
    swap( a._mat, b._mat );
    swap( a._lu, b._lu );
}
