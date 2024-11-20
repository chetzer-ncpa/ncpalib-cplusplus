#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/abstract_linear_system_solver.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/vector.hpp"
#include "NCPA/math.hpp"
#include "NCPA/types.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_GENERIC_TEMPLATE_NO_SUPERCLASS( Solver );
    }
}  // namespace NCPA

NCPA_LINEARALGEBRA_DECLARE_FRIEND_FUNCTIONS( NCPA::linear::Solver,
                                             ELEMENTTYPE );

NCPA_LINEARALGEBRA_DECLARE_FRIEND_BINARY_OPERATORS( NCPA::linear::Solver,
                                                    ELEMENTTYPE )

namespace NCPA {
    namespace linear {
        NCPA_LINEARALGEBRA_DECLARE_SPECIALIZED_TEMPLATE  //
            class Solver<ELEMENTTYPE, _ENABLE_IF_ELEMENTTYPE_IS_NUMERIC> {
            public:
                Solver() {}

                Solver( std::unique_ptr<
                        details::abstract_linear_system_solver<ELEMENTTYPE>>
                            ptr ) :
                    Solver<ELEMENTTYPE>() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                Solver( const Solver<ELEMENTTYPE>& other ) :
                    Solver<ELEMENTTYPE>() {
                    _ptr = std::move( other._ptr->clone() );
                }

                /**
                 * Move constructor.
                 * @param source The vector to assimilate.
                 */
                Solver( Solver<ELEMENTTYPE>&& source ) noexcept :
                    Solver<ELEMENTTYPE>() {
                    ::swap( *this, source );
                }

                ~Solver() {}

                /**
                 * Assignment operator.
                 * @param other The solver to assign to this.
                 */
                Solver<ELEMENTTYPE>& operator=( Solver<ELEMENTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                std::unique_ptr<Solver<ELEMENTTYPE>> clone() const {
                    return std::unique_ptr<Solver<ELEMENTTYPE>>(
                        new Solver<ELEMENTTYPE>( *this ) );
                }

                friend void ::swap<ELEMENTTYPE>(
                    Solver<ELEMENTTYPE>& a, Solver<ELEMENTTYPE>& b ) noexcept;

                explicit operator bool() const {
                    return ( _ptr ? true : false );
                }

                virtual Solver<ELEMENTTYPE>& set_system_matrix(
                    const Matrix<ELEMENTTYPE>& M ) {
                    check_pointer();
                    _ptr->set_system_matrix( M );
                    return *this;
                }

                virtual Solver<ELEMENTTYPE>& clear() {
                    if ( _ptr ) {
                        _ptr.reset();
                    }
                    return *this;
                }

                virtual Matrix<ELEMENTTYPE> solve(
                    const Matrix<ELEMENTTYPE>& RHS ) {
                    check_pointer();
                    return _ptr->solve( RHS );
                }

                virtual Matrix<ELEMENTTYPE> solve(
                    const Vector<ELEMENTTYPE>& RHS ) {
                    check_pointer();
                    return _ptr->solve( RHS );
                }

                virtual void check_pointer() const {
                    if ( !_ptr ) {
                        throw std::logic_error(
                            "Solver: Internal pointer has not been set!" );
                    }
                }

            private:
                std::unique_ptr<
                    details::abstract_linear_system_solver<ELEMENTTYPE>>
                    _ptr;
        };
    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap( NCPA::linear::Solver<T>& a,
                  NCPA::linear::Solver<T>& b ) noexcept {
    // using std::swap;
    a._ptr.swap( b._ptr );
}