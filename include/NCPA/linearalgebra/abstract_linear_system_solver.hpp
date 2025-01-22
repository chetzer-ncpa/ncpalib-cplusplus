#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/linearalgebra/declarations.hpp"
#include "NCPA/linearalgebra/defines.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/vector.hpp"
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
    NCPA::linear::abstract_linear_system_solver, ELEMENTTYPE );

namespace NCPA {
    namespace linear {


        template<typename ELEMENTTYPE>
        class abstract_linear_system_solver {
            public:
                virtual ~abstract_linear_system_solver() {}

                friend void ::swap<ELEMENTTYPE>(
                    abstract_linear_system_solver<ELEMENTTYPE>& a,
                    abstract_linear_system_solver<ELEMENTTYPE>& b ) noexcept;


                virtual abstract_linear_system_solver<ELEMENTTYPE>&
                    set_system_matrix(
                        const NCPA::linear::Matrix<ELEMENTTYPE>& M )
                    = 0;
                virtual abstract_linear_system_solver<ELEMENTTYPE>& clear()
                    = 0;
                virtual NCPA::linear::Vector<ELEMENTTYPE> solve(
                    const NCPA::linear::Matrix<ELEMENTTYPE>& RHS )
                    = 0;
                virtual NCPA::linear::Vector<ELEMENTTYPE> solve(
                    const NCPA::linear::Vector<ELEMENTTYPE>& RHS )
                    = 0;
        };

    }  // namespace linear
}  // namespace NCPA

template<typename T>
static void swap(
    NCPA::linear::abstract_linear_system_solver<T>& a,
    NCPA::linear::abstract_linear_system_solver<T>& b ) noexcept {}
