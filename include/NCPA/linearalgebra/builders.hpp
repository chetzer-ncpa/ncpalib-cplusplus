#pragma once

#include "NCPA/linearalgebra/dense_matrix.hpp"
#include "NCPA/linearalgebra/dense_vector.hpp"
#include "NCPA/linearalgebra/matrix.hpp"
#include "NCPA/linearalgebra/sparse_matrix.hpp"
#include "NCPA/linearalgebra/sparse_vector.hpp"
#include "NCPA/linearalgebra/vector.hpp"

#include <stdexcept>

namespace NCPA {
    namespace linear {
        enum class family_t { INVALID, NCPA_DENSE, NCPA_SPARSE };

        template<typename ELEMENTTYPE>
        class VectorFactory {
            public:
                static Vector<ELEMENTTYPE> build( family_t family ) {
                    switch ( family ) {
                        case family_t::NCPA_DENSE:
                            return Vector<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_vector<ELEMENTTYPE>>(
                                    new details::dense_vector<
                                        ELEMENTTYPE>() ) );
                            break;
                        case family_t::NCPA_SPARSE:
                            return Vector<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_vector<ELEMENTTYPE>>(
                                    new details::sparse_vector<
                                        ELEMENTTYPE>() ) );
                            break;
                        default:
                            throw std::logic_error(
                                "Unknown or unsupported linear algebra family "
                                "requested" );
                    }
                }
        };

        template<typename ELEMENTTYPE>
        class MatrixFactory {
            public:
                static Matrix<ELEMENTTYPE> build( family_t family ) {
                    switch ( family ) {
                        case family_t::NCPA_DENSE:
                            return Matrix<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_matrix<ELEMENTTYPE>>(
                                    new details::dense_matrix<
                                        ELEMENTTYPE>() ) );
                            break;
                        case family_t::NCPA_SPARSE:
                            return Matrix<ELEMENTTYPE>(
                                std::unique_ptr<
                                    details::abstract_matrix<ELEMENTTYPE>>(
                                    new details::dense_matrix<
                                        ELEMENTTYPE>() ) );
                            break;
                        default:
                            throw std::logic_error(
                                "Unknown or unsupported linear algebra family "
                                "requested" );
                    }
                }
        };
    }  // namespace linear
}  // namespace NCPA
