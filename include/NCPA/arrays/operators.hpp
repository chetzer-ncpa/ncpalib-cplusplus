#pragma once

#include "NCPA/arrays/declarations.hpp"

#include <vector>

namespace NCPA {
    namespace arrays {
        template<typename T>
        bool operator==( const std::vector<T>& a, const std::vector<T>& b ) {
            if (a.size() != b.size()) return false;
            auto ait = a.begin();
            auto bit = b.begin();
            while (ait != a.end() && bit != b.end()) {
                if (*ait != *bit) return false;
                ++ait;
                ++bit;
            }
            return true;
        }

        template<typename T>
        bool operator!=( const std::vector<T>& a, const std::vector<T>& b ) {
            return !( a == b );
        }
    }  // namespace arrays
}  // namespace NCPA

/**
 * @brief Operator overload for vector-vector addition.
 */
template<typename T>
std::vector<T> operator+( const std::vector<T>& a, const std::vector<T>& b ) {
    return NCPA::arrays::add_vectors<std::vector<T>>( a, b );
}

/**
 * @brief Operator overload for vector-scalar addition.
 */
template<typename T>
std::vector<T> operator+( const std::vector<T>& a, const T& b ) {
    return NCPA::arrays::offset_vector<std::vector<T>>( a, b );
}

/**
 * @brief Operator overload for vector-vector subtraction.
 */
template<typename T>
std::vector<T> operator-( const std::vector<T>& a, const std::vector<T>& b ) {
    return NCPA::arrays::subtract_vectors<std::vector<T>>( a, b );
}

/**
 * @brief Operator overload for vector-scalar subtraction.
 */
template<typename T>
std::vector<T> operator-( const std::vector<T>& a, const T& b ) {
    return NCPA::arrays::offset_vector<std::vector<T>>( a, -b );
}

/**
 * @brief Operator overload for element-wise vector-vector multiplication.
 */
template<typename T>
std::vector<T> operator*( const std::vector<T>& a, const std::vector<T>& b ) {
    return NCPA::arrays::multiply_vectors<std::vector<T>>( a, b );
}

/**
 * @brief Operator overload for vector scaling.
 */
template<typename T>
std::vector<T> operator*( const std::vector<T>& a, const T& b ) {
    return NCPA::arrays::scale_vector<std::vector<T>>( a, b );
}

/**
 * @brief Operator overload for vector-scalar division.
 */
template<typename T>
std::vector<T> operator/( const std::vector<T>& a, const T& b ) {
    return NCPA::arrays::scale_vector<std::vector<T>>( a,
                                                       std::pow( b, -1.0 ) );
}

/**
 * @brief Operator overload for unary vector negation.
 */
template<typename T>
std::vector<T> operator-( const std::vector<T>& a ) {
    return NCPA::arrays::scale_vector<std::vector<T>>(
        a, -NCPA::constants::one<T>() );
}

/**
 * @brief Stream operator for std::vector.
 */
template<typename T>
std::ostream& operator<<( std::ostream& os, const std::vector<T>& v ) {
    os << "{ ";
    for (auto it = v.cbegin(); it != v.cend(); ++it) {
        if (it != v.cbegin()) {
            os << ", ";
        }
        os << *it;
    }
    os << " }";
    return os;
}
