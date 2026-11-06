#pragma once

#include "NCPA/arrays/declarations.hpp"
#include "NCPA/arrays/exceptions.hpp"
#include <vector>

namespace NCPA {
    namespace arrays {
        namespace details {
            
            template<typename T>
            struct _ndstorage {
                    std::vector<T> data;
            };

            template<typename Shape>
            size_t total_size( const Shape& shape ) {
                size_t n = 1;
                for (size_t d : shape) {
                    n *= d;
                }
                return n;
            }

            template<typename Shape>
            void compute_strides( const Shape& shape, Shape& strides ) {
                strides                     = shape;
                strides[ shape.size() - 1 ] = 1;

                for (int i = static_cast<int>( shape.size() ) - 2; i >= 0;
                     --i) {
                    strides[ i ] = strides[ i + 1 ] * shape[ i + 1 ];
                }
            }

            template<typename ShapeA, typename ShapeB>
            size_t linear_index( const ShapeA& indices, const ShapeB& shape,
                                 const ShapeB& strides, size_t offset = 0 ) {
                size_t result = offset;

                if (indices.size() != shape.size()) {
                    throw InvalidRank( "wrong number of indices");
                }

                for (size_t i = 0; i < shape.size(); ++i) {
                    if (indices[ i ] >= shape[ i ]) {
                        throw IndexOutOfRange();
                    }

                    result += indices[ i ] * strides[ i ];
                }

                return result;
            }

            template<typename Shape>
            size_t storage_index( size_t linear_pos, size_t offset,
                                  const Shape& shape, const Shape& strides ) {
                size_t idx = offset;
                size_t rem = linear_pos;

                for (int d = static_cast<int>( shape.size() ) - 1; d >= 0;
                     --d) {
                    size_t coord  = rem % shape[ d ];
                    rem          /= shape[ d ];

                    idx += coord * strides[ d ];
                }

                return idx;
            }

            template<typename ShapeIn, typename ShapeOut>
            void remove_dimension( const ShapeIn& shape,
                                   const ShapeIn& strides, size_t dim,
                                   ShapeOut& new_shape,
                                   ShapeOut& new_strides ) {
                size_t out = 0;

                for (size_t i = 0; i < shape.size(); ++i) {
                    if (i == dim) continue;

                    new_shape[ out ]   = shape[ i ];
                    new_strides[ out ] = strides[ i ];
                    ++out;
                }
            }

            template<typename ContainerA, typename ContainerB>
            void deal_into( const ContainerA& a, ContainerB& b ) {
                if (a.size() == b.size()) {
                    auto ait = a.begin();
                    auto bit = b.begin();
                    while (ait != a.end() && bit != b.end()) {
                        *bit = *ait;
                        ++ait;
                        ++bit;
                    }
                } else {
                    throw ShapesDoNotMatch(
                        "Cannot deal N-D containers of different sizes" );
                }
            }

            template<typename ContainerA, typename ContainerB>
            void copy_to( const ContainerA& a, ContainerB& b ) {
                if (a.rank() == b.rank() || b.rank_mutable()) {
                    b.resize( a.shape() );
                    deal_into( a, b );
                } else {
                    throw InvalidRank(
                        "Cannot copy N-D containers of different ranks" );
                }
            }

            template<typename ContainerA, typename ContainerB>
            bool same_shape( const ContainerA& a, const ContainerB& b ) {
                if (a.rank() != b.rank()) return false;
                for (auto i = 0; i < a.rank(); ++i) {
                    if (a.dim(i) != b.dim(i)) return false;
                }
                return true;
            }

            template<typename Container>
            std::vector<typename Container::shape_type> indices(const Container& c) {
                auto shape = c.shape();
                auto digits = shape;
                digits.clear();
                std::vector<typename Container::shape_type> out;
                out.reserve( c.size() );
                do {
                    out.push_back( digits );
                }
                while (increment_digits( digits, shape ));
                return out;
            }
        }
    }
}