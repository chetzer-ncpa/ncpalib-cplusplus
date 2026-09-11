#pragma once

#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/units.hpp"

#include <unordered_map>
#include <vector>

namespace NCPA {
    namespace atmos {
        class indexed_axis {
            public:
                indexed_axis( size_t index, const vector_u_t& v ) :
                    _internal { std::make_pair( index, v ) } {}

                size_t index() const { return _internal.first; }

                const vector_u_t& values() const { return _internal.second; }

                vector_u_t& values() { return _internal.second; }

                size_t size() const { return this->values().size(); }

            private:
                std::pair<size_t, vector_u_t> _internal;
        };

        class indexed_axis_set {
            public:
                indexed_axis_set() {}

                indexed_axis_set( std::vector<dimension_t> dims ) {
                    for (size_t i = 0; i < dims.size(); ++i) {
                        _internal.insert( std::make_pair(
                            dims[ i ], indexed_axis { i, {} } ) );
                    }
                }

            private:
                std::unordered_map<dimension_t, indexed_axis> _internal;
        };
    }  // namespace atmos
}  // namespace NCPA
