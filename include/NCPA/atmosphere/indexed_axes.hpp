#pragma once

#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/units.hpp"

#include <unordered_map>
#include <vector>

namespace NCPA {
    namespace atmos {
        class indexed_axis : public std::pair<size_t, vector_u_t> {
            public:
                indexed_axis( size_t index, const vector_u_t& v ) {
                    this->first  = index;
                    this->second = v;
                }

                size_t index() const { return this->first; }

                const vector_u_t& values() const { return this->second; }

                vector_u_t& values() { return this->second; }

                size_t size() const {
                    return this->values().size();
                }
        };

        class indexed_axis_set
            : public std::unordered_map<dimension_t, indexed_axis> {
            public:
                indexed_axis_set() {}

                indexed_axis_set( std::vector<dimension_t> dims ) {
                    for ( size_t i = 0; i < dims.size(); ++i ) {
                        this->insert( std::make_pair(
                            dims[ i ], indexed_axis { i, {} } ) );
                    }
                }
        };
    }  // namespace atmos
}  // namespace NCPA
