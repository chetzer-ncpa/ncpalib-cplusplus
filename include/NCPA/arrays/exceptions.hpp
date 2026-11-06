#pragma once

#include "NCPA/arrays/declarations.hpp"
#include "NCPA/exceptions.hpp"

#include <stdexcept>

namespace NCPA {
    namespace arrays {
        NCPA_DEFINE_EXCEPTION( DimensionOutOfRange, std::out_of_range )
        NCPA_DEFINE_EXCEPTION( InvalidRank, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( IndexOutOfRange, std::out_of_range )
        NCPA_DEFINE_EXCEPTION( ShapesDoNotMatch, std::invalid_argument )
    }
}  // namespace NCPA
