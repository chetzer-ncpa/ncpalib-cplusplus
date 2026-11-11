#pragma once

#include "NCPA/exceptions.hpp"

namespace NCPA {
    namespace atmos {
        NCPA_DEFINE_EXCEPTION( WrongNumberOfDimensions, std::invalid_argument )
        NCPA_DEFINE_EXCEPTION( InvalidAxis, std::out_of_range )
    }
}