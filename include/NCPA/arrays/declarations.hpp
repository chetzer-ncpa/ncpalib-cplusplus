#pragma once

#include <vector>
#include "NCPA/exceptions.hpp"

namespace NCPA {
    namespace arrays {
        
        template<typename>
        class SmartVector;

        template<typename>
        class SmartVectorView;

        template<size_t, typename>
        class SmartArray;

        template<size_t, typename>
        class SmartArrayView;

         template<size_t N, typename T>
        class ndvector;

        template<typename T>
        class vector2d_t;

        template<typename T>
        class vector3d_t;

        // namespace details {
        //     struct _ndstorage;
            
        // }
    }
}