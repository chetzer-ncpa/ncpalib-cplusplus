#pragma once

// #include "NCPA/interpolation/defines.hpp"
// #include "NCPA/interpolation/types.hpp"

namespace NCPA {
    namespace interpolation {
        namespace LANL {
            // template<typename INDEPTYPE, typename DEPTYPE, typename = void,
            //          typename = void, typename = void>
            template<typename FloatType>
            class _lanl_spline_1d;
            // template<typename INDEPTYPE, typename DEPTYPE, typename = void,
            //          typename = void, typename = void>
            template<typename FloatType>
            class LinearSpline1D;
            // template<typename INDEPTYPE, typename DEPTYPE, typename = void,
            //          typename = void, typename = void>
            template<typename FloatType>
            class NaturalCubicSpline1D;
            // template<typename INDEPTYPE, typename DEPTYPE, typename = void,
            //          typename = void, typename = void>
            template<typename FloatType>
            class BicubicSpline2D;
            // template<typename INDEPTYPE, typename DEPTYPE, typename = void,
                    //  typename = void, typename = void>
            template<typename FloatType>
            class NaturalSpline2D;
            // template<typename INDEPTYPE, typename DEPTYPE, typename = void,
            //          typename = void, typename = void>
            template<typename FloatType>
            class HybridSpline3D;
        }  // namespace LANL
    }  // namespace interpolation
}  // namespace NCPA
