#pragma once

#include "NCPA/interpolation/defines.hpp"

#include <memory>
#include <vector>

namespace NCPA {
    namespace interpolation {

        // template<typename T>
        // class vector2d_t : public std::vector<std::vector<T>> {
        //     public:
        //         vector2d_t() : std::vector<std::vector<T>>() {}

        //         vector2d_t( size_t nx1, size_t nx2, const T& val ) :
        //             std::vector<std::vector<T>>() {
        //             this->resize2d( nx1, nx2, val );
        //         }

        //         virtual void resize2d( size_t nx1, size_t nx2,
        //                                const T& val = (T)0 ) {
        //             this->resize( nx1, std::vector<T>( nx2, val ) );
        //         }

        //         virtual void size2d( size_t& nx1, size_t& nx2 ) const {
        //             nx1 = this->size();
        //             if ( nx1 > 0 ) {
        //                 nx2 = this->at( 0 ).size();
        //             } else {
        //                 nx2 = 0;
        //             }
        //         }
        // };

        // template<typename T>
        // class vector3d_t : public std::vector<vector2d_t<T>> {
        //     public:
        //         vector3d_t() : std::vector<vector2d_t<T>>() {}

        //         vector3d_t( size_t nx1, size_t nx2, size_t nx3,
        //                     const T& val ) :
        //             std::vector<vector2d_t<T>>() {
        //             this->resize3d( nx1, nx2, val );
        //         }

        //         virtual void resize3d( size_t nx1, size_t nx2, size_t nx3,
        //                                const T& val = (T)0 ) {
        //             this->resize( nx1, vector2d_t<T>( nx2, nx3, val ) );
        //         }

        //         virtual void size3d( size_t& nx1, size_t& nx2,
        //                              size_t& nx3 ) const {
        //             nx1 = this->size();
        //             if ( nx1 > 0 ) {
        //                 this->at( 0 ).size2d( nx2, nx3 );
        //             } else {
        //                 nx2 = 0;
        //                 nx3 = 0;
        //             }
        //         }
        // };


        enum class interpolator_1d_type_t {
            NEAREST_NEIGHBOR,
            LANL_LINEAR,
            LANL_CUBIC,
            GSL_LINEAR,
            GSL_POLYNOMIAL,
            GSL_CUBIC,
            GSL_AKIMA,
            GSL_STEFFEN,
            GSL_CUBIC_PERIODIC,
            GSL_AKIMA_PERIODIC
        };

        enum class interpolator_2d_type_t { LANL_NATURAL, LANL_BICUBIC };

        enum class interpolator_3d_type_t { LANL_HYBRID };

        enum class extrapolator_type_t {
            FORBIDDEN,
            CONSTANT,
            LINEAR,
            PERIODIC
        };

        template<typename INDEPTYPE, typename DEPTYPE>
        class _abstract_spline_1d;

        template<typename INDEPTYPE, typename DEPTYPE>
        class _abstract_spline_2d;

        template<typename INDEPTYPE, typename DEPTYPE>
        class _abstract_spline_3d;

        DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _spline_1d,
                                               _abstract_spline_1d );
        DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _spline_2d,
                                               _abstract_spline_2d );
        DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( _spline_3d,
                                               _abstract_spline_3d );

        DECLARE_GENERIC_INTERPOLATOR_TEMPLATE( nearest_neighbor_spline_1d,
                                               _spline_1d );

        template<typename INDEPTYPE, typename DEPTYPE>
        using spline_engine_1d_t
            = std::unique_ptr<_spline_1d<INDEPTYPE, DEPTYPE>>;

        template<typename INDEPTYPE, typename DEPTYPE>
        using spline_engine_2d_t
            = std::unique_ptr<_spline_2d<INDEPTYPE, DEPTYPE>>;

        template<typename INDEPTYPE, typename DEPTYPE>
        using spline_engine_3d_t
            = std::unique_ptr<_spline_3d<INDEPTYPE, DEPTYPE>>;

        template<typename INDEPTYPE, typename DEPTYPE>
        class Interpolator1D;

        template<typename INDEPTYPE, typename DEPTYPE>
        class Interpolator2D;

        template<typename INDEPTYPE, typename DEPTYPE>
        class Interpolator3D;
    }  // namespace interpolation
}  // namespace NCPA
