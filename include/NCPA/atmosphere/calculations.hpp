#pragma once

#include "NCPA/atmosphere/constants.hpp"
#include "NCPA/atmosphere/types.hpp"
#include "NCPA/units.hpp"

#include <cmath>

namespace NCPA {
    namespace atmos {

        template<typename T>
        scalar_t<T> t2c(
            const scalar_t<T>& t ) {
            return scalar_t<T>(
                std::sqrt( t.get_as( "K" ) * _GAMMA_ * _R_ ), "m/s" );
        }

        template<typename T>
        scalar_t<T> t2c( T t, const std::string& units ) {
            return t2c<T>( scalar_t<T>( t, units ) );
        }

        template<typename T>
        scalar_t<T> c2t(
            const scalar_t<T>& c ) {
            T Cmps = c.get_as( "m/s" );
            return scalar_t<T>(
                Cmps * Cmps / _GAMMA_ / _R_, "K" );
        }

        template<typename T>
        scalar_t<T> c2t( T t, const std::string& units ) {
            return c2t<T>( scalar_t<T>( t, units ) );
        }

        template<typename T>
        scalar_t<T> pd2c(
            const scalar_t<T>& p,
            const scalar_t<T>& d ) {
            return scalar_t<T>(
                std::sqrt( _GAMMA_ * p.get_as( "Pa" ) / d.get_as( "kg/m3" ) ),
                "m/s" );
        }

        template<typename T>
        scalar_t<T> pd2c( T p, const std::string& p_units,
        T d, const std::string& d_units ) {
            return pd2c<T>( scalar_t<T> )
        }

        

    }  // namespace atmos
}  // namespace NCPA
