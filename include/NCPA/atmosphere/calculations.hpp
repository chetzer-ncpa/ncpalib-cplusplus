#pragma once

#include "NCPA/atmosphere/constants.hpp"
#include "NCPA/atmosphere/types.hpp"
#include "NCPA/units.hpp"

#include <cmath>

namespace NCPA {
    namespace atmos {

        static scalar_t t2c( const scalar_t& t ) {
            return scalar_t( std::sqrt( t.get_as( "K" ) * _GAMMA_ * _R_ ),
                             "m/s" );
        }

        static scalar_t t2c( double t, const std::string& units ) {
            return t2c( scalar_t( t, units ) );
        }

        static scalar_t c2t( const scalar_t& c ) {
            double Cmps = c.get_as( "m/s" );
            return scalar_t( Cmps * Cmps / _GAMMA_ / _R_, "K" );
        }

        static scalar_t c2t( double t, const std::string& units ) {
            return c2t( scalar_t( t, units ) );
        }

        static scalar_t pd2c( const scalar_t& p, const scalar_t& d ) {
            return scalar_t(
                std::sqrt( _GAMMA_ * p.get_as( "Pa" ) / d.get_as( "kg/m3" ) ),
                "m/s" );
        }

        static scalar_t pd2c( double p, const std::string& p_units, double d,
                              const std::string& d_units ) {
            return pd2c( scalar_t( p, p_units ), scalar_t( d, d_units ) );
        }


    }  // namespace atmos
}  // namespace NCPA
