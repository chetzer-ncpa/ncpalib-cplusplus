#pragma once

// #include "NCPA/atmosphere/constants.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/ndvector.hpp"
#include "NCPA/units.hpp"

#include <cmath>

namespace NCPA {
    namespace atmos {

        // temperature to sound speed
        static scalar_u_t t2c( const scalar_u_t& t ) {
            return scalar_u_t( std::sqrt( t.get_as( "K" ) * constants::GAMMA()
                                          * constants::R() ),
                               "m/s" );
        }

        static scalar_u_t t2c( double t, const std::string& units ) {
            return t2c( scalar_u_t( t, units ) );
        }

        static vector_u_t t2c( const vector_u_t& t ) {
            vector_u_t c( t.size(), NCPA::units::Units::from_string( "m/s" ) );
            for (size_t i = 0; i < t.size(); i++) {
                c.set( i, t2c( t.get_scalar( i ) ) );
            }
            return c;
        }

        static vector_u_t t2c( const std::vector<double>& t,
                               const std::string& units ) {
            return t2c( vector_u_t( t, units ) );
        }

        std::vector<double> t2c( const NCPA::arrays::ndvector<1, double>& t,
                                 const units_ptr_t units ) {
            std::vector<double> c( t.size() );
            for (size_t i = 0; i < t.size(); ++i) {
                c[ i ]
                    = t2c( scalar_u_t( *static_cast<const double *>( &t[ i ] ),
                                       units ) )
                          .get();
            }
            return c;
        }

        double t2c( const NCPA::arrays::ndvector<0, double>& t,
                    const units_ptr_t units ) {
            return t2c( scalar_u_t( t, units ) ).get();
        }

        template<size_t N>
        NCPA::arrays::ndvector<N, double> t2c(
            const NCPA::arrays::ndvector<N, double>& t,
            const units_ptr_t units ) {
            NCPA::arrays::ndvector<N, double> c = t;

            for (size_t i = 0; i < t.size(); ++i) {
                c[ i ] = t2c( NCPA::arrays::ndvector<N - 1, double>( t[ i ] ),
                              units );
            }
            return c;
        }

        // sound speed to temperature
        static scalar_u_t c2t( const scalar_u_t& c ) {
            double Cmps = c.get_as( "m/s" );
            return scalar_u_t(
                Cmps * Cmps / constants::GAMMA() / constants::R(), "K" );
        }

        static scalar_u_t c2t( double t, const std::string& units ) {
            return c2t( scalar_u_t( t, units ) );
        }

        static vector_u_t c2t( const vector_u_t& c ) {
            vector_u_t t( c.size(), "K" );
            for (size_t i = 0; i < c.size(); i++) {
                t.set( i, c2t( c.get_scalar( i ) ) );
                // t[ i ] = c2t( c.get_scalar( i ) ).get();
                // double Cmps = c.get_as( i, "m/s" );
                // t[ i ]      = std::sqrt( Cmps * Cmps / constants::GAMMA() /
                // constants::R() );
            }
            return t;
        }

        static vector_u_t c2t( const std::vector<double>& c,
                               const std::string& units ) {
            return c2t( vector_u_t( c, units ) );
        }

        std::vector<double> c2t( const NCPA::arrays::ndvector<1, double>& c,
                                 const units_ptr_t units ) {
            std::vector<double> t( c.size() );
            for (size_t i = 0; i < c.size(); ++i) {
                t[ i ]
                    = c2t( scalar_u_t( *static_cast<const double *>( &c[ i ] ),
                                       units ) )
                          .get();
            }
            return t;
        }

        double c2t( const NCPA::arrays::ndvector<0, double>& c,
                    const units_ptr_t units ) {
            return c2t( scalar_u_t( c, units ) ).get();
        }

        template<size_t N>
        NCPA::arrays::ndvector<N, double> c2t(
            const NCPA::arrays::ndvector<N, double>& c,
            const units_ptr_t units ) {
            NCPA::arrays::ndvector<N, double> t = c;

            for (size_t i = 0; i < c.size(); ++i) {
                t[ i ] = c2t( NCPA::arrays::ndvector<N - 1, double>( c[ i ] ),
                              units );
            }
            return t;
        }

        // pressure and density to sound speed
        static scalar_u_t pd2c( const scalar_u_t& p, const scalar_u_t& d ) {
            return scalar_u_t( std::sqrt( constants::GAMMA() * p.get_as( "Pa" )
                                          / d.get_as( "kg/m3" ) ),
                               "m/s" );
        }

        static scalar_u_t pd2c( double p, const std::string& p_units, double d,
                                const std::string& d_units ) {
            return pd2c( scalar_u_t( p, p_units ), scalar_u_t( d, d_units ) );
        }

        static vector_u_t pd2c( const vector_u_t& p, const vector_u_t& d ) {
            vector_u_t c( p.size(), "m/s" );
            for (size_t i = 0; i < p.size(); i++) {
                c.set( i, pd2c( p.get_scalar( i ), d.get_scalar( i ) ) );
            }
            return c;
        }

        static vector_u_t pd2c( const std::vector<double>& p,
                                const std::string& p_units,
                                const std::vector<double>& d,
                                const std::string& d_units ) {
            return pd2c( vector_u_t( p, p_units ), vector_u_t( d, d_units ) );
        }

        std::vector<double> pd2c( const NCPA::arrays::ndvector<1, double>& p,
                                  const units_ptr_t p_units,
                                  const NCPA::arrays::ndvector<1, double>& d,
                                  const units_ptr_t d_units ) {
            std::vector<double> c( p.size() );
            for (size_t i = 0; i < p.size(); ++i) {
                c[ i ]
                    = pd2c(
                          scalar_u_t( *static_cast<const double *>( &p[ i ] ),
                                      p_units ),
                          scalar_u_t( *static_cast<const double *>( &d[ i ] ),
                                      d_units ) )
                          .get();
            }
            return c;
        }

        double pd2c( const NCPA::arrays::ndvector<0, double>& p,
                     const units_ptr_t p_units,
                     const NCPA::arrays::ndvector<0, double>& d,
                     const units_ptr_t d_units ) {
            return pd2c( scalar_u_t( p, p_units ), scalar_u_t( d, d_units ) )
                .get();
        }

        template<size_t N>
        NCPA::arrays::ndvector<N, double> pd2c(
            const NCPA::arrays::ndvector<N, double>& p,
            const units_ptr_t p_units,
            const NCPA::arrays::ndvector<N, double>& d,
            const units_ptr_t d_units ) {
            NCPA::arrays::ndvector<N, double> c = p;

            for (size_t i = 0; i < p.size(); ++i) {
                c[ i ] = pd2c(
                    NCPA::arrays::ndvector<N - 1, double>( p[ i ] ), p_units,
                    NCPA::arrays::ndvector<N - 1, double>( d[ i ] ), d_units );
            }
            return c;
        }

        // U & V winds to wind speed
        static scalar_u_t uv2ws( const scalar_u_t& u, const scalar_u_t& v ) {
            scalar_u_t ws = u;
            ws            = std::sqrt( ( u.get() * u.get() )
                                       + ( v.get_as( *u.get_units() )
                                * v.get_as( *u.get_units() ) ) );
            return ws;
        }

        static scalar_u_t uv2ws( double u, const std::string& u_units,
                                 double v, const std::string& v_units ) {
            return uv2ws( scalar_u_t( u, u_units ), scalar_u_t( v, v_units ) );
        }

        static vector_u_t uv2ws( const vector_u_t& u, const vector_u_t& v ) {
            vector_u_t ws( u.size(), u.get_units() );
            for (size_t i = 0; i < u.size(); i++) {
                ws.set( i, uv2ws( u.get_scalar( i ), v.get_scalar( i ) ) );
                // ws[ i ] = uv2ws( u.get_scalar(i), v.get_scalar(i) ).get();
            }
            return ws;
        }

        static vector_u_t uv2ws( const std::vector<double>& u,
                                 const std::string& u_units,
                                 const std::vector<double>& v,
                                 const std::string& v_units ) {
            return uv2ws( vector_u_t( u, u_units ), vector_u_t( v, v_units ) );
        }

        std::vector<double> uv2ws( const NCPA::arrays::ndvector<1, double>& u,
                                   const units_ptr_t u_units,
                                   const NCPA::arrays::ndvector<1, double>& v,
                                   const units_ptr_t v_units ) {
            std::vector<double> ws( u.size() );
            for (size_t i = 0; i < u.size(); ++i) {
                ws[ i ]
                    = uv2ws(
                          scalar_u_t( *static_cast<const double *>( &u[ i ] ),
                                      u_units ),
                          scalar_u_t( *static_cast<const double *>( &v[ i ] ),
                                      v_units ) )
                          .get();
            }
            return ws;
        }

        double uv2ws( const NCPA::arrays::ndvector<0, double>& u,
                      const units_ptr_t u_units,
                      const NCPA::arrays::ndvector<0, double>& v,
                      const units_ptr_t v_units ) {
            return uv2ws( scalar_u_t( u, u_units ), scalar_u_t( v, v_units ) )
                .get();
        }

        template<size_t N>
        NCPA::arrays::ndvector<N, double> uv2ws(
            const NCPA::arrays::ndvector<N, double>& u,
            const units_ptr_t u_units,
            const NCPA::arrays::ndvector<N, double>& v,
            const units_ptr_t v_units ) {
            NCPA::arrays::ndvector<N, double> ws = u;

            for (size_t i = 0; i < u.size(); ++i) {
                ws[ i ] = uv2ws(
                    NCPA::arrays::ndvector<N - 1, double>( u[ i ] ), u_units,
                    NCPA::arrays::ndvector<N - 1, double>( v[ i ] ), v_units );
            }
            return ws;
        }

        // U & V winds to wind direction
        static scalar_u_t uv2wd( const scalar_u_t& u, const scalar_u_t& v ) {
            return scalar_u_t(
                NCPA::math::math2az( NCPA::math::rad2deg(
                    // NCPA::math::PI / 2.0 -
                    std::atan2( v.get_as( *u.get_units() ), u.get() ) ) ),
                "deg" );
        }

        static scalar_u_t uv2wd( double u, const std::string& u_units,
                                 double v, const std::string& v_units ) {
            return uv2wd( scalar_u_t( u, u_units ), scalar_u_t( v, v_units ) );
        }

        static vector_u_t uv2wd( const vector_u_t& u, const vector_u_t& v ) {
            vector_u_t wd( u.size(), "deg" );
            for (size_t i = 0; i < u.size(); i++) {
                wd.set( i, uv2wd( u.get_scalar( i ), v.get_scalar( i ) ) );
            }
            return wd;
        }

        static vector_u_t uv2wd( const std::vector<double>& u,
                                 const std::string& u_units,
                                 const std::vector<double>& v,
                                 const std::string& v_units ) {
            return uv2wd( vector_u_t( u, u_units ), vector_u_t( v, v_units ) );
        }

        std::vector<double> uv2wd( const NCPA::arrays::ndvector<1, double>& u,
                                   const units_ptr_t u_units,
                                   const NCPA::arrays::ndvector<1, double>& v,
                                   const units_ptr_t v_units ) {
            std::vector<double> wd( u.size() );
            for (size_t i = 0; i < u.size(); ++i) {
                wd[ i ]
                    = uv2wd(
                          scalar_u_t( *static_cast<const double *>( &u[ i ] ),
                                      u_units ),
                          scalar_u_t( *static_cast<const double *>( &v[ i ] ),
                                      v_units ) )
                          .get();
            }
            return wd;
        }

        double uv2wd( const NCPA::arrays::ndvector<0, double>& u,
                      const units_ptr_t u_units,
                      const NCPA::arrays::ndvector<0, double>& v,
                      const units_ptr_t v_units ) {
            return uv2wd( scalar_u_t( u, u_units ), scalar_u_t( v, v_units ) )
                .get();
        }

        template<size_t N>
        NCPA::arrays::ndvector<N, double> uv2wd(
            const NCPA::arrays::ndvector<N, double>& u,
            const units_ptr_t u_units,
            const NCPA::arrays::ndvector<N, double>& v,
            const units_ptr_t v_units ) {
            NCPA::arrays::ndvector<N, double> wd = u;

            for (size_t i = 0; i < u.size(); ++i) {
                wd[ i ] = uv2wd(
                    NCPA::arrays::ndvector<N - 1, double>( u[ i ] ), u_units,
                    NCPA::arrays::ndvector<N - 1, double>( v[ i ] ), v_units );
            }
            return wd;
        }

        // wind speed and direction to wind component along an azimuth
        static scalar_u_t w2wc( const scalar_u_t& ws, const scalar_u_t& wd,
                                double az ) {
            return scalar_u_t( ws.get()
                                   * std::cos( wd.get_as( "rad" )
                                               - NCPA::math::deg2rad( az ) ),
                               ws.get_units() );
        }

        static scalar_u_t w2wc( double ws, const std::string& ws_units,
                                double wd, const std::string& wd_units,
                                double az ) {
            return w2wc( scalar_u_t( ws, ws_units ),
                         scalar_u_t( wd, wd_units ), az );
        }

        static vector_u_t w2wc( const vector_u_t& ws, const vector_u_t& wd,
                                double az ) {
            vector_u_t wc = ws;
            for (size_t i = 0; i < ws.size(); i++) {
                wc.set( i,
                        w2wc( ws.get_scalar( i ), wd.get_scalar( i ), az ) );
            }
            return wc;
        }

        static vector_u_t w2wc( const std::vector<double>& ws,
                                const std::string& ws_units,
                                const std::vector<double>& wd,
                                const std::string& wd_units, double az ) {
            return w2wc( vector_u_t( ws, ws_units ),
                         vector_u_t( wd, wd_units ), az );
        }

        std::vector<double> w2wc( const NCPA::arrays::ndvector<1, double>& u,
                                  const units_ptr_t u_units,
                                  const NCPA::arrays::ndvector<1, double>& v,
                                  const units_ptr_t v_units, double az ) {
            std::vector<double> wc( u.size() );
            for (size_t i = 0; i < u.size(); ++i) {
                wc[ i ]
                    = w2wc(
                          scalar_u_t( *static_cast<const double *>( &u[ i ] ),
                                      u_units ),
                          scalar_u_t( *static_cast<const double *>( &v[ i ] ),
                                      v_units ),
                          az )
                          .get();
            }
            return wc;
        }

        double w2wc( const NCPA::arrays::ndvector<0, double>& u,
                     const units_ptr_t u_units,
                     const NCPA::arrays::ndvector<0, double>& v,
                     const units_ptr_t v_units, double az ) {
            return w2wc( scalar_u_t( u, u_units ), scalar_u_t( v, v_units ),
                         az )
                .get();
        }

        template<size_t N>
        NCPA::arrays::ndvector<N, double> w2wc(
            const NCPA::arrays::ndvector<N, double>& u,
            const units_ptr_t u_units,
            const NCPA::arrays::ndvector<N, double>& v,
            const units_ptr_t v_units, double az ) {
            NCPA::arrays::ndvector<N, double> wc = u;

            for (size_t i = 0; i < u.size(); ++i) {
                wc[ i ] = w2wc(
                    NCPA::arrays::ndvector<N - 1, double>( u[ i ] ), u_units,
                    NCPA::arrays::ndvector<N - 1, double>( v[ i ] ), v_units,
                    az );
            }
            return wc;
        }

        // attenuation from altitude, temperature, pressure,  density and
        // frequency
        static scalar_u_t attenuation_sutherland_bass( const scalar_u_t& zs,
                                                       const scalar_u_t& ts,
                                                       const scalar_u_t& ps,
                                                       const scalar_u_t& ds,
                                                       double freq ) {
            double P_z = ps.get_as( "Pa" ), D_z = ds.get_as( "kg/m3" ),
                   T_z = ts.get_as( "K" ), z = zs.get_as( "km" );

            constexpr double mu_o
                = 18.192E-6;                // Reference viscosity [kg/(m*s)]
            constexpr double T_o = 293.15;  // Reference temperature [K]
            constexpr double P_o = 101325;  // Reference pressure [Pa]
            constexpr double S   = 117.0;   // Sutherland constant [K]

            // heat capacity|volume for O2, N2, CO2, and O3
            constexpr double Cv_R[]  = { 5.0 / 2.0, 5.0 / 2.0, 3.0, 3.0 };
            // heat capacity|pressure for O2, N2, CO2, and O3
            constexpr double Cp_R[]  = { 7.0 / 2.0, 7.0 / 2.0, 4.0, 4.0 };
            // characteristic temperature for O2, N2, CO2, and O3
            constexpr double theta[] = { 2239.1, 3352.0, 915.0, 1037.0 };

            // calculated parameters
            double c_snd_z
                = std::sqrt( constants::GAMMA() * P_z / D_z );  // m/s
            double mu = mu_o * std::sqrt( T_z / T_o )
                      * ( ( 1.0 + S / T_o )
                          / ( 1.0 + S / T_z ) );  // Viscosity [kg/(m*s)]
            double nu = ( 8.0 * NCPA::math::PI * freq * mu )
                      / ( 3.0 * P_z );            // Nondimensional frequency

                                                  // Gas fractions
            double X[ 7 ];

            //-------- Gas fraction polynomial fits
            //-----------------------------------
            if (z > 90.0) {  // O2 profile
                X[ 0 ] = std::pow( 10.0,
                                   49.296 - ( 1.5524 * z )
                                       + ( 1.8714E-2 * std::pow( z, 2 ) )
                                       - ( 1.1069E-4 * std::pow( z, 3 ) )
                                       + ( 3.199E-7 * std::pow( z, 4 ) )
                                       - ( 3.6211E-10 * std::pow( z, 5 ) ) );
            } else {
                X[ 0 ] = std::pow( 10.0, -0.67887 );
            }

            if (z > 76.0) {  // N2 profile
                X[ 1 ]
                    = std::pow( 10.0, ( 1.3972E-1 ) - ( 5.6269E-3 * z )
                                          + ( 3.9407E-5 * std::pow( z, 2 ) )
                                          - ( 1.0737E-7 * std::pow( z, 3 ) ) );
            } else {
                X[ 1 ] = std::pow( 10.0, -0.10744 );
            }

            X[ 2 ] = std::pow( 10.0, -3.3979 );  // CO2 profile

            if (z > 80.0) {                      // O3 profile
                X[ 3 ] = std::pow( 10.0, -4.234 - ( 3.0975E-2 * z ) );
            } else {
                X[ 3 ]
                    = std::pow( 10.0, -19.027 + ( 1.3093 * z )
                                          - ( 4.6496E-2 * std::pow( z, 2 ) )
                                          + ( 7.8543E-4 * std::pow( z, 3 ) )
                                          - ( 6.5169E-6 * std::pow( z, 4 ) )
                                          + ( 2.1343E-8 * std::pow( z, 5 ) ) );
            }

            if (z > 95.0) {  // O profile
                X[ 4 ]
                    = std::pow( 10.0, -3.2456 + ( 4.6642E-2 * z )
                                          - ( 2.6894E-4 * std::pow( z, 2 ) )
                                          + ( 5.264E-7 * std::pow( z, 3 ) ) );
            } else {
                X[ 4 ]
                    = std::pow( 10.0, -11.195 + ( 1.5408E-1 * z )
                                          - ( 1.4348E-3 * std::pow( z, 2 ) )
                                          + ( 1.0166E-5 * std::pow( z, 3 ) ) );
            }

            // N profile
            X[ 5 ] = std::pow( 10.0, -53.746 + ( 1.5439 * z )
                                         - ( 1.8824E-2 * std::pow( z, 2 ) )
                                         + ( 1.1587E-4 * std::pow( z, 3 ) )
                                         - ( 3.5399E-7 * std::pow( z, 4 ) )
                                         + ( 4.2609E-10 * std::pow( z, 5 ) ) );

            if (z > 30.0) {  // H2O profile
                X[ 6 ] = std::pow( 10.0,
                                   -4.2563 + ( 7.6245E-2 * z )
                                       - ( 2.1824E-3 * std::pow( z, 2 ) )
                                       - ( 2.3010E-6 * std::pow( z, 3 ) )
                                       + ( 2.4265E-7 * std::pow( z, 4 ) )
                                       - ( 1.2500E-09 * std::pow( z, 5 ) ) );
            } else {
                if (z > 100.0) {
                    X[ 6 ] = std::pow( 10.0, -0.62534 - ( 8.3665E-2 * z ) );
                } else {
                    X[ 6 ] = std::pow(
                        10.0, -1.7491 + ( 4.4986E-2 * z )
                                  - ( 6.8549E-2 * std::pow( z, 2 ) )
                                  + ( 5.4639E-3 * std::pow( z, 3 ) )
                                  - ( 1.5539E-4 * std::pow( z, 4 ) )
                                  + ( 1.5063E-06 * std::pow( z, 5 ) ) );
                }
            }
            double X_ON = ( X[ 0 ] + X[ 1 ] ) / 0.9903;

            //-------- Rotational collision
            // number-------------------------------------
            double Z_rot_0
                = 54.1
                * std::exp( -17.3 * ( std::pow( T_z, -1.0 / 3.0 ) ) );  // O2
            double Z_rot_1
                = 63.3
                * std::exp( -16.7 * ( std::pow( T_z, -1.0 / 3.0 ) ) );  // N2
            double Z_rot_
                = 1.0 / ( ( X[ 1 ] / Z_rot_1 ) + ( X[ 0 ] / Z_rot_0 ) );

            //-------- Nondimensional atmospheric
            // quantities---------------------------
            double sigma = 5.0 / std::sqrt( 21.0 );
            double nn    = ( 4.0 / 5.0 ) * std::sqrt( 3.0 / 7.0 ) * Z_rot_;
            double chi   = 3.0 * nn * nu / 4.0;
            double cchi  = 2.36 * chi;

            //---------Classical + rotational
            // loss/dispersion--------------------------
            double a_cl
                = ( 2 * NCPA::math::PI * freq / c_snd_z )
                * std::sqrt( 0.5 * ( std::sqrt( 1 + std::pow( nu, 2 ) ) - 1 )
                             * ( 1 + std::pow( cchi, 2 ) )
                             / ( ( 1 + std::pow( nu, 2 ) )
                                 * ( 1 + std::pow( sigma * cchi, 2 ) ) ) );
            double a_rot
                = ( 2 * NCPA::math::PI * freq / c_snd_z ) * X_ON
                * ( ( pow( sigma, 2 ) - 1 ) * chi / ( 2 * sigma ) )
                * sqrt( 0.5 * ( sqrt( 1 + pow( nu, 2 ) ) + 1 )
                        / ( ( 1 + pow( nu, 2 ) ) * ( 1 + pow( cchi, 2 ) ) ) );
            double a_diff = 0.003 * a_cl;

            //---------Vibrational
            // relaxation-------------------------------------------
            double Tr = std::pow( T_z / T_o, -1.0 / 3.0 ) - 1.0;
            double A1 = ( X[ 0 ] + X[ 1 ] ) * 24.0 * std::exp( -9.16 * Tr );
            double A2 = ( X[ 4 ] + X[ 5 ] ) * 2400.0;
            double B  = 40400.0 * std::exp( 10.0 * Tr );
            double C  = 0.02 * std::exp( -11.2 * Tr );
            double D  = 0.391 * std::exp( 8.41 * Tr );
            double E  = 9.0 * std::exp( -19.9 * Tr );
            double F  = 60000.0;
            double G  = 28000.0 * std::exp( -4.17 * Tr );
            double H  = 22000.0 * std::exp( -7.68 * Tr );
            double I  = 15100.0 * std::exp( -10.4 * Tr );
            double J  = 11500.0 * std::exp( -9.17 * Tr );
            double K  = ( 8.48E08 ) * std::exp( 9.17 * Tr );
            double L  = std::exp( -7.72 * Tr );
            double ZZ = H * X[ 2 ] + I * ( X[ 0 ] + 0.5 * X[ 4 ] )
                      + J * ( X[ 1 ] + 0.5 * X[ 5 ] )
                      + K * ( X[ 6 ] + X[ 3 ] );
            double hu = 100.0 * ( X[ 3 ] + X[ 6 ] );
            double f_vib[ 4 ], a_vib_c[ 4 ];
            f_vib[ 0 ] = ( P_z / P_o ) * ( mu_o / mu )
                       * ( A1 + A2 + B * hu * ( C + hu ) * ( D + hu ) );
            f_vib[ 1 ] = ( P_z / P_o ) * ( mu_o / mu )
                       * ( E + F * X[ 3 ] + G * X[ 6 ] );
            f_vib[ 2 ] = ( P_z / P_o ) * ( mu_o / mu ) * ZZ;
            f_vib[ 3 ] = ( P_z / P_o ) * ( mu_o / mu ) * ( 1.2E5 ) * L;

            double a_vib = 0.0;
            for (size_t m = 0; m < 4; m++) {
                double C_R
                    = ( ( std::pow( theta[ m ] / T_z, 2 ) )
                        * std::exp( -theta[ m ] / T_z ) )
                    / ( std::pow( 1 - std::exp( -theta[ m ] / T_z ), 2 ) );
                double A_max = ( X[ m ] * ( NCPA::math::PI / 2 ) * C_R )
                             / ( Cp_R[ m ] * ( Cv_R[ m ] + C_R ) );
                // A_max      = (X[m]*(PI/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
                a_vib_c[ m ] = ( A_max / c_snd_z )
                             * ( ( 2 * ( pow( freq, 2 ) ) / f_vib[ m ] )
                                 / ( 1 + pow( freq / f_vib[ m ], 2 ) ) );
                a_vib = a_vib + a_vib_c[ m ];
            }

            return scalar_u_t( a_cl + a_rot + a_diff + a_vib, "np/m" );
        }

        static scalar_u_t attenuation_sutherland_bass(
            double z, const std::string& z_units, double t,
            const std::string& t_units, double p, const std::string& p_units,
            double d, const std::string& d_units, double freq ) {
            return attenuation_sutherland_bass(
                scalar_u_t( z, z_units ), scalar_u_t( t, t_units ),
                scalar_u_t( p, p_units ), scalar_u_t( d, d_units ), freq );
        }

        static vector_u_t attenuation_sutherland_bass( const vector_u_t& zs,
                                                       const vector_u_t& ts,
                                                       const vector_u_t& ps,
                                                       const vector_u_t& ds,
                                                       double freq ) {
            vector_u_t alpha( zs.size(), "np/m" );
            for (size_t i = 0; i < zs.size(); i++) {
                alpha.set( i, attenuation_sutherland_bass(
                                  zs.get_scalar( i ), ts.get_scalar( i ),
                                  ps.get_scalar( i ), ds.get_scalar( i ),
                                  freq ) );
            }
            return alpha;
        }

        static vector_u_t attenuation_sutherland_bass(
            std::vector<double>& z, const std::string& z_units,
            std::vector<double>& t, const std::string& t_units,
            std::vector<double>& p, const std::string& p_units,
            std::vector<double>& d, const std::string& d_units, double freq ) {
            return attenuation_sutherland_bass(
                vector_u_t( z, z_units ), vector_u_t( t, t_units ),
                vector_u_t( p, p_units ), vector_u_t( d, d_units ), freq );
        }

        std::vector<scalar_u_t> attenuation_sutherland_bass(
            const NCPA::arrays::ndvector<1, scalar_u_t>& z,
            const NCPA::arrays::ndvector<1, scalar_u_t>& t,
            const NCPA::arrays::ndvector<1, scalar_u_t>& p,
            const NCPA::arrays::ndvector<1, scalar_u_t>& d, double freq ) {
            std::vector<scalar_u_t> alpha( z.size() );
            for (size_t i = 0; i < z.size(); ++i) {
                alpha[ i ] = attenuation_sutherland_bass(
                    z[ i ], t[ i ], p[ i ], d[ i ], freq );
            }
            return alpha;
        }

        // scalar_u_t attenuation_sutherland_bass(
        //     const NCPA::arrays::ndvector<0, scalar_u_t>& z,
        //     const NCPA::arrays::ndvector<0, scalar_u_t>& t,
        //     const NCPA::arrays::ndvector<0, scalar_u_t>& p,
        //     const NCPA::arrays::ndvector<0, scalar_u_t>& d, double freq  ) {
        //     return attenuation_sutherland_bass( z, t, u, p, freq )
        //         .get();
        // }

        template<size_t N>
        NCPA::arrays::ndvector<N, scalar_u_t> attenuation_sutherland_bass(
            const NCPA::arrays::ndvector<N, scalar_u_t>& z,
            const NCPA::arrays::ndvector<N, scalar_u_t>& t,
            const NCPA::arrays::ndvector<N, scalar_u_t>& p,
            const NCPA::arrays::ndvector<N, scalar_u_t>& d, double freq ) {
            NCPA::arrays::ndvector<N, scalar_u_t> alpha = z;

            for (size_t i = 0; i < z.size(); ++i) {
                alpha[ i ] = attenuation_sutherland_bass(
                    NCPA::arrays::ndvector<N - 1, scalar_u_t>( z[ i ] ),
                    NCPA::arrays::ndvector<N - 1, scalar_u_t>( t[ i ] ),
                    NCPA::arrays::ndvector<N - 1, scalar_u_t>( p[ i ] ),
                    NCPA::arrays::ndvector<N - 1, scalar_u_t>( d[ i ] ), freq );
            }
            return alpha;
        }


    }  // namespace atmos
}  // namespace NCPA
