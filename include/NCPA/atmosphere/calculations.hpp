#pragma once

// #include "NCPA/atmosphere/constants.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/ndvector.hpp"
#include "NCPA/units.hpp"

#include <cmath>

// Units for t2c() calculations
#define NCPA_ATMOS_T2C_T_UNITS_STR "K"
#define NCPA_ATMOS_T2C_C_UNITS_STR "m/s"
#define NCPA_ATMOS_T2C_T_UNITS \
    NCPA::units::Units::from_string( NCPA_ATMOS_T2C_T_UNITS_STR )
#define NCPA_ATMOS_T2C_C_UNITS \
    NCPA::units::Units::from_string( NCPA_ATMOS_T2C_C_UNITS_STR )

// Units for c2t() calculations
#define NCPA_ATMOS_C2T_T_UNITS_STR "K"
#define NCPA_ATMOS_C2T_C_UNITS_STR "m/s"
#define NCPA_ATMOS_C2T_T_UNITS \
    NCPA::units::Units::from_string( NCPA_ATMOS_C2T_T_UNITS_STR )
#define NCPA_ATMOS_C2T_C_UNITS \
    NCPA::units::Units::from_string( NCPA_ATMOS_C2T_C_UNITS_STR )

// units for pd2c() calculations
#define NCPA_ATMOS_PD2C_P_UNITS_STR "Pa"
#define NCPA_ATMOS_PD2C_D_UNITS_STR "kg/m3"
#define NCPA_ATMOS_PD2C_C_UNITS_STR "m/s"
#define NCPA_ATMOS_PD2C_P_UNITS \
    NCPA::units::Units::from_string( NCPA_ATMOS_PD2C_P_UNITS_STR )
#define NCPA_ATMOS_PD2C_D_UNITS \
    NCPA::units::Units::from_string( NCPA_ATMOS_PD2C_D_UNITS_STR )
#define NCPA_ATMOS_PD2C_C_UNITS \
    NCPA::units::Units::from_string( NCPA_ATMOS_PD2C_C_UNITS_STR )

// units for uv2wd() calculations
#define NCPA_ATMOS_UV2WD_WD_UNITS_STR "deg"
#define NCPA_ATMOS_UV2WD_WD_UNITS \
    NCPA::units::Units::from_string( NCPA_ATMOS_UV2WD_WD_UNITS_STR )

// units for sutherland_bass_attenuation() calculations
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_P_UNITS_STR     "Pa"
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_D_UNITS_STR     "kg/m3"
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_T_UNITS_STR     "K"
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_Z_UNITS_STR     "km"
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS_STR "np/m"
// #define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_F_UNITS_STR     "Hz"
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_P_UNITS \
    NCPA::units::Units::from_string(                   \
        NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_P_UNITS_STR )
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_D_UNITS \
    NCPA::units::Units::from_string(                   \
        NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_D_UNITS_STR )
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_T_UNITS \
    NCPA::units::Units::from_string(                   \
        NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_T_UNITS_STR )
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_Z_UNITS \
    NCPA::units::Units::from_string(                   \
        NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_Z_UNITS_STR )
#define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS \
    NCPA::units::Units::from_string(                       \
        NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS_STR )

// #define NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_F_UNITS \
//     NCPA::units::Units::from_string(                   \
//         NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_F_UNITS_STR )

namespace NCPA {
    namespace atmos {

        // temperature to sound speed
        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> t2c(
            const NCPA::units::ScalarWithUnits<T>& t ) {
            return NCPA::units::ScalarWithUnits<T>(
                std::sqrt( t.get_as( NCPA_ATMOS_T2C_T_UNITS )
                           * (T)( constants::GAMMA() * constants::R() ) ),
                NCPA_ATMOS_T2C_C_UNITS );
        }

        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> t2c(
            T t, const std::string& units ) {
            return t2c( NCPA::units::ScalarWithUnits<T>( t, units ) );
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> t2c(
            const NCPA::units::VectorWithUnits<T>& t ) {
            NCPA::units::VectorWithUnits<T> c( t.size(),
                                               NCPA_ATMOS_T2C_C_UNITS );
            for (size_t i = 0; i < t.size(); i++) {
                c.set( i, t2c( t.get_scalar( i ) ) );
            }
            return c;
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> t2c(
            const std::vector<T>& t, const std::string& units ) {
            return t2c( NCPA::units::VectorWithUnits<T>( t, units ) );
        }

        template<typename T = double>
        static std::vector<T> t2c( const NCPA::arrays::ndvector<1, T>& t,
                                   const units_ptr_t units ) {
            std::vector<T> c( t.size() );
            for (size_t i = 0; i < t.size(); ++i) {
                c[ i ] = t2c( NCPA::units::ScalarWithUnits<T>(
                                  *static_cast<const T *>( &t[ i ] ), units ) )
                             .get();
            }
            return c;
        }

        template<typename T = double>
        static T t2c( const NCPA::arrays::ndvector<0, T>& t,
                      const units_ptr_t units ) {
            return t2c( NCPA::units::ScalarWithUnits<T>( t, units ) ).get();
        }

        template<size_t N, typename T = double>
        NCPA::arrays::ndvector<N, T> t2c(
            const NCPA::arrays::ndvector<N, T>& t, const units_ptr_t units ) {
            NCPA::arrays::ndvector<N, T> c = t;

            for (size_t i = 0; i < t.size(); ++i) {
                c[ i ]
                    = t2c( NCPA::arrays::ndvector<N - 1, T>( t[ i ] ), units );
            }
            return c;
        }

        // sound speed to temperature
        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> c2t(
            const NCPA::units::ScalarWithUnits<T>& c ) {
            T Cmps = c.get_as( NCPA_ATMOS_C2T_C_UNITS );
            return NCPA::units::ScalarWithUnits<T>(
                Cmps * Cmps / constants::GAMMA() / constants::R(),
                NCPA_ATMOS_C2T_T_UNITS );
        }

        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> c2t(
            T t, const std::string& units ) {
            return c2t( NCPA::units::ScalarWithUnits<T>( t, units ) );
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> c2t(
            const NCPA::units::VectorWithUnits<T>& c ) {
            NCPA::units::VectorWithUnits<T> t( c.size(),
                                               NCPA_ATMOS_C2T_T_UNITS );
            for (size_t i = 0; i < c.size(); i++) {
                t.set( i, c2t( c.get_scalar( i ) ) );
            }
            return t;
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> c2t(
            const std::vector<T>& c, const std::string& units ) {
            return c2t( NCPA::units::VectorWithUnits<T>( c, units ) );
        }

        template<typename T = double>
        static std::vector<T> c2t( const NCPA::arrays::ndvector<1, T>& c,
                                   const units_ptr_t units ) {
            std::vector<T> t( c.size() );
            for (size_t i = 0; i < c.size(); ++i) {
                t[ i ] = c2t( NCPA::units::ScalarWithUnits<T>(
                                  *static_cast<const T *>( &c[ i ] ), units ) )
                             .get();
            }
            return t;
        }

        template<typename T = double>
        static T c2t( const NCPA::arrays::ndvector<0, T>& c,
                      const units_ptr_t units ) {
            return c2t( NCPA::units::ScalarWithUnits<T>( c, units ) ).get();
        }

        template<size_t N, typename T = double>
        NCPA::arrays::ndvector<N, T> c2t(
            const NCPA::arrays::ndvector<N, T>& c, const units_ptr_t units ) {
            NCPA::arrays::ndvector<N, T> t = c;

            for (size_t i = 0; i < c.size(); ++i) {
                t[ i ]
                    = c2t( NCPA::arrays::ndvector<N - 1, T>( c[ i ] ), units );
            }
            return t;
        }

        // pressure and density to sound speed
template<typename T = double>
        static NCPA::units::Vector3DWithUnits<T> pd2c(
            const NCPA::units::Vector3DWithUnits<T>& p,
            const NCPA::units::Vector3DWithUnits<T>& d,
            const units_ptr_t units_out = NCPA_ATMOS_PD2C_C_UNITS ) {
            return NCPA::units::Vector3DWithUnits<T>(
                pd2c( static_cast<const NCPA::arrays::ndvector<3, T>&>( p ),
                      p.get_units(),
                      static_cast<const NCPA::arrays::ndvector<3, T>&>( d ),
                      d.get_units(), units_out ),
                units_out );
            // NCPA::arrays::ndvector<3, T> c = pd2c(
            //     static_cast<const NCPA::arrays::ndvector<3, T>&>( p ),
            //     p.get_units(),
            //     static_cast<const NCPA::arrays::ndvector<3, T>&>( d ),
            //     d.get_units(), units_out );
            // return NCPA::units::Vector3DWithUnits<T>( c, units_out );
        }

        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> pd2c(
            const NCPA::units::ScalarWithUnits<T>& p,
            const NCPA::units::ScalarWithUnits<T>& d,
            const units_ptr_t units_out = NCPA_ATMOS_PD2C_C_UNITS ) {
            return NCPA::units::ScalarWithUnits<T>(
                       std::sqrt( constants::GAMMA()
                                  * p.get_as( NCPA_ATMOS_PD2C_P_UNITS )
                                  / d.get_as( NCPA_ATMOS_PD2C_D_UNITS ) ),
                       NCPA_ATMOS_PD2C_C_UNITS )
                .as( *units_out );
        }

        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> pd2c(
            T p, const std::string& p_units, T d, const std::string& d_units,
            const std::string& units_out = NCPA_ATMOS_PD2C_C_UNITS_STR ) {
            return pd2c( NCPA::units::ScalarWithUnits<T>( p, p_units ),
                         NCPA::units::ScalarWithUnits<T>( d, d_units ),
                         NCPA::units::Units::from_string( units_out ) );
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> pd2c(
            const NCPA::units::VectorWithUnits<T>& p,
            const NCPA::units::VectorWithUnits<T>& d,
            const units_ptr_t units_out = NCPA_ATMOS_PD2C_C_UNITS ) {
            NCPA::units::VectorWithUnits<T> c( p.size(), units_out );
            for (size_t i = 0; i < p.size(); i++) {
                c.set( i, pd2c( p.get_scalar( i ), d.get_scalar( i ) ) );
            }
            return c;
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> pd2c(
            const std::vector<T>& p, const std::string& p_units,
            const std::vector<T>& d, const std::string& d_units,
            const std::string& units_out = NCPA_ATMOS_PD2C_C_UNITS_STR ) {
            return pd2c( NCPA::units::VectorWithUnits<T>( p, p_units ),
                         NCPA::units::VectorWithUnits<T>( d, d_units ),
                         NCPA::units::Units::from_string( units_out ) );
        }

        template<typename T = double>
        static T pd2c( const T& p, const units_ptr_t p_units, const T& d,
                       const units_ptr_t d_units,
                       const units_ptr_t units_out
                       = NCPA_ATMOS_PD2C_C_UNITS ) {
            return pd2c( NCPA::units::ScalarWithUnits<T>( p, p_units ),
                         NCPA::units::ScalarWithUnits<T>( d, d_units ) )
                .get_as( units_out );
        }

        template<typename T = double>
        static T pd2c( const NCPA::arrays::ndvector<0, T>& p,
                       const units_ptr_t p_units,
                       const NCPA::arrays::ndvector<0, T>& d,
                       const units_ptr_t d_units,
                       const units_ptr_t units_out
                       = NCPA_ATMOS_PD2C_C_UNITS ) {
            return pd2c( NCPA::units::ScalarWithUnits<T>( p, p_units ),
                         NCPA::units::ScalarWithUnits<T>( d, d_units ) )
                .get_as( units_out );
        }

        template<typename T = double>
        static std::vector<T> pd2c( const NCPA::arrays::ndvector<1, T>& p,
                                    const units_ptr_t p_units,
                                    const NCPA::arrays::ndvector<1, T>& d,
                                    const units_ptr_t d_units,
                                    const units_ptr_t units_out
                                    = NCPA_ATMOS_PD2C_C_UNITS ) {
            std::vector<T> c( p.size() );
            for (size_t i = 0; i < p.size(); ++i) {
                c[ i ]
                    = pd2c( NCPA::units::ScalarWithUnits<T>(
                                *static_cast<const T *>( &p[ i ] ), p_units ),
                            NCPA::units::ScalarWithUnits<T>(
                                *static_cast<const T *>( &d[ i ] ), d_units ) )
                          .get_as( units_out );
            }
            return c;
        }

        template<size_t N, typename T = double>
        NCPA::arrays::ndvector<N, T> pd2c(
            const NCPA::arrays::ndvector<N, T>& p, const units_ptr_t p_units,
            const NCPA::arrays::ndvector<N, T>& d, const units_ptr_t d_units,
            const units_ptr_t units_out = NCPA_ATMOS_PD2C_C_UNITS ) {
            NCPA::arrays::ndvector<N, T> c = p;

            for (size_t i = 0; i < p.size(); ++i) {
                c[ i ] = pd2c( NCPA::arrays::ndvector<N - 1, T>( p[ i ] ),
                               p_units,
                               NCPA::arrays::ndvector<N - 1, T>( d[ i ] ),
                               d_units, units_out );
            }
            return c;
        }

        

        // U & V winds to wind speed
        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> uv2ws(
            const NCPA::units::ScalarWithUnits<T>& u,
            const NCPA::units::ScalarWithUnits<T>& v,
            const units_ptr_t units_out = nullptr ) {
            units_ptr_t uout
                = ( units_out == nullptr ? u.get_units() : units_out );
            T uu = u.get_as( uout );
            T vv = v.get_as( uout );
            return NCPA::units::ScalarWithUnits<T>(
                std::sqrt( uu * uu + vv * vv ), uout );
            // return ( units_out == nullptr
            //              ? NCPA::units::ScalarWithUnits<T>( std::sqrt(
            //              uu * uu + vv * vv ), uout ) :
            //              NCPA::units::ScalarWithUnits<T>( std::sqrt( uu
            //              * uu + vv * vv ), u.get_units())
            //                    .as( units_out ) );
        }

        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> uv2ws(
            T u, const std::string& u_units, T v, const std::string& v_units,
            const std::string& units_out = "" ) {
            units_ptr_t uout
                = ( units_out.length() == 0
                        ? nullptr
                        : NCPA::units::Units::from_string( units_out ) );
            return uv2ws( NCPA::units::ScalarWithUnits<T>( u, u_units ),
                          NCPA::units::ScalarWithUnits<T>( v, v_units ),
                          uout );
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> uv2ws(
            const NCPA::units::VectorWithUnits<T>& u,
            const NCPA::units::VectorWithUnits<T>& v,
            const units_ptr_t units_out = nullptr ) {
            NCPA::units::VectorWithUnits<T> ws( u.size(), u.get_units() );
            for (size_t i = 0; i < u.size(); i++) {
                ws.set( i, uv2ws( u.get_scalar( i ), v.get_scalar( i ),
                                  units_out ) );
            }
            return ws;
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> uv2ws(
            const std::vector<T>& u, const std::string& u_units,
            const std::vector<T>& v, const std::string& v_units,
            const std::string& units_out = "" ) {
            units_ptr_t uout
                = ( units_out.length() == 0
                        ? nullptr
                        : NCPA::units::Units::from_string( units_out ) );
            return uv2ws( NCPA::units::VectorWithUnits<T>( u, u_units ),
                          NCPA::units::VectorWithUnits<T>( v, v_units ),
                          uout );
        }

        template<typename T = double>
        static T uv2ws( const NCPA::arrays::ndvector<0, T>& u,
                        const units_ptr_t u_units,
                        const NCPA::arrays::ndvector<0, T>& v,
                        const units_ptr_t v_units,
                        const units_ptr_t units_out = nullptr ) {
            return uv2ws( NCPA::units::ScalarWithUnits<T>( u, u_units ),
                          NCPA::units::ScalarWithUnits<T>( v, v_units ),
                          units_out )
                .get();
        }

        template<typename T = double>
        static std::vector<T> uv2ws( const NCPA::arrays::ndvector<1, T>& u,
                                     const units_ptr_t u_units,
                                     const NCPA::arrays::ndvector<1, T>& v,
                                     const units_ptr_t v_units,
                                     const units_ptr_t units_out = nullptr ) {
            std::vector<T> ws( u.size() );
            for (size_t i = 0; i < u.size(); ++i) {
                ws[ i ]
                    = uv2ws( NCPA::units::ScalarWithUnits<T>(
                                 *static_cast<const T *>( &u[ i ] ), u_units ),
                             NCPA::units::ScalarWithUnits<T>(
                                 *static_cast<const T *>( &v[ i ] ), v_units ),
                             units_out )
                          .get();
            }
            return ws;
        }

        template<size_t N, typename T = double>
        NCPA::arrays::ndvector<N, T> uv2ws(
            const NCPA::arrays::ndvector<N, T>& u, const units_ptr_t u_units,
            const NCPA::arrays::ndvector<N, T>& v, const units_ptr_t v_units,
            const units_ptr_t units_out = nullptr ) {
            NCPA::arrays::ndvector<N, T> ws = u;

            for (size_t i = 0; i < u.size(); ++i) {
                ws[ i ] = uv2ws( NCPA::arrays::ndvector<N - 1, T>( u[ i ] ),
                                 u_units,
                                 NCPA::arrays::ndvector<N - 1, T>( v[ i ] ),
                                 v_units, units_out );
            }
            return ws;
        }

        // U & V winds to wind direction
        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> uv2wd(
            const NCPA::units::ScalarWithUnits<T>& u,
            const NCPA::units::ScalarWithUnits<T>& v,
            const units_ptr_t units_out = NCPA_ATMOS_UV2WD_WD_UNITS ) {
            return NCPA::units::ScalarWithUnits<T>(
                       NCPA::math::math2az( NCPA::math::rad2deg( std::atan2(
                           v.get_as( *u.get_units() ), u.get() ) ) ),
                       NCPA_ATMOS_UV2WD_WD_UNITS )
                .as( units_out );
        }

        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> uv2wd(
            T u, const std::string& u_units, T v, const std::string& v_units,
            const std::string& units_out = NCPA_ATMOS_UV2WD_WD_UNITS_STR ) {
            return uv2wd( NCPA::units::ScalarWithUnits<T>( u, u_units ),
                          NCPA::units::ScalarWithUnits<T>( v, v_units ),
                          NCPA::units::Units::from_string( units_out ) );
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> uv2wd(
            const NCPA::units::VectorWithUnits<T>& u,
            const NCPA::units::VectorWithUnits<T>& v,
            const units_ptr_t units_out = NCPA_ATMOS_UV2WD_WD_UNITS ) {
            NCPA::units::VectorWithUnits<T> wd( u.size(), units_out );
            for (size_t i = 0; i < u.size(); i++) {
                wd.set( i, uv2wd( u.get_scalar( i ), v.get_scalar( i ),
                                  units_out ) );
            }
            return wd;
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> uv2wd(
            const std::vector<T>& u, const std::string& u_units,
            const std::vector<T>& v, const std::string& v_units,
            const std::string& units_out = NCPA_ATMOS_UV2WD_WD_UNITS_STR ) {
            return uv2wd( NCPA::units::VectorWithUnits<T>( u, u_units ),
                          NCPA::units::VectorWithUnits<T>( v, v_units ),
                          NCPA::units::Units::from_string( units_out ) );
        }

        template<typename T = double>
        static std::vector<T> uv2wd( const NCPA::arrays::ndvector<1, T>& u,
                                     const units_ptr_t u_units,
                                     const NCPA::arrays::ndvector<1, T>& v,
                                     const units_ptr_t v_units,
                                     const units_ptr_t units_out
                                     = NCPA_ATMOS_UV2WD_WD_UNITS ) {
            std::vector<T> wd( u.size() );
            for (size_t i = 0; i < u.size(); ++i) {
                wd[ i ]
                    = uv2wd( NCPA::units::ScalarWithUnits<T>(
                                 *static_cast<const T *>( &u[ i ] ), u_units ),
                             NCPA::units::ScalarWithUnits<T>(
                                 *static_cast<const T *>( &v[ i ] ), v_units ),
                             units_out )
                          .get();
            }
            return wd;
        }

        template<typename T = double>
        static T uv2wd( const NCPA::arrays::ndvector<0, T>& u,
                        const units_ptr_t u_units,
                        const NCPA::arrays::ndvector<0, T>& v,
                        const units_ptr_t v_units,
                        const units_ptr_t units_out
                        = NCPA_ATMOS_UV2WD_WD_UNITS ) {
            return uv2wd( NCPA::units::ScalarWithUnits<T>( u, u_units ),
                          NCPA::units::ScalarWithUnits<T>( v, v_units ),
                          units_out )
                .get();
        }

        template<size_t N, typename T = double>
        NCPA::arrays::ndvector<N, T> uv2wd(
            const NCPA::arrays::ndvector<N, T>& u, const units_ptr_t u_units,
            const NCPA::arrays::ndvector<N, T>& v, const units_ptr_t v_units,
            const units_ptr_t units_out = NCPA_ATMOS_UV2WD_WD_UNITS ) {
            NCPA::arrays::ndvector<N, T> wd = u;

            for (size_t i = 0; i < u.size(); ++i) {
                wd[ i ] = uv2wd( NCPA::arrays::ndvector<N - 1, T>( u[ i ] ),
                                 u_units,
                                 NCPA::arrays::ndvector<N - 1, T>( v[ i ] ),
                                 v_units, units_out );
            }
            return wd;
        }

        // wind speed and direction to wind component along an azimuth
        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> w2wc(
            const NCPA::units::ScalarWithUnits<T>& ws,
            const NCPA::units::ScalarWithUnits<T>& wd, T az,
            const units_ptr_t units_out = nullptr ) {
            NCPA::units::ScalarWithUnits<T> wc(
                ws.get()
                    * std::cos( wd.get_as( "rad" )
                                - NCPA::math::deg2rad( az ) ),
                ws.get_units() );
            return ( units_out == nullptr ? wc : wc.as( units_out ) );
        }

        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> w2wc(
            T ws, const std::string& ws_units, T wd,
            const std::string& wd_units, T az,
            const std::string& units_out = "" ) {
            return w2wc(
                NCPA::units::ScalarWithUnits<T>( ws, ws_units ),
                NCPA::units::ScalarWithUnits<T>( wd, wd_units ), az,
                ( units_out.length() == 0
                      ? nullptr
                      : NCPA::units::Units::from_string( units_out ) ) );
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> w2wc(
            const NCPA::units::VectorWithUnits<T>& ws,
            const NCPA::units::VectorWithUnits<T>& wd, T az,
            const units_ptr_t units_out = nullptr ) {
            NCPA::units::VectorWithUnits<T> wc = ws;
            for (size_t i = 0; i < ws.size(); i++) {
                wc.set( i, w2wc( ws.get_scalar( i ), wd.get_scalar( i ), az,
                                 units_out ) );
            }
            return wc;
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> w2wc(
            const std::vector<T>& ws, const std::string& ws_units,
            const std::vector<T>& wd, const std::string& wd_units, T az,
            const std::string& units_out = "" ) {
            return w2wc(
                NCPA::units::VectorWithUnits<T>( ws, ws_units ),
                NCPA::units::VectorWithUnits<T>( wd, wd_units ), az,
                ( units_out.length() == 0
                      ? nullptr
                      : NCPA::units::Units::from_string( units_out ) ) );
        }

        template<typename T = double>
        static std::vector<T> w2wc( const NCPA::arrays::ndvector<1, T>& u,
                                    const units_ptr_t u_units,
                                    const NCPA::arrays::ndvector<1, T>& v,
                                    const units_ptr_t v_units, T az,
                                    const units_ptr_t units_out = nullptr ) {
            std::vector<T> wc( u.size() );
            for (size_t i = 0; i < u.size(); ++i) {
                wc[ i ]
                    = w2wc( NCPA::units::ScalarWithUnits<T>(
                                *static_cast<const T *>( &u[ i ] ), u_units ),
                            NCPA::units::ScalarWithUnits<T>(
                                *static_cast<const T *>( &v[ i ] ), v_units ),
                            az, units_out )
                          .get();
            }
            return wc;
        }

        template<typename T = double>
        static T w2wc( const NCPA::arrays::ndvector<0, T>& u,
                       const units_ptr_t u_units,
                       const NCPA::arrays::ndvector<0, T>& v,
                       const units_ptr_t v_units, T az,
                       const units_ptr_t units_out = nullptr ) {
            return w2wc( NCPA::units::ScalarWithUnits<T>( u, u_units ),
                         NCPA::units::ScalarWithUnits<T>( v, v_units ), az,
                         units_out )
                .get();
        }

        template<size_t N, typename T = double>
        NCPA::arrays::ndvector<N, T> w2wc(
            const NCPA::arrays::ndvector<N, T>& u, const units_ptr_t u_units,
            const NCPA::arrays::ndvector<N, T>& v, const units_ptr_t v_units,
            T az, const units_ptr_t units_out = nullptr ) {
            NCPA::arrays::ndvector<N, T> wc = u;
            for (size_t i = 0; i < u.size(); ++i) {
                wc[ i ] = w2wc( NCPA::arrays::ndvector<N - 1, T>( u[ i ] ),
                                u_units,
                                NCPA::arrays::ndvector<N - 1, T>( v[ i ] ),
                                v_units, az, units_out );
            }
            return wc;
        }

        // attenuation from altitude, temperature, pressure,  density and
        // frequency
        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> attenuation_sutherland_bass(
            const NCPA::units::ScalarWithUnits<T>& zs,
            const NCPA::units::ScalarWithUnits<T>& ts,
            const NCPA::units::ScalarWithUnits<T>& ps,
            const NCPA::units::ScalarWithUnits<T>& ds, T freq,
            const units_ptr_t units_out
            = NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS ) {
            T P_z
                = ps.get_as( NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_P_UNITS ),
                D_z
                = ds.get_as( NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_D_UNITS ),
                T_z
                = ts.get_as( NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_T_UNITS ),
                z
                = zs.get_as( NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_Z_UNITS );
            // freq = f.get_as( NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_F_UNITS
            // );

            constexpr T mu_o = 18.192E-6;  // Reference viscosity [kg/(m*s)]
            constexpr T T_o  = 293.15;     // Reference temperature [K]
            constexpr T P_o  = 101325;     // Reference pressure [Pa]
            constexpr T S    = 117.0;      // Sutherland constant [K]

            // heat capacity|volume for O2, N2, CO2, and O3
            constexpr T Cv_R[]  = { 5.0 / 2.0, 5.0 / 2.0, 3.0, 3.0 };
            // heat capacity|pressure for O2, N2, CO2, and O3
            constexpr T Cp_R[]  = { 7.0 / 2.0, 7.0 / 2.0, 4.0, 4.0 };
            // characteristic temperature for O2, N2, CO2, and O3
            constexpr T theta[] = { 2239.1, 3352.0, 915.0, 1037.0 };

            // calculated parameters
            T c_snd_z = std::sqrt( constants::GAMMA() * P_z / D_z );  // m/s
            T mu      = mu_o * std::sqrt( T_z / T_o )
                 * ( ( 1.0 + S / T_o )
                     / ( 1.0 + S / T_z ) );  // Viscosity [kg/(m*s)]
            T nu = ( 8.0 * NCPA::math::PI * freq * mu )
                 / ( 3.0 * P_z );            // Nondimensional frequency

                                             // Gas fractions
            T X[ 7 ];

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
            T X_ON = ( X[ 0 ] + X[ 1 ] ) / 0.9903;

            //-------- Rotational collision
            // number-------------------------------------
            T Z_rot_0
                = 54.1
                * std::exp( -17.3 * ( std::pow( T_z, -1.0 / 3.0 ) ) );  // O2
            T Z_rot_1
                = 63.3
                * std::exp( -16.7 * ( std::pow( T_z, -1.0 / 3.0 ) ) );  // N2
            T Z_rot_ = 1.0 / ( ( X[ 1 ] / Z_rot_1 ) + ( X[ 0 ] / Z_rot_0 ) );

            //-------- Nondimensional atmospheric
            // quantities---------------------------
            T sigma = 5.0 / std::sqrt( 21.0 );
            T nn    = ( 4.0 / 5.0 ) * std::sqrt( 3.0 / 7.0 ) * Z_rot_;
            T chi   = 3.0 * nn * nu / 4.0;
            T cchi  = 2.36 * chi;

            //---------Classical + rotational
            // loss/dispersion--------------------------
            T a_cl
                = ( 2 * NCPA::math::PI * freq / c_snd_z )
                * std::sqrt( 0.5 * ( std::sqrt( 1 + std::pow( nu, 2 ) ) - 1 )
                             * ( 1 + std::pow( cchi, 2 ) )
                             / ( ( 1 + std::pow( nu, 2 ) )
                                 * ( 1 + std::pow( sigma * cchi, 2 ) ) ) );
            T a_rot
                = ( 2 * NCPA::math::PI * freq / c_snd_z ) * X_ON
                * ( ( pow( sigma, 2 ) - 1 ) * chi / ( 2 * sigma ) )
                * sqrt( 0.5 * ( sqrt( 1 + pow( nu, 2 ) ) + 1 )
                        / ( ( 1 + pow( nu, 2 ) ) * ( 1 + pow( cchi, 2 ) ) ) );
            T a_diff = 0.003 * a_cl;

            //---------Vibrational
            // relaxation-------------------------------------------
            T Tr = std::pow( T_z / T_o, -1.0 / 3.0 ) - 1.0;
            T A1 = ( X[ 0 ] + X[ 1 ] ) * 24.0 * std::exp( -9.16 * Tr );
            T A2 = ( X[ 4 ] + X[ 5 ] ) * 2400.0;
            T B  = 40400.0 * std::exp( 10.0 * Tr );
            T C  = 0.02 * std::exp( -11.2 * Tr );
            T D  = 0.391 * std::exp( 8.41 * Tr );
            T E  = 9.0 * std::exp( -19.9 * Tr );
            T F  = 60000.0;
            T G  = 28000.0 * std::exp( -4.17 * Tr );
            T H  = 22000.0 * std::exp( -7.68 * Tr );
            T I  = 15100.0 * std::exp( -10.4 * Tr );
            T J  = 11500.0 * std::exp( -9.17 * Tr );
            T K  = ( 8.48E08 ) * std::exp( 9.17 * Tr );
            T L  = std::exp( -7.72 * Tr );
            T ZZ = H * X[ 2 ] + I * ( X[ 0 ] + 0.5 * X[ 4 ] )
                 + J * ( X[ 1 ] + 0.5 * X[ 5 ] ) + K * ( X[ 6 ] + X[ 3 ] );
            T hu = 100.0 * ( X[ 3 ] + X[ 6 ] );
            T f_vib[ 4 ], a_vib_c[ 4 ];
            f_vib[ 0 ] = ( P_z / P_o ) * ( mu_o / mu )
                       * ( A1 + A2 + B * hu * ( C + hu ) * ( D + hu ) );
            f_vib[ 1 ] = ( P_z / P_o ) * ( mu_o / mu )
                       * ( E + F * X[ 3 ] + G * X[ 6 ] );
            f_vib[ 2 ] = ( P_z / P_o ) * ( mu_o / mu ) * ZZ;
            f_vib[ 3 ] = ( P_z / P_o ) * ( mu_o / mu ) * ( 1.2E5 ) * L;

            T a_vib = 0.0;
            for (size_t m = 0; m < 4; m++) {
                T C_R = ( ( std::pow( theta[ m ] / T_z, 2 ) )
                          * std::exp( -theta[ m ] / T_z ) )
                      / ( std::pow( 1 - std::exp( -theta[ m ] / T_z ), 2 ) );
                T A_max = ( X[ m ] * ( NCPA::math::PI / 2 ) * C_R )
                        / ( Cp_R[ m ] * ( Cv_R[ m ] + C_R ) );
                // A_max      = (X[m]*(PI/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
                a_vib_c[ m ] = ( A_max / c_snd_z )
                             * ( ( 2 * ( pow( freq, 2 ) ) / f_vib[ m ] )
                                 / ( 1 + pow( freq / f_vib[ m ], 2 ) ) );
                a_vib = a_vib + a_vib_c[ m ];
            }

            NCPA::units::ScalarWithUnits<T> alpha(
                a_cl + a_rot + a_diff + a_vib,
                NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS );
            return alpha.as( units_out );
        }

        template<typename T = double>
        static NCPA::units::ScalarWithUnits<T> attenuation_sutherland_bass(
            T z, const std::string& z_units, T t, const std::string& t_units,
            T p, const std::string& p_units, T d, const std::string& d_units,
            T freq,
            const std::string& units_out
            = NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS_STR ) {
            return attenuation_sutherland_bass(
                NCPA::units::ScalarWithUnits<T>( z, z_units ),
                NCPA::units::ScalarWithUnits<T>( t, t_units ),
                NCPA::units::ScalarWithUnits<T>( p, p_units ),
                NCPA::units::ScalarWithUnits<T>( d, d_units ), freq,
                NCPA::units::Units::from_string( units_out ) );
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> attenuation_sutherland_bass(
            const NCPA::units::VectorWithUnits<T>& zs,
            const NCPA::units::VectorWithUnits<T>& ts,
            const NCPA::units::VectorWithUnits<T>& ps,
            const NCPA::units::VectorWithUnits<T>& ds, T freq,
            const units_ptr_t units_out
            = NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS ) {
            NCPA::units::VectorWithUnits<T> alpha( zs.size(), "np/m" );
            for (size_t i = 0; i < zs.size(); i++) {
                alpha.set( i, attenuation_sutherland_bass(
                                  zs.get_scalar( i ), ts.get_scalar( i ),
                                  ps.get_scalar( i ), ds.get_scalar( i ),
                                  freq ) );
            }
            return alpha.as( units_out );
        }

        template<typename T = double>
        static NCPA::units::VectorWithUnits<T> attenuation_sutherland_bass(
            std::vector<T>& z, const std::string& z_units, std::vector<T>& t,
            const std::string& t_units, std::vector<T>& p,
            const std::string& p_units, std::vector<T>& d,
            const std::string& d_units, T freq,
            const std::string& units_out
            = NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS_STR ) {
            return attenuation_sutherland_bass(
                NCPA::units::VectorWithUnits<T>( z, z_units ),
                NCPA::units::VectorWithUnits<T>( t, t_units ),
                NCPA::units::VectorWithUnits<T>( p, p_units ),
                NCPA::units::VectorWithUnits<T>( d, d_units ), freq,
                NCPA::units::Units::from_string( units_out ) );
        }

        template<typename T = double>
        static std::vector<T> attenuation_sutherland_bass(
            const NCPA::arrays::ndvector<1, T>& z, const units_ptr_t z_units,
            const NCPA::arrays::ndvector<1, T>& t, const units_ptr_t t_units,
            const NCPA::arrays::ndvector<1, T>& p, const units_ptr_t p_units,
            const NCPA::arrays::ndvector<1, T>& d, const units_ptr_t d_units,
            T freq,
            const units_ptr_t units_out
            = NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS ) {
            std::vector<T> alpha( z.size() );
            for (size_t i = 0; i < z.size(); ++i) {
                alpha[ i ]
                    = attenuation_sutherland_bass(
                          NCPA::units::ScalarWithUnits<T>( z[ i ], z_units ),
                          NCPA::units::ScalarWithUnits<T>( t[ i ], t_units ),
                          NCPA::units::ScalarWithUnits<T>( p[ i ], p_units ),
                          NCPA::units::ScalarWithUnits<T>( d[ i ], d_units ),
                          freq, units_out )
                          .get();
            }
            return alpha;
        }

        template<typename T = double>
        static T attenuation_sutherland_bass(
            const NCPA::arrays::ndvector<0, T>& z, const units_ptr_t z_units,
            const NCPA::arrays::ndvector<0, T>& t, const units_ptr_t t_units,
            const NCPA::arrays::ndvector<0, T>& p, const units_ptr_t p_units,
            const NCPA::arrays::ndvector<0, T>& d, const units_ptr_t d_units,
            T freq,
            const units_ptr_t units_out
            = NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS ) {
            return attenuation_sutherland_bass(
                       NCPA::units::ScalarWithUnits<T>( z, z_units ),
                       NCPA::units::ScalarWithUnits<T>( t, t_units ),
                       NCPA::units::ScalarWithUnits<T>( p, p_units ),
                       NCPA::units::ScalarWithUnits<T>( d, d_units ), freq,
                       units_out )
                .get();
        }

        template<size_t N, typename T = double>
        NCPA::arrays::ndvector<N, T> attenuation_sutherland_bass(
            const NCPA::arrays::ndvector<N, T>& z, const units_ptr_t z_units,
            const NCPA::arrays::ndvector<N, T>& t, const units_ptr_t t_units,
            const NCPA::arrays::ndvector<N, T>& p, const units_ptr_t p_units,
            const NCPA::arrays::ndvector<N, T>& d, const units_ptr_t d_units,
            T freq,
            const units_ptr_t units_out
            = NCPA_ATMOS_SUTHERLAND_BASS_ATTENUATION_ALPHA_UNITS ) {
            NCPA::arrays::ndvector<N, T> alpha = z;

            for (size_t i = 0; i < z.size(); ++i) {
                alpha[ i ] = attenuation_sutherland_bass(
                    NCPA::arrays::ndvector<N - 1, T>( z[ i ] ), z_units,
                    NCPA::arrays::ndvector<N - 1, T>( t[ i ] ), t_units,
                    NCPA::arrays::ndvector<N - 1, T>( p[ i ] ), p_units,
                    NCPA::arrays::ndvector<N - 1, T>( d[ i ] ), d_units, freq,
                    units_out );
            }
            return alpha;
        }


    }  // namespace atmos
}  // namespace NCPA
