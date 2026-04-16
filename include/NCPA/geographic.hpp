#pragma once

#include "NCPA/constants.hpp"
#include "NCPA/math.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

namespace NCPA {
    namespace geographic {

        // forward declarations
        template<typename T = double>
        class Location;
        double azimuth( double lat1_deg, double lon1_deg, double lat2_deg,
                        double lon2_deg );
        double backazimuth( double lat1_deg, double lon1_deg, double lat2_deg,
                            double lon2_deg );
        double deg2km( double angular );
        std::string dms_str( double coord );
        constexpr double earthradius();
        double earthradius( double lat );
        void great_circle( double startlat, double startlon, double azimuth,
                           double range_km, int length, double *pathlat,
                           double *pathlon );
        void great_circle( double startlat, double startlon, double azimuth,
                           double range_km, int length,
                           std::vector<double>& pathlat,
                           std::vector<double>& pathlon );
        void intersection( double lat1_deg, double lon1_deg,
                           double bearing1_deg, double lat2_deg,
                           double lon2_deg, double bearing2_deg,
                           double& lat3_deg, double& lon3_deg );
        double km2deg( double angular );
        void latlon2xy( size_t npts, double *lat, double *lon, double *x,
                        double *y, double reflat );
        void latlon2xy( const std::vector<double>& lat,
                        const std::vector<double>& lon, std::vector<double>& x,
                        std::vector<double>& y, double reflat );
        template<typename T>
        T normalize_azimuth( T in, bool radians = false );
        double normalizeLon( double lon );
        /*
         * Computes the range along the surface of the earth using the
         * Haversine formula
         */
        double range( double lat1_deg, double lon1_deg, double lat2_deg,
                      double lon2_deg );
        Location<double> xy2latlon( double x, double y, double lat0,
                                    double lon0 );
        void xy2latlon( double x, double y, double lat0, double lon0,
                        double& newlat, double& newlon );

        template<typename T>
        class Location {
            protected:
                T _lat, _lon, _elev;

            public:
                Location() : Location( -999.0, -999.0, 0.0 ) {}

                Location( T lat, T lon, T elev = 0 ) {
                    _lat  = lat;
                    _lon  = lon;
                    _elev = elev;
                }

                bool operator==( const Location& other ) const {
                    if (_lat == other.lat() && _lon == other.lon()
                        && _elev == other.elev())
                        return true;
                    return false;
                }

                bool operator!=( const Location& other ) const {
                    return !( ( *this ) == other );
                }

                bool operator<( const Location& other ) const {
                    if (_lat < other.lat())
                        return true;
                    else if (_lat > other.lat())
                        return false;
                    else if (_lon < other.lon())
                        return true;
                    else
                        return false;
                }

                bool operator>( const Location& other ) const {
                    if (_lat > other.lat())
                        return true;
                    else if (_lat < other.lat())
                        return false;
                    else if (_lon > other.lon())
                        return true;
                    else
                        return false;
                }

                T& lat() { return _lat; }

                const T& lat() const { return _lat; }

                T& lon() { return _lon; }

                const T& lon() const { return _lon; }

                T& elev() { return _elev; }

                const T& elev() const { return _elev; }

                void setLat( T newlat ) {
                    if (std::fabs( newlat ) > 90.0) {
                        throw std::invalid_argument(
                            "Latitude must be in the range [-90,90]!" );
                    }
                    _lat = newlat;
                }

                void setLon( T newlon ) {
                    if (std::fabs( newlon ) > 180.0) {
                        throw std::invalid_argument(
                            "Longitude must be in the range [-180,180]!" );
                    }
                    _lon = newlon;
                }

                void setElev( T newelev ) { _elev = newelev; }
        };

        inline double azimuth( double lat1_deg, double lon1_deg,
                               double lat2_deg, double lon2_deg ) {
            double dlon = NCPA::math::deg2rad( lon2_deg - lon1_deg );
            double lat1 = NCPA::math::deg2rad( lat1_deg );
            double lat2 = NCPA::math::deg2rad( lat2_deg );
            double y    = std::sin( dlon ) * std::cos( lat2 );
            double x = std::cos( lat1 ) * std::sin( lat2 )
                     - std::sin( lat1 ) * std::cos( lat2 ) * std::cos( dlon );
            double bearing = NCPA::math::rad2deg( std::atan2( y, x ) );
            return normalize_azimuth( bearing );
        }

        inline double backazimuth( double lat1_deg, double lon1_deg,
                                   double lat2_deg, double lon2_deg ) {
            return azimuth( lat2_deg, lon2_deg, lat1_deg, lon1_deg );
        }

        inline double deg2km( double angular ) {
            return 2.0 * NCPA::constants::PI * earthradius() / 360.0 * angular;
        }

        inline std::string dms_str( double coord ) {
            int sign  = coord >= 0 ? 1 : -1;
            coord     = std::fabs( coord );
            int degs  = (int)( std::floor( coord ) );
            coord    -= degs;
            coord    *= 60.0;
            int mins  = (int)( std::floor( coord ) );
            coord    -= mins;
            coord    *= 60.0;
            degs     *= sign;

            char buffer[ 256 ];
            std::sprintf( buffer, "%dd%02d'%06.3f\"", degs, mins, coord );
            std::string s( buffer );
            return s;
        }

        inline constexpr double earthradius() {
            return 6371.009;
        }

        inline double earthradius( double lat ) {
            double a = 6378.137;
            double b = 6356.7523;
            lat      = NCPA::math::deg2rad( lat );
            double numerator
                = a * a * a * a * std::cos( lat ) * std::cos( lat )
                + b * b * b * b * std::sin( lat ) * std::sin( lat );
            double denominator
                = a * a * std::cos( lat ) + b * b * std::sin( lat );
            return std::sqrt( numerator / denominator );
        }

        inline void great_circle( double startlat, double startlon,
                                  double azimuth, double range_km, int length,
                                  double *pathlat, double *pathlon ) {
            double lat1  = NCPA::math::deg2rad( startlat );
            double lon1  = NCPA::math::deg2rad( startlon );
            // double d = range_km / 6371.0;
            double brng  = NCPA::math::deg2rad( azimuth );
            pathlat[ 0 ] = startlat;
            pathlon[ 0 ] = startlon;
            for (int i = 1; i < length; i++) {
                double d     = ( range_km / ( length - 1 ) ) * i / 6371.0;
                pathlat[ i ] = NCPA::math::rad2deg( std::asin(
                    std::sin( lat1 ) * std::cos( d )
                    + std::cos( lat1 ) * std::sin( d ) * std::cos( brng ) ) );
                pathlon[ i ] = normalizeLon( NCPA::math::rad2deg(
                    lon1
                    + std::atan2( std::sin( brng ) * std::sin( d )
                                      * std::cos( lat1 ),
                                  std::cos( d )
                                      - std::sin( lat1 )
                                            * std::sin( NCPA::math::deg2rad(
                                                pathlat[ i ] ) ) ) ) );
            }
        }

        inline void great_circle( double startlat, double startlon,
                                  double azimuth, double range_km, int length,
                                  std::vector<double>& pathlat,
                                  std::vector<double>& pathlon ) {
            pathlat.resize( length );
            pathlon.resize( length );
            great_circle( startlat, startlon, azimuth, range_km, length,
                          pathlat.data(), pathlon.data() );
        }

        inline void intersection( double lat1_deg, double lon1_deg,
                                  double bearing1_deg, double lat2_deg,
                                  double lon2_deg, double bearing2_deg,
                                  double& lat3_deg, double& lon3_deg ) {
            double lat1   = NCPA::math::deg2rad( lat1_deg );
            double lat2   = NCPA::math::deg2rad( lat2_deg );
            double lon1   = NCPA::math::deg2rad( lon1_deg );
            double lon2   = NCPA::math::deg2rad( lon2_deg );
            double theta1 = NCPA::math::deg2rad( bearing1_deg );
            double theta2 = NCPA::math::deg2rad( bearing2_deg );
            double dlat   = lat2 - lat1;
            double dlon   = lon2 - lon1;

            double d12
                = 2.0
                * std::asin( std::sqrt(
                    std::sin( dlat / 2 ) * std::sin( dlat / 2 )
                    + std::cos( lat1 ) * std::cos( lat2 )
                          * std::sin( dlon / 2 ) * std::sin( dlon / 2 ) ) );
            double phi1
                = std::acos( std::sin( lat2 )
                             - std::sin( lat1 ) * std::cos( d12 )
                                   / std::sin( d12 ) * std::cos( lat1 ) );
            double phi2
                = std::acos( std::sin( lat1 )
                             - std::sin( lat2 ) * std::cos( d12 )
                                   / std::sin( d12 ) * std::cos( lat2 ) );

            double theta12, theta21;
            if (std::sin( lon2 - lon1 ) > 0) {
                theta12 = phi1;
                theta21 = 2 * NCPA::constants::PI - phi2;
            } else {
                theta12 = 2 * NCPA::constants::PI - phi1;
                theta21 = phi2;
            }

            double alpha1 = std::fmod( theta1 - theta12 + NCPA::constants::PI,
                                       2 * NCPA::constants::PI )
                          - NCPA::constants::PI;
            double alpha2 = std::fmod( theta21 - theta2 + NCPA::constants::PI,
                                       2 * NCPA::constants::PI )
                          - NCPA::constants::PI;
            double alpha3 = std::acos(
                -std::cos( alpha1 ) * std::cos( alpha2 )
                + std::sin( alpha1 ) * std::sin( alpha2 ) * std::cos( d12 ) );
            double d13 = std::atan2(
                std::sin( d12 ) * std::sin( alpha1 ) * std::sin( alpha2 ),
                std::cos( alpha2 ) + std::cos( alpha1 ) * std::cos( alpha3 ) );
            double lat3   = std::asin( std::sin( lat1 ) * std::cos( d13 )
                                       + std::cos( lat1 ) * std::sin( d13 )
                                             * std::cos( theta1 ) );
            lat3_deg      = NCPA::math::rad2deg( lat3 );
            double dlon13 = std::atan2(
                std::sin( theta1 ) * std::sin( d13 ) * std::cos( lat1 ),
                std::cos( d13 ) - std::sin( lat1 ) * std::sin( lat3 ) );
            lon3_deg = NCPA::math::rad2deg(
                std::fmod( lon1 + dlon13 + NCPA::constants::PI,
                           2 * NCPA::constants::PI )
                - NCPA::constants::PI );
        }

        inline double km2deg( double linear ) {
            return linear * 360.0
                 / ( 2.0 * NCPA::constants::PI * earthradius() );
        }

        // convert a series of lat/lon points to absolute x/y points using an
        // equirectangular projection, with reference latitide chosen to be the
        // center of the latitude range
        inline void latlon2xy( size_t npts, double *lat, double *lon,
                               double *x, double *y, double reflat ) {
            if (reflat < -90.0) {
                reflat = 0.5
                       * ( NCPA::math::max( lat, npts )
                           - NCPA::math::min( lat, npts ) );
            }
            double Re = earthradius( reflat );
            for (size_t i = 0; i < npts; i++) {
                x[ i ] = Re * NCPA::math::deg2rad( lon[ i ] )
                       * std::cos( NCPA::math::deg2rad( reflat ) );
                y[ i ] = Re * NCPA::math::deg2rad( lat[ i ] );
            }
        }

        inline void latlon2xy( const std::vector<double>& lat,
                               const std::vector<double>& lon,
                               std::vector<double>& x, std::vector<double>& y,
                               double reflat ) {
            if (lat.size() != lon.size()) {
                throw std::runtime_error(
                    "Mismatched input lat/lon vector lengths!" );
            }

            if (reflat < -90.0) {
                reflat = 0.5
                       * ( NCPA::math::max( lat ) - NCPA::math::min( lat ) );
            }
            double Re = earthradius( reflat );
            x.clear();
            y.clear();
            x.resize( lat.size(), 0.0 );
            y.resize( lat.size(), 0.0 );

            for (size_t i = 0; i < lat.size(); i++) {
                x[ i ] = Re * NCPA::math::deg2rad( lon[ i ] )
                       * std::cos( NCPA::math::deg2rad( reflat ) );
                y[ i ] = Re * NCPA::math::deg2rad( lat[ i ] );
            }
        }

        /**
        @brief Normalizes an azimuth into the range [0,360)
        @param in The azimuth value to normalize, in degrees.
        @returns The normalized azimuth, in degrees.
        */
        template<typename T>
        T normalize_azimuth( T in, bool radians ) {
            double interval = radians ? 2.0 * M_PI : 360.0;
            double out      = (double)in;
            while (out < 0.0) {
                out += interval;
            }
            while (out >= interval) {
                out -= interval;
            }
            return (T)out;
        }

        inline double normalizeLon( double lon ) {
            while (lon > 180.0) {
                lon -= 360.0;
            }
            while (lon <= -180.0) {
                lon += 360.0;
            }
            return lon;
        }

        /*
         * Computes the range along the surface of the earth using the
         * Haversine formula
         */
        inline double range( double lat1_deg, double lon1_deg, double lat2_deg,
                             double lon2_deg ) {
            double dLat = NCPA::math::deg2rad( lat2_deg - lat1_deg );
            double dLon = NCPA::math::deg2rad( lon2_deg - lon1_deg );
            double lat1 = NCPA::math::deg2rad( lat1_deg );
            double lat2 = NCPA::math::deg2rad( lat2_deg );

            double a = std::sin( dLat / 2 ) * std::sin( dLat / 2 )
                     + std::sin( dLon / 2 ) * std::sin( dLon / 2 )
                           * std::cos( lat1 ) * std::cos( lat2 );
            double c = 2.0 * std::atan2( std::sqrt( a ), std::sqrt( 1 - a ) );
            return c * earthradius( 0.5 * lat1_deg + 0.5 * lat2_deg );
        }

        inline void xy2latlon( double x, double y, double lat0, double lon0,
                               double& newlat, double& newlon ) {
            // Simple: go north by y km, then east by x km
            double templat[ 2 ], templon[ 2 ];
            double northaz = 0;
            double eastaz  = 90;
            if (y < 0) {
                northaz = 180;
                y       = -y;
            }
            if (x < 0) {
                eastaz = 270;
                x      = -x;
            }

            great_circle( lat0, lon0, northaz, y, 2, templat, templon );
            double templat2 = templat[ 1 ];
            double templon2 = templon[ 1 ];
            great_circle( templat2, templon2, eastaz, x, 2, templat, templon );
            newlat = templat[ 1 ];
            newlon = templon[ 1 ];
        }

        inline Location<double> xy2latlon( double x, double y, double lat0,
                                           double lon0 ) {
            Location<double> loc;
            xy2latlon( x, y, lat0, lon0, loc.lat(), loc.lon() );
            return loc;
        }
    }  // namespace geographic
}  // namespace NCPA
