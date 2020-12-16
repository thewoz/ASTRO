/*
 * GNU GENERAL PUBLIC LICENSE
 *
 * Copyright (C) 2017
 * Created by Leonardo Parisi (leonardo.parisi[at]gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef _H_ASTRO_LLA2ECEF_H
#define _H_ASTRO_LLA2ECEF_H

#include <cstdlib>
#include <cmath>

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
  
  /* -----------------------------------------------------------------------------
   *
   * LLA2ECEF - convert latitude, longitude, and altitude to
   *            earth-centered, earth-fixed (ECEF) cartesian
   *
   * USAGE:
   * [x,y,z] = lla2ecef(lat,lon,alt)
   *
   * x = ECEF X-coordinate (km)
   * y = ECEF Y-coordinate (km)
   * z = ECEF Z-coordinate (km)
   * lat = geodetic latitude (radians)
   * lon = longitude (radians)
   * alt = height above WGS84 ellipsoid (m)
   *
   * Notes: This function assumes the WGS84 model.
   *        Latitude is customary geodetic (not geocentric).
   *
   * Source: "Department of Defense World Geodetic System 1984"
   *         Page 4-4
   *         National Imagery and Mapping Agency
   *         Last updated June, 2004
   *         NIMA TR8350.2
   *
   * Michael Kleder, July 2005
   * -----------------------------------------------------------------------------*/
  
  
  /*****************************************************************************/
  // lla2ecef
  /*****************************************************************************/
  void lla2ecef(double lat, double lon, double alt, double recef[3]) {
    
    // WGS84 Model constants:
    double a  = 6378137;               // Mean Equatorial radius [m]
    double f  = 1.0 / 298.257223563;   // Flattering factor
    double e  = sqrt(f*(2-f));         // Eccentricity
    double esq = e*e;                  // Powed Eccentricity
    
    // Covert deg to radius
    lat = astro::radians(lat);
    lon = astro::radians(lon);

    // intermediate calculation
    // (prime vertical radius of curvature)
    double N = a / sqrt(1 - esq * pow(sin(lat),2));
    
    // results:
    recef[0] = ((N+alt)  * cos(lat) * cos(lon)) / 1000.0;
    recef[1] = ((N+alt)  * cos(lat) * sin(lon)) / 1000.0;
    recef[2] = (((1-esq) * N + alt) * sin(lat)) / 1000.0;
    
  }
  
  /*****************************************************************************/
  // lla2ecef
  /*****************************************************************************/
  template <class T>
  void lla2ecef(double lat, double lon, double alt, T & recef) {
    
    // WGS84 Model constants:
    double a  = 6378137;               // Mean Equatorial radius [m]
    double f  = 1.0 / 298.257223563;   // Flattering factor
    double e  = sqrt(f*(2-f));         // Eccentricity
    double esq = e*e;                  // Powed Eccentricity
    
    // intermediate calculation
    // (prime vertical radius of curvature)
    double N = a / sqrt(1 - esq * pow(sin(lat),2));
    
    // results:
    recef.x = ((N+alt)  * cos(lat) * cos(lon)) / 1000.0;
    recef.y = ((N+alt)  * cos(lat) * sin(lon)) / 1000.0;
    recef.z = (((1-esq) * N + alt) * sin(lat)) / 1000.0;
    
  }
  
} /* namespace astro */

#endif /* _H_ASTRO_LLA2ECEF_H */
