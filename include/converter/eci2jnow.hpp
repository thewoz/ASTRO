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

#ifndef _H_ASTRO_J2KTOJNOW_H
#define _H_ASTRO_J2KTOJNOW_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vector>


//****************************************************************************/
// namespace astro
//****************************************************************************/
namespace astro {
  
  /* ******************* Function J2K_JNOW ***************************** * \
   *                                                                     *
   * Transform radec in J2000 to radec in JNOW.                          *
   *                                                                     *
   * Input:                                                              *
   * jd (taking into account also jdfrac)                                *
   * ra,dec in J2000                                     [deg]           *
   *                                                                     *
   * Ouput:                                                              *
   * ranow,decnow in JNOW                                [deg]           *
   *                                                                     *
   * ******************************************************************* */
  void radecEci2jnow(double raEci, double decEci, double & raJnow, double & decJnow, double jDay) {
    
    // astronomical algorithms meeus
    double t = (jDay-2451545.0) / 36525.0;
    
    raEci  = astro::radians(raEci);
    decEci = astro::radians(decEci);

    double zita  = astro::radians((2306.2181*t+0.30188*t*t+0.017998*t*t*t)/3600.0);
    double zeta  = astro::radians((2306.2181*t+1.09468*t*t+0.018203*t*t*t)/3600.0);
    double theta = astro::radians((2004.3109*t-0.42665*t*t-0.041833*t*t*t)/3600.0);
    
    double A = cos(decEci)*sin(raEci+zita);
    double B = cos(theta)*cos(decEci)*cos(raEci+zita)-sin(theta)*sin(decEci);
    double C = sin(theta)*cos(decEci)*cos(raEci+zita)+cos(theta)*sin(decEci);
    
    raJnow  = astro::degrees(atan2(A,B) + zeta);
    decJnow = astro::degrees(asin(C));

  }

} /* namespace astro */

#endif /* _H_ASTRO_J2KTOJNOW_H
*/
