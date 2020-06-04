/*
 * MIT License
 *
 * Copyright © 2017 S5Lab
 * Created by Leonardo Parisi (leonardo.parisi[at]gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
