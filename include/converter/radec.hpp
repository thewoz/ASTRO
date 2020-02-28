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

#ifndef _H_ASTRO_RADEC_H
#define _H_ASTRO_RADEC_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
  
  /* ********************* Function radecGeo2Topo ********************** * \
   *                                                                     *
   * Compute the topocentric RA and DEC                                  *
   *                                                                     *
   * Input:                                                              *
   * satCoord satellite coordinate                        [km]           *
   * objCoord object coordinate                           [km]           *
   *                                                                     *
   * NOTE: the satCoord and objCoord muat be in the same system          *
   *                                                                     *
   * Ouput:                                                              *
   * topocentric RA DEC                                  [deg]           *
   *                                                                     *
   * ******************************************************************* */
  void radecGeo2Topo(double satCoord[3], double objCoord[3], double & ra, double $ dec) {
    
//    rho = Stato_sat(:,1:3)-obs_pos;

//    normRHO = sqrt(dot(rho',rho'))'; 

//    Cel_Coord = [atan2(rho(:,2),rho(:,1)), asin(rho(:,3)./normRHO)];

//    Cel_Coord(Cel_Coord(:,1)<0,1) = Cel_Coord(Cel_Coord(:,1)<0,1)+2*pi;

//    Cel_Coord = rad2deg(Cel_Coord);

  }
  

} /* namespace astro */

#endif /* _H_ASTRO_RADEC_H
*/
