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

#ifndef _H_ASTRO_CONVERTER_H
#define _H_ASTRO_CONVERTER_H

#include "./converter/eopc.hpp"

#include "./converter/ecef.hpp"
#include "./converter/teme.hpp"
#include "./converter/eci.hpp"

#include "./converter/lla2ecef.hpp"

#include "./converter/eci2jnow.hpp"

namespace astro {
  
/*****************************************************************************
 %  this function converts the position of on object in right ascension and declination values.
 % uses velocity vector to find the solution of singular cases.
 %
 %  inputs:       description                    range / units
 %    r           -  position vector             km
 %    v           -  velocity vector             km/s
 %
 %  outputs:
 %    ra         - right ascension               rad
 %    dec        - declination                   rad
 %
 *****************************************************************************/
  void rv2radec(const double r[3], const double v[3], double & ra, double & dec) {
    
    const double small = 0.00000001;
    
    double temp = sqrt(r[0]*r[0] + r[1]*r[1]);
    
    if(temp < small) {
      
      double temp1 = sqrt(v[0]*v[0] + v[1]*v[1]);
      
      ra = 0.0;
      
      if(fabs(temp1) > small) ra = atan2(v[1]/temp1, v[0]/temp1);
      
    } else { ra = atan2(r[1]/temp, r[0]/temp); }
    
    dec = asin(r[2] / astMath::mag(r));
    
  }
  
}

/*****************************************************************************/
 // converte da ra e dec geocentrice a ra e dec topocentriche
 // le ra e dec devono essere in gradi
 // site invece le cordinate del sito in TEME
 // sat invece le cordinate del satellite in TEME
 /*****************************************************************************/
 void ra2tradec(double sat[3],  double tra, double tdec, double site[3]) {
   
   double rho[3];
   
   rho[0] = sat[0] - site[0];
   rho[1] = sat[1] - site[1];
   rho[2] = sat[2] - site[2];
   
   double normRho = astMath::mag(rho);
   
   double coord[2];
   
   coord[0] = atan2(rho[1],rho[1]);
   coord[1] = asin(rho[2] / normRho);

   if(coord[0] < 0.1) coord[0] += 2 * M_PI;
   
   tra  = Degrees(coord[0]);
   tdec = Degrees(coord[1]);

 }


#endif /* _H_ASTRO_CONVERTER_H */
