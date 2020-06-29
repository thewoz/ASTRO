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

#ifndef _H_ASTRO_UTILS_H_
#define _H_ASTRO_UTILS_H_

#include <cstdlib>
#include <cstdio>

#include <vector>

//****************************************************************************/
// namespace astro::utils
//****************************************************************************/
namespace astro::utils {
  
  //****************************************************************************/
  // elevazione
  //****************************************************************************/
  double elevazione(double observatory[3], double satellite[3]) {
    
    satellite[0]   *= 6371.0; satellite[1]   *= 6371.0; satellite[2]   *= 6371.0;
    observatory[0] *= 6371.0; observatory[1] *= 6371.0; observatory[2] *= 6371.0;
    
    double dist[3];
    dist[0] = observatory[0] - satellite[0];
    dist[1] = observatory[1] - satellite[1];
    dist[2] = observatory[2] - satellite[2];
    
    return  astro::degrees(acos(astMath::dot(observatory,dist) / (astMath::mag(observatory)*astMath::mag(dist)))) - 90.0f;
    
  }
  
  //****************************************************************************/
  // teta
  //****************************************************************************/
  double teta(double observatory[3], double satellite[3]) {
    
    satellite[0]   *= 6371.0; satellite[1]   *= 6371.0; satellite[2]   *= 6371.0;
    observatory[0] *= 6371.0; observatory[1] *= 6371.0; observatory[2] *= 6371.0;
    
    double dist[3];
    dist[0] = observatory[0] - satellite[0];
    dist[1] = observatory[1] - satellite[1];
    dist[2] = observatory[2] - satellite[2];
    
    float elevazione = astro::degrees(acos(astMath::dot(observatory,dist) / (astMath::mag(observatory)*astMath::mag(dist)))) - 90.0f;
    
    return 90.0 - elevazione;;
    
  }
  
  
  //****************************************************************************/
  // air attenuation factor
  //****************************************************************************/
  double air_attenuation_factor(double observatory[3], double satellite[3]) {
    
    satellite[0]   *= 6371.0; satellite[1]   *= 6371.0; satellite[2]   *= 6371.0;
    observatory[0] *= 6371.0; observatory[1] *= 6371.0; observatory[2] *= 6371.0;
    
    double dist[3];
    dist[0] = observatory[0] - satellite[0];
    dist[1] = observatory[1] - satellite[1];
    dist[2] = observatory[2] - satellite[2];

    float elevazione = astro::degrees(acos(astMath::dot(observatory,dist) / (astMath::mag(observatory)*astMath::mag(dist)))) - 90.0f;
    
    // teta è la distanza dallo zenith
    float teta = 90.0 - elevazione;
    
    float secd = 1.0 / cos(astro::radians(teta));
    
    // valore dell’air mass alla data teta
    float am = secd-0.0018167*(secd-1)-0.0002875*pow((secd-1),2)-0.0008083*pow((secd-1),3);
    
    // assorbività dell'atmosfera
    float alfa = 0.25;
    float Q0 = 1.0;
    
    // la frazione di luce che arriva rispetto ad 1
    float Q = Q0 * exp(-alfa*am);
    
    return Q;
    
  }
  
  
  /*****************************************************************************/
  // Obliquity of ecliptic
  //
  // Obliquety of ecliptic eps (in radians)
  // at time t (in Julian Centuries from J2000.0)
  // Taken from the KDE AStroLib
  /*****************************************************************************/
  double eps(double t) {

    double tp;
    
    tp = 23.43929111 - (46.815+(0.00059-0.001813*t)*t)*t/3600.0;
    tp = 1.74532925199e-2 * tp;
    
    return tp;
    
  }

  /*****************************************************************************/
  //  Ecliptic to Equatorial
  //
  // Convert position vector r1 from ecliptic into equatorial coordinates
  // at t (in Julian Centuries from J2000.0)
  /*****************************************************************************/
  void eclequ(double t, double * rIn, double * rOut) {
    
    std::vector<std::vector<double>> m;
    
    astMath::rot1mat(-eps(t), m); // m = xrot (-eps(t));
    
    astMath::matvecmult(m, rIn, rOut); // r2 = mxvct(m, r1);
    
  }

  /* ****************** Function refraction ******************* *\
   *                                                            *
   * this function compute the correction to true elevation     *
   * due to light refraction                                    *
   * it should be used only for positive elevation              *
   * input:                                                     *
   *       el rad                                               *
   * output:                                                    *
   *     appel rad                                              *
   * we use formula 16.4 of Astronomical Algorithms p.106       *
   * we do not add the correction to recover no refraction for  *
   * elevation=90deg to be consistent with previsat             *
   *                                                            *
   \* ********************************************************** */
  void refraction(double el, double& appel) {
    
    const double rad2deg = 180.0 / M_PI; //rad2deg
    double eltmp;
    
    if(el>0.0) {
        eltmp = el*rad2deg;
        eltmp = (eltmp+10.3/(eltmp+5.11));
	appel = el + 1.02/(60.*tan(eltmp/rad2deg))/rad2deg;
        if(appel<0.0) appel = el;
      }
    
  }

/* ********************* Back-refraction ********************* *\
 *                                                             *
 * this function obtains the true elevation from the apparent  *
 * elevation, inverting the previous formula by means of       *
 * a bisection method.                                         *
 * input: appel (rad)   ->   output: el (rad)                  *  
 *                                                             *
 \* *********************************************************** */
  void back_refraction(double appel, double& el) {
    
    //const double rad2deg = 180.0 / M_PI;
    double el1, el2, elm=0, del, appel1, f;
    //double eltmp;
    
    if(appel>0.0){
      el1=0.0; el2=M_PI/2.;
      del=el2-el1;
      
      while (del > 1.e-8){
        elm=(el2+el1)/2.;
        refraction(elm, appel1);
        f = appel1 - appel;
        if (f < 0.0){
          el1=elm;
        }
        else
        {
          el2=elm;
        }
        del=el2-el1;
      }
      
      el = elm;
      
    }
    
  }

} /* namespace astro */


#endif /* _H_ASTRO_UTILS_H_ */
