/*
 * MIT License
 *
 * Copyright © 2017 S5Lab
 * Created by Leonardo Parisi (leonardo.parisi[at]gmail.com)
 * Modified by Gaetano Zarcone, Lorenzo Mariani
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

#include <cstdio>
#include <cstdlib>
#include <cmath>

//****************************************************************************/
// namespace astro
//****************************************************************************/
namespace astro {
  
  //   % Programed by Darin C. Koblick 01/23/2010
  //   %-------------------------------------------------------------------------
  //   % Function Description:
  //   %--------------------------------------------------------------------------
  //   % RaDec2AzEl will take the Right Ascension and Declination in the topocentric
  //   % reference frame, site latitude and longitude as well as a time in GMT
  //   % and output the Azimuth and Elevation in the local horizon
  //   % reference frame.
  //   %
  //   % Inputs:                                                       Format:
  //   %--------------------------------------------------------------------------
  //   % Topocentric Right Ascension (Degrees)                         [N x 1]
  //   % Topocentric Declination Angle (Degrees)                       [N x 1]
  //   % Lat (Site Latitude in degrees -90:90 -> S(-) N(+))            [N x 1]
  //   % Lon (Site Longitude in degrees -180:180 W(-) E(+))            [N x 1]
  //   % Jday (Modifica Fabio)
  //   %
  //   % Outputs:                                                      Format:
  //   %--------------------------------------------------------------------------
  //   % Local Azimuth Angle   (degrees)                               [N x 1]
  //   % Local Elevation Angle (degrees)                               [N x 1]
  //   %
  //   %
  //   % External Source References:
  //   % Fundamentals of Astrodynamics and Applications
  //   % D. Vallado, Second Edition
  //   % Example 3-5. Finding Local Siderial Time (pg. 192)
  //   % Algorithm 28: AzElToRaDec (pg. 259)
  //   %-------------------------------------------------------------------------
  void tRaDec2AzEl(double ra, double dec, double lat, double lon, double JD, double & Az, double & El) {
    
    
    double lat_rad = astro::radians(lat);
    double dec_rad = astro::radians(dec);
    
    double T_UT1 = (JD-2451545)/36525.0;
    double ThetaGMST = 67310.54841 + (876600.0*3600.0 + 8640184.812866)*T_UT1 + 0.093104*(pow(T_UT1,2)) -(6.2*pow(10.0,-6))*(pow(T_UT1,3)); //deg
    
    // funzione modulo fatta da noi
    double a = ThetaGMST;
    double b = 86400.0*(ThetaGMST/abs(ThetaGMST));
    
    double mod1 = (a - floor(a/b)*b)/240.0;
    ThetaGMST = mod1 - floor(mod1/360.0)*360.0;
    double ThetaLST = ThetaGMST + lon; //deg
    
    //Define Siderial Time LHA)
    double LHA = (ThetaLST-ra) - floor((ThetaLST-ra)/360.0)*360.0; //funzione modulo [deg]
    double LHA_rad = astro::radians(LHA);
    
    //Elevation deg
    double El_rad = asin(sin(lat_rad)*sin(dec_rad)+cos(lat_rad)*cos(dec_rad)*cos(LHA_rad));
    El = astro::degrees(El_rad);
    
    //Azimuth deg
    double c_rad = atan2(-sin(LHA_rad)*cos(dec_rad)/cos(El_rad),(sin(dec_rad)-sin(El_rad)*sin(lat_rad))/(cos(El_rad)*cos(lat_rad)));
    double c     = astro::degrees(c_rad);
    Az = c - floor(c/360.0)*360.0; //deg
    
  }
  
  
  /*------------------------------------------------------------------------------
   *
   *                           tradec_azel
   *
   * this function converts topocentric right ascension declination values with
   * azimuth, and elevation.
   *
   *
   *  inputs          description                                units
   *    rtasc       - topocentric right ascension                rad
   *    decl        - topocentric declination                    rad
   *    latgd       - geodetic latitude                          rad
   *    JD          - julian date
   *    lon         - site longitude                             deg
   *
   *
   *  outputs       :
   *    az          - azimuth                                    deg
   *    el          - elevation                                  deg
   *
   *  use vallado function radec_azel
   -----------------------------------------------------------------------------*/
  void tradec_azel(double & rtasc, double & decl, double & latgd, double & JD, double & lon,double & az, double & el) {
    
    //aggiunta questa parte oer calcolare lst, quindi lst è stato tolto com input e sono stati aggiunti jday e longitudine in gradi
    double T_UT1 = (JD-2451545)/36525.0;
    double ThetaGMST = 67310.54841 + (876600.0*3600.0 + 8640184.812866)*T_UT1 + 0.093104*(pow(T_UT1,2)) -(6.2*pow(10.0,-6))*(pow(T_UT1,3)); //deg
    
    // funzione modulo fatta da noi
    double a = ThetaGMST;
    double b = 86400.0*(ThetaGMST/abs(ThetaGMST));
    
    double mod1 = (a - floor(a/b)*b)/240.0;
    ThetaGMST = mod1 - floor(mod1/360.0)*360.0;
    double lst = astro::radians(ThetaGMST + lon); //rad
    
    astIOD::radec_azel(rtasc, decl, lst, latgd, eTo, az, el);
    
    az = astro::degrees(az);
    el = astro::degrees(el);
    
  } // procedure radec_azel
  
  
  /*------------------------------------------------------------------------------
   *
   *                           rv_tradec
   *
   *  this procedure converts topocentric right-ascension declination with
   *  position and velocity vectors. uses velocity vector to find the
   *  solution of singular cases.
   *
   *  inputs          description                    range / units
   *    rijk        - ijk position vector            km
   *    vijk        - ijk velocity vector            km/s
   *    rsijk       - ijk site position vector       km
   *
   *  outputs       :
   *    trtasc      - top right ascension            rad
   *    tdecl       - top declination                rad
   *
   -----------------------------------------------------------------------------*/
  void rv_tradec(double rijk[3], double vijk[3], double rsijk[3], edirection direct, double & trtasc, double & tdecl) {
    
    double rho, drho, dtrtasc, dtdecl;
    
    astIOD::rv_tradec(rijk,vijk,rsijk, eTo, rho, trtasc, tdecl, drho, dtrtasc, dtdecl);
    
    if(trtasc < 0.0) trtasc += 2.0 * M_PI;
    
    trtasc = astro::degrees(trtasc);
    tdecl  = astro::degrees(tdecl);
    
  }
  
} /* namespace astro */

#endif /* _H_ASTRO_CONVERTER_H */
