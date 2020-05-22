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

#include "./converter/eopc.hpp"

#include "./converter/ecef.hpp"
#include "./converter/teme.hpp"
#include "./converter/eci.hpp"

#include "./converter/lla2ecef.hpp"

#include "./converter/eci2jnow.hpp"

#include <iostream>
#include <vector>

using namespace std;

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

  ////// CI HAI PROVATO LEO, MA PURTROPPO SEI SOLO UN INFORMATICO!!!!!!!!!!!!!!!!!!!
  /*
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
  */
}


void RaDec2AzEl(double ra, double dec, double lat, double lon, double JD, double &Az, double &El)
{
  /*
  % Modifica Fabio: input tempo in jd
  %
  % Programed by Darin C. Koblick 01/23/2010
  %--------------------------------------------------------------------------
  % External Function Call Sequence:
  % [Az El] = RaDec2AzEl(0,0,0,-104,'1992/08/20 12:14:00')
  %
  % Worked Example: pg. 262 Vallado
  %[Az El] = RaDec2AzEl(294.9891115,-20.8235624,39.007,-104.883,'1994/05/14 13:11:20.59856')
  %[210.7514  23.9036] = RaDec2AzEl(294.9891115,-20.8235624,39.007,-104.883,'1994/05/14 13:11:20.59856')
  %
  % Worked Example: http://www.stargazing.net/kepler/altaz.html
  % [Az El] = RaDec2AzEl(344.95,42.71667,52.5,-1.91667,'1997/03/14 19:00:00')
  % [311.92258 22.40100] = RaDec2AzEl(344.95,42.71667,52.5,-1.91667,'1997/03/14 19:00:00')
  %
  % [Beta,el] = RaDec2AzEl(alpha_t,delta_t,phi,lamda,'yyyy/mm/dd hh:mm:ss')
  %
  % Function Description:
  %--------------------------------------------------------------------------
  % RaDec2AzEl will take the Right Ascension and Declination in the topocentric 
  % reference frame, site latitude and longitude as well as a time in GMT
  % and output the Azimuth and Elevation in the local horizon
  % reference frame.
  %
  % Inputs:                                                       Format:
  %--------------------------------------------------------------------------
  % Topocentric Right Ascension (Degrees)                         [N x 1]
  % Topocentric Declination Angle (Degrees)                       [N x 1]
  % Lat (Site Latitude in degrees -90:90 -> S(-) N(+))            [N x 1]
  % Lon (Site Longitude in degrees -180:180 W(-) E(+))            [N x 1]
  % UTC (Coordinated Universal Time YYYY/MM/DD hh:mm:ss)          [N x 1]
  %
  % Outputs:                                                      Format:
  %--------------------------------------------------------------------------
  % Local Azimuth Angle   (degrees)                               [N x 1]
  % Local Elevation Angle (degrees)                               [N x 1]
  %
  %
  % External Source References:
  % Fundamentals of Astrodynamics and Applications 
  % D. Vallado, Second Edition
  % Example 3-5. Finding Local Siderial Time (pg. 192) 
  % Algorithm 28: AzElToRaDec (pg. 259)
  % -------------------------------------------------------------------------

  %Example 3-5
  %[yyyy mm dd HH MM SS] = datevec(datenum(time,'yyyy/mm/dd HH:MM:SS'));
  %JD = juliandate(yyyy,mm,dd,HH,MM,SS);
  */

  double lat_rad = astro::Radians(lat);
  double lon_rad = astro::Radians(lon);
  double ra_rad  = astro::Radians(ra);
  double dec_rad = astro::Radians(dec);

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
  double LHA_rad = astro::Radians(LHA);
  //Elevation deg
  double El_rad = asin(sin(lat_rad)*sin(dec_rad)+cos(lat_rad)*cos(dec_rad)*cos(LHA_rad));
  El = astro::Degrees(El_rad);

  //Azimuth deg
  double c_rad = atan2(-sin(LHA_rad)*cos(dec_rad)/cos(El_rad),(sin(dec_rad)-sin(El_rad)*sin(lat_rad))/(cos(El_rad)*cos(lat_rad)));
  double c     = astro::Degrees(c_rad);
  Az = c - floor(c/360.0)*360.0; //deg

  //vector<double> AzEl;

  //AzEl.push_back(Az);
  //AzEl.push_back(El);

  //return AzEl;
}

/*------------------------------------------------------------------------------
   *
   *                           procedure radec_azel
   *
   * this procedure converts right ascension declination values with
   *   azimuth, and elevation.  notice the range is not defined because
   *   right ascension declination only allows a unit vector to be formed.
   *
   *  author        : david vallado                  719-573-2600   22 jun 2002
   *
   *  inputs          description                    range / units
   *    rtasc       - right ascension                0.0 to 2pi rad
   *    decl        - declination                    -pi/2 to pi/2 rad
   *    lst         - local sidedouble time            -2pi to 2pi rad
   *    latgd       - geodetic latitude              -pi/2 to pi/2 rad
   *    direct      -  direction to convert          eFrom  eTo
   *
   *  outputs       :
   *    az          - azimuth                        0.0 to 2pi rad
   *    el          - elevation                      -pi/2 to pi/2 rad
   *
   *  locals        :
   *    lha         - local hour angle               -2pi to 2pi rad
   *    sinv        - sine value
   *    cosv        - cosine value
   *
   *  coupling      :
   *    arcsin      - arc sine function
   *    atan2       - arc tangent function that resolves quadrant ambiguites
   *
   *  references    :
   *    vallado       2013, 265, alg 27
   -----------------------------------------------------------------------------*/

  void radec_azel
  (
   double& rtasc, double& decl, double& latgd, double &JD, double &lon,
   edirection direct,
   double& az, double& el
   )
  { 

    //aggiunta questa parte oer calcolare lst, quindi lst è stato tolto com input e sono stati aggiunti jday e longitudine in gradi
    double T_UT1 = (JD-2451545)/36525.0;                                                      
    double ThetaGMST = 67310.54841 + (876600.0*3600.0 + 8640184.812866)*T_UT1 + 0.093104*(pow(T_UT1,2)) -(6.2*pow(10.0,-6))*(pow(T_UT1,3)); //deg
    // funzione modulo fatta da noi
    double a = ThetaGMST;
    double b = 86400.0*(ThetaGMST/abs(ThetaGMST));

    double mod1 = (a - floor(a/b)*b)/240.0;
    ThetaGMST = mod1 - floor(mod1/360.0)*360.0;
    double lst = astro::Radians(ThetaGMST + lon); //rad
    double sinv, cosv, lha;
    if (direct == eFrom)
    {
      decl = asin(sin(el) * sin(latgd) + cos(el) * cos(latgd) * cos(az));
      
      sinv = -(sin(az) * cos(el) * cos(latgd)) / (cos(latgd) * cos(decl));
      cosv = (sin(el) - sin(latgd) * sin(decl)) / (cos(latgd) * cos(decl));
      lha = atan2(sinv, cosv);
      rtasc = lst - lha;
      astro::utils::rebox(rtasc);
    }
    else
    {
      lha = lst - rtasc;
      el = asin(sin(decl) * sin(latgd) + cos(decl) * cos(latgd) * cos(lha));
      sinv = -sin(lha) * cos(decl) * cos(latgd) / (cos(el) * cos(latgd));
      cosv = (sin(decl) - sin(el) * sin(latgd)) / (cos(el) * cos(latgd));
      az = atan2(sinv, cosv);
      astro::utils::rebox(az);
    }
    
    //  if (show == 'y')
    //    if (fileout != null)
    //      fprintf(fileout, "%f\n", lha * 180.0 / pi);
  } // procedure radec_azel

  /*------------------------------------------------------------------------------
   *
   *                           procedure rv_tradec
   *
   *  this procedure converts topocentric right-ascension declination with
   *    position and velocity vectors. uses velocity vector to find the
   *    solution of singular cases.
   *
   *  author        : david vallado                  719-573-2600   22 jun 2002
   *
   *  inputs          description                    range / units
   *    rijk        - ijk position vector            er
   *    vijk        - ijk velocity vector            er/tu
   *    rsijk       - ijk site position vector       er
   *    direct      -  direction to convert          eFrom  eTo
   *
   *  outputs       :
   *    rho         - top radius of the sat          er
   *    trtasc      - top right ascension            rad
   *    tdecl       - top declination                rad
   *    drho        - top radius of the sat rate     er/tu
   *    tdrtasc     - top right ascension rate       rad/tu
   *    tddecl      - top declination rate           rad/tu
   *
   *  locals        :
   *    rhov        - ijk range vector from site     er
   *    drhov       - ijk velocity vector from site  er / tu
   *    temp        - temporary extended value
   *    temp1       - temporary extended value
   *    i           - index
   *
   *  coupling      :
   *    astMath::mag         - astMath::magnitude of a vector
   *    atan2       - arc tangent function that resolves the quadrant ambiguities
   *    arcsin      - arc sine function
   *    lncom2      - linear combination of 2 vectors
   *    addvec      - add two vectors
   *    dot         - dot product of two vectors
   *
   *  references    :
   *    vallado       2013, 260, alg 26
   -----------------------------------------------------------------------------*/
  
  void rv_tradec
  (
   double rijk[3], double vijk[3], double rsijk[3],
   edirection direct,
   double& rho, double& trtasc, double& tdecl,
   double& drho, double& dtrtasc, double& dtdecl
   )
  {
    const double small = 0.00000001;
    const double raanearth = 0.05883359221938136;  // earth rot rad/tu
    
    double earthrate[3], rhov[3], drhov[3], vsijk[3];
    double   latgc, temp, temp1;
    
    latgc = asin(rsijk[2] / astMath::mag(rsijk));
    earthrate[0] = 0.0;
    earthrate[1] = 0.0;
    earthrate[2] = raanearth;
    astMath::cross(earthrate, rsijk, vsijk);
    
    if (direct == eFrom)
    {
      /* --------  calculate topocentric vectors ------------------ */
      rhov[0] = (rho * cos(tdecl) * cos(trtasc));
      rhov[1] = (rho * cos(tdecl) * sin(trtasc));
      rhov[2] = (rho * sin(tdecl));
      
      drhov[0] = (drho * cos(tdecl) * cos(trtasc) -
                  rho * sin(tdecl) * cos(trtasc) * dtdecl -
                  rho * cos(tdecl) * sin(trtasc) * dtrtasc);
      drhov[1] = (drho * cos(tdecl) * sin(trtasc) -
                  rho * sin(tdecl) * sin(trtasc) * dtdecl +
                  rho * cos(tdecl) * cos(trtasc) * dtrtasc);
      drhov[2] = (drho * sin(tdecl) + rho * cos(tdecl) * dtdecl);
      
      /* ------ find ijk range vector from site to satellite ------ */
      astMath::addvec(1.0, rhov, 1.0, rsijk, rijk);
      astMath::addvec(1.0, drhov, cos(latgc), vsijk, vijk);
      /*
       if (show == 'y')
       if (fileout != null)
       {
       fprintf(fileout, "rtb %18.7f %18.7f %18.7f %18.7f er\n",
       rhov[1], rhov[2], rhov[3], astMath::mag(rhov));
       fprintf(fileout, "vtb %18.7f %18.7f %18.7f %18.7f\n",
       drhov[1], drhov[2], drhov[3], astMath::mag(drhov));
       }
       */
    }
    else
    {
      /* ------ find ijk range vector from site to satellite ------ */
      astMath::addvec(1.0, rijk, -1.0, rsijk, rhov);
      astMath::addvec(1.0, vijk, -cos(latgc), vsijk, drhov);
      
      /* -------- calculate topocentric angle and rate values ----- */
      rho = astMath::mag(rhov);
      temp = sqrt(rhov[0] * rhov[0] + rhov[1] * rhov[1]);
      if (temp < small)
      {
        temp1 = sqrt(drhov[0] * drhov[0] + drhov[1] * drhov[1]);
        trtasc = atan2(drhov[1] / temp1, drhov[0] / temp1);
      }
      else
        trtasc = atan2(rhov[1] / temp, rhov[0] / temp);
      
      tdecl = asin(rhov[2] / astMath::mag(rhov));

      if(trtasc < 0.0) trtasc += 2.0 * M_PI; 

      temp1 = -rhov[1] * rhov[1] - rhov[0] * rhov[0];
      drho = astMath::dot(rhov, drhov) / rho;
      if (fabs(temp1) > small)
        dtrtasc = (drhov[0] * rhov[1] - drhov[1] * rhov[0]) / temp1;
      else
        dtrtasc = 0.0;
      if (fabs(temp) > small)
        dtdecl = (drhov[2] - drho * sin(tdecl)) / temp;
      else
        dtdecl = 0.0;
      
      astro::utils::rebox(dtrtasc);

      trtasc = astro::Degrees(trtasc);
      tdecl  = astro::Degrees(tdecl);
      
      /*
       if (show == 'y')
       if (fileout != null)
       {
       fprintf(fileout, "rta %18.7f %18.7f %18.7f %18.7f er\n",
       rhov[1], rhov[3], rhov[3], astMath::mag(rhov));
       fprintf(fileout, "vta %18.7f %18.7f %18.7f %18.7f er\n",
       drhov[1], drhov[3], drhov[3], astMath::mag(drhov));
       }
       */
    }
  } // rv_tradec
  

/*------------------------------------------------------------------------------
   *
   *                           procedure rvsez_razel
   *
   *  this procedure converts range, azimuth, and elevation values with slant
   *    range and velocity vectors for a satellite from a radar site in the
   *    topocentric horizon (sez) system.
   *
   *  author        : david vallado                  719-573-2600   22 jun 2002
   *
   *  inputs          description                    range / units
   *    rhovec      - sez satellite range vector     km
   *    drhovec     - sez satellite velocity vector  km/s
   *    direct      -  direction to convert          eFrom  eTo
   *
   *  outputs       :
   *    rho         - satellite range from site      mk
   *    az          - azimuth                        0.0 to 2pi rad
   *    el          - elevation                      -pi/2 to pi/2 rad
   *    drho        - range rate                     km/s
   *    daz         - azimuth rate                   rad/s
   *    del         - elevation rate                 rad/s
   *
   *  locals        :
   *    sinel       - variable for sin( el )
   *    cosel       - variable for cos( el )
   *    sinaz       - variable for sin( az )
   *    cosaz       - variable for cos( az )
   *    temp        -
   *    temp1       -
   *
   *  coupling      :
   *    astMath::mag         - astMath::magnitude of a vector
   *    sgn         - returns the sign of a variable
   *    dot         - dot product
   *    arcsin      - arc sine function
   *    atan2       - arc tangent function that resolves quadrant ambiguites
   *
   *  references    :
   *    vallado       2013, 261, eq 4-4, eq 4-5
   -----------------------------------------------------------------------------*/
  
  void rvsez_razel
  (
   double rhosez[3], double drhosez[3],
   edirection direct,
   double& rho, double& az, double& el, double& drho, double& daz, double& del
   )
  {
    const double small = 0.00000001;
    const double halfpi = M_PI / 2.0;
    
    double temp1, temp, sinel, cosel, sinaz, cosaz;
    
    if (direct == eFrom)
    {
      sinel = sin(el);
      cosel = cos(el);
      sinaz = sin(az);
      cosaz = cos(az);
      
      /* ----------------- form sez range vector ------------------ */
      rhosez[0] = (-rho * cosel * cosaz);
      rhosez[1] = (rho * cosel * sinaz);
      rhosez[2] = (rho * sinel);
      
      /* --------------- form sez velocity vector ----------------- */
      drhosez[0] = (-drho * cosel * cosaz +
                    rhosez[2] * del * cosaz + rhosez[1] * daz);
      drhosez[1] = (drho * cosel * sinaz -
                    rhosez[2] * del * sinaz - rhosez[0] * daz);
      drhosez[2] = (drho * sinel + rho * del * cosel);
    }
    else
    {
      /* ------------ calculate azimuth and elevation ------------- */
      temp = sqrt(rhosez[0] * rhosez[0] + rhosez[1] * rhosez[1]);
      if (fabs(rhosez[1]) < small)
        if (temp < small)
        {
          temp1 = sqrt(drhosez[0] * drhosez[0] +
                       drhosez[1] * drhosez[1]);
          az = atan2(drhosez[1] / temp1, drhosez[0] / temp1);
        }
        else
          if (drhosez[0]  > 0.0)
            az = M_PI;
          else
            az = 0.0;
          else
            az = atan2(rhosez[1] / temp, rhosez[0] / temp);
      
      if (temp < small)   // directly over the north pole
        el = astMath::sgn(rhosez[2]) * halfpi;  // +- 90
      else
        el = asin(rhosez[2] / astMath::mag(rhosez));
      
      /* ------  calculate range, azimuth and elevation rates ----- */
      drho = astMath::dot(rhosez, drhosez) / rho;
      if (fabs(temp * temp) > small)
        daz = (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) /
        (temp * temp);
      else
        daz = 0.0;
      
      if (fabs(temp) > small)
        del = (drhosez[2] - drho * sin(el)) / temp;
      else
        del = 0.0;
    }
  }   // rvsez_razel

/*------------------------------------------------------------------------------
   *
   *                           procedure rv_razel
   *
   *  this procedure converts range, azimuth, and elevation and their rates with
   *    the geocentric equatorial (ecef) position and velocity vectors.  notice the
   *    value of small as it can affect rate term calculations. uses velocity
   *    vector to find the solution of singular cases.
   *
   *  author        : david vallado                  719-573-2600   22 jun 2002
   *
   *  inputs          description                    range / units
   *    recef       - ecef position vector           km
   *    vecef       - ecef velocity vector           km/s
   *    rsecef      - ecef site position vector      km
   *    latgd       - geodetic latitude              -pi/2 to pi/2 rad
   *    lon         - geodetic longitude             -2pi to pi rad
   *    direct      -  direction to convert          eFrom  eTo
   *
   *  outputs       :
   *    rho         - satellite range from site      km
   *    az          - azimuth                        0.0 to 2pi rad
   *    el          - elevation                      -pi/2 to pi/2 rad
   *    drho        - range rate                     km/s
   *    daz         - azimuth rate                   rad/s
   *    del         - elevation rate                 rad/s
   *
   *  locals        :
   *    rhovecef    - ecef range vector from site    km
   *    drhovecef   - ecef velocity vector from site km/s
   *    rhosez      - sez range vector from site     km
   *    drhosez     - sez velocity vector from site  km
   *    tempvec     - temporary vector
   *    temp        - temporary extended value
   *    temp1       - temporary extended value
   *    i           - index
   *
   *  coupling      :
   *    astMath::mag         - astMath::magnitude of a vector
   *    addvec      - add two vectors
   *    rot3        - rotation about the 3rd axis
   *    rot2        - rotation about the 2nd axis
   *    atan2       - arc tangent function which also resloves quadrants
   *    dot         - dot product of two vectors
   *    rvsez_razel - find r2 and v2 from site in topocentric horizon (sez) system
   *    lncom2      - combine two vectors and constants
   *    arcsin      - arc sine function
   *    sgn         - returns the sign of a variable
   *
   *  references    :
   *    vallado       2013, 265, alg 27
   -----------------------------------------------------------------------------*/
  
  void rv_razel
  (
   double recef[3], double vecef[3], double rsecef[3], double latgd, double lon,
   edirection direct,
   double& rho, double& az, double& el, double& drho, double& daz, double& del
   )
  {
    const double halfpi = M_PI / 2.0;
    const double small = 0.0000001;
    
    double temp, temp1;
    double rhoecef[3], drhoecef[3], rhosez[3], drhosez[3], tempvec[3];
    
    if (direct == eFrom)
    {
      /* ---------  find sez range and velocity vectors ----------- */
      rvsez_razel(rhosez, drhosez, direct, rho, az, el, drho, daz, del);
      
      /* ----------  perform sez to ecef transformation ------------ */
      astMath::rot2(rhosez, latgd - halfpi, tempvec);
      astMath::rot3(tempvec, -lon, rhoecef);
      astMath::rot2(drhosez, latgd - halfpi, tempvec);
      astMath::rot3(tempvec, -lon, drhoecef);
      
      /* ---------  find ecef range and velocity vectors -----------*/
      astMath::addvec(1.0, rhoecef, 1.0, rsecef, recef);
      vecef[0] = drhoecef[0];
      vecef[1] = drhoecef[1];
      vecef[2] = drhoecef[2];
    }
    else
    {
      /* ------- find ecef range vector from site to satellite ----- */
      astMath::addvec(1.0, recef, -1.0, rsecef, rhoecef);
      drhoecef[0] = vecef[0];
      drhoecef[1] = vecef[1];
      drhoecef[2] = vecef[2];
      rho = astMath::mag(rhoecef);
      
      /* ------------ convert to sez for calculations ------------- */
      astMath::rot3(rhoecef, lon, tempvec);
      astMath::rot2(tempvec, halfpi - latgd, rhosez);
      astMath::rot3(drhoecef, lon, tempvec);
      astMath::rot2(tempvec, halfpi - latgd, drhosez);
      
      /* ------------ calculate azimuth and elevation ------------- */
      temp = sqrt(rhosez[0] * rhosez[0] + rhosez[1] * rhosez[1]);
      if (fabs(rhosez[1]) < small)
        if (temp < small)
        {
          temp1 = sqrt(drhosez[0] * drhosez[0] +
                       drhosez[1] * drhosez[1]);
          az = atan2(drhosez[1] / temp1, -drhosez[0] / temp1);
        }
        else
          if (rhosez[0] > 0.0)
            az = M_PI;
          else
            az = 0.0;
          else
            az = atan2(rhosez[1] / temp, -rhosez[0] / temp);
      
      if (temp < small)  // directly over the north pole
        el = astMath::sgn(rhosez[2]) * halfpi; // +- 90
      else
        el = asin(rhosez[2] / astMath::mag(rhosez));
      
      /* ----- calculate range, azimuth and elevation rates ------- */
      drho = astMath::dot(rhosez, drhosez) / rho;
      if (fabs(temp * temp) > small)
        daz = (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) /
        (temp * temp);
      else
        daz = 0.0;
      
      if (fabs(temp) > 0.00000001)
        del = (drhosez[2] - drho * sin(el)) / temp;
      else
        del = 0.0;
    }
  }  // rv_razel

































/*****************************************************************************/
 // converte da ra e dec geocentrice a ra e dec topocentriche
 // le ra e dec devono essere in gradi
 // site invece le cordinate del sito in TEME
 // sat invece le cordinate del satellite in TEME
 /*****************************************************************************/
/*
 void rv2tradec(double sat[3],  double & tra, double & tdec, double site[3]) {
   
   double rho[3];
   
   rho[0] = sat[0] - site[0];
   rho[1] = sat[1] - site[1];
   rho[2] = sat[2] - site[2];
   
   double normRho = astMath::mag(rho);
   
   double coord[2];
   
   coord[0] = atan2(rho[1],rho[0]);
   coord[1] = asin(rho[2] / normRho);

   if(coord[0] < 0.1) coord[0] += 2 * M_PI;
   
   tra  = astro::Degrees(coord[0]);
   tdec = astro::Degrees(coord[1]);

 }
  */

#endif /* _H_ASTRO_CONVERTER_H */
