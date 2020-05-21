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

#include <iostream>
#include <vector>

using namespace std;

namespace astro {

//double dot(double vect_A[3], double vect_B[3])
//{
//  double product = 0;
//  for (int i = 0; i < 3; i++){product = product + vect_A[i] * vect_B[i];}
//  return product;
//}


vector<double> rv2radec(double r_sat_ECI[3], double v_sat_ECI[3],double r_obs_ECEF[3], double jday)// double rho_eci[3], double ra, double dec, double drho_eci[3], double dtrtasc, double dtdecl)
{

  // FIXED: Modifica della funzione, da void a vector<double>
  vector<double> coord;
  double rho_eci[3], ra, dec, drho_eci[3], dtrtasc, dtdecl;
  const double small = 0.00000001;
  
  // --------------------- implementation ------------------------
  //----------------- get site vector in ecef -------------------
  //double r_obs_ECEF[3];
  //astIOD::site(latgd,lon,alt,r_obs_ECEF);
  
  double trtasc, tdecl;

  // -------------------- convert ecef to eci -------------------
  //METTERE DIRETTAMENTE LA POSIZIONE DELL'OSSERVATORIO IN ECI
  double r_obs_ECI[3];
  double v_obs_ECI[3];
  double v_obs_ECEF[3] = {0,0,0};
  astro::ecef2eci(r_obs_ECEF,v_obs_ECEF, jday,r_obs_ECI,v_obs_ECI);

  //---------------------range ------------------------
  rho_eci[0]  = r_sat_ECI[0] - r_obs_ECI[0];
  rho_eci[1]  = r_sat_ECI[1] - r_obs_ECI[1];
  rho_eci[2]  = r_sat_ECI[2] - r_obs_ECI[2];
  drho_eci[0] = v_sat_ECI[0] - v_obs_ECI[0];
  drho_eci[1] = v_sat_ECI[1] - v_obs_ECI[1];
  drho_eci[2] = v_sat_ECI[2] - v_obs_ECI[2];
  
  double normRho = astMath::mag(rho_eci);

  double temp = sqrt(pow(rho_eci[0],2) + pow(rho_eci[1],2));
  if(temp < small){trtasc = atan2(drho_eci[1], drho_eci[0]);}
  else{trtasc = atan2(rho_eci[1], rho_eci[0]);} //rad
  
  //FIXED: controllo sulla ra (RIGA AGGIUNTA PER CORREGGERE CODICE VALLADO)
  if(trtasc < 0.1) trtasc += 2.0 * M_PI; 

  if (temp < small)           //directly over the north pole
  {
    if(rho_eci[2] > 0){tdecl= M_PI/2;}    // +90 deg
    else if(rho_eci[2] == 0){tdecl = 0.0;}  // 0
    else {tdecl= -M_PI/2;}                // -90 deg
  }
  else
  {
    double magrhoeci = astMath::mag(rho_eci);
    tdecl= asin(rho_eci[2]/magrhoeci); //rad
  }


  double temp1 = -rho_eci[1]*rho_eci[1] - rho_eci[0]*rho_eci[0];
  double d_rho = astMath::dot(rho_eci,drho_eci)/normRho;
  if (abs(temp1)>small){dtrtasc = (drho_eci[0]*rho_eci[1] - drho_eci[1]*rho_eci[0])/temp1;}
  else {dtrtasc = 0.0;}
  if (abs(temp)>small){dtdecl = (drho_eci[2] - d_rho*sin(tdecl))/temp;}
  else{dtdecl= 0.0;}


  //conversione in gradi ra e dec topocentriche
  ra  = astro::degrees(trtasc);
  dec = astro::degrees(tdecl);

  coord.push_back(ra);
  coord.push_back(dec);

  return coord;
  
}

vector<double> RaDec2AzEl(double ra, double dec, double lat, double lon, double JD)
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

  double lat_rad = astro::radians(lat);
  //double lon_rad = astro::radians(lon);
  //double ra_rad  = astro::radians(ra);
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
  double El = astro::degrees(El_rad);

  //Azimuth deg
  double c_rad = atan2(-sin(LHA_rad)*cos(dec_rad)/cos(El_rad),(sin(dec_rad)-sin(El_rad)*sin(lat_rad))/(cos(El_rad)*cos(lat_rad)));
  double c     = astro::degrees(c_rad);
  double Az    = c - floor(c/360.0)*360.0; //deg

  vector<double> AzEl;

  AzEl.push_back(Az);
  AzEl.push_back(El);

  return AzEl;
}

} /* namespace astro */

#endif /* _H_ASTRO_CONVERTER_H */
