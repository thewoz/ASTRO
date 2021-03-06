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

#ifndef _H_ASTRO_EFFECT_H
#define _H_ASTRO_EFFECT_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vector>


//****************************************************************************/
// namespace astro
//****************************************************************************/
namespace astro {  

  /* ***************** Function Precession ******************* *\
 *                                                           *  
 * Transform Ra,Dec in J2000 to Ra,Dec in JNOW.              * 
 * taking into account precession.                           *
 * Input:                                                    *  
 *    t         - julian centuries                           *
 *   ra         - right ascension in J2000 [rad]             *  
 *   dec        - declination in J2000 [rad]                 * 
 *   ra1        - right ascension in JNOW [rad]              *
 *   dec1       - declination in JNOW [rad]                  *
 *                                                           * 
 * ********************************************************* */
void precession(double ra, double dec, double t, double& ra1, double& dec1)
{
  double zita, zeta, theta;
  double A, B, C;
  const double rad2deg = 180./M_PI;

  zita =  ( 2306.2181*t + 0.30188*t*t + 0.017998*t*t*t ) / rad2deg / 3600.;
  zeta =  ( 2306.2181*t + 1.09468*t*t + 0.018203*t*t*t ) / rad2deg / 3600.;
  theta = ( 2004.3109*t - 0.42665*t*t - 0.041833*t*t*t ) / rad2deg / 3600.;

  A = cos(dec)*sin(ra+zita);
  B = cos(theta)*cos(dec)*cos(ra+zita)-sin(theta)*sin(dec);
  C = sin(theta)*cos(dec)*cos(ra+zita)+cos(theta)*sin(dec);

  ra1  = atan2(A,B)+zeta;
  dec1 = asin(C);
}

/* ***************** Function Precession2 ******************* *\
 *                                                            * 
 * General transform taking into account precession.          *
 * from t0=(jd_start-jd_j2k)/36525 to t=(jd_end-jd_j2k)/36525 *
 *                                                            *
 * ********************************************************** */
void precession2(double ra, double dec, double t0, double t, double& ra1, double& dec1)
{
  double zita, zeta, theta;
  double A, B, C;
  const double rad2deg = 180./M_PI;

  zita =  (2306.2181+1.39656*t0-0.000139*t0*t0)*t + (0.30188-0.000344*t0)*t*t + 0.017998*t*t*t;
  zita = zita/rad2deg/3600.;
  zeta =  (2306.2181+1.39656*t0-0.000139*t0*t0)*t + (1.09468+0.000066*t0)*t*t + 0.018203*t*t*t;
  zeta = zeta/rad2deg/3600.;
  theta = (2004.3109-0.85330*t0-0.000217*t0*t0)*t - (0.42665+0.000217*t0)*t*t - 0.041833*t*t*t;
  theta = theta/rad2deg/3600.;

  A = cos(dec)*sin(ra+zita);
  B = cos(theta)*cos(dec)*cos(ra+zita)-sin(theta)*sin(dec);
  C = sin(theta)*cos(dec)*cos(ra+zita)+cos(theta)*sin(dec);

  ra1  = atan2(A,B)+zeta;
  dec1 = asin(C);
}

  

/* ****************** Function aberration ******************* *\
 *                                                            *
 * this function compute the correction to ra/dec (in jnow)   *
 * due to annual aberration                                   *
 * input:                                                     *
 *       ra  [rad]                                            *
 *       dec [rad]                                            *
 *       t   julian centuries                                 *
 * output:                                                    *
 *       dra   [rad]                                          *
 *       ddec  [rad]                                          *
 * we use formula 23.3 of Astronomical Algorithms p.151       *
 *                                                            *
 \* ********************************************************* */
  void aberration(double Ra, double Dec, double t, double& deltaRa, double& deltaDec)
  {
  const double rad2deg = 180.0 / M_PI;
  double cosRa, sinRa, cosDec, sinDec;
  double SunMeanLon, SunMeanAnomaly, SunEquationOfTheCenter, SunTrueLon;
  double MoonMeanLon, MoonMeanAnomaly, MoonAscendingNodeLon;
  double EarthPerihelionLon, eEarthOrbit;
  double MeanObliquityEcliptic, nutationLon, nutationObliquity, eclipticObliquity;
  double k, cosSunTrueLon, sinSunTrueLon;
  double cosEclipticObliquity, tanEclipticObliquity;
  double cosEarthPerihelionLon, sinEarthPerihelionLon;

  // aberration constant
  k = 20.49552 / 3600.;

  // Earth orbit eccentricity and longitude of perihelion
  eEarthOrbit = 0.016708634 - 0.000042037*t - 0.0000001267*t*t;
  EarthPerihelionLon = (102.93735 + 1.71946*t + 0.00046*t*t) / rad2deg;
  astUtils::rebox(EarthPerihelionLon);

  // Moon Longitude
  MoonMeanLon = (218.3165 + 481267.8813*t) / rad2deg;
  astUtils::rebox(MoonMeanLon);
  MoonMeanAnomaly = (134.96298 + 477198.867398*t - 0.0086972*t*t + t*t*t/56250)  / rad2deg;
  astUtils::rebox(MoonMeanAnomaly);
  MoonAscendingNodeLon = (125.04452 - 1934.136261*t + 0.0020708*t*t + t*t*t / 450000)  / rad2deg;
  astUtils::rebox(MoonAscendingNodeLon);

  // Sun Longitude, Anomaly
  SunMeanLon = (280.46646 + 36000.76983*t) / rad2deg;
  astUtils::rebox(SunMeanLon);
  SunMeanAnomaly = (357.52772 + 35999.050340*t - 0.0001603*t*t - t*t*t/300000)  / rad2deg;
  astUtils::rebox(SunMeanAnomaly);
  SunEquationOfTheCenter = (1.914602 - 0.004817*t - 0.000014*t*t) * sin(SunMeanAnomaly)
    + (0.019993 - 0.000101*t) * sin(2*SunMeanAnomaly) + 0.000289 * sin(3 * SunMeanAnomaly);
  SunEquationOfTheCenter /= rad2deg;
  SunTrueLon = SunMeanLon + SunEquationOfTheCenter;
  astUtils::rebox(SunTrueLon);

  // Ecliptic Obliquity
  MeanObliquityEcliptic = (23.43929 - 46.815/3600.*t - 0.00059/3600.*t*t + 0.001813/3600.*t*t*t) / rad2deg;
  nutationLon = - 17.2*sin(MoonAscendingNodeLon) - 1.32*sin(2.0*SunMeanLon) - 0.23*sin(2.0*MoonMeanLon) + 0.21*sin(2*MoonAscendingNodeLon);
  nutationLon = nutationLon/3600./rad2deg;
  nutationObliquity = 9.2*cos(MoonAscendingNodeLon) + 0.57*cos(2*SunMeanLon) + 0.1*cos(2*MoonMeanLon) - 0.09*cos(2*MoonAscendingNodeLon);
  nutationObliquity = nutationObliquity/3600./rad2deg;
  eclipticObliquity = MeanObliquityEcliptic + nutationObliquity;

  cosRa = cos(Ra);
  sinRa = sin(Ra);
  cosDec = cos(Dec);
  sinDec = sin(Dec);
  cosSunTrueLon = cos(SunTrueLon);
  sinSunTrueLon = sin(SunTrueLon);
  cosEclipticObliquity = cos(eclipticObliquity);
  tanEclipticObliquity = tan(eclipticObliquity);
  cosEarthPerihelionLon = cos(EarthPerihelionLon);
  sinEarthPerihelionLon = sin(EarthPerihelionLon);

  deltaRa = - cosRa * cosSunTrueLon * cosEclipticObliquity - sinRa * sinSunTrueLon
    + eEarthOrbit * cosRa * cosEarthPerihelionLon * cosEclipticObliquity
    + eEarthOrbit * sinRa * sinEarthPerihelionLon;
  deltaRa = deltaRa * k / (rad2deg * cosDec);

  deltaDec = -cosSunTrueLon * cosEclipticObliquity * (tanEclipticObliquity * cosDec - sinRa * sinDec)
    - cosRa * sinDec * sinSunTrueLon
    + eEarthOrbit * cosEarthPerihelionLon * cosEclipticObliquity * (tanEclipticObliquity * cosDec - sinRa * sinDec)
    + eEarthOrbit * cosRa * sinDec * sinEarthPerihelionLon;
  deltaDec = deltaDec * k / (rad2deg);

}

/* ******************* Function nutation ******************* *\
 *                                                           *
 * this function compute the correction in terms of right    *
 * ascension and declination due to nutation                 *
 *                                                           *
 \* ******************************************************** */
void nutation_radec(double Ra, double Dec, double t, double& deltaRa, double& deltaDec)
{  
  const double rad2deg = 180.0 / M_PI;
  double cosRa, sinRa, tanDec;
  double SunMeanLon; //, SunMeanAnomaly;//, SunTrueLon;
  double MoonMeanLon, MoonAscendingNodeLon; //MoonMeanAnomaly;
  double MeanObliquityEcliptic, nutationLon, nutationObliquity, eclipticObliquity;
  double cosEclipticObliquity, sinEclipticObliquity;
  
  // Moon Longitude
  MoonMeanLon = (218.3165 + 481267.8813*t) / rad2deg;
  astUtils::rebox(MoonMeanLon);  
  MoonAscendingNodeLon = (125.04452 - 1934.136261*t + 0.0020708*t*t + t*t*t / 450000)  / rad2deg;
  astUtils::rebox(MoonAscendingNodeLon);

  // Sun Longitude, Anomaly
  SunMeanLon = (280.46646 + 36000.76983*t) / rad2deg;
  astUtils::rebox(SunMeanLon);
  
  // Nutation & Ecliptic Obliquity
  MeanObliquityEcliptic = (23.43929 - 46.815/3600.*t - 0.00059/3600.*t*t + 0.001813/3600.*t*t*t) / rad2deg;
  nutationLon = - 17.2*sin(MoonAscendingNodeLon) - 1.32*sin(2.0*SunMeanLon) - 0.23*sin(2.0*MoonMeanLon) + 0.21*sin(2*MoonAscendingNodeLon);
  nutationLon =nutationLon/3600./rad2deg;
  nutationObliquity = 9.2*cos(MoonAscendingNodeLon) + 0.57*cos(2*SunMeanLon) + 0.1*cos(2*MoonMeanLon) - 0.09*cos(2*MoonAscendingNodeLon);
  nutationObliquity = nutationObliquity/3600./rad2deg;
  eclipticObliquity = MeanObliquityEcliptic + nutationObliquity;  
  
  cosRa = cos(Ra);
  sinRa = sin(Ra);
  tanDec = tan(Dec);
  cosEclipticObliquity = cos(eclipticObliquity);
  sinEclipticObliquity = sin(eclipticObliquity);
  
  deltaRa = nutationLon*(cosEclipticObliquity + sinEclipticObliquity*sinRa*tanDec) - nutationObliquity*cosRa*tanDec;
  deltaDec = nutationLon*sinEclipticObliquity*cosRa + nutationObliquity*sinRa;
  
}

/* ******************* Function J2K_JNOW ******************* *\
 * Transformations between radec J2K and radec JNOW,         *
 * including aberration, precession and nutation.            *
\* ********************************************************* */
void j2k_jnow(double& raj2k, double& decj2k, double jd, edirection direct, double& rajnow, double& decjnow)
{
  
  double t, t0, t1;
  double ra1, dec1, deltara, deltadec, dra_prec, ddec_prec;
  double dra_nut, ddec_nut;//, dra_abe, ddec_abe;
  
  if (direct == eTo)
    {    
      t = (jd - 2451545.)/36525.;
      astro::precession(raj2k, decj2k, t, ra1, dec1);
      astUtils::rebox(ra1);
      dra_prec = ra1-raj2k; ddec_prec = dec1-decj2k;
      astro::nutation_radec(raj2k, decj2k, t, deltara, deltadec);
      dra_nut = deltara; ddec_nut = deltadec;
      //astro::aberration(raj2k, decj2k, t, deltara, deltadec);
      //dra_abe = deltara; ddec_abe = deltadec;
      rajnow = raj2k + dra_prec + dra_nut;// + dra_abe;
      decjnow = decj2k + ddec_prec + ddec_nut;// + ddec_abe;
    }
  
  if (direct == eFrom)
    {
      t0 = (jd - 2451545.)/36525.;
      t1 = (2451545. - jd)/36525.;
      astro::precession2(rajnow, decjnow, t0, t1, raj2k, decj2k);
      astUtils::rebox(raj2k);
      dra_prec = rajnow-raj2k; ddec_prec = decjnow-decj2k;
      astro::nutation_radec(raj2k, decj2k, t0, deltara, deltadec);
      dra_nut = deltara; ddec_nut = deltadec;
      //astro::aberration(raj2k, decj2k, t0, deltara, deltadec);
      //dra_abe = deltara; ddec_abe = deltadec;
      raj2k = rajnow - dra_prec - dra_nut;// - dra_abe;
      decj2k = decjnow - ddec_prec - ddec_nut;// - ddec_abe;
    }
  
}

// convert az/el (ECEF/JNOW) in ra/dec (ECI/J2K) and viceversa
void j2k_azel(double& ff, double& ra0, double& dec0, double latobs, double lonobs, double altobs,
              double jd, edirection edirect, double& gg, double& az0, double& el0)
{
  
  double rsat_ecef[3], vsat_ecef[3];
  double rsat_eci[3], vsat_eci[3];
  double robs_ecef[3], vobs_ecef[3];
  double robs_eci[3];//, vobs_eci[3];
  double drange, daz0, del0, dra0, ddec0;
  
  // site in ecef/eci
  astIOD::site(latobs, lonobs, altobs, robs_ecef, vobs_ecef);
  astro::ecef2eci(robs_ecef, jd, robs_eci);

  for(int k=0;k<3;k++){
    vsat_ecef[k]=vsat_eci[k]=0.0;
  }
  
  drange=0.0;
  daz0=0.0; del0=0.0;
  dra0=0.0; ddec0=0.0;
  
  if (edirect == eFrom)
    {
      astIOD::rv_razel(rsat_ecef, vsat_ecef, robs_ecef, latobs, lonobs, eFrom, gg, az0, el0, drange, daz0, del0);
      astro::ecef2eci(rsat_ecef, jd, rsat_eci);
      astIOD::rv_tradec(rsat_eci, vsat_eci, robs_eci, eTo, ff, ra0, dec0, drange, dra0, ddec0);
      astUtils::rebox(ra0);
    }
  
  if (edirect == eTo)
    {
      astIOD::rv_tradec(rsat_eci, vsat_eci, robs_eci, eFrom, ff, ra0, dec0, drange, dra0, ddec0);
      astro::eci2ecef(rsat_eci, jd, rsat_ecef);
      astIOD::rv_razel(rsat_ecef, vsat_ecef, robs_ecef, latobs, lonobs, eTo, gg, az0, el0, drange, daz0, del0);
      astUtils::rebox(az0);
    }
  
}  

  
} /* namespace astro */

#endif /* _H_ASTRO_EFFECT_H
*/
