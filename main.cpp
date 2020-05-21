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

//#define TESTEOPC
//#define TESTDATE
//#define TESTSATELLITE
//#define TESTSUN
//#define TESTSTATION
//#define TESTATTITUDE
//#define TESTCONVERT
#define TESTRADEC

#include <cstdio>
#include <cstdlib>

#include <cerrno>
#include <cstring>

#include <string>
#include <iostream>

#include "astro/astro.hpp"


using namespace std;

//****************************************************************************/
// main
//****************************************************************************/
int main(int argc, char *argv[]) {

//****************************************************************************/
// Test Convert
//****************************************************************************/
#ifdef TESTCONVERT
  
  double jday_prova = astro::Date(19,04,2020,11,00,00).getJDay();
  double a[3] = {1000, 2000, 3000}; double b[3]; double c[3];

  astro::teme2ecef(a, jday_prova, b); 
  astro::ecef2teme(b, jday_prova, c); 
  printf("teme2ecef %e\n", sqrt((((a[0]-c[0])*(a[0]-c[0]))+((a[1]-c[1])*(a[1]-c[1]))+((a[2]-c[2])*(a[2]-c[2])))));
  astro::teme2eci(a,  jday_prova, b);  
  astro::eci2teme(b, jday_prova, c); 
  printf("eci2teme  %e\n",  sqrt((((a[0]-c[0])*(a[0]-c[0]))+((a[1]-c[1])*(a[1]-c[1]))+((a[2]-c[2])*(a[2]-c[2])))));

  astro::ecef2teme(a, jday_prova, b); 
  astro::teme2ecef(b, jday_prova, c); 
  printf("ecef2teme %e\n", sqrt((((a[0]-c[0])*(a[0]-c[0]))+((a[1]-c[1])*(a[1]-c[1]))+((a[2]-c[2])*(a[2]-c[2])))));
  astro::ecef2eci(a,  jday_prova, b);  
  astro::eci2ecef(b, jday_prova, c); 
  printf("ecef2eci  %e\n",  sqrt((((a[0]-c[0])*(a[0]-c[0]))+((a[1]-c[1])*(a[1]-c[1]))+((a[2]-c[2])*(a[2]-c[2])))));

  astro::eci2ecef(a, jday_prova, b); 
  astro::ecef2eci(b, jday_prova, c); 
  printf("eci2ecef %e\n", sqrt((((a[0]-c[0])*(a[0]-c[0]))+((a[1]-c[1])*(a[1]-c[1]))+((a[2]-c[2])*(a[2]-c[2])))));
  astro::eci2teme(a, jday_prova, b); 
  astro::teme2eci(b, jday_prova, c); 
  printf("eci2teme %e\n", sqrt((((a[0]-c[0])*(a[0]-c[0]))+((a[1]-c[1])*(a[1]-c[1]))+((a[2]-c[2])*(a[2]-c[2])))));
  

  //ATLAS CENTAUR 2         
  //1 00694U 63047A   20140.15326762  .00000182  00000-0  82935-5 0  9991
  //2 00694  30.3610  52.1406 0585972 292.3446  61.6075 14.02595276832603

  char TLE_line1[] = "1 00694U 63047A   20140.15326762  .00000182  00000-0  82935-5 0  9991";
  char TLE_line2[] = "2 00694  30.3610  52.1406 0585972 292.3446  61.6075 14.02595276832603";

  double jday = astro::Date(19,05,2020,11,00,0.0).getJDay();

  astro::Satellite(TLE_line1, TLE_line2).position(jday, a, CRS::TEME);

  double vijk[3];

  double rr; double ra; double dec; double drr; double drtasc; double ddecl;

  astIOD::rv_radec(a, vijk, edirection::eTo, rr, ra, dec, drr, drtasc, ddecl);

  printf("ra %f dec %f\n", astro::degrees(ra), astro::degrees(dec));

  astro::teme2ecef(a, jday, b); 
  astro::ecef2teme(b, jday, c); 
  printf("teme2ecef %e\n", sqrt((((a[0]-c[0])*(a[0]-c[0]))+((a[1]-c[1])*(a[1]-c[1]))+((a[2]-c[2])*(a[2]-c[2])))));
  astro::teme2eci(a,  jday, b);  
  astro::eci2teme(b, jday, c); 
  printf("eci2teme  %e\n",  sqrt((((a[0]-c[0])*(a[0]-c[0]))+((a[1]-c[1])*(a[1]-c[1]))+((a[2]-c[2])*(a[2]-c[2])))));

#endif


//****************************************************************************/
// Test Attitude
//****************************************************************************/
#ifdef TESTATTITUDE
  
  // angoli in radienati
  double pitch = 0;
  double yaw   = 0;
  double roll  = 0;

  //pitch = astro::radians(1.0);
  //yaw   = astro::radians(1.0);
  //roll  = astro::radians(1.0);
  
  astro::quaternion_t q(pitch, yaw, roll);

  // velocita in radianti al secondo
  double Wx = 0.0;
  double Wy = 0.0;
  double Wz = 0.0;
  
  //Wx = astro::radians(1.0);
  //Wy = astro::radians(1.0);
  //Wz = astro::radians(1.0);
  
  double Jx = 1;
  double Jy = 1;
  double Jz = 1;
  
  astro::attitude_t start(Wx,Wy,Wz,q.q4,q.q1,q.q2,q.q3,Jx,Jy,Jz);

  std::vector<astro::attitude_t> result;
  
  // numero passi
  size_t steps = 100;
  
  // dt integrazione in secondi
  float dt = 1;
  
  astro::attitude::compute(start, steps, dt, result);
  
  for(int i=0; i<steps; ++i) {

    result[i].getAngles(pitch, yaw, roll);
    
    //printf("%d) %f %f %f - %f %f %f\n", i, result[i].Wx(), result[i].Wy(), result[i].Wz(), Degrees(pitch), Degrees(yaw), Degrees(roll));
    
  //  q.println();
    
  }

#endif
  
  
//****************************************************************************/
// Test Date class
//****************************************************************************/
#ifdef TESTDATE

  fprintf(stderr, "Test Date class:\n\n");
  
  astro::Date T1 = astro::Date(19, 11, 2017, 4, 8, 35); T1.println();
  astro::Date(T1.jDay, T1.jDayFrac).println();
  astro::Date(T1.jDay+T1.jDayFrac).println();
  
  astro::Date T2 = astro::Date(19, 11, 2017, 20, 8, 35); T2.println();
  astro::Date(T2.jDay, T2.jDayFrac).println();
  astro::Date(T2.jDay+T2.jDayFrac).println();
  
  T1('+',  10, astro::Date::MINUTES).println();
  T1('+', -10, astro::Date::MINUTES).println();

#endif
  
  
//****************************************************************************/
// Test Epoc class
//****************************************************************************/
#ifdef TESTEOPC

  fprintf(stderr, "Test Epoc class:\n\n");
  
  double xp, yp, lod, dpsi, deps, jdut1, jdut1Frac, ttt;
  
  astro::eopc::getParameters(astro::Date(29, 1, 2017, 13, 34, 1).getJDay(), 'l', 'f', xp, yp, lod, dpsi, deps, jdut1, jdut1Frac, ttt);
  
  printf("\n xp %e\n yp %e\n lod %e\n dpsi %e\n deps %e\n jdut1 %e\n jdut1Frac %e\n ttt %e\n\n\n", xp, yp, lod, dpsi, deps, jdut1, jdut1Frac, ttt);
  
#endif
  
  
//****************************************************************************/
// Test Orbit Propagator from Satellite
//****************************************************************************/
#ifdef TESTSATELLITE
    
  fprintf(stderr, "Test Satellite:\n\n");

  std::string outStrSatellite = "/home/gaetano/Desktop/satellite.dat";

  FILE * outputSatellite = fopen(outStrSatellite.c_str(), "w");

  if(outputSatellite==NULL){
    fprintf(stderr, "error in opening the file '%s': %s\n", outStrSatellite.c_str(), strerror(errno));
    abort();
  }

  char TLE_line1[] = "1 00694U 63047A   20140.15326762  .00000182  00000-0  82935-5 0  9991";
  char TLE_line2[] = "2 00694  30.3610  52.1406 0585972 292.3446  61.6075 14.02595276832603";
  
  std::vector<astro::SatelliteState> states;

  astro::Satellite(TLE_line1, TLE_line2).orbit(astro::Date(19,05,2020,11,00,00), astro::Date(19,05,2020,12,00,00), 60, states, CRS::ECI);

  for(std::size_t i=0; i<states.size(); ++i){
    states[i].println(outputSatellite);
    states[i].println(stdout);
  }
  
  fclose(outputSatellite);

  fprintf(stderr, "\n\n\n");
    
#endif
  
  
//****************************************************************************/
// Test Orbit Propagator for Sun
//****************************************************************************/
#ifdef TESTSUN

  std::string outStrSun = "/Users/thewoz/Desktop/sun_kde_precise.dat";
  
  FILE * outputSun = fopen(outStrSun.c_str(), "w");
  
  if(outputSun==NULL){
    fprintf(stderr, "error in opening the file '%s': %s\n", outStrSun.c_str(), strerror(errno));
    abort();
  }
  
  std::vector<astro::SunState> states;
  
  astro::Sun::orbit(astro::Date(20,4,2017,17,55,00).getJDay(), astro::Date(20,5,2017,17,55,00).getJDay(), 60, states, CRS::ECI);
  
  for(std::size_t i=0; i<states.size(); ++i){
    fprintf(outputSun, "%s ", astro::Date(states[i].jDay).toGregorianString());
    states[i].println(outputSun);
    //fprintf(stdout, "%s ", astro::Date(states[i].jDay).toGregorianString());
    //states[i].println(stdout);
  }
  
  fclose(outputSun);
  
  fprintf(stderr, "\n\n\n");

#endif

  
//****************************************************************************/
// Test Orbit Propagator for Station
//****************************************************************************/
#ifdef TESTSTATION
    
  fprintf(stderr, "Test Station:\n\n");

  std::string outStrStation = "/home/gaetano/Desktop/urbe.dat";
  
  FILE * outputStation = fopen(outStrStation.c_str(), "w");
  
  if(outputStation==NULL){
    fprintf(stderr, "error in opening the file '%s': %s\n", outStrStation.c_str(), strerror(errno));
    abort();
  }

  double latitude  = 41.9558333333333;
  double longitude = 12.5055555555556;
  double altitude  = 76.0;
  
  double dt = 60;

  //lat e lon sono in gradi. La conversione in rad la fa dentro Observatory
  astro::Observatory urbe(latitude, longitude, altitude);
  
  std::vector<astro::ObservatoryState> states;
  
  urbe.orbit(astro::Date(19,05,2020,11,00,00).getJDay(), astro::Date(19,05,2020,12,00,00).getJDay(), dt, states, CRS::ECI);

  for(std::size_t i=0; i<states.size(); ++i){
    fprintf(stdout, "%s ", astro::Date(states[i].jDay).toGregorianString());
    states[i].println(stdout);
    fprintf(outputStation, "%s ", astro::Date(states[i].jDay).toGregorianString());
    states[i].println(outputStation);
  }
  
  fclose(outputStation);
  
  fprintf(stderr, "\n\n\n");

#endif

  
//****************************************************************************/
// Test finale
//****************************************************************************/
#ifdef TESTRADEC

  //OSSERVATORIO (URBE)
  double latitude  = 41.9558333333333; //deg
  double longitude = 12.5055555555556; //deg
  double altitude  = 76.0;             // m

  //tempo di propagazione (secondi)
  double dt_minuti = 3;
  //step di propagazione (secondi)
  double step = 60.0;
  //Data inizio propagazione
  int day = 19;
  int month = 05;
  int year = 2020;
  int hour = 11;
  int min = 00;
  int sec = 00;

  // TLE
  char TLE_line1[] = "1 00694U 63047A   20141.86162056  .00000196  00000-0  10739-4 0  9991";
  char TLE_line2[] = "2 00694  30.3607  42.7351 0585922 307.1727  47.6741 14.02595553833226";

  ////////////////////////////////
  ////    START  ////////////////
  ///////////////////////////////

  //conversione data inizio propagazione in data giuliana
  double jday_start = astro::Date(day,month,year,hour,min,sec).getJDay();

  //conversione step di propagazione in data giuliana
  double dt2jd = astro::Date::convert(dt_minuti*60,astro::Date::FROM_SECOND_TO_JD);

  // intervalli di propagazione
  double jday_end   = jday_start + dt2jd;

  //inizializzazione variabile data giuliana per il ciclo while
  double jday = jday_start;

  //lat e lon sono in gradi. La conversione in rad la fa dentro Observatory
  astro::Observatory urbe_ecef(latitude, longitude, altitude);
  astro::Observatory urbe_eci(latitude, longitude, altitude);
  std::vector<astro::ObservatoryState> obs_ecef_states;
  std::vector<astro::ObservatoryState> obs_eci_states;

  //propagazione posizione osservatorio in eci
  urbe_ecef.orbit(jday_start, jday_end, step, obs_ecef_states, CRS::ECEF);
  urbe_eci.orbit(jday_start, jday_end, step, obs_eci_states, CRS::ECI);

  // propagazione stato satellite in eci
  std::vector<astro::SatelliteState> sat_states;
  astro::Satellite(TLE_line1, TLE_line2).orbit(astro::Date(19,05,2020,11,00,00), astro::Date(19,05,2020,11,03,00), step, sat_states, CRS::ECI);

  //contatore di ausilio
  int j = 0;

  //ciclo sul tempo
  while(jday_end-jday > 0.0 )
  {
    //inizializzazione variabili
    double r_sat_eci[3], v_sat_eci[3], r_obs_ecef[3], r_obs_eci[3];

    //poszione satellite e osservatorio in ECI e ECEF
    for(int i = 0; i < 3; i++)
    {
      r_sat_eci[i]  = sat_states[j].position[i];
      v_sat_eci[i]  = sat_states[j].velocity[i];
      r_obs_ecef[i] = obs_ecef_states[j].position[i];
      r_obs_eci[i]  = obs_eci_states[j].position[i];
    }


    double ra,dec;
    //conversione da r [km],v[km/s] a ra,dec topocentriche [rad]
    astro::rv2radec(r_sat_eci, v_sat_eci, r_obs_ecef, jday, ra, dec);
    //printf("RA  = %f, DEC\n",ra);
    //printf("DEC = %f\n",dec);

    double Az,El;
    //conversione ra,dec [deg,deg] topocentriche in azimut ed elevazione [deg,deg]
    astro::RaDec2AzEl(ra,dec,latitude,longitude,jday,Az,El);
    //printf("Az  = %f\n",Az);
    //printf("El  = %f\n",El);

    //condizioni di visibilità
    int jday_int = jday;
    double jday_frac = jday - jday_int;
    double El_min = 10; //deg
    bool ill_sat = astUtils::light(r_sat_eci,(double)jday_int,jday_frac,'e');
    bool ill_obs = astUtils::light(r_obs_eci,(double)jday_int,jday_frac,'e');
    bool el = El > El_min;
    
    /*
    if(ill_sat == false){cout << "satellite non illuminato" << endl;}
    else{cout << "satellite illuminato" << endl;}
    if(ill_obs == false){cout << "stazione in ombra" << endl;}
    else{cout << "stazione illuminata" << endl;}
    if(el == false){cout << "satellite troppo basso" << endl;}
    else{cout << "satellite in linea di vista" << endl;}
    */

    if(ill_sat == true && ill_obs == false && el == true)
    {cout << "Il satellite è in visibilità: " << "RA = " << ra << "°, " << "DEC = " << dec << "°, " << "Az = " << Az << "°, "  << "El = " << El << "°, " << endl;}
    else
    {cout << "Il satellite NON è in visibilità: " << "RA = " << ra << "°, " << "DEC = " << dec << "°, " << "Az = " << Az << "°, "  << "El = " << El << "°, " << endl;}

    //aggiornamento tempo di propagazione e contatore
    jday += step/(24.0*3600.0);
    j++;
    
  }

    
#endif

  return 0;
  
}
