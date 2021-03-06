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
#include <iostream>
#include <fstream>

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
  
  //autovalori matrice di inerzia
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

  astro::Satellite(TLE_line1, TLE_line2).orbit(astro::Date(19,05,2020,11,00,00).getJDay(), astro::Date(19,05,2020,12,00,00).getJDay(), 60, states, CRS::ECI);

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

////////////////////////////////
////    INPUT  ////////////////
///////////////////////////////

double Rt = 6378.137;

//OSSERVATORIO (URBE)
/*
double latitude  = 41.9558333333333; //deg
double longitude = 12.5055555555556; //deg
double altitude  = 76.0;             // m
*/

//OSSERVATORIO (COLLEPARDO)
double latitude  = 41.765242; //deg
double longitude = 13.375038; //deg
double altitude  = 555.0;     // m


//tempo di propagazione (minuti) 
double dt_minuti = 3;//0.066666667*1000/60;
//step di propagazione (secondi)
double step = 1;//0.066666667;
//Data inizio propagazione
int day = 26;
int month = 9;
int year = 2020;
int hour = 19;
int min = 0;
int sec = 0;

// TLE 
char TLE_line1[] = "1 13068U 82013B   20269.92133849 -.00000280 +00000-0 -19928-4 0  9995";
char TLE_line2[] = "2 13068 081.2009 295.6177 0020273 327.1233 032.8728 15.07472599103500";


/////////////////////////////////
////    START    ////////////////
/////////////////////////////////

//conversione data inizio propagazione in data giuliana
double jday_start = astro::Date(day,month,year,hour,min,sec).getJDay();

//conversione tempo di propagazione in data giuliana
double dt2jd = astro::Date::convert(dt_minuti*60,astro::Date::FROM_SECOND_TO_JD);

// intervalli di propagazione
double jday_end   = jday_start + dt2jd;

//inizializzazione variabile data giuliana per il ciclo while
double jday = jday_start;

//lat e lon sono in gradi. La conversione in rad la fa dentro Observatory
astro::Observatory obs_ecef(latitude, longitude, altitude);
astro::Observatory obs_eci(latitude, longitude, altitude);
std::vector<astro::ObservatoryState> obs_ecef_states;
std::vector<astro::ObservatoryState> obs_eci_states;

//conversione latitudine geocentrica a latitudine geodetica
double eesqrd = 0.006694385000;
double latgd = atan(tan(astro::radians(latitude))/(1.0 - eesqrd));

//propagazione posizione osservatorio in eci
obs_ecef.orbit(jday_start, jday_end, step, obs_ecef_states, CRS::ECEF);
obs_eci.orbit(jday_start, jday_end, step, obs_eci_states, CRS::ECI);

// propagazione stato satellite in eci
std::vector<astro::SatelliteState> sat_eci_states;
std::vector<astro::SatelliteState> sat_ecef_states;
std::vector<astro::SatelliteState> sat_teme_states;
astro::Satellite(TLE_line1, TLE_line2).orbit(jday_start,jday_end, step, sat_eci_states, CRS::ECI);
astro::Satellite(TLE_line1, TLE_line2).orbit(jday_start,jday_end, step, sat_ecef_states, CRS::ECEF);
astro::Satellite(TLE_line1, TLE_line2).orbit(jday_start,jday_end, step, sat_teme_states, CRS::TEME);

//contatore di ausilio
int j = 0;

// Scrittura file di testo con RA e DEC
ofstream myfile;
myfile.open ("./RADEC.txt");
//ciclo sul tempo
while(jday_end-jday > 0.0 )
{
  //inizializzazione variabili
  double r_sat_eci[3], v_sat_eci[3],r_sat_ecef[3],v_sat_ecef[3], r_obs_ecef[3], r_obs_eci[3];
  double ra,dec, Az, El, ra_rad, dec_rad, Az2, El2;;

  //posizione satellite e osservatorio in ECI e ECEF
  for(int i = 0; i < 3; i++)
  {
    r_sat_eci[i]  = sat_eci_states[j].position[i];
    v_sat_eci[i]  = sat_eci_states[j].velocity[i];
    r_sat_ecef[i]  = sat_ecef_states[j].position[i];
    v_sat_ecef[i]  = sat_ecef_states[j].velocity[i];
    r_obs_ecef[i] = obs_ecef_states[j].position[i];
    r_obs_eci[i]  = obs_eci_states[j].position[i];
  }

  //conversione da r [km],v[km/s] a ra,dec topocentriche [rad]
  astro::rv_tradec(r_sat_eci, v_sat_eci, r_obs_eci, eTo, ra, dec); //vallado
  printf("RA  = %f\n",ra);
  printf("DEC = %f\n",dec);

  char stringa[50];
  sprintf(stringa,"%f %f\n", ra, dec);
  myfile << stringa;
  myfile.flush();

  //PROVA CODICE NOSTRO (CORRETTO)
  astro::tRaDec2AzEl(ra,dec,latitude,longitude,jday, Az, El);
  printf("Az  = %f\n",Az);
  printf("El  = %f\n",El);

  //PROVA VALLADO radec_azel (CORRETTO, piccole differenze con il nostro codice)
  ra_rad = astro::radians(ra);
  dec_rad = astro::radians(dec);
  astro::radec2azel(ra_rad, dec_rad, latgd, jday, longitude, Az2, El2);
  printf("Az2_vall  = %f\n",Az2);
  printf("El2_vall  = %f\n",El2);


  //condizioni di visibilità
  int jday_int = jday;
  double jday_frac = jday - jday_int;
  double El_min = 10; //deg

  // Posizione satellite e osservatorio in raggi terrestri
  double r_sat[3];
  r_sat[0] = r_sat_eci[0]/Rt;
  r_sat[1] = r_sat_eci[1]/Rt;
  r_sat[2] = r_sat_eci[2]/Rt;
  double r_obs[3];
  r_obs[0] = r_obs_eci[0]/Rt;
  r_obs[1] = r_obs_eci[1]/Rt;
  r_obs[2] = r_obs_eci[2]/Rt;
  
  bool ill_sat = astUtils::light(r_sat,(double)jday_int,jday_frac,'e');
  bool ill_obs = astUtils::light(r_obs,(double)jday_int,jday_frac,'e');
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

myfile.close();

    
#endif

return 0;
  
}
