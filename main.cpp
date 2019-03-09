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

//#define TESTEOPC
//#define TESTDATE
//#define TESTSATELLITE
#define TESTSUN
//#define TESTSTATION
//#define TESTATTITUDE

#include <cstdio>
#include <cstdlib>

#include <cerrno>
#include <cstring>

#include <string>

#include "astro.hpp"


/*****************************************************************************/
// main
/*****************************************************************************/
int main(int argc, char *argv[]) {

#ifdef TESTATTITUDE
  /*****************************************************************************/
  // Test Attitude
  /*****************************************************************************/
  
  //astro::attitude::test(); exit(0);
  
  // angoli in radienati
  double pitch = 0;
  double yaw   = 0;
  double roll  = 0;

  //pitch = mpl::geometry::Radians(1.0);
  //yaw   = mpl::geometry::Radians(1.0);
  //roll  = mpl::geometry::Radians(1.0);
  
  quaternion_t q(pitch, yaw, roll);

  // velocita in radianti al secondo
  double Wx = 0.0;
  double Wy = 0.0;
  double Wz = 0.0;
  
  //Wx = mpl::geometry::Radians(1.0);
  //Wy = mpl::geometry::Radians(1.0);
  //Wz = mpl::geometry::Radians(1.0);
  
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
  
#ifdef TESTDATE
  /*****************************************************************************/
  // Test Date class
  /*****************************************************************************/
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
  
#ifdef TESTEOPC
  /*****************************************************************************/
  // Test Epoc class
  /*****************************************************************************/
  fprintf(stderr, "Test Epoc class:\n\n");

  astro::eopc::init(1, true);
  
  double xp, yp, lod, dpsi, deps, jdut1, jdut1Frac, ttt;
  
  astro::eopc::getParameters(astro::Date(29, 1, 2017, 13, 34, 1).getJDay(), 11, xp, yp, lod, dpsi, deps, jdut1, jdut1Frac, ttt);
  
  printf("\n xp %e\n yp %e\n lod %e\n dpsi %e\n deps %e\n jdut1 %e\n jdut1Frac %e\n ttt %e\n\n\n", xp, yp, lod, dpsi, deps, jdut1, jdut1Frac, ttt);
  
#endif
  
  
#ifdef TESTSATELLITE
  /*****************************************************************************/
  // Test Orbit Propagator from Satellite
  /*****************************************************************************/
  //while(1)
  {
    
    fprintf(stderr, "Test Satellite:\n\n");
  
    std::string outStrSatellite = "/Users/thewoz/Desktop/satellite.dat";
  
    FILE * outputSatellite = fopen(outStrSatellite.c_str(), "w");
  
    if(outputSatellite==NULL){
      fprintf(stderr, "error in opening the file '%s': %s\n", outStrSatellite.c_str(), strerror(errno));
      abort();
    }

    char TLE_line1[] = "1 25544U 98067A   19047.93785069  .00001237  00000-0  26703-4 0  9993";
    char TLE_line2[] = "2 25544  51.6407 239.1790 0005867  36.2033  74.9248 15.53268833156571";
    
    std::vector<astro::SatelliteState> states;

    astro::Satellite(TLE_line1, TLE_line2).orbit(astro::Date(17,2,2019,11,00,00), astro::Date(17,2,2019,12,00,00), 60, states, CRS::ECI);

    for(std::size_t i=0; i<states.size(); ++i){
      states[i].println(outputSatellite);
      states[i].println(stdout);
    }
    
    fclose(outputSatellite);
  
    fprintf(stderr, "\n\n\n");
  
  }
  
#endif
  
  
#ifdef TESTSUN
  /*****************************************************************************/
  // Test Orbit Propagator for Sun
  /*****************************************************************************/
  {
    
    std::string outStrSun = "/Users/thewoz/Desktop/sun.dat";
    
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
      fprintf(stdout, "%s ", astro::Date(states[i].jDay).toGregorianString());
      states[i].println(stdout);
    }
    
    fclose(outputSun);
    
    fprintf(stderr, "\n\n\n");
    
  }

#endif

  
#ifdef TESTSTATION
  /*****************************************************************************/
  // Test Orbit Propagator for Station
  /*****************************************************************************/
  {
    
    fprintf(stderr, "Test Station:\n\n");

    std::string outStrStation = "/Users/thewoz/Desktop/station.dat";
    
    FILE * outputStation = fopen(outStrStation.c_str(), "w");
    
    if(outputStation==NULL){
      fprintf(stderr, "error in opening the file '%s': %s\n", outStrStation.c_str(), strerror(errno));
      abort();
    }
    
    double latitude  = Radians(41.9577777777);
    double longitude = Radians(12.505555555);
    double altitude  = 16.0;
    
    double dt = 60;

    astro::Observatory collepardo(latitude, longitude, altitude);
    
    std::vector<astro::ObservatoryState> states;
    
    collepardo.positions(astro::Date(12,2,2018,4,19,59).getJDay(), astro::Date(12,2,2018,5,20,59).getJDay(), dt, states, CRS::TEME);
    
    for(std::size_t i=0; i<states.size(); ++i){
      fprintf(stdout, "%s ", astro::Date(states[i].jDay).toGregorianString());
      states[i].println(stdout);
      fprintf(outputStation, "%s ", astro::Date(states[i].jDay).toGregorianString());
      states[i].println(outputStation);
    }
    
    fclose(outputStation);
    
    fprintf(stderr, "\n\n\n");

  }
    
#endif
  
  return 0;
  
}



