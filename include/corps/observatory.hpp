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

#ifndef _H_ASTRO_OBSERVATORY_H_
#define _H_ASTRO_OBSERVATORY_H_

#include <cstdlib>
#include <cstdio>

#include <vector>
#include <string>

#include <vallado/astIOD.h>

#include "../utils.h"

#include "../converter.h"

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
  
  class sun;

  /*****************************************************************************/
  // class ObservatoryState
  /*****************************************************************************/
  class ObservatoryState {
    
  public:
    
    double jDay;
    
    double latitude  = 0.0;
    double longitude = 0.0;
    double height    = 0.0;
    
    double position[3];
    
    bool isInShadow = true;
    
    void finalize() { computeInShadow(); }
    
    void computeInShadow() {
      
      
    }
    
    inline void print(FILE * output = stdout)   { fprintf(output, "%f %f %f %f",   jDay, position[0], position[1], position[2]); }
    inline void println(FILE * output = stdout) { fprintf(output, "%f %f %f %f\n", jDay, position[0], position[1], position[2]); }
    
    inline bool operator < (const ObservatoryState & state) { return (jDay<state.jDay); }
    
  }; /* class ObservatoryState */
  
  /*****************************************************************************/
  // class Observatory
  /*****************************************************************************/
  class Observatory {

  private:

    // IN ECEF
    double coord[3];
    
  public:
    
    // in radianti
    double latitude  = 0.0;
    double longitude = 0.0;
    double height    = 0.0;

    
    /*****************************************************************************/
    // constructor
    /*****************************************************************************/
    //   latitude   [deg]
    //   longitude  [deg]
    //   height     [m]
    /*****************************************************************************/
    Observatory(double _latitude, double _longitude, double _height = 0.0, int crs = CRS::ECEF) {
      
      latitude  = _latitude;
      longitude = _longitude;
      height    = _height;
      
      //_convert(states[i].jDay, states[i].position);
      astro::lla2ecef(latitude, longitude, height, coord);

    }

    /*****************************************************************************/
    // position
    /*****************************************************************************/
    void positions(double jDayStart, double jDayEnd, double integrationTimeSec, std::vector<astro::ObservatoryState> & states, int crs = CRS::ECEF) {
      
      double integrationTimeJD = astro::Date::convert(integrationTimeSec, astro::Date::FROM_SECOND_TO_JD);
      
      std::size_t samples = ((jDayEnd - jDayStart) / integrationTimeJD) + 1;

      states.resize(samples);
      
      for(std::size_t i=0; i<samples; ++i) {
        
        states[i].jDay = jDayStart + (integrationTimeJD*i);
        
        states[i].position[0] = coord[0]; states[i].position[1] = coord[1]; states[i].position[2] = coord[2];
        
        //TOGLIERE SOLO SE NON HO FATTO INIT
        // METTERE INIT
        //_convert(states[i].jDay, states[i].position);
        astro::lla2ecef(latitude, longitude, height, coord);

      }
      
      if(crs != CRS::ECEF) {
        
        double dummy[3];
        double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;

        for(std::size_t i=0; i<samples; ++i) {

          astro::eopc::getParameters(states[i].jDay, 'l', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
          
          if(crs == CRS::TEME)
            astro::ecef2teme(states[i].position, dummy, dummy, states[i].position, dummy, dummy, ttt, jdut1+jdut1Frac, lod, xp, yp);
          
          if(crs == CRS::ECI)
            astro::ecef2eci(states[i].position, dummy, dummy, states[i].position, dummy, dummy, ttt, jdut1+jdut1Frac, lod, xp, yp);
          
        }
        
      }
      
    }

    // position
    /*****************************************************************************/
    void position(double jDay, double _coord[3], int crs = CRS::ECEF) {
      
      //VEDI SOPRA
      astro::lla2ecef(latitude, longitude, height, _coord);
      
      if(crs != CRS::ECEF) {
        
        double dummy[3];
        double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;
      
        astro::eopc::getParameters(jDay, 'l', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
        
        if(crs == CRS::TEME)
          astro::ecef2teme(_coord, dummy, dummy, _coord, dummy, dummy, ttt, jdut1+jdut1Frac, lod, xp, yp);
        
        if(crs == CRS::ECI)
          astro::ecef2eci(_coord, dummy, dummy, _coord, dummy, dummy, ttt, jdut1+jdut1Frac, lod, xp, yp);
        
      }
      
    }

    
  private:

    /*****************************************************************************/
    // _convert
    /*****************************************************************************/
    void _convert(double jDay, double r[3]) {

      double gst = astTime::gstime(jDay);

      double a = 6378.137;
      double b = 6356.7523142;

      double f = (a-b)/a;
      double equad = 2.0*f-f*f;

      double teta = gst + longitude*M_PI/180.0;
      double L = latitude * M_PI/180.0;

      double N = a / sqrt(1-equad*pow(sin(L),2));

      r[0] = N*cos(L)*cos(teta);
      r[1] = N*cos(L)*sin(teta);
      r[2] = (height/1000+(1-equad)*N)*sin(L);

    }
    
  }; /* class Observatory */
  
} /* namespace astro */

#endif /* _H_ASTRO_OBSERVATORY_H_ */
