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


//****************************************************************************/
// namespace astro
//****************************************************************************/
namespace astro {
  
  class sun;

  //****************************************************************************/
  // class ObservatoryState
  //****************************************************************************/
  class ObservatoryState {
    
  public:
    
    double jDay;
    
    double latitude  = 0.0;
    double longitude = 0.0;
    double height    = 0.0;
    
    double position[3];
    
    //bool isInShadow = true;
    
    //void finalize() { computeInShadow(); }
    
    inline void print(FILE * output = stdout)   { fprintf(output, "%f %f %f %f",   jDay, position[0], position[1], position[2]); }
    inline void println(FILE * output = stdout) { fprintf(output, "%f %f %f %f\n", jDay, position[0], position[1], position[2]); }
    
    inline bool operator < (const ObservatoryState & state) { return (jDay<state.jDay); }
    
  }; /* class ObservatoryState */
  
  //****************************************************************************/
  // class Observatory
  //****************************************************************************/
  class Observatory {

  private:

    // IN ECEF
    double coord[3];
    
  public:
    
    // in radianti
    double latitude  = 0.0;
    double longitude = 0.0;
    
    // in km
    double height    = 0.0;
    
    //****************************************************************************/
    // constructor
    //****************************************************************************/
    //   latitude   [deg]
    //   longitude  [deg]
    //   height     [m]
    /*****************************************************************************/
    Observatory(double _latitude, double _longitude, double _height = 0.0) {
      
      // Passo da gradi a radianti
      latitude  = astro::radians(_latitude);
      longitude = astro::radians(_longitude);
      
      // Passo da metri a chilometri
      height    = _height / 1000.0;
            
      astIOD::site(latitude, longitude, height, coord);
      
    }

    //****************************************************************************/
    // orbit
    //****************************************************************************/
    void orbit(double jDayStart, double jDayStop, double integrationTimeSec, std::vector<ObservatoryState> & states, int crs = CRS::ECEF) {
      
      double integrationTimeJD = astro::Date::convert(integrationTimeSec, astro::Date::FROM_SECOND_TO_JD);
      
      std::size_t samples = round((jDayStop - jDayStart) / integrationTimeJD);

      states.resize(samples);
      
      for(std::size_t i=0; i<samples; ++i) {
        
        states[i].jDay = jDayStart + (integrationTimeJD*i);
        
        _position(states[i].jDay, states[i].position, crs);

      }
      
    }

    //****************************************************************************/
    // position
    //****************************************************************************/
    void position(double jDay, double _coord[3], int crs = CRS::ECEF) {
      
      _position(jDay, _coord, crs);
      
    }

  private:
    
    //****************************************************************************/
    // _position
    //****************************************************************************/
    void _position(double jDay, double _coord[3], int crs = CRS::ECEF) {
      
      if(crs == CRS::ECEF) {
        _coord[0] = coord[0];
        _coord[1] = coord[1];
        _coord[2] = coord[2];
      }
                              
      if(crs == CRS::TEME)
        astro::ecef2teme(coord, jDay, _coord);
        
      if(crs == CRS::ECI)
        astro::ecef2eci(coord, jDay, _coord);
              
    }
    
  }; /* class Observatory */
  
} /* namespace astro */

#endif /* _H_ASTRO_OBSERVATORY_H_ */
