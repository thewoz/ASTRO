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
