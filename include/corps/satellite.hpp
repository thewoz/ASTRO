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

#ifndef _H_ASTRO_SATELLITE_H_
#define _H_ASTRO_SATELLITE_H_

#include <cstdio>
#include <cstdlib>

#include <vector>
#include <string>

//****************************************************************************/
// namespace astro
//****************************************************************************/
namespace astro {
  
  class observatory;
  class sun;
  
  //****************************************************************************/
  // class SatelliteState
  //****************************************************************************/
  class SatelliteState {
    
  public:
    
    double jDay;
    
    double position[3];
    double velocity[3];
    
    //double azimuth  = 0.0;
    //double altitude = 0.0;
    
    //double rightAscension = 0.0;
    //double declination    = 0.0;
    
    //bool isInShadow = true;
    
    inline void print(FILE * output = stdout)   { fprintf(output, "%f %f %f %f",   jDay, position[0], position[1], position[2]); }
    inline void println(FILE * output = stdout) { fprintf(output, "%f %f %f %f\n", jDay, position[0], position[1], position[2]); }
    
    inline bool operator < (const SatelliteState & state) { return (jDay<state.jDay); }
    
  }; /* class SatelliteState */
  
  
  
  //****************************************************************************/
  // class Satellite
  //****************************************************************************/
  class Satellite {
    
  public:
    
    // satellite name
    std::string name;
    
    Satellite() : isInit(false) { }
    
    Satellite(FILE * input) : isInit(true) { tle.init(input); name = tle.name; }
    
    Satellite(const char * input) : isInit(true) {
      FILE * tleFile = fopen(input, "r");
      if(tleFile == NULL) {
        fprintf(stderr, "error in opening file \"%s\" \n", input);
        abort();
      }
      tle.init(tleFile); name = tle.name;
      fclose(tleFile);
    }

    Satellite(const char * TLE_line1, const char * TLE_line2) : isInit(true) {
      tle.init(TLE_line1, TLE_line2);
      name = tle.name;
    }

    //****************************************************************************/
    // orbit
    //****************************************************************************/
    void orbit(double jDayStart, double jDayStop, double integrationTimeSec, std::vector<astro::SatelliteState> & states, int crs = CRS::TEME)
    //void orbit(astro::Date startDate, astro::Date stopDate, double integrationTimeSec, std::vector<astro::SatelliteState> & states, int crs = CRS::TEME) 
    {
      
      if(!isInit){
        fprintf(stderr, "satellite must init before\n");
        abort();
      }
      
      // tempo di propagazione (in minuti)
      //double propagationTimeMin = astro::Date::difference(stopDate,startDate,astro::Date::MINUTES);
      
      // converto il passo di integrazione (in minuti)
      //double integrationTimeMin = astro::Date::convert(integrationTimeSec, astro::Date::FROM_SECOND_TO_MINUTE);
      
      // tempo d'inzio della propagazione dal rilascio (in minuti)
      //double startTimeMin = astro::Date::difference(startDate, tle.releaseDate, astro::Date::MINUTES);
      

      double propagationTimeMin = (jDayStop - jDayStart)*(24.0*60.0);
      
      double integrationTimeMin = astro::Date::convert(integrationTimeSec, astro::Date::FROM_SECOND_TO_MINUTE);
      
      double startTimeMin = (jDayStart - tle.releaseDate)*(24.0*60.0);

      // numero di campionamenti da fare
      // FIXME: sto piu uno non mi torna
      int samples = round(propagationTimeMin / integrationTimeMin);
      
      // alloco lo spazio
      states.resize(samples);
      
      // converto il tempo di inizio in (jDay)
      //double startTimeJD = startDate.getJDay();
      double startTimeJD = jDayStart;

      // converto il tempo di integrazione (jDay)
      double integrationTimeJD = astro::Date::convert(integrationTimeSec, astro::Date::FROM_SECOND_TO_JD);
      
      for(std::size_t step=0; step<samples; ++step) {
        
        double sinceTimeMin = startTimeMin + (integrationTimeMin * step);
        
        states[step].jDay = startTimeJD + (integrationTimeJD * step);
        
        _position(states[step].jDay, sinceTimeMin, states[step].position, states[step].velocity, crs);
        
      } // for(step)
      
    }
    
    //****************************************************************************/
    // position
    //****************************************************************************/
    // NOTE: si potrebbe chiamare
    void position(double jDay, double coord[3], int crs = CRS::TEME) {
      
      if(!isInit){
        fprintf(stderr, "satellite must init before\n");
        abort();
      }
      
      double sinceTimeMin = astro::Date::difference(jDay, tle.releaseDate, astro::Date::MINUTES);
      
      double dummy[3];
      
      _position(jDay, sinceTimeMin, coord, dummy, crs);
      
    }
    
  private:
    
    // satellite tle
    astro::TLE tle;
    
    // true if the satellite in inited
    bool isInit;
    
    //****************************************************************************/
    // _position
    //****************************************************************************/
    void _position(double jDay, double sinceTimeMin, double coord[3], double vel[3], int crs = CRS::TEME) {
      
      double tmpCoord[3]; double tmpVel[3];
      
      astro::sgp4(tle.satrec, sinceTimeMin, tmpCoord, tmpVel);
      
      if(crs != CRS::TEME) {
        
        if(crs == CRS::ECI)
          astro::teme2eci(tmpCoord, tmpVel, jDay, coord, vel);
        
        if(crs == CRS::ECEF)
          astro::teme2ecef(tmpCoord, tmpVel, jDay, coord, vel);
        
      } else {
        
        coord[0] = tmpCoord[0];
        coord[1] = tmpCoord[1];
        coord[2] = tmpCoord[2];
        
        vel[0] = tmpVel[0];
        vel[1] = tmpVel[1];
        vel[2] = tmpVel[2];
        
      }
      
      
    }
    
  }; /* class satellite */
  
} /* namespace astro */

#endif /* _H_ASTRO_SATELLITE_H_ */
