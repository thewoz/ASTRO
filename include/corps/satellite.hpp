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

#ifndef _H_ASTRO_SATELLITE_H_
#define _H_ASTRO_SATELLITE_H_

#include <cstdio>
#include <cstdlib>

#include <vector>
#include <string>

#include "../utils.h"
#include "../propagator.h"
#include "../converter.h"

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
  
  class observatory;
  class sun;
  
  /*****************************************************************************/
  // class SatelliteState
  /*****************************************************************************/
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
  
  
  
  /*****************************************************************************/
  // class Satellite
  /*****************************************************************************/
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

    /*****************************************************************************/
    // orbit
    /*****************************************************************************/
    void orbit(astro::Date startDate, astro::Date stopDate, double integrationTimeSec, std::vector<SatelliteState> & states, int crs = CRS::TEME) {
      _orbit(startDate, stopDate, integrationTimeSec, states, crs);
    }
    
    /*****************************************************************************/
    // orbit
    /*****************************************************************************/
    std::vector<SatelliteState> orbit(astro::Date startDate, astro::Date stopDate, double integrationTimeSec, int crs = CRS::TEME) {
      
      std::vector<SatelliteState> states;
      
      _orbit(startDate, stopDate, integrationTimeSec, states, crs);
    
      return states;
      
    }
    
    /*****************************************************************************/
    // position
    /*****************************************************************************/
    void position(double jDay, double coord[3], int crs = CRS::TEME) {
      
      if(!isInit){
        fprintf(stderr, "satellite must init before\n");
        abort();
      }
      
      astro::eopc::init();

      double startTimeMin = Date::difference(jDay, tle.releaseDate, Date::MINUTES);
      
      double dummy[3];
      
      astro::sgp4(tle.satrec, startTimeMin, coord, dummy);

      if(crs != CRS::TEME) {
        
        double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;
        
        astro::eopc::getParameters(jDay, 'l', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
        
        if(crs == CRS::ECI)
          astro::teme2eci(coord, dummy, coord, dummy, ttt, ddpsi, ddeps);
        
        if(crs == CRS::ECEF)
          astro::teme2ecef(coord, dummy, coord, dummy, ttt, jdut1+jdut1Frac, lod, xp, yp);
        
      } // crs != CRS::TEME
      
    }
    
  private:
    
    // satellite tle
    astro::TLE tle;
    
    // true if the satellite in inited
    bool isInit;
    
    /*****************************************************************************/
    // _orbit
    /*****************************************************************************/
    void _orbit(astro::Date & startDate, astro::Date & stopDate, double integrationTimeSec, std::vector<SatelliteState> & states, int crs) {
      
      if(!isInit){
        fprintf(stderr, "satellite must init before\n");
        abort();
      }
      
      astro::eopc::init();
      
      // tempo di propagazione (in minuti)
      double propagationTimeMin = Date::difference(stopDate, startDate, Date::MINUTES);
      
      // converto il tempo di integrazione (in minuti)
      double integrationTimeMin = astro::Date::convert(integrationTimeSec, astro::Date::FROM_SECOND_TO_MINUTE);
      
      // tempo d'inzio della propagazione dal rilascio (in minuti)
      double startTimeMin = Date::difference(startDate, tle.releaseDate, Date::MINUTES);
      
      // numero di campionamenti da fare
      // FIXME: sto piu uno non mi torna
      int samples = (propagationTimeMin / integrationTimeMin) + 1;
           
      // alloco lo spazio
      states.resize(samples);
      
      // converto il tempo di inizio in (jDay)
      double startTimeJD = startDate.getJDay();
      
      // converto il tempo di integrazione (jDay)
      double integrationTimeJD = astro::Date::convert(integrationTimeSec, astro::Date::FROM_SECOND_TO_JD);
      
      for(std::size_t step=0; step<samples; ++step) {
        
        double sinceTimeMin = startTimeMin + (integrationTimeMin * step);
                
        states[step].jDay = startTimeJD + (integrationTimeJD * step);

        astro::sgp4(tle.satrec, sinceTimeMin, &states[step].position[0], &states[step].velocity[0]);
        
        if(crs != CRS::TEME) {
        
          double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;

          astro::eopc::getParameters(states[step].jDay, 'l', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);

          if(crs == CRS::ECI)
            astro::teme2eci(&states[step].position[0], &states[step].velocity[0], &states[step].position[0], &states[step].velocity[0], ttt, ddpsi, ddeps);
        
          if(crs == CRS::ECEF)
            astro::teme2ecef(&states[step].position[0], &states[step].velocity[0], &states[step].position[0], &states[step].velocity[0], ttt, jdut1+jdut1Frac, lod, xp, yp);
        
        } // crs != CRS::TEME
      
      } // for(step)

    }
    
    /*****************************************************************************/
    // _position
    /*****************************************************************************/
    inline void _position(double jDay, astro::SatelliteState & state, int crs) {
      
      if(!isInit){
        fprintf(stderr, "satellite must init before\n");
        abort();
      }
      
      // tempo d'inzio della propagazione dal rilascio (in minuti)
      double sinceTimeMin = Date::difference(jDay, tle.releaseDate, Date::MINUTES);
      
      state.jDay = jDay;
      
      astro::sgp4(tle.satrec, sinceTimeMin, &state.position[0], &state.velocity[0]);
      
      if(crs != CRS::TEME) {
        
        double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;
        
        astro::eopc::getParameters(state.jDay, 'l', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
        
        if(crs == CRS::ECI)
          astro::teme2eci(&state.position[0], &state.velocity[0], &state.position[0], &state.velocity[0], ttt, ddpsi, ddeps);
        
        if(crs == CRS::ECEF)
          astro::teme2ecef(&state.position[0], &state.velocity[0], &state.position[0], &state.velocity[0], ttt, jdut1+jdut1Frac, lod, xp, yp);
        
      } // crs != CRS::TEME
      
    }
    
  }; /* class satellite */
  
} /* namespace astro */

#endif /* _H_ASTRO_SATELLITE_H_ */
