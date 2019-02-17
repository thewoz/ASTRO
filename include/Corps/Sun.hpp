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

#ifndef _H_ASTRO_SUN_H_
#define _H_ASTRO_SUN_H_

#include <cstdlib>
#include <cstdio>

#include <vector>
#include <string>

#include <vallado/ast2Body.h>

#include "../utils.h"
#include "../converter.h"

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
  
  class observatory;

  /*****************************************************************************/
  // class SunState
  /*****************************************************************************/
  class SunState {
    
  public:
    
    double jDay;
    
    double position[3];
    
    double azimuth  = 0.0;
    double altitude = 0.0;
    
    double rightAscension = 0.0;
    double declination    = 0.0;
    
    inline void finalize(const astro::observatory & station) { computeRaDec(station); }
    
    inline void computeRaDec(const astro::observatory & station) { }
    
    inline void print(FILE * output = stdout)   { fprintf(output, "%f %f %f %f",   jDay, position[0], position[1], position[2]); }
    inline void println(FILE * output = stdout) { fprintf(output, "%f %f %f %f\n", jDay, position[0], position[1], position[2]); }
        
    inline bool operator < (const SunState & state) { return (jDay<state.jDay); }
    
  }; /* class SunState */
  
  /*****************************************************************************/
  // class sun
  /*****************************************************************************/
  class Sun {
    
    // NOTE:
    // la posizione del sole e' in J2000
    // dt di integrazione e' in secondi
    
  private:
    
    Sun() { }
    
  public:
    
    /*****************************************************************************/
    // operator ()
    /*****************************************************************************/
    inline astro::SunState operator () (double jDay, int crs = CRS::ECI) const { return position(jDay,crs); }
    
    /*****************************************************************************/
    // position
    /*****************************************************************************/
    inline static astro::SunState position(double jDay, int crs = CRS::ECI) { astro::SunState state; _position(jDay, state, crs); return state; }
    inline static void position(double jDay, astro::SunState & state, int crs = CRS::ECI) { _position(jDay, state, crs); }

    /*****************************************************************************/
    // orbit
    /*****************************************************************************/
    static void orbit(double jDayStart, double jDayStop, double integrationTime, std::vector<astro::SunState> & states, int crs = CRS::ECI) {
      _orbit(jDayStart, jDayStop, integrationTime, states, crs);
    }
    
    static std::vector<astro::SunState> orbit(double jDayStart, double jDayStop, double integrationTime, int crs = CRS::ECI) {
      std::vector<astro::SunState> states;
      _orbit(jDayStart, jDayStop, integrationTime, states, crs);
      return states;
    }
    
  private:

    /*****************************************************************************/
    // _orbit
    /*****************************************************************************/
    static void _orbit(double jDayStart, double jDayStop, double integrationTime, std::vector<astro::SunState> & states, int crs = CRS::ECI) {
      
      astro::eopc::init();
      
      double dummy[3];

      // converto il tempo di integrazione da secondi in jDay
      double integrationTimeJD = astro::Date::convert(integrationTime, astro::Date::FROM_SECOND_TO_JD);
      
      // mi calcolo il numero di samples che dovro' fare
      std::size_t samples = ((jDayStop - jDayStart) / integrationTimeJD) + 1;
      
      // alloco lo spazio
      states.resize(samples);
      
      double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;

      // propago
      for(std::size_t step=0; step<samples; ++step) {
        
        states[step].jDay = jDayStart + (integrationTimeJD * step);

        _position(states[step].jDay, states[step]);
        
        if(crs == CRS::ECEF) {
          
          astro::eopc::getParameters(states[step].jDay, 'l', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
 
          astro::eci2ecef(&states[step].position[0], dummy, dummy, &states[step].position[0], dummy, dummy, ttt, jdut1+jdut1Frac, lod, xp, yp);
          
        }
        
        if(crs == CRS::TEME) {
          
          astro::eopc::getParameters(states[step].jDay, 11, xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
          
          astro::eci2teme(&states[step].position[0], dummy, dummy, &states[step].position[0], dummy, dummy, ttt, ddpsi, ddeps);
          
        }
        
      }
      
    }
    
    /*****************************************************************************/
    // _position
    /*****************************************************************************/
    inline static void _position(double jDay, astro::SunState & state, int crs = CRS::ECI) {
      
      double rtasc, decl;
      
      // chiamo la funzione di vallado
      ast2Body::sun(jDay, 0.0, state.position, rtasc, decl);
      
      state.position[0] *= 149597870.7;
      state.position[1] *= 149597870.7;
      state.position[2] *= 149597870.7;
      
//      if(crs == CRS::TEME) ;
//      if(crs == CRS::ECEF) ;
      
    }
    
  }; /* class sun */
  
} /* namespace astro */

#endif /* _H_ASTRO_SUN_H_ */
