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

#define SUN_VALLADO
//#define SUN_KDE_QUICK
//define SUN_KDE_PRECISE

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
    
    //double azimuth  = 0.0;
    //double altitude = 0.0;
    
    //double rightAscension = 0.0;
    //double declination    = 0.0;
    
    //inline void finalize(const astro::observatory & station) { computeRaDec(station); }
    
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
    // la posizione del sole dovrebbe essere in J2000
    // dt di integrazione e' in secondi
    
  private:
    
    Sun() { }
    
#ifdef SUN_KDE_PRECISE
    // Used in _position taken from the KDE AstroLib
    static double t, a, d, uu, dl, dr, db, m2, m3, m4, m5, m6;
    static double c[9], c3[9], s[9], s3[9];
#endif
    
  public:
    
    /*****************************************************************************/
    // position
    /*****************************************************************************/
    inline static void position(double jDay, astro::SunState & state, int crs = CRS::ECI) {
      
      _position(jDay, state, crs);
      
    }

    /*****************************************************************************/
    // orbit
    /*****************************************************************************/
    static void orbit(double jDayStart, double jDayStop, double integrationTimeSec, std::vector<astro::SunState> & states, int crs = CRS::ECI) {
      
      // converto il tempo di integrazione da secondi in jDay
      double integrationTimeJD = astro::Date::convert(integrationTimeSec, astro::Date::FROM_SECOND_TO_JD);
      
      // mi calcolo il numero di samples che dovro' fare
      std::size_t samples = round((jDayStop - jDayStart) / integrationTimeJD);
      
      // alloco lo spazio
      states.resize(samples);
      
      // propago
      for(std::size_t step=0; step<samples; ++step) {
        
        states[step].jDay = jDayStart + (integrationTimeJD * step);
        
        _position(states[step].jDay, states[step], crs);
        
      }
      
    }
    
  private:
    
#ifdef SUN_VALLADO
    /*****************************************************************************/
    // _position
    /*****************************************************************************/
    inline static void _position(double jDay, astro::SunState & state, int crs = CRS::ECI) {
      
      double rtasc, decl;
      
      double tmpCoord[3];
      
      // chiamo la funzione di vallado
      ast2Body::sun(jDay, 0.0, tmpCoord, rtasc, decl);
      
      state.position[0] *= 149597870.7;
      state.position[1] *= 149597870.7;
      state.position[2] *= 149597870.7;
      
      if(crs != CRS::TEME) {

        if(crs == CRS::ECEF)
          astro::eci2ecef(tmpCoord, jDay, state.position);
        
        if(crs == CRS::TEME)
          astro::eci2teme(tmpCoord, jDay, state.position);
     
      } else {
        
        state.position[0] = tmpCoord[0];
        state.position[1] = tmpCoord[1];
        state.position[2] = tmpCoord[2];

      }
      
      
    }
#endif

    
#ifdef SUN_KDE_QUICK
    /*****************************************************************************/
    // _position
    // (Taken fro KDE AstroLib)
    // Low precision position valid only between 1950 and 2050
    // Ecliptic coordinates (in A.U.) of the sun for
    // Equinox of Date given in Julian centuries since J2000
    // NOTE: jDay must be in standard Julian date the position is in ECI (J2000)
    /*****************************************************************************/
    inline static void _position(double jDay, astro::SunState & state, int crs = CRS::ECI) {

      // converto from jDay to Julian centuries since J2000
      // ovvero jDay - jDay(01/01/2000 12:00 PM) / 36525.0
      double t = (jDay - 2451545.0) / 36525.0;
      
      double n, g, l;
      
      n = 36525.0 * t;
      l = 280.460 + 0.9856474 * n; // mean longitude
      g = 357.528 + 0.9856003 * n; // mean anomaly
      
      l = 2.0 * M_PI * fmod(l/360.0,1.0);
      g = 2.0 * M_PI * fmod(g/360.0,1.0);
      l = l + (1.915 * sin(g) + 0.020 * sin(2.0*g)) * 1.74532925199e-2; // ecl.long.
      g = 1.00014 - 0.01671 * cos(g) - 0.00014 * cos(2.0*g);            // radius
      
      double rEcli[3];
      
      rEcli[0] = g * cos(l) * 149597870.7;
      rEcli[1] = g * sin(l) * 149597870.7;
      rEcli[2] = 0;
      
      astUtils::eclequ(t, rEcli, state.position);
      
      //state.position[0] *= 149597870.7;
      //state.position[1] *= 149597870.7;
      //state.position[2] *= 149597870.7;
      
      //FIXME: not used crs
      
    }
#endif
    
#ifdef SUN_KDE_PRECISE
    /*****************************************************************************/
    // _position
    // (Taken fro KDE AstroLib)
    // Ecliptic coordinates (in A.U.) and velocity (in A.U./day)
    // of the sun for Equinox of Date given in Julian centuries since J2000
    // NOTE: jDay must be in standard Julian date the position is in ECI (J2000)
    /*****************************************************************************/
    inline static void _position(double jDay, astro::SunState & state, int crs = CRS::ECI) {
      
      // converto from jDay to Julian centuries since J2000
      // ovvero jDay - jDay(01/01/2000 12:00 PM) / 36525.0
      t = (jDay - 2451545.0) / 36525.0;
      
      const double  p2 = 2.0 * M_PI;
      
      dl = 0; dr = 0.0; db = 0.0;
      
      m2 = p2 * fmod(0.1387306 + 162.5485917 * t, 1.0);
      m3 = p2 * fmod(0.9931266 + 99.9973604  * t, 1.0);
      m4 = p2 * fmod(0.0543250 + 53.1666028  * t, 1.0);
      m5 = p2 * fmod(0.0551750 + 8.4293972   * t, 1.0);
      m6 = p2 * fmod(0.8816500 + 3.3938722   * t, 1.0);
      d  = p2 * fmod(0.8274    + 1236.8531   * t, 1.0);
      a  = p2 * fmod(0.3749    + 1325.5524   * t, 1.0);
      uu = p2 * fmod(0.2591    + 1342.2278   * t, 1.0);
      
      c3[1] = 1.0;     s3[1] = 0.0;
      c3[2] = cos(m3); s3[2] = sin(m3);
      c3[0] = c3[2];   s3[0] = -s3[2];
      
      for(int i=3; i<9; i++)
        addthe(c3[i-1], s3[i-1], c3[2], s3[2], c3[i], s3[i]);
      
      pertven(); pertmar(); pertjup(); pertsat(); pertmoo();
      
      dl = dl +  6.4 * sin(p2*(0.6983 + 0.0561*t))
              + 1.87 * sin(p2*(0.5764 + 0.4174*t))
              + 0.27 * sin(p2*(0.4189 + 0.3306*t))
              + 0.20 * sin(p2*(0.3581 + 2.4814*t));
      
      double l = p2 * fmod(0.7859453 + m3/p2 + ((6191.2+1.1*t)*t+dl)/1296.0E3, 1.0);
      double r = 1.0001398 - 0.0000007 * t + dr * 1E-6;
      double b = db * 4.8481368111E-6;
      
      double cl = cos(l);
      double sl = sin(l);
      double cb = cos(b);
      double sb = sin(b);
      
      double rEcli[3];
      
      rEcli[0] = r * cl * cb;
      rEcli[1] = r * sl * cb;
      rEcli[2] = r * sb;
      
      astUtils::eclequ(t, rEcli, state.position);
      
      state.position[0] *= 149597870.7;
      state.position[1] *= 149597870.7;
      state.position[2] *= 149597870.7;
      
      //FIXME: not used crs
      
    }

  private:
    
    inline static void addthe(double c1, double s1, double c2, double s2, double & cc, double & ss) {
      
      cc = c1 * c2 - s1 * s2;
      ss = s1 * c2 + c1 * s2;
      
    }
    
    void static term(int i1, int i, int it, double dlc, double dls, double drc, double drs, double dbc, double dbs) {
      
      double u, v;
      
      if(it == 0) addthe(c3[i1+1], s3[i1+1], c[i+8], s[i+8], u, v);
      else {
        u = u * t;
        v = v * t;
      }
      
      dl = dl + dlc*u + dls*v;
      dr = dr + drc*u + drs*v;
      db = db + dbc*u + dbs*v;
      
    }
    
    // Kepler terms and perturbations by Venus
    static void pertven() {
      
      c[8] = 1.0; s[8] = 0.0; c[7] = cos(m2); s[7] = -sin(m2);
      
      for(int i=7; i>2; i--)
        addthe(c[i], s[i], c[7], s[7], c[i-1], s[i-1]);
      
      term (1, 0, 0, -0.22, 6892.76, -16707.37, -0.54, 0.0, 0.0);
      term (1, 0, 1, -0.06, -17.35, 42.04, -0.15, 0.0, 0.0);
      term (1, 0, 2, -0.01, -0.05, 0.13, -0.02, 0.0, 0.0);
      term (2, 0, 0, 0.0, 71.98, -139.57, 0.0, 0.0, 0.0);
      term (2, 0, 1, 0.0, -0.36, 0.7, 0.0, 0.0, 0.0);
      term (3, 0, 0, 0.0, 1.04, -1.75, 0.0, 0.0, 0.0);
      term (0, -1, 0, 0.03, -0.07, -0.16, -0.07, 0.02, -0.02);
      term (1, -1, 0, 2.35, -4.23, -4.75, -2.64, 0.0, 0.0);
      term (1, -2, 0, -0.1, 0.06, 0.12, 0.2, 0.02, 0.0);
      term (2, -1, 0, -0.06, -0.03, 0.2, -0.01, 0.01, -0.09);
      term (2, -2, 0, -4.7, 2.9, 8.28, 13.42, 0.01, -0.01);
      term (3, -2, 0, 1.8, -1.74, -1.44, -1.57, 0.04, -0.06);
      term (3, -3, 0, -0.67, 0.03, 0.11, 2.43, 0.01, 0.0);
      term (4, -2, 0, 0.03, -0.03, 0.1, 0.09, 0.01, -0.01);
      term (4, -3, 0, 1.51, -0.4, -0.88, -3.36, 0.18, -0.1);
      term (4, -4, 0, -0.19, -0.09, -0.38, 0.77, 0.0, 0.0);
      term (5, -3, 0, 0.76, -0.68, 0.3, 0.37, 0.01, 0.0);
      term (5, -4, 0, -0.14, -0.04, -0.11, 0.43, -0.03, 0.0);
      term (5, -5, 0, -0.05, -0.07, -0.31, 0.21, 0.0, 0.0);
      term (6, -4, 0, 0.15, -0.04, -0.06, -0.21, 0.01, 0.0);
      term (6, -5, 0, -0.03, -0.03, -0.09, 0.09, -0.01, 0.0);
      term (6, -6, 0, 0.0, -0.04, -0.18, 0.02, 0.0, 0.0);
      term (7, -5, 0, -0.12, -0.03, -0.08, 0.31, -0.02, -0.01);
      
    }
    
    // Kepler terms and perturbations by Mars
    static void pertmar() {
      
      c[7] = cos(m4); s[7] = -sin(m4);
      
      for(int i=7; i>0; i--)
        addthe(c[i], s[i], c[7], s[7], c[i-1], s[i-1]);
      
      term (1, -1, 0, -0.22, 0.17, -0.21, -0.27, 0.0, 0.0);
      term (1, -2, 0, -1.66, 0.62, 0.16, 0.28, 0.0, 0.0);
      term (2, -2, 0, 1.96, 0.57, -1.32, 4.55, 0.0, 0.01);
      term (2, -3, 0, 0.4, 0.15, -0.17, 0.46, 0.0, 0.0);
      term (2, -4, 0, 0.53, 0.26, 0.09, -0.22, 0.0, 0.0);
      term (3, -3, 0, 0.05, 0.12, -0.35, 0.15, 0.0, 0.0);
      term (3, -4, 0, -0.13, -0.48, 1.06, -0.29, 0.01, 0.0);
      term (3, -5, 0, -0.04, -0.2, 0.2, -0.04, 0.0, 0.0);
      term (4, -4, 0, 0.0, -0.03, 0.1, 0.04, 0.0, 0.0);
      term (4, -5, 0, 0.05, -0.07, 0.2, 0.14, 0.0, 00);
      term (4, -6, 0, -0.1, 0.11, -0.23, -0.22, 0.0, 0.0);
      term (5, -7, 0, -0.05, 0.0, 0.01, -0.14,  0.0, 0.0);
      term (5, -8, 0, 0.05, 0.01, -0.02, 0.1, 0.0, 0.0);
      
    }
    
    // Kepler terms and perturbations by Jupiter
    static void pertjup() {
      
      c[7] = cos(m5); s[7] = -sin(m5);
      
      for(int i=7; i>4; i--)
        addthe(c[i], s[i], c[7], s[7], c[i-1], s[i-1]);
      
      term (1, -1, 0, 0.01, 0.07, 0.18, -0.02, 0.0, -0.02);
      term (0, -1, 0, -0.31, 2.58, 0.52, 0.34, 0.02, 0.0);
      term (1, -1, 0, -7.21, -0.06, 0.13, -16.27, 0.0, -0.02);
      term (1, -2, 0, -0.54, -1.52, 3.09, -1.12, 0.01, -0.17);
      term (1, -3, 0, -0.03, -0.21, 0.38, -0.06, 0.0, -0.02);
      term (2, -1, 0, -0.16, 0.05, -0.18, -0.31, 0.01, 0.0);
      term (2, -2, 0, 0.14, -2.73, 9.23, 0.48, 0.0, 0.0);
      term (2, -3, 0, 0.07, -0.55, 1.83, 0.25, 0.01,  0.0);
      term (2, -4, 0, 0.02, -0.08, 0.25, 0.06, 0.0, 0.0);
      term (3, -2, 0, 0.01, -0.07, 0.16, 0.04, 0.0, 0.0);
      term (3, -3, 0, -0.16, -0.03, 0.08, -0.64, 0.0, 0.0);
      term (3, -4, 0, -0.04, -0.01, 0.03, -0.17, 0.0, 0.0);
      
    }
    
    // Kepler terms and perturbations by Saturn
    static void pertsat() {
      
      c[7] = cos(m6); s[7] = -sin(m6);
      
      addthe(c[7], s[7], c[7], s[7], c[6], s[6]);
      
      term (0, -1, 0, 0.0, 0.32, 0.01, 0.0, 0.0, 0.0);
      term (1, -1, 0, -0.08, -0.41, 0.97, -0.18, 0.0, -0.01);
      term (1, -2, 0, 0.04, 0.1, -0.23, 0.1, 0.0, 0.0);
      term (2, -2, 0, 0.04, 0.1, -0.35, 0.13, 0.0, 0.0);
      
    }
    
    // corrections for Earth-Moon center of gravity
    static void pertmoo() {
      
      dl = dl + 6.45  * sin(d) - 0.42 * sin(d-a) + 0.18 * sin(d+a) + 0.17 * sin(d-m3) - 0.06 * sin(d+m3);
      dr = dr + 30.76 * cos(d) - 3.06 * cos(d-a) + 0.85 * cos(d+a) - 0.58 * cos(d+m3) + 0.57 * cos(d-m3);
      db = db + 0.576 * sin(uu);
      
    }
    
#endif

  }; /* class sun */
  
#ifdef SUN_KDE_PRECISE
  // Used in _position taken from the KDE AstroLib
  double Sun::t  = 0.0;
  double Sun::a  = 0.0;
  double Sun::d  = 0.0;
  double Sun::uu = 0.0;
  double Sun::dl = 0.0;
  double Sun::dr = 0.0;
  double Sun::db = 0.0;
  double Sun::m2 = 0.0;
  double Sun::m3 = 0.0;
  double Sun::m4 = 0.0;
  double Sun::m5 = 0.0;
  double Sun::m6 = 0.0;
  double Sun::c[9]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double Sun::c3[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double Sun::s[9]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double Sun::s3[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif
  
  
  
} /* namespace astro */

#endif /* _H_ASTRO_SUN_H_ */
