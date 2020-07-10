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

#ifndef _H_ASTRO_TEME2ECEF_H
#define _H_ASTRO_TEME2ECEF_H

#include <cstdio>
#include <cstdlib>

#include <vector>


//****************************************************************************/
// namespace astro
//****************************************************************************/
namespace astro {

  //****************************************************************************/
  // teme2ecef
  //****************************************************************************/
  //
  // This function trsnforms a vector from the true equator mean equniox frame (teme),
  // to an earth fixed (itrf) frame.
  // The results take into account the effects of sidereal time, and polar motion.
  //
  //****************************************************************************/
  void teme2ecef(double rteme[3], double vteme[3], double ateme[3],
                 double jDay,
                 double recef[3], double vecef[3], double aecef[3]) {
    
    double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;
    
    astro::eopc::getParameters(jDay, 'l', 'f', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
    
    coordFK5::teme_ecef(rteme, vteme, ateme, edirection::eTo, recef, vecef, aecef, ttt, jdut1+jdut1Frac, lod, xp, yp);
    
  }
  
  //****************************************************************************/
  // teme2ecef
  //****************************************************************************/
  void teme2ecef(double rteme[3], double vteme[3],
                 double jDay,
                 double recef[3], double vecef[3]) {
    
    double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;
    
    astro::eopc::getParameters(jDay, 'l', 'f', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
    
    double dummy[3];

    coordFK5::teme_ecef(rteme, vteme, dummy, edirection::eTo, recef, vecef, dummy, ttt, jdut1+jdut1Frac, lod, xp, yp);
    
  }
  
  //****************************************************************************/
  // teme2ecef
  //****************************************************************************/
  void teme2ecef(double rteme[3],
                 double jDay,
                 double recef[3]) {
    
    double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;
    
    astro::eopc::getParameters(jDay, 'l', 'f', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
    
    double dummy[3];

    coordFK5::teme_ecef(rteme, dummy, dummy, edirection::eTo, recef, dummy, dummy, ttt, jdut1+jdut1Frac, lod, xp, yp);
    
  }
  
  
  
  //****************************************************************************/
  // teme2eci
  //****************************************************************************/
  // This function transforms a vector from the true equator mean equinox system,
  // (teme) to the mean equator mean equinox (j2000) system.
  //
  //****************************************************************************/
  void teme2eci(double rteme[3], double vteme[3], double ateme[3],
                double jDay,
                double reci[3],  double veci[3],  double aeci[3]) {
    
    double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;
    
    astro::eopc::getParameters(jDay, 'l', 'f', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
        
    coordFK5::teme_eci(rteme, vteme, ateme, edirection::eTo, reci, veci, aeci, astro::iau80::get(), ttt, ddpsi, ddeps);
    
  }
  
  //****************************************************************************/
  // teme2eci
  //****************************************************************************/
  void teme2eci(double rteme[3], double vteme[3],
                double jDay,
                double reci[3],  double veci[3]) {
    
    double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;
    
    astro::eopc::getParameters(jDay, 'l', 'f', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
    
    double dummy[3];
    
    coordFK5::teme_eci(rteme, vteme, dummy, edirection::eTo, reci, veci, dummy, astro::iau80::get(), ttt, ddpsi, ddeps);
    
  }
  
  //****************************************************************************/
  // teme2eci
  //****************************************************************************/
  void teme2eci(double rteme[3],
                double jDay,
                double reci[3]) {
    
    double xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt;
    
    astro::eopc::getParameters(jDay, 'l', 'f', xp, yp, lod, ddpsi, ddeps, jdut1, jdut1Frac, ttt);
    
    double dummy[3];
    
    coordFK5::teme_eci(rteme, dummy, dummy, edirection::eTo, reci, dummy, dummy, astro::iau80::get(), ttt, ddpsi, ddeps);
    
  }
  
} /* namespace astro */

#endif /* _H_ASTRO_TEME2ECEF_H */

