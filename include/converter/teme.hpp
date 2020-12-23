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

