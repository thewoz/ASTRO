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

#include <vallado/coordFK5.h>

#include "iau80.hpp"

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {

  /********************************************************************************
   * teme2ecef
   ********************************************************************************
   *
   *  This function trsnforms a vector from the true equator mean equniox frame (teme),
   *  to an earth fixed (itrf) frame.
   *  The results take into account the effects of sidereal time, and polar motion.
   *
   *  inputs:         description                     range / units
   *    rteme       - position vector teme            km
   *    vteme       - velocity vector teme            km/s
   *    ateme       - acceleration vector teme        km/s2
   *    ttt         - julian centuries of tt          centuries
   *    jdut1       - julian date of ut1              days from 4713 bc
   *    lod         - excess length of day            sec
   *    xp          - polar motion coefficient        rad
   *    yp          - polar motion coefficient        rad
   *    eqeterms    - terms for ast calculation       0,2 use extra two terms (kinematic) after 1997  0, 2
   *
   *  outputs:
   *    recef       - position vector earth fixed     km
   *    vecef       - velocity vector earth fixed     km/s
   *    aecef       - acceleration vector earth fixed km/s2
   *
   **********************************************************************************/
  void teme2ecef(double rteme[3], double vteme[3], double ateme[3],
                 double recef[3], double vecef[3], double aecef[3],
                 double ttt,      double jdut1,    double lod,
                 double xp,       double yp,       int eqeterms = 2) {
    
    coordFK5::teme_ecef(rteme, vteme, ateme, edirection::eTo, recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms);
    
  }
  
  /********************************************************************************
   * teme2ecef
   ********************************************************************************/
  void teme2ecef(std::vector<double> & rteme, std::vector<double> & vteme, std::vector<double> & ateme,
                 std::vector<double> & recef, std::vector<double> & vecef, std::vector<double> & aecef,
                 const std::vector<double> & ttt, const std::vector<double> & jdut1, const std::vector<double> & lod,
                 const std::vector<double> & xp,  const std::vector<double> & yp, int eqeterms = 2) {
    
    //TODO: controlla che tutte le size sono uguali
    
    recef.resize(rteme.size());
    vecef.resize(rteme.size());
    aecef.resize(rteme.size());
    
    std::size_t size = rteme.size();
    
    for(std::size_t i=0; i<size; ++i)
      coordFK5::teme_ecef(&rteme[i], &vteme[i], &ateme[i], edirection::eTo, &recef[i], &vecef[i], &aecef[i], ttt[i], jdut1[i], lod[i], xp[i], yp[i], eqeterms);
    
  }
  
  /********************************************************************************
   * teme2ecef
   ********************************************************************************/
  void teme2ecef(double rteme[3], double vteme[3], double recef[3], double vecef[3],
                 double ttt,      double jdut1,    double lod,
                 double xp,       double yp,       int eqeterms = 2) {
    
    double dummy[3];

    coordFK5::teme_ecef(rteme, vteme, dummy, edirection::eTo, recef, vecef, dummy, ttt, jdut1, lod, xp, yp, eqeterms);
    
  }
  
  /********************************************************************************
   * teme2ecef
   ********************************************************************************/
  void teme2ecef(std::vector<double> & rteme, std::vector<double> & vteme, std::vector<double> & recef, std::vector<double> & vecef,
                 const std::vector<double> & ttt, const std::vector<double> & jdut1, const std::vector<double> & lod,
                 const std::vector<double> & xp,  const std::vector<double> & yp, int eqeterms = 2) {
    
    //TODO: controlla che tutte le size sono uguali
    
    double dummy[3];

    recef.resize(rteme.size());
    vecef.resize(rteme.size());
    
    std::size_t size = rteme.size();
    
    for(std::size_t i=0; i<size; ++i)
      coordFK5::teme_ecef(&rteme[i], &vteme[i], dummy, edirection::eTo, &recef[i], &vecef[i], dummy, ttt[i], jdut1[i], lod[i], xp[i], yp[i], eqeterms);
    
  }
  
  /********************************************************************************
   * teme2eci
   ********************************************************************************
   *  This function transforms a vector from the true equator mean equinox system,
   *  (teme) to the mean equator mean equinox (j2000) system.
   *
   *  inputs:         description                     range / units
   *    rteme       - position vector of date         km
   *    vteme       - velocity vector of date         km / s
   *    ateme       - acceleration vector of date     km / s2
   *    ttt         - julian centuries of tt          centuries
   *    ddpsi       - delta psi correction to gcrf    rad
   *    ddeps       - delta eps correction to gcrf    rad
   *
   *  outputs :
   *    reci        - position vector eci             km
   *    veci        - velocity vector eci             km / s
   *    aeci        - acceleration vector eci         km / s2
   *
   **********************************************************************************/
  void teme2eci(double rteme[3], double vteme[3], double ateme[3],
                double reci[3],  double veci[3],  double aeci[3],
                double ttt,      double ddpsi,    double ddeps) {
    
    coordFK5::teme_eci(rteme, vteme, ateme, edirection::eTo, reci, veci, aeci, astro::iau80::get(), ttt, ddpsi, ddeps);
    
  }
  
  /********************************************************************************
   * teme2eci
   ********************************************************************************/
  void teme2eci(double rteme[3], double vteme[3], double reci[3],  double veci[3],
                double ttt,      double ddpsi,    double ddeps) {
    
    double dummy[3];
    
    coordFK5::teme_eci(rteme, vteme, dummy, edirection::eTo, reci, veci, dummy, astro::iau80::get(), ttt, ddpsi, ddeps);
    
  }
  
  /********************************************************************************
   * teme2eci
   ********************************************************************************/
  void teme2eci(std::vector<double> & rteme, std::vector<double> & vteme, std::vector<double> & ateme,
                std::vector<double> & reci,  std::vector<double> & veci,  std::vector<double> & aeci,
                const std::vector<double> & ttt, const std::vector<double> & ddpsi, const std::vector<double> & ddeps) {

    //TODO: controlla che tutte le size sono uguali
    
    reci.resize(rteme.size());
    veci.resize(rteme.size());
    aeci.resize(rteme.size());

    std::size_t size = rteme.size();
    
    for(std::size_t i=0; i<size; ++i)
      coordFK5::teme_eci(&rteme[i], &vteme[i], &ateme[i], edirection::eTo, &reci[i], &veci[i], &aeci[i], astro::iau80::get(), ttt[i], ddpsi[i], ddeps[i]);

  }
  
  /********************************************************************************
   * teme2eci
   ********************************************************************************/
  void teme2eci(std::vector<double> & rteme, std::vector<double> & vteme, std::vector<double> & reci,  std::vector<double> & veci,
                const std::vector<double> & ttt, const std::vector<double> & ddpsi, const std::vector<double> & ddeps) {
    
    
    //TODO: controlla che tutte le size sono uguali
    
    double dummy[3];

    reci.resize(rteme.size());
    veci.resize(rteme.size());
    
    std::size_t size = rteme.size();
    
    for(std::size_t i=0; i<size; ++i)
      coordFK5::teme_eci(&rteme[i], &vteme[i], dummy, edirection::eTo, &reci[i], &veci[i], dummy, astro::iau80::get(), ttt[i], ddpsi[i], ddeps[i]);
    
  }
  
} /* namespace astro */

#endif /* _H_ASTRO_TEME2ECEF_H */

