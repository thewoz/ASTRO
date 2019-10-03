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

#ifndef _H_ASTRO_ECEF2TEME_H
#define _H_ASTRO_ECEF2TEME_H

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
   * ecef2teme
   ********************************************************************************
   *
   *  This function trsnforms a vector from earth fixed (itrf) frame
   *  to the true equator mean equniox frame (teme).
   *  The results take into account the effects of sidereal time, and polar motion.
   *
   *  inputs:         description                     range / units
   *    recef       - position vector earth fixed     km
   *    vecef       - velocity vector earth fixed     km/s
   *    aecef       - acceleration vector earth fixed km/s2
   *    ttt         - julian centuries of tt          centuries
   *    jdut1       - julian date of ut1              days from 4713 bc
   *    lod         - excess length of day            sec
   *    xp          - polar motion coefficient        rad
   *    yp          - polar motion coefficient        rad
   *    eqeterms    - terms for ast calculation       0,2 use extra two terms (kinematic) after 1997  0, 2
   *
   *  outputs:
   *    rteme       - position vector teme            km
   *    vteme       - velocity vector teme            km/s
   *    ateme       - acceleration vector teme        km/s2
   *
   **********************************************************************************/
  void ecef2teme(double recef[3], double vecef[3], double aecef[3],
                 double rteme[3], double vteme[3], double ateme[3],
                 double ttt,      double jdut1,    double lod,
                 double xp,       double yp,       int eqeterms = 2) {
    
    coordFK5::teme_ecef(rteme, vteme, ateme, edirection::eFrom, recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms);
    
  }
  
  /********************************************************************************
   * ecef2teme
   ********************************************************************************/
  void ecef2teme(std::vector<double> & recef, std::vector<double> & vecef, std::vector<double> & aecef,
                 std::vector<double> & rteme, std::vector<double> & vteme, std::vector<double> & ateme,
                 const std::vector<double> & ttt, const std::vector<double> & jdut1, const std::vector<double> & lod,
                 const std::vector<double> & xp,  const std::vector<double> & yp, int eqeterms = 2) {
    
    //TODO: controlla che tutte le size sono uguali
    
    std::size_t size = recef.size();
    
    rteme.resize(size);
    vteme.resize(size);
    ateme.resize(size);
    
    for(std::size_t i=0; i<size; ++i)
      coordFK5::teme_ecef(&rteme[i], &vteme[i], &ateme[i], edirection::eFrom, &recef[i], &vecef[i], &aecef[i], ttt[i], jdut1[i], lod[i], xp[i], yp[i], eqeterms);
    
  }
  
  /********************************************************************************
   * ecef2teme
   ********************************************************************************/
  void ecef2teme(double recef[3], double vecef[3], double rteme[3], double vteme[3],
                 double ttt,      double jdut1,    double lod,
                 double xp,       double yp,       int eqeterms = 2) {
    
    double dummy[3];
    
    coordFK5::teme_ecef(rteme, vteme, dummy, edirection::eFrom, recef, vecef, dummy, ttt, jdut1, lod, xp, yp, eqeterms);
    
  }
  
  /********************************************************************************
   * ecef2teme
   ********************************************************************************/
  void ecef2teme(std::vector<double> & recef, std::vector<double> & vecef, std::vector<double> & rteme, std::vector<double> & vteme,
                 const std::vector<double> & ttt, const std::vector<double> & jdut1, const std::vector<double> & lod,
                 const std::vector<double> & xp,  const std::vector<double> & yp, int eqeterms = 2) {
    
    //TODO: controlla che tutte le size sono uguali
    
    double dummy[3];
    
    std::size_t size = recef.size();

    rteme.resize(size);
    vteme.resize(size);
    
    for(std::size_t i=0; i<size; ++i)
      coordFK5::teme_ecef(&rteme[i], &vteme[i], dummy, edirection::eFrom, &recef[i], &vecef[i], dummy, ttt[i], jdut1[i], lod[i], xp[i], yp[i], eqeterms);
    
  }
  
  /********************************************************************************
   * ecef2eci
   ********************************************************************************
   *
   *  This function transforms a vector from the ECI mean Equator Mean Equinox (j2000)
   *  to the ECEF Earth Fixed (itrf) frame.
   *
   *  inputs:         description                     range / units
   *    recef       - position vector earth fixed     km
   *    vecef       - velocity vector earth fixed     km/s
   *    aecef       - acceleration vector earth fixed km/s2
   *    ttt         - julian centuries of tt          centuries
   *    jdut1       - julian date of ut1              days from 4713 bc
   *    lod         - excess length of day            sec
   *    xp          - polar motion coefficient        rad
   *    yp          - polar motion coefficient        rad
   *    eqeterms    - terms for ast calculation       0,2
   *
   *  outputs:
   *    reci        - position vector eci             km
   *    veci        - velocity vector eci             km/s
   *    aeci        - acceleration vector eci         km/s2
   *
   **********************************************************************************/
  void ecef2eci(double recef[3], double vecef[3], double aecef[3], double reci[3], double veci[3], double aeci[3],
                double ttt,      double jdut1,    double lod,      double xp,      double yp,      int eqeterms = 2){
    
    std::vector< std::vector<double> > trans;
    
    // NOTE: ultimo argomento non serve
    coordFK5::itrf_j2k(recef, vecef, aecef, edirection::eTo, reci, veci, aeci, astro::iau80::get(), ttt, jdut1, lod, xp, yp, eqeterms, trans);
    
  }
  
  /********************************************************************************
   * ecef2eci
   ********************************************************************************/
  void ecef2eci(double recef[3], double vecef[3], double reci[3], double veci[3],
                double ttt,      double jdut1,    double lod,     double xp,      double yp, int eqeterms = 2){
    
    std::vector< std::vector<double> > trans;
    
    double dummy[3];
    
    // NOTE: ultimo argomento non serve
    coordFK5::itrf_j2k(recef, vecef, dummy, edirection::eTo, reci, veci, dummy, astro::iau80::get(), ttt, jdut1, lod, xp, yp, eqeterms, trans);
    
  }
  
  /********************************************************************************
   * ecef2eci
   ********************************************************************************/
  void ecef2eci(std::vector<double> & recef, std::vector<double> & vecef, std::vector<double> & aecef,
                std::vector<double> & reci,  std::vector<double> & veci,  std::vector<double> & aeci,
                std::vector<double> & ttt,   std::vector<double> & jdut1, std::vector<double> & lod,
                std::vector<double> & xp,    std::vector<double> & yp,    int eqeterms = 2) {
    
    //TODO: controlla che tutte le size sono uguali
    
    std::vector< std::vector<double> > trans;
    
    std::size_t size = recef.size();
    
    reci.resize(size);
    veci.resize(size);
    aeci.resize(size);
    
    for(std::size_t i=0; i<size; ++i)
      coordFK5::itrf_j2k(&recef[i], &vecef[i], &aecef[i], edirection::eTo, &reci[i], &veci[i], &aeci[i], astro::iau80::get(), ttt[i], jdut1[i], lod[i], xp[i], yp[i], eqeterms, trans);
    
  }
  
  /********************************************************************************
   * ecef2eci
   ********************************************************************************/
  void ecef2eci(std::vector<double> & recef, std::vector<double> & vecef, std::vector<double> & reci, std::vector<double> & veci,
                std::vector<double> & ttt,   std::vector<double> & jdut1, std::vector<double> & lod,
                std::vector<double> & xp,    std::vector<double> & yp,    int eqeterms = 2) {
    
    //TODO: controlla che tutte le size sono uguali
    
    std::vector< std::vector<double> > trans;
    
    double dummy[3];
    
    std::size_t size = recef.size();

    reci.resize(size);
    veci.resize(size);
    
    for(std::size_t i=0; i<size; ++i)
      coordFK5::itrf_j2k(&recef[i], &vecef[i], dummy, edirection::eTo, &reci[i], &veci[i], dummy, astro::iau80::get(), ttt[i], jdut1[i], lod[i], xp[i], yp[i], eqeterms, trans);
    
  }
  
} /* namespace astro */

#endif /* _H_ASTRO_ECEF2TEME_H */


