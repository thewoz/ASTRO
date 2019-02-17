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

#ifndef _H_ASTRO_ECI2ECEF_H
#define _H_ASTRO_ECI2ECEF_H

#include <cstdio>
#include <cstdlib>

#include <cmath>

#include <vector>

#include <vallado/coordFK5.h>

#include "iau80.hpp"


/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
  
  /********************************************************************************
   * eci2ecef
   ********************************************************************************
   *
   *  This function transforms a vector from the eci mean equator mean equinox (j2000) to the earth fixed (itrf) frame.
   *
   *  inputs:         description                     range / units
   *    reci        - position vector eci             km
   *    veci        - velocity vector eci             km/s
   *    aeci        - acceleration vector eci         km/s2
   *    ttt         - julian centuries of tt          centuries
   *    jdut1       - julian date of ut1              days from 4713 bc
   *    lod         - excess length of day            sec
   *    xp          - polar motion coefficient        rad
   *    yp          - polar motion coefficient        rad
   *    eqeterms    - terms for ast calculation       0,2
   *
   *  outputs:
   *    recef       - position vector earth fixed     km
   *    vecef       - velocity vector earth fixed     km/s
   *    aecef       - acceleration vector earth fixed km/s2

   *
   **********************************************************************************/
  void eci2ecef(double recef[3], double vecef[3], double aecef[3], double reci[3], double veci[3], double aeci[3],
                double ttt,      double jdut1,    double lod,      double xp,      double yp,  int eqeterms = 2){
    
    astro::iau80::init();
    
    std::vector< std::vector<double> > trans;
    
    trans.resize(3);  // rows
    
    for(std::vector< std::vector<double> >::iterator it = trans.begin(); it != trans.end(); ++it)
      it->resize(3);
    
    std::vector< std::vector<double> > prec, nut, st, stdot, pm, temp,
    pmp, stp, nutp, precp;
    
    double psia, wa, epsa, chia, deltapsi, deltaeps, trueeps, meaneps,
    omega, thetasa, omegaearth[3], rpef[3], vpef[3], apef[3], omgxr[3], omgxomgxr[3],
    omgxv[3], tempvec1[3], tempvec[3];
    
    // ---- find matrices
    coordFK5::precess(ttt, e80, psia, wa, epsa, chia, prec);
    coordFK5::nutation(ttt, 0.0, 0.0, astro::iau80::get(), deltapsi, deltaeps, trueeps, meaneps, omega, nut);
    coordFK5::sidereal(jdut1, deltapsi, meaneps, omega, lod, 2, st, stdot);
    coordFK5::polarm(xp, yp, ttt, e80, pm);
    
    // ---- perform transformations
    thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
    
    omegaearth[0] = 0.0;
    omegaearth[1] = 0.0;
    omegaearth[2] = thetasa;
    
    astMath::mattrans(pm, pmp, 3, 3);
    astMath::mattrans(st, stp, 3, 3);
    astMath::mattrans(nut, nutp, 3, 3);
    astMath::mattrans(prec, precp, 3, 3);
    
    astMath::matmult(stp, nutp, temp, 3, 3, 3);
    astMath::matmult(temp, precp, trans, 3, 3, 3);
    astMath::matvecmult(trans, reci, rpef);
    astMath::matvecmult(pmp, rpef, recef);
    
    astMath::cross(omegaearth, rpef, omgxr);
    astMath::matvecmult(trans, veci, tempvec1);
    astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
    astMath::matvecmult(pmp, vpef, vecef);
    
    astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
    astMath::cross(omegaearth, vpef, omgxv);
    astMath::cross(omegaearth, omgxr, omgxomgxr);
    astMath::matvecmult(trans, aeci, tempvec1);
    astMath::addvec(1.0, tempvec1, -1.0, omgxomgxr, tempvec);
    astMath::addvec(1.0, tempvec, -2.0, omgxv, apef);
    astMath::matvecmult(pmp, apef, aecef);
    
  }
  
  
  /********************************************************************************
   * eci2ecef
   ********************************************************************************/
  template <typename T>
  void eci2ecef(std::vector<T> & recef,      std::vector<T> & vecef,      std::vector<T> & aecef,
                std::vector<T> & reci,       std::vector<T> & veci,       std::vector<T> & aeci,
                std::vector<double> & ttt,   std::vector<double> & jdut1, std::vector<double> & lod,
                std::vector<double> & xp,    std::vector<double> & yp, int eqeterms = 2) {
    
    //TODO: controlla che tutte le size sono uguali
    
    std::size_t size = recef.size();
    
    for(std::size_t i=0; i<size; ++i)
      eci2ecef(&reci[i][0], &veci[i][0], &aeci[i][0], &recef[i][0], &vecef[i][0], &aecef[i][0], ttt[i], jdut1[i], lod[i], xp[i], yp[i], eqeterms);
    
  }
  
} /* namespace astro */

#endif /* _H_ASTRO_ECI2ECEF_H */


