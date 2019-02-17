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
#include <cmath>

#include <vallado/astMath.h>
#include <vallado/coordFK5.h>

#include "eopc.hpp"

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
  void teme2ecef(const double rteme[3], const double vteme[3], const double ateme[3],
                 double recef[3], double vecef[3], double aecef[3],
                 double ttt,      double jdut1,    double lod,
                 double xp,       double yp,       int eqeterms = 2) {
    
    std::vector< std::vector<double> > st      = std::vector< std::vector<double> >(3, std::vector<double>(3));
    std::vector< std::vector<double> > stdot   = std::vector< std::vector<double> >(3, std::vector<double>(3));
    std::vector< std::vector<double> > temp    = std::vector< std::vector<double> >(3, std::vector<double>(3));
    std::vector< std::vector<double> > tempmat = std::vector< std::vector<double> >(3, std::vector<double>(3));
    
    std::vector< std::vector<double> > pm, pmp, stp;
    
    double thetasa, omegaearth[3], rpef[3], vpef[3], apef[3], omgxr[3], omgxomgxr[3], omgxv[3], tempvec1[3], tempvec[3], gmst, deg2rad, omega, gmstg;
    
    deg2rad = M_PI / 180.0;
    
    // find omeage from nutation theory
    omega = 125.04452222 + (-6962890.5390 *ttt + 7.455 *ttt*ttt + 0.008 *ttt*ttt*ttt) / 3600.0;
    omega = fmod(omega, 360.0) * deg2rad;
    
    // ------------------------find gmst--------------------------
    gmst = astTime::gstime(jdut1);
    
    // teme does not include the geometric terms here after 1997, kinematic terms apply
    if((jdut1 > 2450449.5) && (eqeterms > 0)) {
      
      gmstg = gmst + 0.00264*M_PI / (3600 * 180)*sin(omega) + 0.000063*M_PI / (3600 * 180)*sin(2.0 *omega);
      
    } else gmstg = gmst;
    
    gmstg = fmod(gmstg, 2.0*M_PI);
    
    thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
    
    omegaearth[0] = 0.0;
    omegaearth[1] = 0.0;
    omegaearth[2] = thetasa;
    
    st[0][0] = cos(gmstg);
    st[0][1] = -sin(gmstg);
    st[0][2] = 0.0;
    st[1][0] = sin(gmstg);
    st[1][1] = cos(gmstg);
    st[1][2] = 0.0;
    st[2][0] = 0.0;
    st[2][1] = 0.0;
    st[2][2] = 1.0;
    
    // compute sidereal time rate matrix
    stdot[0][0] = -omegaearth[2] * sin(gmstg);
    stdot[0][1] = -omegaearth[2] * cos(gmstg);
    stdot[0][2] = 0.0;
    stdot[1][0] = omegaearth[2] * cos(gmstg);
    stdot[1][1] = -omegaearth[2] * sin(gmstg);
    stdot[1][2] = 0.0;
    stdot[2][0] = 0.0;
    stdot[2][1] = 0.0;
    stdot[2][2] = 0.0;
    
    coordFK5::polarm(xp, yp, ttt, e80, pm);
    
    astMath::mattrans(pm, pmp, 3, 3);
    astMath::mattrans(st, stp, 3, 3);
    
    astMath::matvecmult(stp, rteme, rpef);
    astMath::matvecmult(pmp, rpef, recef);
    
    astMath::cross(omegaearth, rpef, omgxr);
    astMath::matvecmult(stp, vteme, tempvec1);
    astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
    astMath::matvecmult(pmp, vpef, vecef);
    
    astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
    astMath::cross(omegaearth, vpef, omgxv);
    astMath::cross(omegaearth, omgxr, omgxomgxr);
    astMath::matvecmult(stp, ateme, tempvec1);
    astMath::addvec(1.0, tempvec1, -1.0, omgxomgxr, tempvec);
    astMath::addvec(1.0, tempvec, -2.0, omgxv, apef);
    astMath::matvecmult(pmp, apef, aecef);
    
  }
  
  
} /* namespace astro */

#endif /* _H_ASTRO_TEME2ECEF_H */

