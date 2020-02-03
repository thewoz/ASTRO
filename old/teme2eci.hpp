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

#if(0)

#ifndef _H_ASTRO_TEME2ECI_H
#define _H_ASTRO_TEME2ECI_H

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
              
    std::vector< std::vector<double> > prec, nut, temp, tempmat, nutp, precp, eqep;
    
    double psia, wa, epsa, chia, deltapsi, deltaeps, trueeps, meaneps, omega, eqeg;
    
    std::vector< std::vector<double> > eqe = std::vector< std::vector<double> >(3, std::vector<double>(3));

    //astro::iau80::init();
    
    // mi serve: prec
    coordFK5::precess(ttt, e80, psia, wa, epsa, chia, prec);
    
    // mi serve: deltapsi, meaneps, nut
    coordFK5::nutation(ttt, ddpsi, ddeps, astro::iau80::get(), deltapsi, deltaeps, trueeps, meaneps, omega, nut);

    // ------------------------find eqeg----------------------
    // rotate teme through just geometric terms
    eqeg = deltapsi * cos(meaneps);
    eqeg = fmod(eqeg, 2.0*M_PI);
    
    //printf("%f %f %f %f %f %f %f %f\n", psia, wa, epsa, chia, deltapsi, deltaeps, trueeps, meaneps, omega, eqeg);
    //printf("psia %f wa %f epsa %f chia %f deltapsi %f trueeps %f meaneps %f omega %f eqeg %f\n\n", psia, wa, epsa, chia, deltapsi, trueeps, meaneps, omega, eqeg);

    eqe[0][0] = cos(eqeg);
    eqe[0][1] = sin(eqeg);
    eqe[0][2] = 0.0;
    eqe[1][0] = -sin(eqeg);
    eqe[1][1] = cos(eqeg);
    eqe[1][2] = 0.0;
    eqe[2][0] = 0.0;
    eqe[2][1] = 0.0;
    eqe[2][2] = 1.0;
    
    // tempmat = prec * nut * eqe';
    astMath::mattrans(eqe, eqep, 3, 3);
    astMath::matmult(nut,  eqep, temp, 3, 3, 3);
    astMath::matmult(prec, temp, tempmat, 3, 3, 3);
    
    astMath::matvecmult(tempmat, rteme, reci);
    astMath::matvecmult(tempmat, vteme, veci);
    astMath::matvecmult(tempmat, ateme, aeci);
                
  }

  /********************************************************************************
   * teme2eci
   ********************************************************************************/
  template <typename T>
  void teme2eci(std::vector<T> & rteme, std::vector<T> & vteme, std::vector<T> & ateme,
                std::vector<T> & reci,  std::vector<T> & veci,  std::vector<T> & aeci,
                std::vector<double> & ttt, std::vector<double> & ddpsi, std::vector<double> & ddeps) {
   
    //TODO: controlla che tutte le size sono uguali
   
    std::size_t size = rteme.size();
   
    for(std::size_t i=0; i<size; ++i)
      teme2eci(&rteme[i][0], &vteme[i][0], &ateme[i][0], &reci[i][0], &veci[i][0], &aeci[i][0], ttt[i], ddpsi[i], ddeps[i]);
   
  }
   
} /* namespace astro */

#endif /* _H_ASTRO_TEME2ECI_H */

#endif
