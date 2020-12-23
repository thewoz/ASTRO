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

#ifndef _H_ASTRO_SGP4_H_
#define _H_ASTRO_SGP4_H_

#include <cstdlib>
#include <cstdio>

#include <cmath>

//**********************************************************************************/
// namespace astro
//**********************************************************************************/
namespace astro {
  
  //**********************************************************************************/
  // namespace util
  //**********************************************************************************/
  namespace util {
    
    //**********************************************************************************/
    // getElesetrecErrorString
    //**********************************************************************************/
    char elesetrecErrorString[7][1024] = { "no error",
      "mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er",
      "mean motion less than 0.0",
      "pert elements, ecc < 0.0  or  ecc > 1.0",
      "semi-latus rectum < 0.0",
      "epoch elements are sub-orbital",
      "satellite has decayed" };
    
    const char * getElesetrecErrorString(int error){ return elesetrecErrorString[error]; }
    
  }  /* namespace util */
  
  //**********************************************************************************/
  // sgp4
  //**********************************************************************************/
  void sgp4(elsetrec & satrec, double tsince, double r[3], double v[3]) {
    
    SGP4Funcs::sgp4(satrec, tsince, r, v);
    
    if(satrec.error > 0) {
      fprintf(stderr, "SGP4 error at time %f: %s\n",  satrec.t, util::getElesetrecErrorString(satrec.error));
      abort();
    }
    
  }
  
  //**********************************************************************************/
  // sgp4
  //**********************************************************************************/
  void sgp4(elsetrec & satrec, double tsince, double r[3]) {
    
    SGP4Funcs::sgp4(satrec, tsince, r);
    
    if(satrec.error > 0) {
      fprintf(stderr, "SGP4 error at time %f: %s\n",  satrec.t, util::getElesetrecErrorString(satrec.error));
      abort();
    }
    
  }
  
} /* namespace astro */

#endif /* _H_ASTRO_SGP4_H_ */
