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

#ifndef _H_ASTRO_SGP4_H_
#define _H_ASTRO_SGP4_H_

#include <cstdlib>
#include <cstdio>

#include <cmath>

#include <vallado/SGP4.h>

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
  
  /*****************************************************************************/
  // namespace util
  /*****************************************************************************/
  namespace util {
    
    /*****************************************************************************/
    // getElesetrecErrorString
    /*****************************************************************************/
    char elesetrecErrorString[7][1024] = { "no error",
      "mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er",
      "mean motion less than 0.0",
      "pert elements, ecc < 0.0  or  ecc > 1.0",
      "semi-latus rectum < 0.0",
      "epoch elements are sub-orbital",
      "satellite has decayed" };
    
    const char * getElesetrecErrorString(int error){ return elesetrecErrorString[error]; }
    
  }  /* namespace util */
  
  /*****************************************************************************/
  // sgp4
  /*****************************************************************************/
  void sgp4(elsetrec & satrec, double tsince, double r[3], double v[3]) {
    
    SGP4Funcs::sgp4(satrec, tsince, r, v);
    
    if(satrec.error > 0) {
      fprintf(stderr, "SGP4 error at time %f: %s\n",  satrec.t, util::getElesetrecErrorString(satrec.error));
      abort();
    }
    
  }
  
  /*****************************************************************************/
  // sgp4
  /*****************************************************************************/
  void sgp4(elsetrec & satrec, double tsince, double r[3]) {
    
    SGP4Funcs::sgp4(satrec, tsince, r);
    
    if(satrec.error > 0) {
      fprintf(stderr, "SGP4 error at time %f: %s\n",  satrec.t, util::getElesetrecErrorString(satrec.error));
      abort();
    }
    
  }
  
} /* namespace astro */

#endif /* _H_ASTRO_SGP4_H_ */
