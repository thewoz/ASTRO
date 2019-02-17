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

#ifndef _H_ASTRO_IAU80IN_H
#define _H_ASTRO_IAU80IN_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vallado/coordFK5.h>

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
  
  /*****************************************************************************/
  // class iau80
  /*****************************************************************************/
  class iau80 {
    
  private:
    
    iau80();
    
    static bool isInited;
    
    //struct iau80data {
    //  int    iar80[107][6];
    //  double rar80[107][5];
    //}
    
    static iau80data iau80rec;
    
  public:

    /*****************************************************************************/
    // init
    /*****************************************************************************/
    static void init() {
      
      if(!isInited) {

        coordFK5::iau80in(iau80rec, "/usr/local/include/vallado/data/nut80.dat");
        
        isInited = true;
        
        //print();
        
      }
      
    }
    
    /*****************************************************************************/
    // get
    /*****************************************************************************/
    static iau80data & get() { return iau80rec; }
    
    /*****************************************************************************/
    // print
    /*****************************************************************************/
    static void print(FILE * output = stdout) {
      
      //double convrt = 3600.0 / 0.0001; // 0.0001" to deg
      double convrt = 0.0001;//0.0001 * pi / (180*3600.0); // 0.0001" to rad

      for(int i=0; i<107; ++i)
        fprintf(output, "%d %d %d %d %d %d %.16lf %.16lf %.16lf %.16lf %.1lf\n",
                iau80rec.iar80[i][0], iau80rec.iar80[i][1], iau80rec.iar80[i][2], iau80rec.iar80[i][3], iau80rec.iar80[i][4], iau80rec.iar80[i][5],
                iau80rec.rar80[i][0] / convrt, iau80rec.rar80[i][1] / convrt, iau80rec.rar80[i][2] / convrt, iau80rec.rar80[i][3] / convrt, iau80rec.rar80[i][4] / convrt);
      
    }
  
  }; /* class iau80 */
  
  bool iau80::isInited = false;
  
  iau80data iau80::iau80rec = iau80data();
  
} /* namespace astro */

#endif /* _H_ASTRO_IAU80IN_H */
