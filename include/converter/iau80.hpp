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

#ifndef _H_ASTRO_IAU80IN_H
#define _H_ASTRO_IAU80IN_H

#include <cstdio>
#include <cstdlib>

#include <cmath>

//****************************************************************************/
// namespace astro
//****************************************************************************/
namespace astro {
  
  //****************************************************************************/
  // class iau80
  //****************************************************************************/
  class iau80 {
    
  private:
    
    iau80();
    
    static bool isInited;
    
    static iau80data iau80rec;

    //****************************************************************************/
    // init
    //****************************************************************************/
    static void init() {
      
      if(!isInited) {

        coordFK5::iau80in(iau80rec, "/usr/local/include/vallado/data/nut80.dat");
                
        isInited = true;
      
      }
      
    }
    
  public:

    //****************************************************************************/
    // get
    //****************************************************************************/
    static iau80data & get() { init();  return iau80rec; }
    
    //****************************************************************************/
    // print
    //****************************************************************************/
    static void print(FILE * output = stdout) {
      
      //double convrt = 3600.0 / 0.0001; // 0.0001" to deg
      double convrt = 0.0001 * M_PI / (180 * 3600.0);   // 0.0001" to rad

      for(int i=1; i<=106; ++i)
        fprintf(output, "%d %d %d %d %d %lf %lf %lf %lf\n",
                iau80rec.iar80[i][1], iau80rec.iar80[i][2], iau80rec.iar80[i][3],
                iau80rec.iar80[i][4], iau80rec.iar80[i][5],
                iau80rec.rar80[i][1]/ convrt, iau80rec.rar80[i][2]/ convrt, iau80rec.rar80[i][3]/ convrt,
                iau80rec.rar80[i][4]/ convrt);
    }
  
  }; /* class iau80 */
  
  bool iau80::isInited = false;
  
  iau80data iau80::iau80rec = iau80data();
  
} /* namespace astro */

#endif /* _H_ASTRO_IAU80IN_H */
