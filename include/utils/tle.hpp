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

#ifndef _H_ASTRO_TLE_H_
#define _H_ASTRO_TLE_H_

#include <cstdlib>
#include <cstdio>

#include <string>


//**********************************************************************************/
// namespace astro
//**********************************************************************************/
namespace astro {
  
  //**********************************************************************************/
  // class TLE
  //**********************************************************************************/
  class TLE {
    
  public:
    
    std::string name;
    
    double releaseDate;

    elsetrec satrec;
    
    TLE() { }
    TLE(FILE * input) { init(input); }
    TLE(const char * TLE_line1, const char * TLE_line2) { init(TLE_line1, TLE_line2); }
    TLE(const char * TLE_line1, const char * TLE_line2, const char * TLE_line3) { init(TLE_line1, TLE_line2, TLE_line3); }
    
    void init(FILE * input) {
      
      if(input == NULL) {
        fprintf(stderr, "the file is not good\n");
        abort();
      }
      
      char * TLE_line1 = NULL;
      char * TLE_line2 = NULL;
      size_t len = 0;
      if(getline(&TLE_line1, &len, input)==-1){
        fprintf(stderr, "error in read the file\n");
        //free(TLE_line1);
        abort();
      }
      
      if(getline(&TLE_line2, &len, input)==-1){
        fprintf(stderr, "error in read the file\n");
        free(TLE_line1);
        free(TLE_line2);
        abort();
      }

      init(TLE_line1, TLE_line2);
      
      free(TLE_line1);
      free(TLE_line2);
      
    }
    
    void init(const char * TLE_line1, const char * TLE_line2) { _init("none", TLE_line1, TLE_line2); }
    
    void init(const char * TLE_line1, const char * TLE_line2, const char * TLE_line3) { _init(TLE_line1, TLE_line2, TLE_line3); }
    
    
  private:
    
    void _init(const char * TLE_line1, const char * TLE_line2, const char * TLE_line3) {
      
      //TODO:
      //name = ;
      
      // sgp4fix addiional parameters to store from the TLE
      satrec.classification = 'U';
      strncpy(satrec.intldesg, "          ", 11 * sizeof(char));
      satrec.ephtype = 0;
      satrec.elnum   = 0;
      satrec.revnum  = 0;
      
      // since we perform a complete catalog evaluation 'c', -+ 1 day
      double startmfe = 0.0;
      double stopmfe;
      double deltamin;
      
      SGP4Funcs::twoline2rv(TLE_line2, TLE_line3, 'c', '*', 'i', wgs84, startmfe, stopmfe, deltamin, satrec);
      
      releaseDate = satrec.jdsatepoch + satrec.jdsatepochF;
      
    }
    
  }; /* class TLE */
  
} /* namespace astro */

#endif /* _H_ASTRO_TLE_H_ */
