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
    
    astro::Date releaseDate;
    
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
      
      releaseDate = astro::Date(satrec.jdsatepoch, satrec.jdsatepochF);
      
      //printf("TLE release time:       %s\n",  TLErelease.toString());
      
    }
    
  }; /* class TLE */
  
} /* namespace astro */

#endif /* _H_ASTRO_TLE_H_ */
