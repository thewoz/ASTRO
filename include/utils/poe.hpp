/*
 * GNU GENERAL PUBLIC LICENSE
 *
 * Copyright (C) 2021
 * Created by Stefania Melillo (stefania.melillo79[at]gmail.com)
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

#ifndef _H_ASTRO_POE_H
#define _H_ASTRO_POE_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <vector>

//****************************************************************************/
// namespace astro
//****************************************************************************/
namespace astro {
  
  //****************************************************************************/
  // class poe
  //****************************************************************************/
  class poe {
    
  private:
    
    //****************************************************************************/
    // struct poe_t
    //****************************************************************************/
    struct poe_t {
      
      // data del POE in UNIX time
      double time;
      
      // posizione in 3D del satellite
      double position[3];
      
      bool operator <  (const poe_t & arg) const { return time <  arg.time; }
      bool operator == (const poe_t & arg) const { return time == arg.time; }
      
    };
    
    // vettore di appoggio
    std::vector<poe_t> data;
    
  public:
    
    //****************************************************************************/
    // costructor
    //****************************************************************************/
    poe(const char * path) { load(path); }
    
    //****************************************************************************/
    // load
    //****************************************************************************/
    // carica data, posizione e velocita' del satellite EOF
    // formato del file anno mese giorno ora minuti secondi x y z vx vy vz
    //****************************************************************************/
    void load(const char * path){
            
      char str[PATH_MAX];
      
      if(path[0] == '~'){
        snprintf(str, PATH_MAX, "%s%s", getenv("HOME"), &path[1]);
      } else {
        strcpy(str, path);
      }
      
      FILE * inFile = fopen(str, "r");
      
      if(inFile == NULL) {
        fprintf(stderr, "Oh dear, something went wrong with open file \"%s\": %s\n", str, strerror(errno));
        exit(1);
      }
      
      data.clear();
      
      int Y, M, D, h, m, s;
      double dummy;
      
      char line[PATH_MAX];
            
      poe_t tmp;
      
      time_t rawtime;
      std::time(&rawtime);
      
      struct tm * timeinfo;

      while(fgets(line, PATH_MAX, inFile)){
        
        if(strlen(line) <= 1) continue;
        
        if(line[0] == '#') continue;
        
        sscanf(line, "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf",
               &Y, &M, &D, &h, &m, &s, &tmp.position[0], &tmp.position[1], &tmp.position[2], &dummy, &dummy, &dummy);
        
        timeinfo = localtime(&rawtime);
        
        timeinfo->tm_year = Y - 1900;
        timeinfo->tm_mon  = M - 1;
        timeinfo->tm_mday = D;
        timeinfo->tm_hour = h;
        timeinfo->tm_min  = m;
        timeinfo->tm_sec  = s;
                
        tmp.time = std::mktime(timeinfo);
        
        data.push_back(tmp);
        
      }
      
      fclose(inFile);
      
    }
    
    //****************************************************************************/
    // get
    //****************************************************************************/
    void get(double time, double * position) {
      
      for(int i=0; i<data.size(); ++i) {
        
        if(time == data[i].time) {
          position[0] = data[i].position[0];
          position[1] = data[i].position[1];
          position[2] = data[i].position[2];
          break;
        }
        
        if(time < data[i].time) {
          
          // distanza in secondi tra i due punti che uso per interpolare
          double deltaT = data[i].time - data[i-1].time;
          
          // distanza in secondi tra il primo punto e il punto che voglio
          double deltaInterpol = time - data[i-1].time;
          
//          printf("%f %f %f %f %f %f %f %f\n", deltaT, deltaInterpol,
//                 data[i-1].position[0], data[i-1].position[1], data[i-1].position[2],
//                 data[i].position[0],   data[i].position[1],   data[i].position[2]);
          
          position[0] = data[i-1].position[0] + (data[i].position[0] - data[i-1].position[0])*deltaInterpol/(double)deltaT;
          position[1] = data[i-1].position[1] + (data[i].position[1] - data[i-1].position[1])*deltaInterpol/(double)deltaT;
          position[2] = data[i-1].position[2] + (data[i].position[2] - data[i-1].position[2])*deltaInterpol/(double)deltaT;
                 
          break;

        }
        
      }
              
    }
    
  }; /* class Poe */
  
} /* namespace astro */

#endif /* _H_ASTRO_POE_H */
