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

#ifndef _H_ASTRO_EOPC_H
#define _H_ASTRO_EOPC_H

#include <cstdlib>
#include <cstdio>

#include <sys/types.h>
#include <sys/stat.h>

#include <cmath>
#include <cfloat>

#include <cstring>
#include <cerrno>

#include <vector>
#include <array>

#include <string>

//****************************************************************************/
// namespace astro
//****************************************************************************/
namespace astro {

  //****************************************************************************/
  // class eopc
  //****************************************************************************/
  class eopc {
    
  private:
    
    // variabile che mi dice se devo fare l'update dei file
    static bool toUpdate;
    
    // data referred to the precession-nutation model IAU 1980
    static const char  urlEopc[PATH_MAX];
    static const char fileEopc[PATH_MAX];

    //
    static std::vector<eopdata> data;
    
    static double jdeopstart;
    
    static bool isInited;
    
    // Hide default constructor
    eopc() {}
    
  public:
    
    //****************************************************************************/
    // getParameters
    //****************************************************************************/
    static void getParameters(const std::vector<double> & jday, char interp, char whichm, std::vector<double> & xp, std::vector<double> & yp, std::vector<double> & lod, std::vector<double> & ddpsi, std::vector<double> & ddeps, std::vector<double> & jdut1, std::vector<double> & jdut1Frac, std::vector<double> & ttt) {
     
      if(!isInited) init();

      xp.resize(jday.size());
      yp.resize(jday.size());
      lod.resize(jday.size());
      ddpsi.resize(jday.size());
      ddeps.resize(jday.size());
      jdut1.resize(jday.size());
      jdut1Frac.resize(jday.size());
      ttt.resize(jday.size());
      
      for(int i=0; i<jday.size(); ++i) getParameters(jday[i], interp, whichm, xp[i], yp[i], lod[i], ddpsi[i], ddeps[i], jdut1[i], jdut1Frac[i], ttt[i]);
      
    }
    
    //****************************************************************************/
    // getParameters
    //****************************************************************************/
    static void getParameters(double jDay, char interp, char whichm, double & xp, double & yp, double & lod, double & ddpsi, double & ddeps, double & jdut1, double & jdut1Frac, double & ttt) {
      
      if(!isInited) init();
      
      int dat;
      
      double x, y, dx, dy, s, deltapsi, deltaeps, dut1;
      
      EopSpw::findeopparam(jDay, 0, interp, whichm, data, jdeopstart, dut1, dat, lod, xp, yp, ddpsi, ddeps, dx, dy, x, y, s, deltapsi, deltaeps);

      int year; int mon; int day; int hr; int min; double sec;
      
      astro::Date::invjday(jDay, 0, year, mon, day, hr, min, sec);
      
      astTime::convtime(year, mon, day, hr, min, sec, 0, dut1, dat, jdut1, jdut1Frac, ttt);
      
    }
    
    //****************************************************************************/
    // init
    //****************************************************************************/
    static void init(int days = 1, bool verbose = false) {
      
      // aggiorno i file
      update(days, verbose);
      
      data = std::vector<eopdata>();
      
      // carico i dati (l'ordine e' importante)
      load();
      
      // mi segno che mi sono aggiornato i file
      toUpdate = false;
      
      isInited = true;
      
    }
    
  private:
        
    
    
    //****************************************************************************/
    // update() - Scarico se serve il file nuovo
    //****************************************************************************/
    static void update(int days, bool verbose) {
      
      if(verbose) fprintf(stderr, "updating eopc data file...\n");
      
      // current time in second since January 1, 1970
      time_t currentTime = time(NULL);
      
      // get creation time of the EOP data file
      struct stat st; time_t time = 0; if(stat(fileEopc, &st) == 0) { time = st.st_mtime; }
      
      // days in second
      time_t maxAge = 60 * 60 * 24 * days;
      
      if(currentTime-time > maxAge) {
        
        if(curl::get(urlEopc, fileEopc, false, true)) fprintf(stderr, "warning epoc not able to update EOP data file\n");
       
        toUpdate = true;
        
      } else { toUpdate = false; }
      
    }
    
    //****************************************************************************/
    // load() - leggo il file dei parametri
    //****************************************************************************/
    static void load() {

      double jdeopstartFrac;

      EopSpw::initeop(data, fileEopc, jdeopstart, jdeopstartFrac);

    }

  };

  //****************************************************************************/
  // static class variable definitions
  //****************************************************************************/
  std::vector<eopdata> eopc::data = std::vector<eopdata>();

  // data referred to the precession-nutation model IAU 1980
  const char eopc::urlEopc[PATH_MAX]  = "https://celestrak.com/SpaceData/eop19620101.txt";
  const char eopc::fileEopc[PATH_MAX] = "/usr/local/include/vallado/data/eopc.dat";

  // jDay del primo parametro eopc
  double eopc::jdeopstart = 0;

  // variablie che mi dice se devo fare l'update dei file
  bool eopc::toUpdate = true;
  
  bool eopc::isInited = false;
  
} /* namespace astro */

#endif /* _H_ASTRO_EOPC_H */














