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

#ifndef _H_ASTRO_DATE_H
#define _H_ASTRO_DATE_H

#include <cstdio>
#include <cstdlib>

#include <cmath>

#include <string>
#include <vector>

#include <vallado/astTime.h>

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
  
  /*****************************************************************************/
  // class Date
  /*****************************************************************************/
  class Date {
    
  private:
    
    static char unitString[6][12];

    mutable char printString[1024];
    
  public:
    
    enum { YEARS, MONTHS, DAYS, HOURS,  MINUTES, SECONDS };
    enum { FROM_SECOND_TO_JD, FROM_SECOND_TO_MINUTE};
    
    double jDay;
    double jDayFrac;
    
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    
    /*****************************************************************************/
    // constructor
    /*****************************************************************************/
    Date() { }
    Date(int _day, int _month, int _year, int _hour, int _minute, double _second){ init(_day, _month, _year, _hour, _minute, _second); }
    Date(double _jDay, double _jDayFrac) { init(_jDay, _jDayFrac); }
    Date(double _jDay) { init(_jDay); }
    Date(const std::string & strDate, const std::string & strTime) { init(strDate, strTime); }

    /*****************************************************************************/
    // init
    /*****************************************************************************/
    void init(int _day, int _month, int _year, int _hour, int _minute, double _second) {
      
      //TODO: non controllo che l'input e' corretto
      year = _year; month = _month; day = _day; hour = _hour; minute = _minute; second = _second;
    
      astTime::jday(year, month, day, hour, minute, second, jDay, jDayFrac);
      
    }
    
    /*****************************************************************************/
    // init
    /*****************************************************************************/
    void init(const std::string & strDate, const std::string & strTime) {

      //TODO: non controllo che l'input e' corretto
      std::string delims = "/";
      std::vector<std::string> tokens;
      size_t beg, pos = 0;
      
      while((beg = strDate.find_first_not_of(delims, pos)) != std::string::npos) {
        pos = strDate.find_first_of(delims, beg + 1);
        tokens.push_back(strDate.substr(beg, pos - beg));
      }

      day = atoi(tokens[0].c_str());
      month = atoi(tokens[1].c_str());
      year = atoi(tokens[2].c_str());
      
      delims = ":";
      tokens.clear();
      pos = 0;
      
      while((beg = strTime.find_first_not_of(delims, pos)) != std::string::npos) {
        pos = strTime.find_first_of(delims, beg + 1);
        tokens.push_back(strTime.substr(beg, pos - beg));
      }
      
      hour = atoi(tokens[0].c_str());
      minute = atoi(tokens[1].c_str());
      second = atof(tokens[2].c_str());
      
      astTime::jday(year, month, day, hour, minute, second, jDay, jDayFrac);
            
    }
    
    void init(double _jDay, double _jDayFrac = 0.0) {
      
      jDay = _jDay; jDayFrac = _jDayFrac;
      
      astTime::invjday(jDay, jDayFrac, year, month, day, hour, minute, second);
      
    }
    
    
    /*****************************************************************************/
    // getJDay
    /*****************************************************************************/
    double getJDay() const { return jDay + jDayFrac; }
    
    
    /*****************************************************************************/
    // convert
    /*****************************************************************************/
    double convert(int type) const {
      
      if(type != MINUTES && type != SECONDS) {
        fprintf(stderr, "The conversion in %s is not yet implemented\n", unitString[type]);
        abort();
      }
      
      if(type == MINUTES) {
        
        printf("%f %f\n", jDay, jDayFrac);
        // converto in secondi e divido per i secondi
        return (jDay+jDayFrac)*(86400/60.0);//((jDay-0.5) * 1440) + ((hour+12) * 60) + minute + (second / 60.0);
        
      }
      
      if(type == SECONDS) {
      
        return (jDay+jDayFrac)*86400;
      }
      
      fprintf(stderr, "The conversion type not reconized\n");
      
      abort();
      
    }
    
    /*****************************************************************************/
    // convert
    /*****************************************************************************/
    static double convert(float value, int mode) {
      
      if(mode == FROM_SECOND_TO_JD) return (value/86400.0);
      if(mode == FROM_SECOND_TO_MINUTE) return (value / 60.0);
      abort();
      
    }
    

    /*****************************************************************************/
    // difference
    /*****************************************************************************/
    static double difference(const Date & dataA, const Date & dataB, int unit) {
      
      return dataA.convert(unit) - dataB.convert(unit);
      
    }

    /*****************************************************************************/
    // operator ()
    /*****************************************************************************/
    template <typename T>
    Date operator () (char mode, T value, int unit) {
      
      if(unit != MINUTES) {
        fprintf(stderr, "The unit '%s' is not reconized\n", unitString[unit]);
        abort();
      }
      
      if(unit == MINUTES && (mode != '+' && mode != '-')) {
        fprintf(stderr, "The operation '%c' is not reconized\n", mode);
        abort();
      }
      
      if(unit == MINUTES) {
        
        double time = convert(MINUTES);
        
        if(mode == '+') time += value;
        if(mode == '-') time -= value;

        // gli tolgo 720 minuti che corrispondo a 0.5 jd
        double _jDay = trunc((time-720) / 1440.0) + 0.5;
        
        // mi prendo i minuti restanti
        double _left = time - (_jDay * 1440);
        
        // converto i minuti mancanti in ore minuti e secondi
        double _hours   = trunc(_left / 60.0);
        double _minutes = trunc(_left - (_hours * 60));
        double _second  = (_left - (_hours * 60) - _minutes) * 60.0;
        
        double _jDayFrac = ((_hours) / 24.0) + (_minutes / 1440.0) + (_second / 86400.0);
        
        return Date(_jDay, _jDayFrac);
        
      }

      fprintf(stderr, "The operation is not valid\n");
      
      abort();
        
    }
    
    /*****************************************************************************/
    // print functions
    /*****************************************************************************/
    void printGregorianDate(FILE * output = stdout) const {
      fprintf(output, "%02d/%02d/%d", day, month, year);
    }
    void printGregorian(FILE * output = stdout) const {
      fprintf(output, "%02d/%02d/%d at %02d:%02d:%05.2f", day, month, year, hour, minute, second);
    }
    void printJulian(FILE * output = stdout) const {
      fprintf(output, "%d %05.3f", (int)trunc(jDay), jDayFrac);
    }
    void print(FILE * output = stdout) const {
      fprintf(output, "%02d/%02d/%d at %02d:%02d:%05.2f JD %.5f", day, month, year, hour, minute, second, jDay+jDayFrac);
    }
    void printGregorianDateln(FILE * output = stdout) const {
      fprintf(output, "%02d/%02d/%d\n", day, month, year);
    }
    void printGregorianln(FILE * output = stdout) const {
      fprintf(output, "%02d/%02d/%d at %02d:%02d:%05.2f\n", day, month, year, hour, minute, second);
    }
    void printJulianln(FILE * output = stdout) const {
      fprintf(output, "%d %05.3f\n", (int)trunc(jDay), jDayFrac);
    }
    void println(FILE * output = stdout) const {
      fprintf(output, "%02d/%02d/%d at %02d:%02d:%05.2f JD %.5f [%f %f]\n", day, month, year, hour, minute, second, jDay+jDayFrac, jDay, jDayFrac);
    }
    
    
    
    /*****************************************************************************/
    // toString functions
    /*****************************************************************************/
    const char * toGregorianDateString(FILE * output = stdout) const {
      sprintf(printString, "%02d/%02d/%d", day, month, year);
      return printString;
    }
    const char * toGregorianString() const {
      sprintf(printString, "%02d/%02d/%d at %02d:%02d:%05.2f", day, month, year, hour, minute, second);
      return printString;
    }
    const char * toJulianString() const {
      sprintf(printString, "%d %05.3f", (int)trunc(jDay), jDayFrac);
      return printString;
    }
    const char * toString() const {
      sprintf(printString, "%02d/%02d/%d at %02d:%02d:%05.2f JD %.5f  [%f %f]", day, month, year, hour, minute, second, jDay+jDayFrac, jDay, jDayFrac);
      return printString;
    }
    
    /*****************************************************************************/
    // invjday
    /*****************************************************************************/
    static void invjday(double _jDay, double _jDayFrac, int & _year, int & _month, int & _day, int & _hour, int & _minute, double &  _second) {
      
      astTime::invjday(_jDay, _jDayFrac, _year, _month, _day, _hour, _minute, _second);

    }
    
    
  }; /* class Date */
  
  char Date::unitString[6][12] = { "years", "months", "days", "hours", "minutes", "seconds" };
  
  
  
  
} /* namespace astro */

#endif /* _H_ASTRO_DATE_H */
