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

#ifndef _H_ASTRO_ATTITUDE_H_
#define _H_ASTRO_ATTITUDE_H_

#include <cstdlib>
#include <cstdio>

#include <cmath>

//**********************************************************************************/
// namespace astro
//**********************************************************************************/
namespace astro {
 
  enum { wx, wy, wz, q4, q1, q2, q3, ix, iy, iz };
  
  //**********************************************************************************/
  // struct attitude_t
  //**********************************************************************************/
  struct attitude_t {
    
    attitude_t() {}
    
    attitude_t(double Wx, double Wy, double Wz, double Q4, double Q1, double Q2, double Q3, double Jx, double Jy, double Jz) {
      (*this)(Wx, Wy, Wz, Q4, Q1, Q2, Q3, Jx, Jy, Jz);
    }
    
    double data[10];

    //**********************************************************************************/
    // get omega function
    //**********************************************************************************/
    double Wx() const { return data[wx]; }
    double Wy() const { return data[wy]; }
    double Wz() const { return data[wz]; }
    
    //**********************************************************************************/
    // get quaternion
    //**********************************************************************************/
    astro::quaternion_t q() const { return astro::quaternion_t(data[q4], data[q1], data[q2], data[q3]); }

    //**********************************************************************************/
    // normalize quaternion
    //**********************************************************************************/
    inline void normalize() {
      
      double norm = sqrt((data[q1] * data[q1]) + (data[q2] * data[q2]) + (data[q3] * data[q3]) + (data[q4] * data[q4]));
      
      data[q1] /= norm;
      data[q2] /= norm;
      data[q3] /= norm;
      data[q4] /= norm;
    
    }
    
    //**********************************************************************************/
    // operator ()
    //**********************************************************************************/
    inline void operator () (double Wx, double Wy, double Wz, double Q4, double Q1, double Q2, double Q3, double Jx, double Jy, double Jz) {
      
      data[wx] = Wx; data[wy] = Wy; data[wz] = Wz;
      
      data[q4] = Q4; data[q1] = Q1; data[q2] = Q2; data[q3] = Q3;

      double Ix = (Jy-Jz) / Jx;
      double Iy = (Jz-Jx) / Jy;
      double Iz = (Jx-Jy) / Jz;
      
      data[ix] = Ix; data[iy] = Iy; data[iz] = Iz;
      
    }
    
    //**********************************************************************************/
    // println
    //**********************************************************************************/
    inline void println(FILE * output = stdout) {
      
      fprintf(output, "%f %f %f %f %f %f %f %f %f %f\n", data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9]);
      
    }
    
    //**********************************************************************************/
    // getAngles
    //**********************************************************************************/
    template <typename T>
    inline void getAngles(T & pitch, T & yaw, T & roll) const {
    
      astro::quaternion_t(data[q4], data[q1], data[q2], data[q3]).getAngles(pitch, yaw, roll);
      
    }
    
  };
  
  
  //**********************************************************************************/
  // namespace attitude
  //**********************************************************************************/
  namespace attitude {
    
    void integration(const double * const input, double * output) {
      
      // data[10] = { Wx Wy Wz q4 q1 q2 q3 Ix Iy Iz }
      //               0  1  2  3  4  5  6  7  8  9
      
      output[wx] = (input[wy]*input[wz]*input[ix]);
      output[wy] = (input[wx]*input[wz]*input[iy]);
      output[wz] = (input[wx]*input[wy]*input[iz]);
      
// QUELLO BUONO ERA QUESTO
      output[q1] =  ( input[wx]*input[q4] - input[wy]*input[q3] + input[wz]*input[q2]) * 0.5;
      output[q2] =  ( input[wx]*input[q3] + input[wy]*input[q4] - input[wz]*input[q1]) * 0.5;
      output[q3] =  (-input[wx]*input[q2] + input[wy]*input[q1] + input[wz]*input[q4]) * 0.5;
      output[q4] =  (-input[wx]*input[q1] - input[wy]*input[q2] - input[wz]*input[q3]) * 0.5;
      
//      output[q4] =  ( inuput[wx]*inuput[q1] - inuput[wy]*inuput[q3] + inuput[wz]*inuput[q2]) * 0.5;
//      output[q1] =  (-inuput[wx]*inuput[q4] - inuput[wy]*inuput[q2] - inuput[wz]*inuput[q3]) * 0.5;
//      output[q2] =  ( inuput[wx]*inuput[q3] + inuput[wy]*inuput[q1] - inuput[wz]*inuput[q4]) * 0.5;
//      output[q3] =  (-inuput[wx]*inuput[q2] + inuput[wy]*inuput[q4] + inuput[wz]*inuput[q1]) * 0.5;

//      output[q1] =  (-inuput[wx]*inuput[q2] - inuput[wy]*inuput[q3] - inuput[wz]*inuput[q4]) * 0.5;
//      output[q2] =  ( inuput[wx]*inuput[q1] - inuput[wy]*inuput[q4] + inuput[wz]*inuput[q3]) * 0.5;
//      output[q3] =  ( inuput[wx]*inuput[q4] + inuput[wy]*inuput[q1] - inuput[wz]*inuput[q2]) * 0.5;
//      output[q4] =  (-inuput[wx]*inuput[q3] + inuput[wy]*inuput[q2] + inuput[wz]*inuput[q1]) * 0.5;

//      output[ix] = input[ix];
//      output[iy] = input[iy];
//      output[iz] = input[iz];
      
      output[ix] = 0;
      output[iy] = 0;
      output[iz] = 0;

      //printf("W %e %e %e\n", input[wx], input[wy], input[wz]);
      //printf("I %e %e %e\n", input[ix], input[iy], input[iz]);
      
    }
  
    //**********************************************************************************/
    // Compute attitude for just multiple dt
    //**********************************************************************************/
    void compute(const attitude_t & attitude0, std::size_t steps, double dt, std::vector<attitude_t> & attitude) {
      
      attitude.resize(steps, attitude0);
      
      //attitude[0] = attitude0;
      
      for(std::size_t i=1; i<steps; ++i) {
        
        astro::rk4<10>(attitude[i-1].data, dt, attitude[i].data, attitude::integration);
        
        attitude[i].normalize();
        
      }
      
    }
    
  } /* namespace attitude */
  
} /* namespace astro */

#endif /* _H_ASTRO_ATTITUDE_H_ */

