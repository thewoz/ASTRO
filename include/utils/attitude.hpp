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

//#include <glm/glm.hpp>
//#include <glm/gtc/quaternion.hpp>

#include "quaternion.hpp"
#include "rk4.hpp"

/*****************************************************************************/
// namespace astro
/*****************************************************************************/
namespace astro {
 
  enum { wx, wy, wz, q4, q1, q2, q3, ix, iy, iz };
  
  /*****************************************************************************/
  // struct attitude_t
  /*****************************************************************************/
  struct attitude_t {
    
    attitude_t() {}
    
    attitude_t(double Wx, double Wy, double Wz, double Q4, double Q1, double Q2, double Q3, double Jx, double Jy, double Jz) {
      (*this)(Wx, Wy, Wz, Q4, Q1, Q2, Q3, Jx, Jy, Jz);
    }
    
    double data[10];

    /*****************************************************************************/
    // get omega function
    /*****************************************************************************/
    double Wx() const { return data[wx]; }
    double Wy() const { return data[wy]; }
    double Wz() const { return data[wz]; }
    
    /*****************************************************************************/
    // get quaternion
    /*****************************************************************************/
    astro::quaternion_t q() const { return astro::quaternion_t(data[q4], data[q1], data[q2], data[q3]); }

    /*****************************************************************************/
    // normalize quaternion
    /*****************************************************************************/
    inline void normalize() {
      
      double norm = sqrt((data[q1] * data[q1]) + (data[q2] * data[q2]) + (data[q3] * data[q3]) + (data[q4] * data[q4]));
      
      //printf("%f %f %f %f - %f\n", data[q1], data[q2], data[q3], data[q4], norm);

      data[q1] /= norm;
      data[q2] /= norm;
      data[q3] /= norm;
      data[q4] /= norm;
      
      //norm = sqrt((data[q1] * data[q1]) + (data[q2] * data[q2]) + (data[q3] * data[q3]) + (data[q4] * data[q4]));
      
      //printf("%f %f %f %f - %f\n\n", data[q1], data[q2], data[q3], data[q4], norm);

    }
    
    /*****************************************************************************/
    // operator ()
    /*****************************************************************************/
    inline void operator () (double Wx, double Wy, double Wz, double Q4, double Q1, double Q2, double Q3, double Jx, double Jy, double Jz) {
      
      data[wx] = Wx; data[wy] = Wy; data[wz] = Wz;
      
      data[q4] = Q4; data[q1] = Q1; data[q2] = Q2; data[q3] = Q3;

      double Ix = (Jy-Jz) / Jx;
      double Iy = (Jz-Jx) / Jy;
      double Iz = (Jx-Jy) / Jz;
      
      data[ix] = Ix; data[iy] = Iy; data[iz] = Iz;
      
    }
    
    /*****************************************************************************/
    // println
    /*****************************************************************************/
    inline void println(FILE * output = stdout) {
      
      fprintf(output, "%f %f %f %f %f %f %f %f %f %f\n", data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9]);
      
    }
    
    /*****************************************************************************/
    // getAngles
    /*****************************************************************************/
    template <typename T>
    inline void getAngles(T & pitch, T & yaw, T & roll) const {
    
//      glm::dvec3 angles = glm::eulerAngles(glm::dquat(data[q4], data[q1], data[q2], data[q3]));
//      pitch = angles.x; yaw = angles.y; roll = angles.z;

      astro::quaternion_t(data[q4], data[q1], data[q2], data[q3]).getAngles(pitch, yaw, roll);
      
    }
    
  };
  
  /*****************************************************************************/
  // namespace attitude
  /*****************************************************************************/
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
    
//    /*****************************************************************************/
//    // Compute attitude for just a dt
//    /*****************************************************************************/
//    void compute(const attitude_t & attitude0, double dt, attitude_t & attitude) {
//      rk4<10>(attitude0.data, dt, attitude.data, attitude::integration);
//    }
//
//    /*****************************************************************************/
//    // Compute attitude for just a dt
//    /*****************************************************************************/
//    attitude_t compute(const attitude_t & attitude0, double dt) {
//      attitude_t attitude;
//      compute(attitude0, dt, attitude);
//      return attitude;
//    }
//
    /*****************************************************************************/
    // Compute attitude for just multiple dt
    /*****************************************************************************/
    void compute(const attitude_t & attitude0, std::size_t steps, double dt, std::vector<attitude_t> & attitude) {
      
      attitude.resize(steps, attitude0);
      
      //attitude[0] = attitude0;
      
      for(std::size_t i=1; i<steps; ++i) {
        
        astro::rk4<10>(attitude[i-1].data, dt, attitude[i].data, attitude::integration);
        
        attitude[i].normalize();
        
      }
      
    }
    
//    /*****************************************************************************/
//    // Compute attitude for just multiple dt
//    /*****************************************************************************/
//    std::vector<attitude_t> compute(const attitude_t & attitude0, std::size_t steps, double dt) {
//      std::vector<attitude_t> attitude;
//      compute(attitude0, steps, dt, attitude);
//      return attitude;
//    }
 
    /*****************************************************************************/
    // test
    /*****************************************************************************/
//    void test() {
//      
//      float pitch = glm::radians(45.0f);
//      float yaw = glm::radians(30.0f);
//      float roll = glm::radians(20.0f);
//      
//      glm::quat qGlm = glm::vec3(pitch, yaw, roll);
//      quaternion_t qMpl(pitch, yaw, roll);
//
//      glm::vec3 euAnglesGlm = glm::eulerAngles(qGlm);
//      glm::vec3 euAnglesMpl; qMpl.getAngles(euAnglesMpl.x, euAnglesMpl.y, euAnglesMpl.z);
//      
//      printf("glm %f %f %f\n", glm::degrees(euAnglesGlm.x), glm::degrees(euAnglesGlm.y), glm::degrees(euAnglesGlm.z));
//      printf("mpl %f %f %f\n", glm::degrees(euAnglesMpl.x), glm::degrees(euAnglesMpl.y), glm::degrees(euAnglesMpl.z));
//      
//      printf("glm q %f %f %f %f\n", qGlm.x, qGlm.y, qGlm.z, qGlm.w);
//      printf("mpl q %f %f %f %f\n", qMpl.x, qMpl.y, qMpl.z, qMpl.w);
//      
//      glm::quat qGlm2 =  glm::quat(qMpl.w, qMpl.x, qMpl.y, qMpl.z);
//      glm::vec3 euAnglesGlm2 = glm::eulerAngles(qGlm2);
//      printf("glm %f %f %f\n", glm::degrees(euAnglesGlm2.x), glm::degrees(euAnglesGlm2.y), glm::degrees(euAnglesGlm2.z));
//
//      attitude_t att = attitude_t(0,0,0, qMpl.w, qMpl.x, qMpl.y, qMpl.z, 1, 1, 1);
//      glm::vec3 euAnglesAtt; att.getAngles(euAnglesAtt.x, euAnglesAtt.y, euAnglesAtt.z);
//      printf("att %f %f %f\n", glm::degrees(euAnglesAtt.x), glm::degrees(euAnglesAtt.y), glm::degrees(euAnglesAtt.z));
//      quaternion_t attQ = att.q();
//      printf("att q %f %f %f %f\n", attQ.x, attQ.y, attQ.z, attQ.w);
//      
//      glm::mat3 matrixGlm = glm::mat3_cast(qGlm);
//      
//      printf("%f %f %f\n", matrixGlm[0][0], matrixGlm[0][1], matrixGlm[0][2]);
//      printf("%f %f %f\n", matrixGlm[1][0], matrixGlm[1][1], matrixGlm[1][2]);
//      printf("%f %f %f\n", matrixGlm[2][0], matrixGlm[2][1], matrixGlm[2][2]);
//      
//      double matrixMpl[3][3];
//      qMpl.getRotationMatrix(matrixMpl);
//      
//      printf("%f %f %f\n", matrixMpl[0][0], matrixMpl[0][1], matrixMpl[0][2]);
//      printf("%f %f %f\n", matrixMpl[1][0], matrixMpl[1][1], matrixMpl[1][2]);
//      printf("%f %f %f\n", matrixMpl[2][0], matrixMpl[2][1], matrixMpl[2][2]);
//      
//    }
    
    //  /*****************************************************************************/
    //  // struct inertia_t
    //  /*****************************************************************************/
    //  struct inertia_t {
    //
    //    // Principal axes of inertia
    //    double axes[3];
    //
    //    inertia_t() {}
    //    inertia_t(double axesX, double axesY, double axesZ) { (*this)(axesX, axesY, axesZ);   };
    //
    //    inline void operator () (double axesX, double axesY, double axesZ) {
    //      axes[0] = axesX; axes[1] = axesY; axes[2] = axesZ;
    //    }
    //
    //    inline double   operator [] (size_t index) const { return axes[index]; }
    //    inline double & operator [] (size_t index)       { return axes[index]; }
    //
    //  };
    //
    //  /*****************************************************************************/
    //  // struct rotation_t
    //  /*****************************************************************************/
    //  struct rotation_t {
    //
    //  public:
    //
    //    // Quaterion
    //    mpl::math::quaternion_t q;
    //
    //    rotation_t() {}
    //    rotation_t(double q4, double q1, double q2, double q3) { (*this)(q4, q1, q2, q3);   };
    //
    //    inline void operator () (double q4, double q1, double q2, double q3) {
    //
    //      q[0] = q4; q[1] = q1; q[2] = q2; q[3] = q3;
    //
    //      double check = fabs((q4*q4)+(q1*q1)+(q2*q2)+(q3*q3)-1);
    //
    //      if(check > 1.0e-08) printf("MERDA %f\n", check);
    //
    //    }
    //
    //    inline double   operator [] (size_t index) const { return q[index]; }
    //    inline double & operator [] (size_t index)       { return q[index]; }
    //
    //  };
    //
    //  /*****************************************************************************/
    //  // struct velocity_t
    //  /*****************************************************************************/
    //  struct velocity_t {
    //
    //    double omega[3];
    //
    //    velocity_t() {}
    //    velocity_t(double omegaX, double omegaY, double omegaZ) { (*this)(omegaX, omegaY, omegaZ);   };
    //
    //    inline void operator () (double omegaX, double omegaY, double omegaZ) {
    //      omega[0] = omegaX; omega[1] = omegaY; omega[2] = omegaZ;
    //    }
    //
    //    inline double   operator [] (size_t index) const { return omega[index]; }
    //    inline double & operator [] (size_t index)       { return omega[index]; }
    //
    //  };
    
    
    
  } /* namespace attitude */
  
} /* namespace astro */

#endif /* _H_ASTRO_ATTITUDE_H_ */

