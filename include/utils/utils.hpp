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

#ifndef _H_ASTRO_UTILS_H_
#define _H_ASTRO_UTILS_H_

#include <cstdlib>
#include <cstdio>

#include <glm/glm.hpp>

#include "angle.hpp"

/*****************************************************************************/
// namespace astro::utils
/*****************************************************************************/
namespace astro::utils {

  
  /*****************************************************************************/
  // elevazione
  /*****************************************************************************/
  double elevazione(glm::vec3 observatory, glm::vec3 satellite) {
    
    satellite.x   *= 6371.0; satellite.y   *= 6371.0; satellite.z   *= 6371.0;
    observatory.x *= 6371.0; observatory.y *= 6371.0; observatory.z *= 6371.0;
    
    glm::vec3 dist = observatory - satellite;
    
    return  glm::degrees(glm::acos(glm::dot(observatory,dist) / (glm::length(observatory)*glm::length(dist)))) - 90.0f;
    
  }
  
  /*****************************************************************************/
  // teta
  /*****************************************************************************/
  double teta(glm::vec3 observatory, glm::vec3 satellite) {
    
    satellite.x   *= 6371.0; satellite.y   *= 6371.0; satellite.z   *= 6371.0;
    observatory.x *= 6371.0; observatory.y *= 6371.0; observatory.z *= 6371.0;
    
    glm::vec3 dist = observatory - satellite;
    
    float elevazione = glm::degrees(glm::acos(glm::dot(observatory,dist) / (glm::length(observatory)*glm::length(dist)))) - 90.0f;
    
    return 90.0 - elevazione;;
    
  }
  
  
  /*****************************************************************************/
  // air attenuation factor
  /*****************************************************************************/
  double air_attenuation_factor(glm::vec3 observatory, glm::vec3 satellite) {
    
    satellite.x   *= 6371.0; satellite.y   *= 6371.0; satellite.z   *= 6371.0;
    observatory.x *= 6371.0; observatory.y *= 6371.0; observatory.z *= 6371.0;
    
    glm::vec3 dist = observatory - satellite;
    
    float elevazione = glm::degrees(glm::acos(glm::dot(observatory,dist) / (glm::length(observatory)*glm::length(dist)))) - 90.0f;
    
    // teta è la distanza dallo zenith
    float teta = 90.0 - elevazione;
    
    float secd = 1.0 / cos(glm::radians(teta));
    
    // valore dell’air mass alla data teta
    float am = secd-0.0018167*(secd-1)-0.0002875*pow((secd-1),2)-0.0008083*pow((secd-1),3);
    
    // assorbività dell'atmosfera
    float alfa = 0.25;
    float Q0 = 1.0;
    
    // la frazione di luce che arriva rispetto ad 1
    float Q = Q0 * exp(-alfa*am);
    
    return Q;
    
  }
  
  /*****************************************************************************/
  // converte da ra e dec geocentrice a ra e dec topocentriche
  // le ra e dec devono essere in gradi
  // site invece le cordinate del sito in TEME
  // sat invece le cordinate del satellite in TEME
  /*****************************************************************************/
  void ra2tradec(double sat[3],  double tra, double tdec, double site[3]) {
    
    double rho[3];
    
    rho[0] = sat[0] - site[0];
    rho[1] = sat[1] - site[1];
    rho[2] = sat[2] - site[2];
    
    double normRho = astMath::mag(rho);
    
    double coord[2];
    
    coord[0] = atan2(rho[1],rho[1]);
    coord[1] = asin(rho[2] / normRho);

    if(coord[0] < 0.1) coord[0] += 2 * M_PI;
    
    tra  = Degrees(coord[0]);
    tdec = Degrees(coord[1]);

  }
  
  
} /* namespace astro */


#endif /* _H_ASTRO_UTILS_H_ */
