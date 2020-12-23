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


#ifndef _H_ASTRO_RK4_H_
#define _H_ASTRO_RK4_H_

#include <cstdlib>
#include <cstdio>

#include <vector>
#include <array>
    
namespace astro {

    //****************************************************************************
    //
    //  Purpose:
    //
    //    RK4VEC takes one Runge-Kutta step for a vector ODE.
    //
    //  Discussion:
    //
    //    It is assumed that an initial value problem, of the form
    //
    //      du/dt = f(t, u)
    //      u(t0) = u0
    //
    //    is being solved.
    //
    //    If the user can supply current values of t, u, a stepsize dt, and a
    //    function to evaluate the derivative, this function can compute the
    //    fourth-order Runge Kutta estimate to the solution at time t+dt.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    09 October 2013
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input:
    //
    //      t0, the current time.
    //      u0, the solution estimate at the current time.
    //      dt, the time step.
    //
    //     void f(in, out)
    //      a function which evaluates the derivative, or right hand side of the problem.
    //
    //    Output: the fourth-order Runge-Kutta solution estimate at time t0+dt.
    //
    
    template<size_t N, size_t M = N>
    //****************************************************************************/
    // RK4 function via double *
    //****************************************************************************/
    void rk4(const double * u0, double dt, double * ut, void f(const double * const input, double * output)) {
      
      double u1[N];
      double u2[N];
      double u3[N];
      
      double k1[N];
      double k2[N];
      double k3[N];
      double k4[N];
      
      //  Get four sample values of the derivative
      f(u0, k1);
      
      for(int i=0; i<M; ++i)
        u1[i] = u0[i] + dt * k1[i] / 2.0;
      
      f(u1, k2);
      
      for(int i=0; i<M; ++i)
        u2[i] = u0[i] + dt * k2[i] / 2.0;
      
      f(u2, k3);
      
      for(int i=0; i<M; ++i)
        u3[i] = u0[i] + dt * k3[i];
      
      f(u3, k4);
      
      //  Combine them to estimate the solution
      for(int i=0; i<M; ++i)
        ut[i] = u0[i] + (dt/6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
            
    }
  
} /* namespace astro */

#endif /* _H_ASTRO_RK4_H_ */
