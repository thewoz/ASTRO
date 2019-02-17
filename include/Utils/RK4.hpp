/*
 * MIT License
 *
 * Copyright © 2017
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


#ifndef _H_RK4_H_
#define _H_RK4_H_

#include <cstdlib>
#include <cstdio>

#include <vector>
#include <array>
    
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
    /*****************************************************************************/
    // RK4 function via double *
    /*****************************************************************************/
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

#endif /* _H_RK4_H_ */
