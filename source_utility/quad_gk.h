/*
Definition of the class QUAD_GK for the adaptive Gauss-Kronrod quadrature.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark.

SDE_quark is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SDE_quark is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SDE_quark. If not, see <https://www.gnu.org/licenses/>.
*/
// class QUAD_GK<typename T>: private GK
using namespace std;
template<typename T>
class QUAD_GK: private GK
{
public:
  // Pure virtual kernl function that takes a real argument and returns T-type, \
  to be overriden by derived classes.
  virtual T ker( double arg ) = 0;
  /* Reset the condtion before running the quadrature.
  fmax: maximum number of function evaluations.
  abs_tol: absolute tolerance on a segment. */
  void reset( unsigned fmax = MAX_FUN_EVAL, double abs_tol = ABS_QUAD_TOL )
  {
    num_fun_eval = 0;
    max_fun_eval = fmax;
    quad_tol = abs_tol;
  }
  // Evaluate the integral on segment [a, b].
  T quad_gk( double a, double b )
  {
    if( num_fun_eval < max_fun_eval )
    {
      num_fun_eval += GK::nk;
      double coe = 0.5 * ( b - a );// middle point of the integral interval
      T y;// Stores the kernel function at a given quadrature point.
      T quad_g = 0.0;// quadrature result from Gauss rule
      T quad_k = 0.0;// quadrature result from Kronrod rule
      for(int i=0; i<nk; i++)
      {
        y = ker( coe * ( 1.0 + root_k[i] ) + a );
        quad_k += y * weight_k[i];
        if( i%2 == 1 )
        { quad_g += y * weight_g[ ( i - 1 ) / 2 ]; }
      }
      quad_g *= coe;
      quad_k *= coe;
      // Function recursion for the adaptive quadrature. Test if tolerance \
      condition is satisfied. Division of the segment is made if not.
      if( abs( quad_g - quad_k ) <= quad_tol )
      { return quad_k; }
      else
      {
        return quad_gk( a, a + coe ) + quad_gk( a + coe, b );
      }
    }
    else
    {
      cout<<"Maximium number of function evaluations reached by quad_gk(). "<<endl;
      return 0.0;
    }
  }
private:
  double quad_tol;
  unsigned num_fun_eval;
  unsigned max_fun_eval;
};
