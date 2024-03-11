/*
Definition of the class QUAD_GCH for fixed-grid Gauss-Chebyshev quadrature
that evaluates half of the integration segment [2a-b, b].
The integrand is supposed to be symmetric with respect to x = a such that the
full integration is the double of the half contribution on x \in [a, b].

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
// class QUAD_GCH<typename T>: private GC; \
#include<cmath>
// Gauss-Chebyshev quadrature, output is multiplied by 2.
using namespace std;
template<typename T>
class QUAD_GCH: private GC
{
public:
  /* full constructor
  type_I: choose if Gauss-Chebyshev quadrture of the first kind is applied.
  n_grid: number of grid points. */
  QUAD_GCH( bool type_I, unsigned n_grid ): GC( type_I, n_grid ) {}
  // Default constructor, applies in conjunction with resize_grid( bool, unsigned ).
  QUAD_GCH( void ) = default;
  // Modify the kind and number of grid points for the quadrature.
  void resize_grid( bool is_type_I_in, unsigned n_grid )
  {
    is_type_I = is_type_I_in;
    n = n_grid;
    init_grid();
  }
  // Modify the number of grid points for the quadrature.
  void resize_grid( unsigned n_grid )
  { n = n_grid; init_grid(); }
  // pure virtual integrand
  virtual T ker_gch( double arg ) = 0;
  // Evaluate the integral on [a,b], applying variable transform \
  u = 2 * ( x - a ) / ( b - a ) - 1; and x = ( b - a ) * ( u + 1 ) / 2 + a;
  T quad_gch( double a, double b )
  {
    double coe = 0.5 * ( b - a );
    double theta = 0.0;
    double rx, x;
    T fun_return = 0.0;
    if( is_type_I )
    { theta = - 0.5 * phi; }
    for( unsigned i = 0; i < n; i ++ )
    {
      theta += phi;
      rx = cos( theta );
      x = coe * ( rx + 1.0 ) + a;
      fun_return += ker_gch( x ) * gc_w( rx );
    }
    fun_return *= 2.0 * coe * phi;
    return fun_return;
  }
};
