/*
Definition of the class GC_GRID for the grid of Gauss-Chebyshev
quadrature rules.

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
// class GC_GRID: private GC
using namespace std;
class GC_GRID: private GC
{
public:
  /* constructor
  is_type_I: if Gauss-Chebyshev quadrature of the first kind is applied.
  n_grid: number of grid points.
  a: lower limit of integration.
  b: upper limit of integration. */
  GC_GRID( bool is_type_I, size_t n_grid, double a, double b ):
  GC( is_type_I, n_grid ), n_grid( n_grid )
  { // applying variable transform u = 2 * ( x - a ) / ( b - a ) - 1; \
    and x = ( b - a ) * ( u + 1 ) / 2 + a;
    double coe = 0.5 * ( b - a );
    double theta = 0.0;
    double rx;
    if( is_type_I )
    { theta = - 0.5 * phi; }
    u_grid = new double[ n_grid ];
    w_grid = new double[ n_grid ];
    is_allocated_u_grid = true;
    is_allocated_w_grid = true;
    for( size_t i = 0; i < n_grid; i++ )
    {
      theta += phi;
      rx = cos( theta );
      u_grid[ n_grid - i - 1 ] = coe * ( rx + 1.0 ) + a;
      w_grid[ n_grid - i - 1 ] = coe * GC::phi * GC::gc_w( rx );
    }
  }
  ~GC_GRID( void )
  {
    if( is_allocated_u_grid ){ delete[] u_grid; u_grid = NULL; }
    if( is_allocated_w_grid ){ delete[] w_grid; w_grid = NULL; }
  }
  bool is_allocated_u_grid = false;
  bool is_allocated_w_grid = false;
  size_t n_grid;
  double *u_grid;// quadrature roots
  double *w_grid;// quadrature weights
};
