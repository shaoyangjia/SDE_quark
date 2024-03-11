/*
Definition of the class GC for the step in the grid of Gauss-Chebyshev
quadrature.

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
// class GC
using namespace std;
class GC
{
public:
  /* full constructor
  is_type_I: if Gauss-Chebyshev quadrature of the first kind is applied.
  n_grid: number of grid points. */
  GC( bool type_I, size_t n_grid ): is_type_I( type_I ), n( n_grid )
  { init_grid(); }
  // Default constructor, applies together with init_grid() after setting \
  is_type_I and n.
  GC( void ) = default;
  /* members */
  bool is_type_I;
  size_t n;// number of grid points
  double phi;// step change in the angular variable
  /* methods */
  void init_grid( void )
  {
    if( is_type_I ){ phi = PI / double( n ); }
    else{ phi = PI / double( n + 1 ); }
  }
  double gc_w( double rx ){ return sqrt( 1.0 - rx * rx ); }// profile factor
};
