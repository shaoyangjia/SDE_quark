/*
Definition of the class TRAPZ for the trapezoidal quadrature.

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
// class TRAPZ<typename T>: public ITP<T> \
#include<vector> #include"itp_1d.h"
/* Uniform one-dimensional trapezoidal quadrature with the last element forced
to be zero. Linear extrapolation based on the quadrature grid can be accessed by
calling TRAPZ::ITP::itp( double arg ) for the integrand. */
using namespace std;
template<typename T>
class TRAPZ: public ITP<T>
{
public:
  // Full constructor, n: number of grid points.
  TRAPZ( unsigned n ): n_grid( n ) {}
  // Default constructor, applies in conjunction with set_n_trapz( unsigned ).
  TRAPZ( void ) = default;
  // Set the number of quadrature points.
  void set_n_trapz( unsigned n ){ n_grid = n; }
  // kernel function of the integral
  virtual T ker_trapz( double arg ) = 0;
  // Evaluate the integral while saving the kernel on the grid to the linear interpolator.
  T quad_trapz( double a, double b )
  {
    vector<double> var_grid;
    vector<T> val_grid;
    var_grid.resize( n_grid - 1 );
    val_grid.resize( n_grid - 1 );
    var_grid[ 0 ] = a;
    val_grid[ 0 ] = ker_trapz( a );
    double h = ( b - a ) / double( n_grid - 1 );
    T fun_return = 0.5 * val_grid[ 0 ];
    for( unsigned i = 1; i < n_grid - 1; i++ )
    {
      var_grid[ i ] = var_grid[ i - 1 ] + h;
      val_grid[ i ] = ker_trapz( var_grid[i] );
      fun_return += val_grid[ i ];
    }
    var_grid.push_back( b );
    val_grid.push_back( 0.0 );
    fun_return *= h;
    ITP<T>::init_grid( var_grid, val_grid, n_grid );
    return fun_return;
  }
private:
  unsigned n_grid;
};
