/*
Definition of the class ITP for linear interpolation.

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
// class ITP<typename T> ( #include<algorithm> )
#include<algorithm>
using namespace std;
template<typename T>
class ITP
{
public:
  // Default constructor, applies in conjunction with init_grid().
  ITP( void ) = default;
  /* Full constructor
  var: grid for the vriable.
  val: values of the function on the variable grid.
  num: number of grid points. */
  ITP( vector<double> var, vector<T> val, unsigned num )
  { init_grid( var, val, num ); }
  void init_grid( vector<double> var, vector<T> val, unsigned num )
  {
    n_grid = num;
    var_0 = var[0]; val_0 = val[0];
    var_n = var[num - 1]; val_n = val[num -1];
    var_grid.resize( n_grid );
    val_grid.resize( n_grid );
    for( unsigned i = 0; i < n_grid; i++ )
    {
      var_grid[i] = var[i];
      val_grid[i] = val[i];
    }
  }
  vector<double> get_var_grid( void ){ return var_grid; }
  vector<T> get_val_grid( void ){ return val_grid; }
  // main function of the interpolation with extrapolations to end points
  T itp( double x )
  {
    if( x <= var_0 ){ return val_0; }
    else if( x >= var_n )
    { return val_n; }
    else
    {
      vector<double>::iterator up;
      up = upper_bound( var_grid.begin(), var_grid.end(), x );
      unsigned ind_up = up - var_grid.begin();
      if( ind_up == 0 ){ return val_grid[0]; }
      return litp_seg(x, var_grid[ind_up - 1], var_grid[ind_up],
        val_grid[ind_up - 1], val_grid[ind_up] );
    }
  }
private:
  /* variables */
  unsigned n_grid;
  double var_0, var_n;
  T val_0, val_n;
  vector<double> var_grid;
  vector<T> val_grid;
  /* methods */
  // interpolator on a single segment
  T litp_seg(double x, double a, double b, T fa, T fb)
  {
    T fun_return;
    fun_return = ( x - a ) / ( b - a );
    fun_return *= fb - fa;
    fun_return += fa;
    return fun_return;
  }
};
