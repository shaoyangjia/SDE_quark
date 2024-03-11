/*
Definition of the class Grid for the Euclidean-space solution of the SDE
for the quark propagator. Applies natural cubic spline interpolation when
accessing quark momenta not at the grid of the solution.

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
// class Grid ( NCS, sde_solution, sde_solution_array<> ) \
( #include"natural_cubic_spline.h" )
// #include<vector> #include<complex> #include<cmath>
// If IS_FREE_NCS ture, the following Fatal Python error could occure with \
  multiprocessing or OpenMP: Segmentation fault. Not present when firstprivate() \
  is not used in OpenMP for the quark propagator.
#define IS_FREE_NCS false
#include"natural_cubic_spline.h"
using namespace std;
class Grid
{
public:
  /* full constructor
  n: number of grid points.
  c: variable transformation scale in GeV^2.
  vector<double> u: grid for the mapped varabiel u.
  vector<double> a: Dirac vector component of the inverse propagator on the u grid.
  vector<double> b: Dirac scalar component of the inverse propagator on the u grid.
  uv_cutoff: cutoff of the radial integrals in GeV^2. */
  Grid( unsigned n, double c, vector<double> u, vector<double> a,
    vector<double> b, double uv_cutoff ): a_itp( u, a, n ), b_itp( u, b, n ),
  n_grid( n ), c_scale( c ), cutoff_p2( uv_cutoff ) {}
  // constructor from struct sde_solution
  Grid( sde_solution sf ): a_itp( sf.u, sf.a, sf.n ), b_itp( sf.u, sf.b, sf.n ),
  n_grid( sf.n ), c_scale( sf.c ), cutoff_p2( sf.uv_cutoff ) {}
  // constructor from struct sde_solution_array
  template<unsigned num_grid>
  Grid( sde_solution_array<num_grid> sf ): n_grid( sf.n ), c_scale( sf.c ),
    cutoff_p2( sf.uv_cutoff )
  {
    vector<double> u_grid, a_grid, b_grid;
    u_grid.resize( sf.n );
    a_grid.resize( sf.n );
    b_grid.resize( sf.n );
    for( unsigned ind = 0; ind < sf.n; ind++ )
    {
      u_grid[ind] = sf.u[ind];
      a_grid[ind] = sf.a[ind];
      b_grid[ind] = sf.b[ind];
    }
    a_itp( u_grid, a_grid, sf.n );
    b_itp( u_grid, b_grid, sf.n );
  }
  // Reduced constructor, applies in conjunction with either of the following \
  overloaded function reset_grid().
  Grid( void ) = default;
  // Update the Euclidean-space solution.
  void reset_grid( unsigned n, double c, vector<double> u, vector<double> a,
    vector<double> b, double uv_cutoff )
  {
    n_grid = n;
    c_scale = c;
    a_itp.init_grid( u, a, n );
    b_itp.init_grid( u, b, n );
    cutoff_p2 = uv_cutoff;
  }
  // Update the Euclidean-space solution by struct sde_solution.
  void reset_grid( sde_solution sf )
  {
    n_grid = sf.n;
    c_scale = sf.c;
    a_itp.init_grid( sf.u, sf.a, sf.n );
    b_itp.init_grid( sf.u, sf.b, sf.n );
    cutoff_p2 = sf.uv_cutoff;
  }
  unsigned sf_n_grid( void ){ return n_grid; }// Get number of grid points.
  double sf_c_scale( void ){ return c_scale; }// Get variable transformation scale.
  // Get the Euclidean-space solution.
  vector<double> sf_u_grid( void ){ return a_itp.get_var_grid(); }
  vector<double> sf_a_grid( void ){ return a_itp.get_val_grid(); }
  vector<double> sf_b_grid( void ){ return b_itp.get_val_grid(); }
  // Get the UV cutoff in GeV^2.
  double sf_uv_cutoff( void ){ return cutoff_p2; }
  // natural cubic spline interpolators for invere propagator in mapped variable
  double itp_au( double uq ){ return a_itp.itp( uq ); }
  double itp_bu( double uq ){ return b_itp.itp( uq ); }
  // natural cubic spline interpolators for invere propagator in momentum variable
  double itp_a( double p2 ){ return a_itp.itp( vt_u( p2 ) ); }
  double itp_b( double p2 ){ return b_itp.itp( vt_u( p2 ) ); }
  // radial variable transformation function
  double vt_p2( double u ){ return c_scale * u / ( 1.0 - u ); }
  // integral measure due to the radial variable transformation
  double vt_w( double u )
  {
    double fun_return = 1.0 / ( 1.0 - u );
    fun_return *= c_scale * fun_return;
    return fun_return;
  }
  // inverse of the radial variable transformation
  double vt_u( double p2 ){ return p2 / ( p2 + c_scale ); }
  // public member
  double cutoff_p2;
private:
  unsigned n_grid;
  double c_scale;
  NCS a_itp, b_itp;
};
