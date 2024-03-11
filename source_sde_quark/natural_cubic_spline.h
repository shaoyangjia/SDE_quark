/* Definition of the class NCS used in SDE_quark as an interface to the
  natural cubic spline interpolator from GNU scientific library.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3). */
// class NCS; ( #include<gsl/gsl_errno.h> )
#include<gsl/gsl_spline.h>
using namespace std;
class NCS
{
public:
  // default constructor
  NCS( void ) = default;
  /* full constructor
  vector<double> x: grid of the variable.
  vector<double> y: value of the function on the grid.
  n: number of grid points. */
  NCS( vector<double> x, vector<double> y, unsigned n )
  { init_grid( x, y, n, false ); }
  // destructor
  #if IS_FREE_NCS
    ~NCS( void ){ free_interp(); }
  #endif
  // Reset the grid.
  void init_grid( vector<double> x, vector<double> y, unsigned n,
    bool is_free_iterp = false )
  {
    // Store the grid and function values.
    var_grid = x;
    val_grid = y;
    // Store the end-point values.
    x_0 = x[0]; x_n = x[n-1];
    y_0 = y[0]; y_n = y[n-1];
    // Convert std::vector into array.
    bool is_allocated_xy = false;
    double* x_grid = NULL;
    double* y_grid = NULL;
    x_grid = new double[n];
    y_grid = new double[n];
    is_allocated_xy = true;
    for( unsigned i = 0; i < n; i++ )
    { x_grid[i] = x[i]; y_grid[i] = y[i]; }
    // Release the interpolator, which could cause an error in certain \
    implementation of OpenMP.
    if( is_free_iterp ){ free_interp(); }
    // Allocate and compute the spline.
    acc = gsl_interp_accel_alloc();// Allocate the accelerator.
    spline = gsl_spline_alloc( gsl_interp_cspline, n );// Allocate the spline.
    gsl_spline_init( spline, x_grid, y_grid, n );// Compute the spline.
    // Free the arrays.
    if( is_allocated_xy )
    {
      delete[] x_grid; x_grid = NULL;
      delete[] y_grid; y_grid = NULL;
      is_allocated_xy = false;
    }
  }
  vector<double> get_var_grid( void ){ return var_grid; }
  vector<double> get_val_grid( void ){ return val_grid; }
  // spline interpolater function with end-point extrapolation
  double itp( double xq )
  {
    if( xq <= x_0 ){ return y_0; }
    else if( xq >= x_n ){ return y_n; }
    else{ return gsl_spline_eval( spline, xq, acc ); }
  }
private:
  /* stored values for the grid */
  vector<double> var_grid, val_grid;
  /* variables */
  double x_0, x_n, y_0, y_n;
  gsl_interp_accel *acc = NULL;
  gsl_spline *spline = NULL;
  /* methods */
  void free_interp( void )
  {
    gsl_spline_free( spline );
    gsl_interp_accel_free( acc );
  }
};
