/*
Definition of the class SDE_quark for the iterative solver of the quark
propagator in the Euclidean space from its SDE.

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
// class SDE_quark ( sde_solution ): public SDE_quark_iter, private TOL_TEST
using namespace std;
class SDE_quark: public SDE_quark_iter, private TOL_TEST
{
public:
  /* full constructor
  d_scale: IR strength of the Maris-Tandy model in GeV^2.
  omega: IR scale of the Maris-Tandy model in GeV.
  is_pv_mt: if Pauli-Villars regularization is applied, default false.
  lambda_pv: mass of the Pauli-Villars regulartor in GeV.
  c_scale: scale of the radial variable transformation in GeV^2.
  n_grid: number of grid point in the Euclidean-space solution.
  m_rnm: renormalized quark mass in GeV.
  mu2: renormalization scale in GeV^2.
  uv_cutoff: UV cutoff of the radial integration in GeV^2.
  iter_rel_tol, iter_abs_tol: relative and absolute tolerances of the iteartive
    solver.
  n_grid_gc: number of fixed-grid Gauss-Chebyshev quadrature points for angular
    integrals.
  max_n_iter: maximum steps of iterations in the iteartive solver.
  quark_flavor: label for the quark flavor.
  is_linear_ren_cond: if linear renormalization condition is applied. */
  SDE_quark( double d_scale, double omega, bool is_pv_mt, double lambda_pv,
    double c_scale, unsigned n_grid, double m_rnm, double mu2, double uv_cutoff,
    double iter_rel_tol, double iter_abs_tol, unsigned max_n_iter,
    unsigned quark_flavor, bool is_linear_ren_cond ): SDE_quark_iter( d_scale,
    omega, is_pv_mt, lambda_pv, c_scale, n_grid, m_rnm, mu2, uv_cutoff,
    quark_flavor, is_linear_ren_cond ), TOL_TEST( iter_rel_tol, iter_abs_tol ),
  n_grid( n_grid ), max_n_iter( max_n_iter )
  {
    a_grid_old = new double[ n_grid ];
    b_grid_old = new double[ n_grid ];
  }
  ~SDE_quark( void ){ delete[] a_grid_old; delete[] b_grid_old; }
  // Run the iterative solver of the SDE with num_processes OpenMP processes.
  bool iter_solve_sde( unsigned num_processes )
  {
    SDE_quark_iter::init_sde_ker( num_processes );
    unsigned iter = 0;
    bool iter_cond = false;
    while( iter < max_n_iter )
    {
      iter += 1;
      for( unsigned i = 0; i < n_grid; i++ )
      { a_grid_old[i] = a_grid[i];
        b_grid_old[i] = b_grid[i]; }
      z_2_old = z_2;
      z_m_old = z_m;
      SDE_quark_iter::run_sde_iter();
      iter_cond = tol_test();
      if( iter_cond )
      { SDE_quark_iter::SF::iter_converge = true;
        cout<<"SDE for the quark propagator SDE converged at step "<<iter<<endl;
        break; }
      else
      { cout<<"number of iterations in the quark propagator SDE = "<<iter<<endl; }
    }
    return iter_cond;
  }
  // Display the solution.
  void print_solutions( void )
  {
    cout<<"z_2 = "<<z_2<<" z_m = "<<z_m<<endl;
    for( unsigned i = 0; i < n_grid; i++ )
    { cout<<"u = "<<u_grid[i]<<" a = "<<a_grid[i]<<" b = "<<b_grid[i]<<endl; }
  }
  sde_solution get_solution_struct( void )
  {
    sde_solution sf;
    sf.n = n_grid;
    sf.c = get_c_scale();
    sf.u.resize( n_grid );
    sf.a.resize( n_grid );
    sf.b.resize( n_grid );
    for( unsigned ind = 0; ind < n_grid; ind++ )
    {
      sf.u[ind] = SDE_quark_iter::SF::GC_GRID::u_grid[ind];
      sf.a[ind] = SDE_quark_iter::SF::a_grid[ind];
      sf.b[ind] = SDE_quark_iter::SF::b_grid[ind];
    }
    sf.mu2 = SDE_quark_iter::SF::mu2;
    sf.m_rnm = SDE_quark_iter::SF::m_rnm;// m_bare = m_rnm * z_m
    sf.uv_cutoff = SDE_quark_iter::SF::uv_cutoff;
    sf.z_2 = SDE_quark_iter::SF::z_2;
    sf.z_m = SDE_quark_iter::SF::z_m;// z_4 = z_2 * z_m
    sf.is_linear_ren_cond = SDE_quark_iter::get_is_linear_ren_cond();
    sf.is_pv_mt = SDE_quark_iter::get_is_pv_mt();
    sf.lambda_pv = SDE_quark_iter::get_lambda_pv();
    return sf;
  }
  // Save the solution in binary file.
  void save_solution( string filename )
  {
    SDE_Solution_IO sf_io( get_solution_struct(),
      SDE_quark_iter::SF::quark_flavor );
    sf_io.save_sf_bin( filename );
  }
  // Load solution from a binary file.
  void load_solution( string filename )
  {
    SDE_Solution_IO sf_io( filename );
    sde_solution sf = sf_io.get_sde_solution();
    n_grid = sf.n;
    set_c_scale( sf.c );
    SDE_quark_iter::SF::GC_GRID::u_grid = new double[ n_grid ];
    SDE_quark_iter::SF::GC_GRID::is_allocated_u_grid = true;
    SDE_quark_iter::SF::a_grid = new double[ n_grid ];
    SDE_quark_iter::SF::is_allocated_a_grid = true;
    SDE_quark_iter::SF::b_grid = new double[ n_grid ];
    SDE_quark_iter::SF::is_allocated_b_grid = true;
    for( unsigned ind = 0; ind < n_grid; ind++ )
    {
      SDE_quark_iter::SF::GC_GRID::u_grid[ind] = sf.u[ind];
      SDE_quark_iter::SF::a_grid[ind] = sf.a[ind];
      SDE_quark_iter::SF::b_grid[ind] = sf.b[ind];
    }
    SDE_quark_iter::SF::mu2 = sf.mu2;
    SDE_quark_iter::SF::m_rnm = sf.m_rnm;// m_bare = m_rnm * z_m
    SDE_quark_iter::SF::uv_cutoff = sf.uv_cutoff;
    SDE_quark_iter::SF::z_2 = sf.z_2;
    SDE_quark_iter::SF::z_m = sf.z_m;// z_4 = z_2 * z_m
    SDE_quark_iter::SF::iter_converge = true;
    SDE_quark_iter::set_is_linear_ren_cond( sf.is_linear_ren_cond );
    SDE_quark_iter::set_is_pv_mt( sf.is_pv_mt );
    SDE_quark_iter::set_lambda_pv( sf.lambda_pv );
  }
  // functions to get model parameters and solution
  double get_iter_rel_tol( void ){ return TOL_TEST::get_rel_tol(); }
  double get_iter_abs_tol( void ){ return TOL_TEST::get_abs_tol(); }
  double get_m_rnm( void ){ return SDE_quark_iter::SF::m_rnm; }
  double get_mu2( void ){ return SDE_quark_iter::SF::mu2; }
  unsigned get_max_n_iter( void ){ return max_n_iter; }
  unsigned get_quark_flavor( void ){ return SDE_quark_iter::SF::quark_flavor; }
  unsigned get_n_grid( void ){ return n_grid; }
  double get_c_scale( void ){ return SDE_quark_iter::SF::RVT::get_c_scale(); }
  vector<double> get_u_grid( void )
  {
    vector<double> fun_return;
    fun_return.resize( n_grid );
    for( unsigned i = 0; i < n_grid; i++ )
    { fun_return[i] = u_grid[i]; }
    return fun_return;
  }
  vector<double> get_a_grid( void )
  {
    vector<double> fun_return;
    fun_return.resize( n_grid );
    for( unsigned i = 0; i < n_grid; i++ )
    { fun_return[i] = a_grid[i]; }
    return fun_return;
  }
  vector<double> get_b_grid( void )
  {
    vector<double> fun_return;
    fun_return.resize( n_grid );
    for( unsigned i = 0; i < n_grid; i++ )
    { fun_return[i] = b_grid[i]; }
    return fun_return;
  }
  double get_z_2( void ){ return SDE_quark_iter::SF::z_2; }
  double get_z_m( void ){ return SDE_quark_iter::SF::z_m; }
  double get_uv_cutoff( void ){ return SDE_quark_iter::SF::uv_cutoff; }
private:
  unsigned n_grid;
  unsigned max_n_iter;
  double* a_grid_old;
  double* b_grid_old;
  double z_2_old, z_m_old;
  // function that test the tolerance conditions
  bool tol_test( void )
  {
    bool fun_return = true;
    fun_return &= test_tol_cond( z_2_old, z_2 );
    if( fun_return )
    {
      fun_return &= test_tol_cond( z_m_old, z_m );
      if( fun_return )
      {
        fun_return &= test_tol_cond( a_grid_old, a_grid, n_grid );
        if( fun_return )
        {
          fun_return &= test_tol_cond( b_grid_old, b_grid, n_grid );
        }
      }
    }
    return fun_return;
  }
};
