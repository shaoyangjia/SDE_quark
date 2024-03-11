/*
Definition of the class SDE_quark_iter for the iteration operations
in the iterative solver of the SDE for the quark propagator in the
Euclidean space.

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
// class SDE_quark_iter ( Zker ): public SF
using namespace std;
class SDE_quark_iter: public SF
{
public:
  /* constructor
  d_scale: IR strength of the Maris-Tandy model in GeV^2.
  omega: IR scale of the Maris-Tandy model in GeV.
  is_pv_mt: if Pauli-Villars regularization is applied, default false.
  lambda_pv: mass of the Pauli-Villars regulartor in GeV.
  c_scale: scale of the radial variable transformation in GeV^2.
  n_grid: number of grid points for the Euclidean-space solution.
  m_rnm: renormalized mass in GeV.
  mu2: renormalization scale in GeV^2.
  uv_cutoff: UV cutoff of the radial integrals in GeV^2.
  quark_flavor: label of the quark flavor. 0 = light, 1 = up, 2 = down,
    3 = charm, 4 = strange, 5 = top, 6 = bottom.
  is_linear_ren_cond: select if linear renormalization condition is applied. */
  SDE_quark_iter( double d_scale, double omega, bool is_pv_mt, double lambda_pv,
    double c_scale, unsigned n_grid, double m_rnm, double mu2, double uv_cutoff,
    unsigned quark_flavor, bool is_linear_ren_cond ): SF( c_scale, n_grid, m_rnm,
      mu2, uv_cutoff, quark_flavor ), omega( omega ), d_scale( d_scale ),
  is_pv_mt( is_pv_mt ), lambda_pv( lambda_pv ), is_linear_ren_cond(
    is_linear_ren_cond )
  { // Allocate for the kernels of the iterative solver.
    sde_ker_mu2[0] = new double[ n_grid ];
    sde_ker_mu2[1] = new double[ n_grid ];
    sde_ker_grid[0] = new double*[ n_grid ];
    sde_ker_grid[1] = new double*[ n_grid ];
    for( unsigned ind_c = 0; ind_c < 2; ind_c ++ )
    {
      for( unsigned i = 0; i < n_grid; i++ )
      { sde_ker_grid[ind_c][i] = new double[ n_grid ]; }
    }
  }
  ~SDE_quark_iter( void )
  {
    delete[] sde_ker_mu2[0]; delete[] sde_ker_mu2[1];
    for( unsigned i = 0; i < n_grid; i++ )
    { delete[] sde_ker_grid[0][i]; delete[] sde_ker_grid[1][i]; }
    delete[] sde_ker_grid[0]; delete[] sde_ker_grid[1];
  }
  // Compute the kernel arrays of the iterative solver. num_processes: number \
  of processes in OpenMP.
  void init_sde_ker( unsigned num_processes )
  {
    double x, y;
    unsigned ind_x, ind_y;
    Zker<double> zker( omega, d_scale, is_pv_mt, lambda_pv );
    #pragma omp parallel for if( num_processes > 0 ) \
      default( shared ) private( ind_x, ind_y, x, y ) firstprivate( zker ) \
      schedule( dynamic ) num_threads( num_processes )
    for( ind_x = 0; ind_x < n_grid; ind_x++ )
    {
      x = SF::p2_grid[ ind_x ];
      for( ind_y = 0; ind_y < n_grid; ind_y++ )
      {
        y = SF::p2_grid[ ind_y ];
        sde_ker_grid[0][ind_x][ind_y] = zker.cz_itg( x, y, true );
        sde_ker_grid[1][ind_x][ind_y] = zker.cz_itg( x, y, false );
      }
      sde_ker_mu2[0][ind_x] = zker.cz_itg( mu2, x, true );
      sde_ker_mu2[1][ind_x] = zker.cz_itg( mu2, x, false );
    }
  }
  // Take a single step of iteration.
  void run_sde_iter( void )
  { iter_AB( true ); iter_AB( false ); ren_cond(); }
  // Change if the linear renormalization condition is applied.
  void set_is_linear_ren_cond( bool is_linear_ren_cond_in )
  { is_linear_ren_cond = is_linear_ren_cond_in; }
  // Change if the Pauli-Villars regularization is applied.
  void set_is_pv_mt( bool is_pv_mt_in ){ is_pv_mt = is_pv_mt_in; }
  // Set the mass of the Pauli-Villars regulator.
  void set_lambda_pv( double lambda_pv_in ){ lambda_pv = lambda_pv_in; }
  // functions to get model paramters
  double get_omega( void ){ return omega; }
  double get_d_scale( void ){ return d_scale; }
  bool get_is_linear_ren_cond( void ){ return is_linear_ren_cond; }
  bool get_is_pv_mt( void ){ return is_pv_mt; }
  double get_lambda_pv( void ){ return lambda_pv; }
private:
  double omega, d_scale, lambda_pv;
  bool is_pv_mt, is_linear_ren_cond;
  double* sde_ker_mu2[2];// kernels for the renormalization point
  double** sde_ker_grid[2];// kernels for the loop integrals
  double coe_sf = 6.0 * pow( PI, 3 );// constant
  /* Single step of iteration for either of the component. is_vector: ture for
  the Dirac vector component, false for the scalar component. */
  void iter_AB( bool is_vector )
  {
    unsigned ind_c = unsigned( not( is_vector ) );
    double self_energy_rdd;
    double *fac;
    fac = new double[n_grid];
    for( unsigned i = 0; i < n_grid; i++ )
    { fac[i] = weighted_prop_ker( i, is_vector ); }
    for( unsigned i = 0; i < n_grid; i++ )
    {
      self_energy_rdd = 0.0;
      for( unsigned j = 0; j < n_grid; j++ )
      { self_energy_rdd += sde_ker_grid[ind_c][i][j] * fac[j]; }
      if( is_vector )
      { SF::a_grid[i] = inv_prop_sde( self_energy_rdd, true ); }
      else
      { SF::b_grid[i] = inv_prop_sde( self_energy_rdd, false ); }
    }
    delete[] fac;
    fac = NULL;
  }
  // Compute the self-energy at the renormalization point followed by the \
  renormalization constants.
  void ren_cond( void )
  {
    double sigma_v_mu2 = 0.0;
    double sigma_s_mu2 = 0.0;
    for( unsigned i = 0; i < n_grid; i++ )
    {
      sigma_v_mu2 += sde_ker_mu2[0][i] * weighted_prop_ker( i, true );
      sigma_s_mu2 += sde_ker_mu2[1][i] * weighted_prop_ker( i, false );
    }
    sigma_v_mu2 /= coe_sf;
    sigma_s_mu2 /= coe_sf;
    ren_cond_elem( sigma_v_mu2, sigma_s_mu2 );
  }
  // factor of the quark propagator in the kernel function
  double prop_ker( double a, double b, double p2, bool is_vector )
  {
    double fun_return = 1.0 / ( a * a + b * b / p2 );
    if( is_vector ){ fun_return *= a; }
    else{ fun_return *= b; }
    return fun_return;
  }
  // Factors of the quark propagator, quadrature weight, and variable \
  transformation weight.
  double weighted_prop_ker( unsigned i, bool is_vector )
  {
    return prop_ker( a_grid[i], b_grid[i], p2_grid[i], is_vector )
      * w_grid[i] * w_vt[i];
  }
  // Applies the renormalization condition to compute the renormalization \
  constants from self-energy.
  void ren_cond_elem( double sigma_v_mu2, double sigma_s_mu2 )
  {
    double a_mu2 = 1.0;
    double b_mu2 = SF::m_rnm;
    double z_4;
    if( is_linear_ren_cond )
    {
      SF::z_2 = a_mu2 - sigma_v_mu2;
      z_4 = ( b_mu2 - sigma_s_mu2 ) / SF::m_rnm;
    }
    else
    {
      double dlt = 4.0 * sigma_v_mu2 * a_mu2 + 1.0;
      if( dlt >= 0.0 )
      {
        SF::z_2 = ( sqrt( dlt ) - 1.0 ) / ( 2.0 * sigma_v_mu2 );
        z_4 = ( b_mu2 - z_2 * z_2 * sigma_s_mu2 ) / SF::m_rnm;
      }
      else
      {
        SF::z_2 = ( a_mu2 + sigma_v_mu2 ) / ( 1.0 + 2.0 * sigma_v_mu2 );
        z_4 = ( b_mu2 - sigma_s_mu2 / sigma_v_mu2 * ( a_mu2 - z_2 ) ) / SF::m_rnm;
      }
    }
    SF::z_m = z_4 / SF::z_2;
  }
  // Compute the inverse of the propagator from self-energy.
  double inv_prop_sde( double self_energy_rdd, bool is_vector )
  {
    double m_bare = SF::z_m * SF::m_rnm;
    if( not( is_linear_ren_cond ) )
    { self_energy_rdd *= SF::z_2 * SF::z_2; }
    if( is_vector )
    { return z_2 + self_energy_rdd / coe_sf; }
    else
    { return z_2 * m_bare + self_energy_rdd / coe_sf; }
  }
};
