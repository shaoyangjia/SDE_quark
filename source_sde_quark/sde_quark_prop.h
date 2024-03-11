/*
Definitions of class SDE_SE for the fixed-grid quadrature computation of the
quark self-energy from its SDE. The template parameter T specifies the shared
type of the momentum variable and of the self-energy.

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
// class SDE_SE<typename T>: public SDE_ker<T>, public Grid, public QUAD_GCH<T>, \
public QUAD_GC<T>, public TRAPZ<T>, public QUAD_GH<T>, public SDE_SE_ASYMP, \
public VT_PMT_Z \
#include<complex> #include<vector> #include<cmath>
using namespace std;
template<typename T>
class SDE_SE: public SDE_ker<T>, public Grid,
  public QUAD_GCH<T>, public QUAD_GC<T>, public TRAPZ<T>, public QUAD_GH<T>,
  public SDE_SE_ASYMP, public VT_PMT_Z
{
public:
  // Default constructor, apply in conjunction with init_se() to initialize \
  the quark propagator.
  SDE_SE( void ): SDE_ker<T>(), Grid(), QUAD_GCH<T>(), QUAD_GC<T>(), TRAPZ<T>(),
  SDE_SE_ASYMP() {}
  // Initialize the self-energy using struct sde_solution.
  void init_se( double omega, double d_scale, unsigned n_grid_gc,
    sde_solution sf )
  {
    omega_scale = omega;
    std_n_grid_gc = n_grid_gc;
    SDE_ker<T>::Gfun::set_gfun_pmts( omega, d_scale, sf.is_pv_mt, sf.lambda_pv );
    Grid::reset_grid( sf );
    QUAD_GCH<T>::resize_grid( false, n_grid_gc );
    QUAD_GC<T>::resize_gc_grid( false, n_grid_gc );
    TRAPZ<T>::set_n_trapz( n_grid_gc );
    SDE_SE_ASYMP::reset_pmts( sf.a[sf.n-1], sf.b[sf.n-1], sf.uv_cutoff,
      sf.is_pv_mt, sf.lambda_pv );
  }
  // Initialize the self-energy using class SDE_Solution_IO.
  void init_se( double omega, double d_scale, unsigned n_grid_gc,
    SDE_Solution_IO sf_io )
  {
    omega_scale = omega;
    std_n_grid_gc = n_grid_gc;
    SDE_ker<T>::Gfun::set_gfun_pmts( omega, d_scale, sf_io.get_is_pv_mt(),
      sf_io.get_lambda_pv() );
    Grid::reset_grid( sf_io.get_sde_solution() );
    QUAD_GCH<T>::resize_grid( false, n_grid_gc );
    QUAD_GC<T>::resize_gc_grid( false, n_grid_gc );
    TRAPZ<T>::set_n_trapz( n_grid_gc );
    SDE_SE_ASYMP::reset_pmts( sf_io.get_a_last(), sf_io.get_b_last(),
      sf_io.get_uv_cutoff(), sf_io.get_is_pv_mt(), sf_io.get_lambda_pv() );
  }
  /* Initialize the self-energy.
  omega: IR scale of the Maris-Tandy model in GeV.
  d_scale: IR strength of the Maris-Tandy model in GeV^2.
  n: number of grid point in the Euclidean-space solution.
  c: scale of the radial variable transformation in GeV^2.
  vector<double> u: grid of the mapped radial variable u.
  vector<double> a: Dirac vector part of the inverse propagator at the u grid.
  vector<double> b: Dirac scalar part of the inverse propagator at the u grid.
  uv_cutoff: UV cutoff of the radial integration in GeV^2.
  n_grid_gc: number of fixed-grid Gauss-Chebyshev quadrature points for angular
    integrals.
  is_pv_mt: if Pauli-Villars regularization is applied, default false.
  lambda_pv: mass of the Pauli-Villars regulartor in GeV. */
  void init_se( double omega, double d_scale, unsigned n, double c,
    vector<double> u, vector<double> a, vector<double> b,
    double uv_cutoff, unsigned n_grid_gc, bool is_pv_mt, double lambda_pv )
  {
    omega_scale = omega;
    std_n_grid_gc = n_grid_gc;
    SDE_ker<T>::Gfun::set_gfun_pmts( omega, d_scale, is_pv_mt, lambda_pv );
    Grid::reset_grid( n, c, u, a, b, uv_cutoff );
    QUAD_GCH<T>::resize_grid( false, n_grid_gc );
    QUAD_GC<T>::resize_gc_grid( false, n_grid_gc );
    TRAPZ<T>::set_n_trapz( n_grid_gc );
    SDE_SE_ASYMP::reset_pmts( a[n-1], b[n-1], uv_cutoff, is_pv_mt, lambda_pv );
  }
  // full constructor
  SDE_SE( double omega, double d_scale, unsigned n, double c, vector<double> u,
    vector<double> a, vector<double> b, double uv_cutoff, unsigned n_grid_gc,
    bool is_pv_mt, double lambda_pv ): omega_scale( omega ), std_n_grid_gc( n_grid_gc ),
  SDE_ker<T>( omega, d_scale, is_pv_mt, lambda_pv ), Grid( n, c, u, a, b, uv_cutoff ),
  QUAD_GCH<T>( false, n_grid_gc ), QUAD_GC<T>( false, n_grid_gc ),
  TRAPZ<T>( n_grid_gc ), SDE_SE_ASYMP( a[n-1], b[n-1], uv_cutoff, is_pv_mt,
    lambda_pv ) {}
  // full constructor using struct sde_solution
  SDE_SE( double omega, double d_scale, unsigned n_grid_gc, sde_solution sf ):
  omega_scale( omega ), std_n_grid_gc( n_grid_gc ), SDE_ker<T>( omega, d_scale,
    sf.is_pv_mt, sf.lambda_pv ), Grid( sf ), QUAD_GCH<T>( false, n_grid_gc ),
  QUAD_GC<T>( false, n_grid_gc ),TRAPZ<T>( n_grid_gc ), SDE_SE_ASYMP( sf.a[sf.n-1],
    sf.b[sf.n-1], sf.uv_cutoff, sf.is_pv_mt, sf.lambda_pv ) {}
  // full constructor using class SDE_Solution_IO
  SDE_SE( double omega, double d_scale, unsigned n_grid_gc, SDE_Solution_IO sf_io ):
  omega_scale( omega ), std_n_grid_gc( n_grid_gc ), SDE_ker<T>( omega, d_scale,
    sf_io.get_is_pv_mt(), sf_io.get_lambda_pv() ), Grid( sf_io.get_sde_solution() ),
  QUAD_GCH<T>( false, n_grid_gc ), QUAD_GC<T>( false, n_grid_gc ),
  TRAPZ<T>( n_grid_gc ), SDE_SE_ASYMP( sf_io.get_a_last(), sf_io.get_b_last(),
  sf_io.get_uv_cutoff(), sf_io.get_is_pv_mt(), sf_io.get_lambda_pv() ) {}
  // functions to get model parameters and Euclidean-space solution
  double se_omega( void ){ return SDE_ker<T>::Gfun::val_omega(); }
  double se_d_scale( void ){ return SDE_ker<T>::Gfun::val_d_scale(); }
  unsigned se_n_grid( void ){ return Grid::sf_n_grid(); }
  double se_c_scale( void ){ return Grid::sf_c_scale(); }
  vector<double> se_u_grid( void ){ return Grid::sf_u_grid(); }
  vector<double> se_a_grid( void ){ return Grid::sf_a_grid(); }
  vector<double> se_b_grid( void ){ return Grid::sf_b_grid(); }
  double se_uv_cutoff( void ){ return Grid::sf_uv_cutoff(); }
  unsigned get_n_grid_gc( void ){ return std_n_grid_gc; }
  bool se_is_pv_mt( void ){ return SDE_ker<T>::Gfun::val_is_pv_mt(); }
  double se_lambda_pv( void ){ return SDE_ker<T>::Gfun::val_lambda_pv(); }
  // wrapper function for test only
  T self_energy_z( T p2, double val_z, bool bool_is_vector )
  { x = p2; is_vector = bool_is_vector; return self_energy_ker_z( val_z ); }
  // integrand for the trapezoidal quadrature of the angular integral
  T ker_trapz( double zq ) override
  { return self_energy_ker_z( zq ); }
  // Self-energy of the quark propagator from the SDE, p2: external momentum. \
  is_vector_bool: true for Dirac vector component, false for scalar.
  T self_energy( T p2, bool is_vector_bool )
  {
    T trapz_offset = 0.0;
    x = p2;
    QUAD_GC<T>::resize_gc_grid( std_n_grid_gc );
    is_vector = is_vector_bool;
    if( real( x ) >= real_p2_threshold )
    {
      trapz_offset = 2.0 * TRAPZ<T>::quad_trapz( 0.0, 1.0 );
      double z_0 = 1.0 - val_det_z_0 / abs( p2 );
      b = VT_PMT_Z::b_of_z0( z_0 );
      inv_atan_b = 1.0 / atan( b );
    }
    T fun_return = QUAD_GCH<T>::quad_gch( 0.0, 1.0 );
    if( is_vector ){ fun_return = fun_return + SDE_SE_ASYMP::ct_v; }
    else{ fun_return = fun_return + SDE_SE_ASYMP::ct_s; }
    return fun_return + trapz_offset;
  }
  // integrand for the Gauss-Chebyshev quadrature of the angular variable
  T ker_gch( double zp ) override
  {
    if( real( x ) >= real_p2_threshold )
    {
      return ( self_energy_ker_z( vt_atan( zp ) ) -
        TRAPZ<T>::ITP::itp( vt_atan( zp ) ) ) * w_vt_atan( zp );
    }
    else{ return self_energy_ker_z( zp ); }
  }
  // radial integral for a given angular variable val_z
  T self_energy_ker_z( double val_z )
  {
    double theta_min = SDE_SE_ASYMP::vt_theta( Grid::cutoff_p2, is_vector, x );
    z = val_z;
    T fun_return = 0.0;
    re_sqrt_p2 = real( sqrt( complex<double>( x ) ) );
    if( ( real( x ) >= separation_threshold ) &&
      ( z * re_sqrt_p2 >= omega_scale * gh_threshold ) )
    {
      is_include_IR = false;// in the radial integrals with quad_gc
      fun_return = QUAD_GH<T>::quad_gh() + QUAD_GC<T>::quad_gc( PI, theta_min );
    }
    else
    {
      is_include_IR = true;// in the radial integrals with quad_gc
      fun_return = QUAD_GC<T>::quad_gc( PI, theta_min );
    }
    return fun_return;
  }
  // IR contribution to the radial integral when conditions are met
  T ker_gh( double arg ) override
  {
    double temp = omega_scale * arg + z * re_sqrt_p2;// sqrt( y )
    if( temp <= 0.0 ){ return 0.0; }// corresponds to the lower limit of the y-integral
    else
    {
      double yq = temp * temp;
      double prop_rdd = fun_prop_rdd( yq );
      T fun_return = SDE_ker<T>::sde_ker_ir( x, yq, z, is_vector )
        * prop_rdd;
      prop_rdd = fun_prop_rdd( yq + Grid::cutoff_p2 );
      fun_return -= SDE_ker<T>::sde_ker_ir( x, yq + Grid::cutoff_p2, z, is_vector )
        * prop_rdd;
      fun_return *= 2.0 * omega_scale * temp;
      return fun_return;
    }
  }
  // kernel function for the Gauss-Chebyshev quadrature of the radial integral
  T ker_gc( double theta ) override
  {
    double yq = SDE_SE_ASYMP::vt_y( theta, is_vector, x );
    double prop_rdd = fun_prop_rdd( yq );
    T fun_return = SDE_ker<T>::sde_ker( x, yq, z, is_vector, is_include_IR ) * prop_rdd;
    fun_return = fun_return - SDE_SE_ASYMP::asp( yq, z, is_vector );
    fun_return *= SDE_SE_ASYMP::w_vt( theta, is_vector, x );
    return fun_return;
  }
private:
  double omega_scale;
  // parameters that control the sepration of the IR and UV terms
  double separation_threshold = 7.5;
  double gh_threshold = 2.0;
  double re_sqrt_p2;
  bool is_include_IR;
  // working variables
  T x; // x = p2
  double z;
  bool is_vector;
  unsigned std_n_grid_gc;
  // quark propagator function in the intergrand
  double fun_prop_rdd( double yq )
  {
    double uq = Grid::vt_u( yq );
    double aq = Grid::itp_au( uq );
    double bq = Grid::itp_bu( uq );
    double prop_rdd = 0.0;
    double dq = aq * aq;
    if( yq != 0.0 )
    {
      dq += bq * bq / yq;
      if( is_vector ){ prop_rdd = aq / dq; }
      else{ prop_rdd = bq / dq; }
    }
    return prop_rdd;
  }
  // parameters and functions for the arctangent variable transform
  double real_p2_threshold = 0.21;
  double val_det_z_0 = 0.1;
  double b, inv_atan_b;
  double vt_atan( double zeta ){ return inv_atan_b * atan( b * zeta ); }
  double w_vt_atan( double zeta )
  { return inv_atan_b * b / ( 1.0 + b * b * zeta * zeta ); }
  double ivt_atan( double z ){ return tan( z / inv_atan_b ) / b; }
};
