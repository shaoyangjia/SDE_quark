/*
Definition of the class Complex_prop that applies fixed-grid quadrature
rules to compute the quark propagator with real or complex-valued momentum
from the Schwinger-Dyson equation (SDE) and Euclidean-space solution. The
template parameter T specifies the shared type of the momentum variable and
of the inverse propagator.

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
// class Complex_prop<typename T> ( sde_solution, SDE_Solution_IO ): public SDE_SE<T>
// #include<vector> #include<complex> #include"sde_quark_prop.h"
using namespace std;
template<typename T>
class Complex_prop: public SDE_SE<T>
{
public:
  // Default constructor, applies in conjunction with set_complex_prop().
  Complex_prop( void ): SDE_SE<T>() {}
  /* Initialize parameters and Euclidean space solution.
  omega: IR scale of the Maris-Tandy model in GeV.
  d_scale: IR strength of the Maris-Tandy model in GeV^2.
  n_grid_gc: number of fixed-grid Gauss-Chebyshev quadrature points for angular
    integrals. */
  // Euclidean-space solution given by sf as a sde_solution struct.
  void set_complex_prop( double omega, double d_scale, unsigned n_grid_gc,
    sde_solution sf )
  {
    SDE_SE<T>::init_se( omega, d_scale, n_grid_gc, sf );
    z_2 = sf.z_2;
    m_bare = sf.z_m * sf.m_rnm;
    is_linear_ren_cond = sf.is_linear_ren_cond;
  }
  // Euclidean-space solution given by sf_io as an instance of class SDE_Solution_IO.
  void set_complex_prop( double omega, double d_scale, unsigned n_grid_gc,
    SDE_Solution_IO sf_io )
  {
    SDE_SE<T>::init_se( omega, d_scale, n_grid_gc, sf_io );
    z_2 = sf_io.get_z_2();
    m_bare = sf_io.get_m_bare();
    is_linear_ren_cond = sf_io.get_is_linear_ren_cond();
  }
  /* Original function for setting parameters and Euclidean-space solution.
  n: number of grid point in the Euclidean-space solution.
  c: scale of the radial variable transformation in GeV^2.
  vector<double> u: grid of the mapped radial variable u.
  vector<double> a: Dirac vector part of the inverse propagator at the u grid.
  vector<double> b: Dirac scalar part of the inverse propagator at the u grid.
  uv_cutoff: UV cutoff of the radial integration in GeV^2.
  val_z_2: wave function renormalization constant.
  val_m_bare: bare mass in GeV.
  bool_is_linear_ren_cond: if linear renormalization condition is applied,
    default true.
  is_pv_mt: if Pauli-Villars regularization is applied, default false.
  lambda_pv: mass of the Pauli-Villars regulartor in GeV. */
  void set_complex_prop( double omega, double d_scale, unsigned n, double c,
    vector<double> u, vector<double> a, vector<double> b, double uv_cutoff,
    double val_z_2, double val_m_bare, bool bool_is_linear_ren_cond,
    unsigned n_grid_gc, bool is_pv_mt, double lambda_pv )
  {
    SDE_SE<T>::init_se( omega, d_scale, n, c, u, a, b, uv_cutoff, n_grid_gc,
      is_pv_mt, lambda_pv );
    z_2 = val_z_2;
    m_bare = val_m_bare;
    is_linear_ren_cond = bool_is_linear_ren_cond;
  }
  // full constructor
  Complex_prop( double omega, double d_scale, unsigned n, double c,
    vector<double> u, vector<double> a, vector<double> b, double uv_cutoff,
    double val_z_2, double val_m_bare, bool is_linear_ren_cond,
    unsigned n_grid_gc, bool is_pv_mt, double lambda_pv ): SDE_SE<T>( omega,
    d_scale, n, c, u, a, b, uv_cutoff, n_grid_gc, is_pv_mt, lambda_pv ),
  z_2( val_z_2 ), m_bare( val_m_bare ), is_linear_ren_cond( is_linear_ren_cond )
  {}
  // constructor using struct sde_solution
  Complex_prop( double omega, double d_scale, unsigned n_grid_gc,
    sde_solution sf ): SDE_SE<T>( omega, d_scale, n_grid_gc, sf ), z_2( sf.z_2 ),
  m_bare( sf.z_m * sf.m_rnm ), is_linear_ren_cond( sf.is_linear_ren_cond ) {}
  // constructor using class SDE_Solution_IO
  Complex_prop( double omega, double d_scale, unsigned n_grid_gc,
    SDE_Solution_IO sf_io ): SDE_SE<T>( omega, d_scale, n_grid_gc, sf_io ),
  z_2( sf_io.get_z_2() ), m_bare( sf_io.get_m_bare() ),
  is_linear_ren_cond( sf_io.get_is_linear_ren_cond() ) {}
  // constructor loading solution from sf_filename named by sf_filename
  Complex_prop( double omega, double d_scale, unsigned n_grid_gc,
    string sf_filename ): SDE_SE<T>()
  {// Load Euclidean-space solution from binary file.
    SDE_Solution_IO sf_io( sf_filename );
    SDE_SE<T>::init_se( omega, d_scale, n_grid_gc, sf_io );
    reset_rnm_consts( sf_io.get_z_2(), sf_io.get_m_bare(),
      sf_io.get_is_linear_ren_cond() );
  }
  // Change the renormalization constants and conditions.
  void reset_rnm_consts( double val_z_2, double val_m_bare,
    bool bool_is_linear_ren_cond )
  {
    z_2 = val_z_2;
    m_bare = val_m_bare;
    is_linear_ren_cond = bool_is_linear_ren_cond;
  }
  // Change the Euclidean-space solution.
  void reset_grid( unsigned n, double c, vector<double> u, vector<double> a_grid,
    vector<double> b_grid, double uv_cutoff )
  { SDE_SE<T>::Grid::reset_grid( n, c, u, a_grid, b_grid, uv_cutoff ); }
  // functions to get model parameters and Euclidean-space solution
  double get_omega( void ){ return SDE_SE<T>::SDE_ker::Gfun::val_omega(); }
  double get_d_scale( void ){ return SDE_SE<T>::SDE_ker::Gfun::val_d_scale(); }
  unsigned get_n_grid( void ){ return SDE_SE<T>::Grid::sf_n_grid(); }
  double get_c_scale( void ){ return SDE_SE<T>::Grid::sf_c_scale(); }
  vector<double> get_u_grid( void ){ return SDE_SE<T>::Grid::sf_u_grid(); }
  vector<double> get_a_grid( void ){ return SDE_SE<T>::Grid::sf_a_grid(); }
  vector<double> get_b_grid( void ){ return SDE_SE<T>::Grid::sf_b_grid(); }
  double get_uv_cutoff( void ){ return SDE_SE<T>::Grid::sf_uv_cutoff(); }
  double get_z_2( void ){ return z_2; }
  double get_m_bare( void ){ return m_bare; }
  bool get_is_linear_ren_cond( void ){ return is_linear_ren_cond; }
  unsigned get_n_grid_gc( void ){ return SDE_SE<T>::get_n_grid_gc(); }
  bool get_is_pv_mt( void ){ return SDE_SE<T>::se_is_pv_mt(); }
  double get_lambda_pv( void ){ return SDE_SE<T>::se_lambda_pv(); }
  /* Function to compute the inverse of the quark propagator, returns scalar
  components as vector<T> = { A(p^2), B(p^2) }. */
  vector<T> inv_prop_complex( T p2 )
  {
    vector<T> fun_return (2, 0.0);
    fun_return[0] = inv_prop_ab( p2, true );
    fun_return[1] = inv_prop_ab( p2, false );
    return fun_return;
  }
  // function to compute the Dirac scalar component of the inverse propagator
  T inv_prop_b( T p2 ){ return inv_prop_ab( p2, false ); }
  /* Function to compute the inverse of the quark propagator.
  p2: momentum square. is_vector: true returns the Dirac vector component A(p^2),
    false for the Dirac scalar component B(p^2). */
  T inv_prop_ab( T p2, bool is_vector )
  {
    T fun_return = SDE_SE<T>::self_energy( p2, is_vector );
    fun_return = inv_prop_sde( fun_return, is_vector );
    return fun_return;
  }
  /* Function to compute the quark propagator. p2: momentum square. is_vector:
  true returns the Dirac vector component, false for the scalar component. */
  T prop_fun( T p2, bool is_vector )
  { // Dirac vector component of the inverse propagator
    T a_p2 = inv_prop_ab( p2, true );
    // Dirac scalar component of the inverse propagator
    T b_p2 = inv_prop_ab( p2, false );
    T fun_return = p2 * pow( a_p2, 2 ) + pow( b_p2 , 2 );// denominator
    if( is_vector )
    { fun_return = a_p2 / fun_return; }
    else
    { fun_return = b_p2 / fun_return; }
    return fun_return;
  }
  /* Overloaded prop_fun() function to compute the quark propagator, returns
  the vector and scalar components as array<T,2> = { vector, scalar }. */
  array<T,2> prop_fun( T p2 )
  {
    T a_p2 = inv_prop_ab( p2, true );
    T b_p2 = inv_prop_ab( p2, false );
    T denominator = p2 * pow( a_p2, 2 ) + pow( b_p2 , 2 );
    array<T,2> fun_return = { a_p2 / denominator, b_p2 / denominator };
    return fun_return;
  }
  // self-energy of the quark propagator modulo constant 1.0 / coe
  T self_energy( T p2, bool is_vector )
  { return SDE_SE<T>::self_energy( p2, is_vector ); }
private:
  /* variables */
  bool is_linear_ren_cond;
  double z_2, m_bare;
  double coe = 6.0 * pow( PI, 3 );// constant for the self-energy
  /* methods */
  // Computes the inverse propagator from the self-energy using the SDE.
  T inv_prop_sde( T self_energy_rdd, bool is_vector )
  {
    if( not( is_linear_ren_cond ) ){ self_energy_rdd *= z_2 * z_2; }
    if( is_vector ){ return z_2 + self_energy_rdd / coe; }
    else{ return z_2 * m_bare + self_energy_rdd / coe; }
  }
};
