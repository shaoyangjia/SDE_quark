/* Definition of the class SDE_SE_ker for the kernel function of the quark
  self-energy used in SDE_quark, with the interpolation of the Euclidean-space
  solution. The template parameter T specifies the shared type of the momentum
  variable and of the kernel function return.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3). */
// class SDE_SE_ker: public SDE_ker<T>, public Grid
// #include<complex> #include<vector> #include<cmath>
using namespace std;
template<typename T>
class SDE_SE_ker: public SDE_ker<T>, public Grid
{
public:
  /* constructor
  omega: IR scale of the Maris-Tandy model in GeV.
  d_scale: IR strength of the Maris-Tandy model in GeV^2.
  n: number of grid point in the Euclidean-space solution.
  c: scale of the radial variable transformation in GeV^2.
  vector<double> u: grid of the mapped radial variable u.
  vector<double> a: Dirac vector part of the inverse propagator at the u grid.
  vector<double> b: Dirac scalar part of the inverse propagator at the u grid.
  uv_cutoff: UV cutoff of the radial integration.
  is_pv_mt: if Pauli-Villars regularization is applied, default false.
  lambda_pv: mass of the Pauli-Villars regulartor in GeV. */
  SDE_SE_ker( double omega, double d_scale, unsigned n, double c,
    vector<double> u, vector<double> a, vector<double> b, double uv_cutoff,
    bool is_pv_mt, double lambda_pv ): SDE_ker<T>( omega, d_scale, is_pv_mt,
    lambda_pv ), Grid( n, c, u, a, b, uv_cutoff ) {}
  // functions to get model paramters
  double ker_omega( void ){ return SDE_ker<T>::Gfun::val_omega(); }
  double ker_d_scale( void ){ return SDE_ker<T>::Gfun::val_d_scale(); }
  unsigned ker_n_grid( void ){ return Grid::sf_n_grid(); }
  double ker_c_scale( void ){ return Grid::sf_c_scale(); }
  // functions to get Euclidean space solution
  vector<double> ker_u_grid( void ){ return Grid::sf_u_grid(); }
  vector<double> ker_a_grid( void ){ return Grid::sf_a_grid(); }
  vector<double> ker_b_grid( void ){ return Grid::sf_b_grid(); }
  // functions to get regularization parameters
  double ker_uv_cutoff( void ){ return Grid::sf_uv_cutoff(); }
  bool ker_is_pv_mt( void ){ return SDE_ker<T>::Gfun::val_is_pv_mt(); }
  double ker_lambda_pv( void ){ return SDE_ker<T>::Gfun::val_lambda_pv(); }
  // functions to comput the self-energy of quark propagator from SDE
  /* wrapper function
  p2: external momentum.
  val_z: angular integration variable.
  y: radial inegration variable.
  bool_is_vector: true returns the Dirac vector component of the kernel
    function, false for the Dirac scalar component. */
  T self_energy_ker( T p2, double val_z, double y, bool bool_is_vector )
  {
    x = p2;
    z = val_z;
    is_vector = bool_is_vector;
    return ker( y );
  }
  // kernel function for the y integral
  T ker( double yq )
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
    T fun_return = SDE_ker<T>::sde_ker(x, yq, z, is_vector) * prop_rdd;
    return fun_return;
  }
private:
  // variables
  T x;
  double z;
  bool is_vector;
};
