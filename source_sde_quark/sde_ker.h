/*
Definition of the class Zker SDE_ker to compute kernel functions of the
self-energy from the SDE for the quark propagator. The template parameter T
specifies the shared type of the momentum variable and of the self-energy.

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
// class SDE_ker<typename T>: public Gfun
template<typename T>
class SDE_ker: public Gfun
{
public:
  /* full constructor
  omega: IR scale of the Maris-Tandy model in GeV.
  d_scale: IR strength of the Maris-Tandy model in GeV^2.
  is_pv_mt: if Pauli-Villars regularization is applied, default false.
  lambda_pv: mass of the Pauli-Villars regulartor in GeV. */
  SDE_ker( double omega, double d_scale, bool is_pv_mt, double lambda_pv ):
    Gfun( omega, d_scale, is_pv_mt, lambda_pv ) {}
  // Reduced constructor, applies in conjunction with \
  Gfun::set_gfun_pmts( bool is_pv_mt_in, double lambda_pv ).
  SDE_ker( double omega, double d_scale ): Gfun( omega, d_scale ) {}
  // Default constructor, applies in conjunction with \
  Gfun::set_gfun_pmts( double, double, bool, double )
  SDE_ker( void ): Gfun() {}
  /* Compute the kernel function applying symmetrization with respect to z = 0.
  x: external momentum.
  y: radial variable of the internal loop momentum.
  z: angular variable of the internal loop momentum.
  is_vector: ture if Dirac vector component is to be returned, false for Dirac
    scalar component. */
  T sde_ker( T x, double y, double z, bool is_vector )
  {
    T temp = z;
    if( z != 0.0 ){ temp *= sqrt( x * y ); }
    T term_p = Gfun::gfun_mt( x + y - 2.0 * temp );
    T term_m = Gfun::gfun_mt( x + y + 2.0 * temp );
    T factor_1s = term_p + term_m;
    T factor_1d = term_p - term_m;
    T factor_2s = 6.0;
    T factor_2d = 0.0;
    if( is_vector )
    {
      term_p = 1.0 / ( x + y - 2.0 * temp );
      term_m = 1.0 / ( x + y + 2.0 * temp );
      temp *= ( 1.0 + y / x );
      term_p *= temp - 2.0 * y;
      term_m *= - temp - 2.0 * y;
      if( z != 0.0 )
      {
        temp = 2.0 * sqrt( y / x ) * z;
        term_p += temp;
        term_m -= temp;
      }
      factor_2s = term_p + term_m;
      factor_2d = term_p - term_m;
    }
    return 0.25 * sqrt( 1.0 - pow( z, 2 ) ) *
      ( factor_1s * factor_2s + factor_1d * factor_2d );
  }
  // Overloaded version of the kernel function with the option to select if \
  the IR term of the Maris-Model is to be included.
  T sde_ker( T x, double y, double z, bool is_vector, bool is_include_IR )
  {
    T temp = z;
    if( z != 0.0 ){ temp *= sqrt( x * y ); }
    T term_p = Gfun::gfun_mt( x + y - 2.0 * temp, is_include_IR );
    T term_m = Gfun::gfun_mt( x + y + 2.0 * temp, is_include_IR );
    T factor_1s = term_p + term_m;
    T factor_1d = term_p - term_m;
    T factor_2s = 6.0;
    T factor_2d = 0.0;
    if( is_vector )
    {
      term_p = 1.0 / ( x + y - 2.0 * temp );
      term_m = 1.0 / ( x + y + 2.0 * temp );
      temp *= ( 1.0 + y / x );
      term_p *= temp - 2.0 * y;
      term_m *= - temp - 2.0 * y;
      if( z != 0.0 )
      {
        temp = 2.0 * sqrt( y / x ) * z;
        term_p += temp;
        term_m -= temp;
      }
      factor_2s = term_p + term_m;
      factor_2d = term_p - term_m;
    }
    return 0.25 * sqrt( 1.0 - pow( z, 2 ) ) *
      ( factor_1s * factor_2s + factor_1d * factor_2d );
  }
  // kernel function with the IR term of the Maris-Model only
  T sde_ker_ir( T x, double y, double z, bool is_vector )
  {
    T temp = z;
    if( z != 0.0 ){ temp *= sqrt( x * y ); }
    T term_p = Gfun::g_mt_ir( x + y - 2.0 * temp );
    T term_m = Gfun::g_mt_ir( x + y + 2.0 * temp );
    T factor_1s = term_p + term_m;
    T factor_1d = term_p - term_m;
    T factor_2s = 6.0;
    T factor_2d = 0.0;
    if( is_vector )
    {
      term_p = 1.0 / ( x + y - 2.0 * temp );
      term_m = 1.0 / ( x + y + 2.0 * temp );
      temp *= ( 1.0 + y / x );
      term_p *= temp - 2.0 * y;
      term_m *= - temp - 2.0 * y;
      if( z != 0.0 )
      {
        temp = 2.0 * sqrt( y / x ) * z;
        term_p += temp;
        term_m -= temp;
      }
      factor_2s = term_p + term_m;
      factor_2d = term_p - term_m;
    }
    return 0.25 * sqrt( 1.0 - pow( z, 2 ) ) *
      ( factor_1s * factor_2s + factor_1d * factor_2d );
  }
};
