/*
Definitions of class CT_ELEM and class SDE_SE_ASYMP for the separation
of the asymptotic contribution to integrals in the quark self-energy.

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
// class CT_ELEM: public SDE_SE_ASYMP_ELEM, public QUAD_GK<double>; \
class SDE_SE_ASYMP: public SDE_SE_ASYMP_ELEM; \
#include"sde_quark_prop_ker_asymp_elem.h"
using namespace std;
class CT_ELEM: public SDE_SE_ASYMP_ELEM, public QUAD_GK<double>
{
public:
  /* Full constructor for elements of counter terms.
  a0: UV value for the Dirac vector component of the inverse propagator.
  b0: UV value for the Dirac scalar component of the inverse propagator.
  p2_cutoff: UV cutoff of the radial integral in GeV^2.
  is_type_0: selection of the profile functions.
  is_vector: true if counter term is for the Dirac vector component, false for
    the scalar component.
  is_pv_mt: if Pauli-Villars regularization is applied.
  lambda_pv: Pauli-Villars regulator mass in GeV.
  */
  CT_ELEM( double a0, double b0, double p2_cutoff, bool is_type_0,
    bool is_vector, bool is_pv_mt, double lambda_pv ): is_type_0( is_type_0 ),
  is_vector( is_vector ), SDE_SE_ASYMP_ELEM( a0, b0, is_pv_mt, lambda_pv )
  {
    QUAD_GK<double>::reset();
    ct = QUAD_GK<double>::quad_gk( 0.0, p2_cutoff );
  }
  // Reduced constructor, applies in conjunction with reset_pmts().
  CT_ELEM( bool is_type_0, bool is_vector ): is_type_0( is_type_0 ),
  is_vector( is_vector ), SDE_SE_ASYMP_ELEM() {}
  // Change parameters for the counter term and update the integral.
  void reset_pmts( double a0, double b0, double p2_cutoff, bool is_pv_mt,
    double lambda_pv )
  {
    SDE_SE_ASYMP_ELEM::reset_asymp_elem( a0, b0, is_pv_mt, lambda_pv );
    QUAD_GK<double>::reset();
    ct = QUAD_GK<double>::quad_gk( 0.0, p2_cutoff );
  }
  // Get the counter term.
  double get_ct( void ){ return ct; }
  // kernel function for the Gauss-Kronrod integral
  double ker( double arg ) override
  {
    double asp, pfl;
    if( is_type_0 ){ asp = asp_0( arg ); }
    else{ asp = asp_1( arg ); }
    if( is_vector ){ pfl = profile_v( arg ); }
    else{ pfl = profile_s( arg ); }
    return asp * pfl;
  }
private:
  bool is_type_0, is_vector;
  double ct;
};

class SDE_SE_ASYMP: public SDE_SE_ASYMP_ELEM
{
public:
  /* Full constructor for the asymptotic subtraction and related variable
    transformations.
  a0: UV value for the Dirac vector component of the inverse propagator.
  b0: UV value for the Dirac scalar component of the inverse propagator.
  p2_cutoff: UV cutoff of the radial integral in GeV^2.
  is_pv_mt: if Pauli-Villars regularization is applied.
  lambda_pv: Pauli-Villars regulator mass in GeV. */
  SDE_SE_ASYMP( double a0, double b0, double p2_cutoff, bool is_pv_mt,
    double lambda_pv ):
  ct0s( a0, b0, p2_cutoff, true, false, is_pv_mt, lambda_pv ),
  ct0v( a0, b0, p2_cutoff, true, true, is_pv_mt, lambda_pv ),
  ct1v( a0, b0, p2_cutoff, false, true, is_pv_mt, lambda_pv ),
  SDE_SE_ASYMP_ELEM( a0, b0, is_pv_mt, lambda_pv )
  {
    ct_v = SDE_SE_ASYMP_ELEM::coe / a0 * (
      - 2.0 * SDE_SE_ASYMP_ELEM::coe_1 * ct0v.get_ct()
      + 6.0 * SDE_SE_ASYMP_ELEM::coe_2 * ct1v.get_ct() );
    ct_s = SDE_SE_ASYMP_ELEM::coe * b0 / ( a0 * a0 ) *
      3.0 * SDE_SE_ASYMP_ELEM::coe_0 * ct0s.get_ct();
  }
  // Reduced constructor, applies in conjunction with reset_pmts().
  SDE_SE_ASYMP( void ): ct0s( true, false ), ct0v( true, true ),
  ct1v( false, true ), SDE_SE_ASYMP_ELEM() {}
  // Change parameters and update counter terms.
  void reset_pmts( double a0, double b0, double p2_cutoff, bool is_pv_mt,
    double lambda_pv )
  {
    ct0s.reset_pmts( a0, b0, p2_cutoff, is_pv_mt, lambda_pv );
    ct0v.reset_pmts( a0, b0, p2_cutoff, is_pv_mt, lambda_pv );
    ct1v.reset_pmts( a0, b0, p2_cutoff, is_pv_mt, lambda_pv );
    SDE_SE_ASYMP_ELEM::reset_asymp_elem( a0, b0, is_pv_mt, lambda_pv );
    ct_v = SDE_SE_ASYMP_ELEM::coe / a0 * (
      - 2.0 * SDE_SE_ASYMP_ELEM::coe_1 * ct0v.get_ct()
      + 6.0 * SDE_SE_ASYMP_ELEM::coe_2 * ct1v.get_ct() );
    ct_s = SDE_SE_ASYMP_ELEM::coe * b0 / ( a0 * a0 ) *
      3.0 * SDE_SE_ASYMP_ELEM::coe_0 * ct0s.get_ct();
  }
  // counter terms as public members
  double ct_v, ct_s;
  /* asymptotic profile to be subtracted
  y: radial variable. z: angular variable.
  is_vecotr: ture for Dirac vector component, false for scalar. */
  double asp( double y, double z, bool is_vector )
  {
    if( is_vector ){ return SDE_SE_ASYMP_ELEM::asp_v( y, z ); }
    else{ return SDE_SE_ASYMP_ELEM::asp_s( y, z ); }
  }
  // inverse transformation for the radial variable y from \theta \
  at a govem externel momentum p2
  template<typename T>
  double vt_y( double theta, bool is_vector, T p2 )
  {
    if( is_vector ){ return ivartran( theta, cc_v ); }
    else{ return ivartran( theta, cc_s( p2 ) ); }
  }
  // integral measure change due to radial variable transformation from y to \theta.
  template<typename T>
  double w_vt( double theta, bool is_vector, T p2 )
  {
    if( is_vector ){ return wvartran( theta, cc_v ); }
    else{ return wvartran( theta, cc_s( p2 ) ); }
  }
  // radial variable transformation from y to \theta \
  at a given externel momentum p2
  template<typename T>
  double vt_theta( double y, bool is_vector, T p2 )
  {
    if( is_vector ){ return vartran( y, cc_v ); }
    else{ return vartran( y, cc_s( p2 ) ); }
  }
  // inverse transforamtion for the radial variable y from u \
  at a given externel momentum p2
  template<typename T>
  double uvt_y( double u, bool is_vector, T p2 )
  {
    if( is_vector ){ return ivt( u, cc_v ); }
    else{ return ivt( u, cc_s( p2 ) ); }
  }
  // integral measure change due to variable transformation from y to u.
  template<typename T>
  double w_uvt( double u, bool is_vector, T p2 )
  {
    if( is_vector ){ return wvt( u, cc_v ); }
    else{ return wvt( u, cc_s( p2 ) ); }
  }
  // radial variable transforamtion fro y to u at a given externel momentum p2
  template<typename T>
  double uvt_u( double y, bool is_vector, T p2 )
  {
    if( is_vector ){ return vt( y, cc_v ); }
    else{ return vt( y, cc_s( p2 ) ); }
  }
private:
  CT_ELEM ct0s, ct0v, ct1v;
  // scale of the radial variable transformation for the Dirac vector component
  double cc_v = 1.0;
  // scale of the radial variable transformation for the Dirac scalar component
  template<typename T>
  double cc_s( T p2 ){ return 1.0; }
  // common functions of the radial variable transformation y to \theta
  double vartran( double y, double c )
  { return acos( ( y - c ) / ( y + c ) ); }
  double ivartran( double theta, double c )
  {
    double cos_theta = cos( theta );
    return c * ( 1.0 + cos_theta ) / ( 1.0 - cos_theta );
  }
  double wvartran( double theta, double c )
  { return - 2.0 * c * sin( theta ) / pow( 1.0 - cos( theta ), 2 ); }
  // common functions of the radial variable transformation y to u
  double vt( double y, double c )
  { return ( y - c ) / ( y + c ); }
  double ivt( double u, double c )
  { return c * ( 1.0 + u ) / ( 1.0 - u ); }
  double wvt( double u, double c )
  { return 2.0 * c / pow( 1.0 - u, 2 ); }
};
