/* Definition of the class SDE_SE_ASYMP_ELEM as the elementary blocks for the
  asymptotic behavior of the quark self-energy integrand used in SDE_quark.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3). */
// class SDE_SE_ASYMP_ELEM;
using namespace std;
class SDE_SE_ASYMP_ELEM
{
public:
  /* full constructor
  a0: UV value for the Dirac vector component of the inverse propagator.
  b0: UV value for the Dirac scalar component of the inverse propagator.
  is_pv_mt: if Pauli-Villars regularization is applied.
  lambda_pv: Pauli-Villars regulator mass in GeV. */
  SDE_SE_ASYMP_ELEM( double a0, double b0, bool is_pv_mt, double lambda_pv ):
  a0( a0 ), b0( b0 ), is_pv_mt( is_pv_mt )
  { lambda_pv_2 = pow( lambda_pv, 2 ); }
  // Default constructor, applies in conjunction with reset_asymp_elem().
  SDE_SE_ASYMP_ELEM( void ): a0( 0.0 ), b0( 0.0 ), is_pv_mt( true ),
  lambda_pv_2( 0.0 ) {}
  // Modify members for updates in paramters.
  void reset_asymp_elem( double a0_in, double b0_in, bool is_pv_mt_in,
    double lambda_pv )
  {
    a0 = a0_in;
    b0 = b0_in;
    is_pv_mt = is_pv_mt_in;
    lambda_pv_2 = pow( lambda_pv, 2 );
  }
  // public members
  double a0, b0;
  double coe = coe_nf( 4 );
  // coe_0 = integrate( sqrt( 1 - zeta ** 2 ), ( zeta, -1, 1 ) )
  double coe_0 = 0.5 * PI;
  // coe_1 = integrate( sqrt( 1 - zeta ** 2 ) ** 3, ( zeta, -1, 1 ) )
  double coe_1 = 3.0 * PI / 8.0;
  // coe_2 = integrate( zeta ** 2 * sqrt( 1 - zeta ** 2 ), (zeta, -1, 1 ) )
  double coe_2 = PI / 8.0;
  // base asymptotic functions
  double asp_0( double y )
  {
    double fun_return = y / cc_0;
    fun_return = ( fun_return + 1.0 ) * log( fun_return + EBN );
    fun_return *= cc_0;
    fun_return = 1.0 / fun_return;
    return fun_return;
  }
  double asp_1( double y )
  {
    double fun_return = y / cc_0;
    fun_return = 1.0 + 1.0 / log( fun_return + EBN );
    fun_return *= asp_0( y );
    return fun_return;
  }
  // profile function for the Dirac vector component
  double profile_v( double y )
  {
    double fun_return = pow( profile_b( y ), 2 );
    if( is_pv_mt )
    { fun_return *= profile_pv( y ); }
    return fun_return;
  }
  // profile function for the Dirac scalar component
  double profile_s( double y )
  {
    double fun_return = profile_b( y );
    if( is_pv_mt )
    { fun_return *= profile_pv( y ); }
    return fun_return;
  }
  // asymptotic function for the Dirac vector component
  double asp_v( double y, double z )
  {
    double z2 = z * z;
    double fun_return = coe / a0 * sqrt( 1.0 - z2 ) * profile_v( y );
    fun_return *= -2.0 * ( 1.0 - z2 ) * asp_0( y ) + 6.0 * z2 * asp_1( y );
    return fun_return;
  }
  // asymptotic function for the Dirac scalar component
  double asp_s( double y, double z )
  {
    double fun_return = coe * b0 / ( a0 * a0 ) * sqrt( 1.0 - z * z );
    fun_return *= 3.0 * asp_0( y ) * profile_s( y );
    return fun_return;
  }
private:
  // members
  bool is_pv_mt;
  double lambda_pv_2;
  double cc_0 = pow( LAMBDA_QCD, 2 );
  double coe_nf( int nf ){ return 48.0 * pow( PI, 2 ) / double( 33 - 2 * nf ); }
  double profile_b( double y ){ return y / ( y + cc_0 ); }
  double profile_pv( double y ){ return 1.0 / ( y / lambda_pv_2 + 1.0 ); }
};
