/* Definition of class VT_PMT_Z that specifies the transformation for the
  angular variable in the fixed-grid quadrature used in SDE_quark.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3). */
// class VT_PMT_Z
using namespace std;
class VT_PMT_Z
{
public:
  // constructor with val_zeta_0 the desired center location for the peak
  VT_PMT_Z( double val_zeta_0 = 0.5 ): zeta_0( val_zeta_0 ) {}
  double b_of_z0( double z_0 )
  {
    if( z_0 <= zeta_0 + 1.0e-9 )
    { return 8.944271768357949e-05; }
    else if( z_0 >= 0.9999 )
    { return alpha / ( 1.0 - z_0 ); }
    else
    { return find_b( rf_eps, 6367.0, z_0); }
  }
private:
  /* members */
  double rf_eps = 1.0e-12;
  double zeta_0 = 0.5;
  double alpha = 0.63678181;
  /* methods */
  double fzero( double b, double z_0 )
  { return atan( b * zeta_0 ) - z_0 * atan( b ); }
  double find_b_elem( double b0, double bn, double f0, double fn, double z_0 )
  {
    double bm = 0.5 * ( b0 + bn );
    if( bn - b0 < 1.0e-8 )
    { return bm; }
    else
    {
      double fm = fzero( bm, z_0 );
      if( f0 * fm < 0.0 )
      { return find_b_elem( b0, bm, f0, fm, z_0 ); }
      else
      { return find_b_elem( bm, bn, fm, fn, z_0 ); }
    }
  }
  // bisection root finder
  double find_b( double b0, double bn, double z_0 )
  {
    double f0 = fzero( b0, z_0 );
    double fn = fzero( bn, z_0 );
    if( f0 * fn > 0.0 )
    { cout<<"incorrect end points"<<endl; return 0.0; }
    else
    { return find_b_elem( b0, bn, f0, fn, z_0 ); }
  }
};
