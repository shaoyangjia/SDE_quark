/* Definition of the class RVT as the radial variable transformation used
  in SDE_quark.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3). */
// class RVT
using namespace std;
class RVT
{
public:
  // Default constructor, applies in conjunction with set_c_scale().
  RVT( void ) = default;
  // Set the variable transformation scale.
  void set_c_scale( double c_scale_in ){ c_scale = c_scale_in; }
  // Full constructor, scale: variable transformation scale in GeV^2.
  RVT( double val_c_scale ): c_scale( val_c_scale ) {}
  // Get the variable transformation scale.
  double get_c_scale( void ){ return c_scale; }
  // inverse of the transformation
  double fun_vt_p2( double u ){ return c_scale * u * vmap( u ); }
  // change in the integral measure due to the transformation
  double fun_vt_w( double u ){ return c_scale * pow( vmap( u ), 2 ); }
  // function of the radial variable transformation
  double fun_vt_u( double p2 ){ return ivmap( p2, c_scale ); }
private:
  double c_scale;
  double vmap( double u ){ return 1.0 / ( 1.0 - u ); }
  double ivmap( double p2, double c ){ return p2 / ( p2 + c ); }
};
