/*
Definition of the class QUAD_GH for 51-point Gauss-Hermite quadrature.

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
// class QUAD_GH<typename T>: private GH51 ( #include"gh_rw.h" );
#include"gh_rw.h"
template<typename T>
class QUAD_GH: private GH51
{
public:
  // pure virtual integrand
  virtual T ker_gh( double x ) = 0;
  // Evaluate the integral.
  T quad_gh( void )
  {
    T fun_return = 0.0;
    for( unsigned i = 0; i < n_grid; i++ )
    { fun_return += w[ i ] * ker_gh( r[i] ); }
    return fun_return;
  }
};
