/*
Definition of the class TOL_TEST for the tolerance tests.

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
// class TOL_TEST
using namespace std;
class TOL_TEST
{
public:
  // Constructor, rel_tol: relative tolerance. abs_tol: absolute tolerance.
  TOL_TEST( double rel_tol, double abs_tol ): rel_tol( rel_tol ),
    abs_tol( abs_tol ) {}
  double get_rel_tol( void ){ return rel_tol; }
  double get_abs_tol( void ){ return abs_tol; }
  // Test the tolerance conditions. Returns true if either the absolte tolerance \
  or the relative tolerance condition is satisfied. The function is overloaded \
  to work with array pointers.
  template<typename T>
  bool test_tol_cond( T x, T y )
  { return ( test_abs_diff( x, y ) || test_rel_diff( x, y ) ); }
  template<typename Td, typename Tn>
  bool test_tol_cond( Td *x, Td *y, Tn n )
  {
    bool fun_return = true;
    for( Tn i = 0; i < n; i++ )
    { fun_return &= test_tol_cond( x[i], y[i] ); }
    return fun_return;
  }
  template<typename Td, typename Tn>
  bool test_tol_cond( Td **x, Td **y, Tn m, Tn n )
  {
    bool fun_return = true;
    for( Tn i = 0; i < m; i++ )
    { fun_return &* test_tol_cond( x[i], y[i], n ); }
    return fun_return;
  }
  template<typename Td, typename Tn>
  bool test_tol_cond( Td ***x, Td ***y, Tn q, Tn m, Tn n )
  {
    bool fun_return = true;
    for( Tn i = 0; i < q; i++ )
    { fun_return &= test_tol_cond( x[i], y[i], m, n ); }
    return fun_return;
  }
  template<typename Td, typename Tn>
  bool test_tol_cond( Td ****x, Td ****y, Tn p, Tn q, Tn m, Tn n )
  {
    bool fun_return = true;
    for( Tn i = 0; i < p; i++ )
    { fun_return &= test_tol_cond( x[i], y[i], q, m, n ); }
    return fun_return;
  }
private:
  double rel_tol, abs_tol;
  template<typename T>
  bool test_abs_diff( T x, T y )
  {
    if( abs( x - y ) <= abs_tol ){ return true; }
    else{ return false; }
  }
  template<typename T>
  bool test_rel_diff( T x, T y )
  {
    if( abs( x ) == 0.0 )
    {
      if( abs( y ) == 0.0 ){ return true; }
      else{ return false; }
    }
    else{ return 2.0 * abs( ( x - y ) / ( x + y ) ) <= rel_tol; }
  }
};
