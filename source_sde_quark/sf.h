/* Definition of class SF that contains the inverse quark propagator in the
  Euclidean space used in SDE_quark.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3). */
// class SF: public RVT, public GC_GRID
class SF: public RVT, public GC_GRID
{
public:
  /* constructor
  c_scale: variable transformation scale in GeV^2.
  n_grid: number of grid points.
  m_rnm: renormalized mass in GeV.
  mu2: renormalization scale in GeV^2.
  uv_cutoff: Uv cutoff of the radial integral in GeV^2.
  quark_flavor: label of the quark flavor. 0 = light, 1 = up, 2 = down,
    3 = charm, 4 = strange, 5 = top, 6 = bottom. */
  SF( double c_scale, unsigned n_grid, double m_rnm, double mu2,
    double uv_cutoff, unsigned quark_flavor ): RVT( c_scale ),
  GC_GRID( IS_SF_TYPE_I, n_grid, 0.0, RVT::fun_vt_u( uv_cutoff ) ),
  m_rnm( m_rnm ), mu2( mu2 ), uv_cutoff( uv_cutoff ), quark_flavor( quark_flavor )
  {
    p2_grid = new double[ n_grid ]; is_allocated_p2_grid = true;
    w_vt = new double[ n_grid ]; is_allocated_w_vt = true;
    a_grid = new double[ n_grid ]; is_allocated_a_grid = true;
    b_grid = new double[ n_grid ]; is_allocated_b_grid = true;
    for( unsigned i = 0; i < n_grid; i++ )
    {
      p2_grid[i] = RVT::fun_vt_p2( GC_GRID::u_grid[i] );
      w_vt[i] = RVT::fun_vt_w( GC_GRID::u_grid[i] );
      a_grid[i] = 1.0;
      b_grid[i] = m_rnm;
    }
  }
  ~SF( void )
  {
    if( is_allocated_p2_grid ){ delete[] p2_grid; }
    if( is_allocated_w_vt ){ delete[] w_vt; }
    if( is_allocated_a_grid ){ delete[] a_grid; }
    if( is_allocated_b_grid ){ delete[] b_grid; }
  }
  bool is_allocated_p2_grid = false;
  bool is_allocated_w_vt = false;
  bool is_allocated_a_grid = false;
  bool is_allocated_b_grid = false;
  double m_rnm, mu2, uv_cutoff;
  unsigned quark_flavor;
  double *p2_grid, *w_vt, *a_grid, *b_grid;
  bool iter_converge = false;
  double z_2 = 1.0;
  double z_m = 1.0;
};
