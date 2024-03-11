/*
"filename_complex_prop.h" contains several file names for class definitions
associated with Complex_prop<typename T> used in SDE_quark.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3).
*/
#include"source_utility/gk_rw.h"// class GK
#include"source_utility/quad_gk.h"// class QUAD_GK<typename T>: private GK
#include"source_gker/g_mt.h"// class Gfun;
// dependencies for the quark kernel arrays in the rest frame
#include"source_utility/gc_rw.h"// class GC
#include"source_utility/quad_gch.h"// class QUAD_GCH<typename T>: private GC;
#include"source_utility/quad_gc.h"// class QUAD_GC<typename T>: private GC;
#include"source_utility/quad_gh.h"// class QUAD_GH<typename T>: private GH51 ( #include"gh_rw.h" );
#include"source_utility/itp_1d.h"// class ITP<typename T> ( #include<algorithm> )
#include"source_utility/quad_trapz.h"// class TRAPZ<typename T>: public ITP<T>;
#include"source_sde_quark/sde_ker.h"// class SDE_ker<typename T>: public Gfun
#include"source_sde_quark/sde_solution.h"// struct sde_solution; \
struct sde_solution_array<unsigned num_grid> class SDE_Solution_IO ( sde_solution )
#include"source_sde_quark/quark_prop.h"// class Grid ( NCS, sde_solution, \
  sde_solution_array<> ) ( #include"natural_cubic_spline.h" )
#include"source_sde_quark/sde_quark_prop_ker_asymp_elem.h"// class SDE_SE_ASYMP_ELEM;
#include"source_sde_quark/sde_quark_prop_ker_asymp.h"// class CT_ELEM: \
public SDE_SE_ASYMP_ELEM, public QUAD_GK<double>; \
class SDE_SE_ASYMP: public SDE_SE_ASYMP_ELEM;
#include"source_sde_quark/var_tran_pf_z.h"// class VT_PMT_Z
#include"source_sde_quark/sde_quark_prop.h"// class SDE_SE<typename T>: \
public SDE_ker<T>, public Grid, public QUAD_GCH<T>, public QUAD_GC<T>, \
public TRAPZ<T>, public QUAD_GH<T>, public SDE_SE_ASYMP, public VT_PMT_Z;
#include"source_sde_quark/complex_prop.h"// class Complex_prop<typename T> \
( sde_solution, SDE_Solution_IO ): public SDE_SE<T>;
