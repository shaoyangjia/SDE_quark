# -*- coding: utf-8 -*-
"""
"test_sf_compiled.py" is the python script that constains the accuracy and
performance tests of the fixed-grid quadrature algorithm in SDE_quark.

Instructions:
(1) Modify the global variable mpp_num_processes to the desired number of
    processes suitable to the hardware running this program.
(2) Set the desired tolerance of quad_vec() by the global variable tol_quad_vec.
    Suggest keeping the default value of 1.0e-3.
(3) Set other paramters in the "running parameters" and "physical paramters"
    section of the main program.
(4) The default path for the binary file of the Euclidean-space solution is the
    parent folder of "/python/". Modify the argument of the line with
    "self.sde_light.load_solution( '..' )# change to the parent directory" if
    the binary file is located in another path.
(5) Run the script in Python environment with packages numpy, scipy, matplotlib,
    and pybind11.

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
"""
# management of working directories
from os import chdir, getcwd
# support of multiprocessing
from multiprocessing import Pool
from itertools import product as itp
from functools import partial
# numpy data structures and functions
import numpy as np
# nquad and quad_vec functions
from scipy import integrate
# configuration of pyplot
import matplotlib.pyplot as plt
font_size = 10
# Set the default text font size
plt.rc('font', size=font_size)
# Set the axes title font size
plt.rc('axes', titlesize=font_size)
# Set the axes labels font size
plt.rc('axes', labelsize=font_size)
# Set the font size for x tick labels
plt.rc('xtick', labelsize=font_size)
# Set the font size for y tick labels
plt.rc('ytick', labelsize=font_size)
# Set the legend font size
plt.rc('legend', fontsize=font_size)
# Set the font size of the figure title
plt.rc('figure', titlesize=font_size)
#
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica",
    #"font.size": 14
    "legend.fontsize": font_size,
})

# supporting functions
import tfc
# pybind11 compiled module for fixed-grid quadrature to compute the inverse propagator
from ker_itg import Complex_prop as ker_Complex_prop
# pybind11 compiled module for the iterative solver of the SDE in the Euclidean space
import sde_iter_solver
# pybind11 compiled module for kernel functions of the SDE
import sde_ker

# number of processes
mpp_num_processes = 36
# Desired tolerances of quad_vec(), default 1.0e-3.
tol_quad_vec = 1.0e-3

# Choose if the valence quark is the same with the anti-quark, True by default.
is_equal_cst = True
# Number of grid points for the fixed-grid angular integral, 201 by default.
sde_n_grid_gc = 201
# Choose if Pauli-Villars regularization is applied, False by default.
is_pv_mt_sde = False
# Mass of the Pauli-Villars regulartor in GeV, irrelevant if is_pv_mt_sde == False.
lambda_pv = 1.0e+4
# Selection of the renormalization condition, True by defult.
is_linear_ren_cond = True
# maximium number of iterations in the iterative solver
max_iter_sde = 75

'''
Class with functions operating the compiled module of the iterative solver of
the SDE.
d_scale: IR strength of the Maris-Tandy model GeV^2.
omega: IR scale of the Maris-Tandy model GeV.
c_scale: variable transformation scale of the solution grid in GeV^2.
n_grid: number of grid points for the solution.
m_rnm: renormalized mass in GeV.
mu2: renormalization scale in GeV^2.
uv_cutoff: UV cutoff in GeV^2.
iter_rel_tol iter_abs_tol: relative and absolute tolerances of the iterative solver
quark_flavor = 'light': flavor label of the quark propagator solved, 'light'
    for light quarks and 'strange' for strange quarks.
solution_file = []: name of the binary file that saves the solution.
'''
class SDE_quark_cpl:
    def __init__(self, d_scale, omega, c_scale, n_grid, m_rnm, \
                 mu2, uv_cutoff, iter_rel_tol, iter_abs_tol, \
                 quark_flavor = 'light', solution_file = [] ):
        self.init_sde_quark( d_scale, omega, c_scale, n_grid, m_rnm, \
                            mu2, uv_cutoff, iter_rel_tol, iter_abs_tol, \
                            quark_flavor = quark_flavor, \
                            solution_file = solution_file )
    # constructor for class SDE_quark_cpl
    def init_sde_quark(self, d_scale, omega, c_scale, n_grid, m_rnm, \
                       mu2, uv_cutoff, iter_rel_tol, iter_abs_tol, \
                       quark_flavor = 'light', solution_file = [] ):
        self.omega = omega
        self.d_scale = d_scale
        self.n_grid = n_grid
        self.c_scale = c_scale
        self.uv_cutoff = uv_cutoff
        self.m_rnm = m_rnm
        self.solution_file = solution_file
        if quark_flavor == 'light':
            unsigned_quark_flavor = 0
        elif quark_flavor == 'strange':
            unsigned_quark_flavor = 4
        else:
            unsigned_quark_flavor = 7
        pmts_sde_cst_quarks = ( d_scale, omega, is_pv_mt_sde, lambda_pv, c_scale, \
                               n_grid, m_rnm, mu2, uv_cutoff, iter_rel_tol, \
                               iter_abs_tol, max_iter_sde, unsigned_quark_flavor, \
                               is_linear_ren_cond )
        # Pass the module and running parameters to the compiled module.
        self.sde_cst_quarks = sde_iter_solver.SDE_quark( *pmts_sde_cst_quarks )
    # iterative solver of the SDE
    def iter_solve_sde(self):
        self.sde_cst_quarks.iter_solve_sde( mpp_num_processes )
        self.sde_cst_quarks.save_solution( self.solution_file )
    # Load the Euclidean-space solution from binary file named by self.solution_file
    # in the folder specified by file_dir.
    def load_solution(self, file_dir = None ):
        cw_dir = getcwd()
        if file_dir != None:
            chdir( file_dir )
        self.sde_cst_quarks.load_solution( self.solution_file )
        chdir( cw_dir )
    # Get parameters needed by the compiled module for the quark propagator
    # with complex-valued momentum with fixed-grid quadrature.
    def pmts_sf_compiled(self):
        self.u_grid = np.array( self.sde_cst_quarks.u_grid() )
        self.a_grid = np.array( self.sde_cst_quarks.a_grid() )
        self.b_grid = np.array( self.sde_cst_quarks.b_grid() )
        self.z_2 = self.sde_cst_quarks.z_2()
        self.z_m = self.sde_cst_quarks.z_m()
        pmts = ( self.omega, self.d_scale, self.n_grid, self.c_scale, \
                self.u_grid, self.a_grid, self.b_grid, self.uv_cutoff, \
                self.z_2, self.z_m * self.m_rnm, is_linear_ren_cond, \
                sde_n_grid_gc, is_pv_mt_sde, lambda_pv )
        return pmts

'''
Class with functions operating the compiled module of the fixed-grid quadrature
    computation of the quark propagator, with support for the iterative solver
    of the SDE in the Euclidean space.
m_rnm_light: renormalized mass for light quarks in GeV.
m_rnm_strange: renormalized mass for strange quarks in GeV.
d_scale: IR strength of the Maris-Tandy model.
omega: IR scale of the Maris-Tandy model.
c_scale_sde: variable transformation scale of the radial momentum.
n_grid: number of grid points for the solution.
mu2: renormalization scale in GeV^2.
uv_cutoff_sde: UV cutoff in GeV^2.
iter_rel_tol iter_abs_tol: relative and absolute tolerances of the iterative solver
quark_flavor = 'light': flavor label of the quark propagator solved, 'light'
    for light quarks and 'strange' for strange quarks.
'''
class SDE_cst_quarks_cpl:
    def __init__(self, m_rnm_light, m_rnm_strange, d_scale, omega, \
                 c_scale_sde, n_grid, mu2, uv_cutoff_sde, iter_rel_tol, iter_abs_tol ):
        self.init_sde_cst_quarks( m_rnm_light, m_rnm_strange, d_scale, omega, \
                                 c_scale_sde, n_grid, mu2, uv_cutoff_sde, \
                                 iter_rel_tol, iter_abs_tol )
    # constructor for class SDE_cst_quarks_cpl
    def init_sde_cst_quarks(self, m_rnm_light, m_rnm_strange, d_scale, omega, \
                            c_scale_sde, n_grid, mu2, uv_cutoff_sde, \
                            iter_rel_tol, iter_abs_tol ):
        self.m_rnm_light = m_rnm_light
        self.m_rnm_strange = m_rnm_strange
        self.d_scale = d_scale
        self.omega = omega
        self.c_scale_sde = c_scale_sde
        self.n_grid = n_grid
        self.mu2 = mu2
        self.uv_cutoff_sde = uv_cutoff_sde
        self.iter_rel_tol = iter_rel_tol
        self.iter_abs_tol = iter_abs_tol
    # Initialize the Euclidean-space solutions, either by the iterative solver
    # ( is_load_sde == False ) or load from binary files ( is_load_sde == True ).
    def iter_solve_sde(self, is_load_sde ):
        self.sde_light = SDE_quark_cpl( self.d_scale, self.omega, self.c_scale_sde, \
                                       self.n_grid, self.m_rnm_light, self.mu2, \
                                       self.uv_cutoff_sde, self.iter_rel_tol, \
                                       self.iter_abs_tol, quark_flavor = 'light', \
                                       solution_file = "sde_solution_light_cpl.bin" )
        if is_load_sde:
            self.sde_light.load_solution( '..' )# change to the parent directory
        else:
            self.sde_light.iter_solve_sde()
        if not is_equal_cst:
            self.sde_strange = SDE_quark_cpl( self.d_scale, self.omega, self.c_scale_sde, \
                                             self.n_grid, self.m_rnm_light, self.mu2, \
                                             self.uv_cutoff_sde, self.iter_rel_tol, \
                                             self.iter_abs_tol, quark_flavor = 'strange', \
                                             solution_file = "sde_solution_strange_cpl.bin" )
            if is_load_sde:
                self.sde_strange.load_solution( '..' )
            else:
                self.sde_strange.iter_solve_sde()
    # Generate members based on the compiled module of the quark propagator with
    # complex-valued momentum.
    def init_prop_pm(self):
        pmts_light = self.sde_light.pmts_sf_compiled()
        self.prop_p = ker_Complex_prop( *pmts_light )
        if is_equal_cst:
            self.prop_m = self.prop_p
            return pmts_light
        else:
            pmts_strange = self.sde_strange.pmts_sf_compiled()
            self.prop_m = ker_Complex_prop( *pmts_strange )
            return ( pmts_light, pmts_strange )
'''
The function that computes the reference value of the self-energy using nquad(),
    not accurate for small modulus or complex-valued momenta.
p2: momentum square of the quark in GeV^2
is_vector: True returns Dirac vector component. False returns Dirac scalar component.
ker_fun( p2, z, y, is_vector ): kernel function of the integrals.
uv_cutoff: UV cutoff in GeV^2 for radial integrals.
'''
def sde_ker_itg( p2, is_vector, ker_fun, uv_cutoff ):
    func_r = lambda y, z: np.real( ker_fun( p2, z, y, is_vector ) )
    func_i = lambda y, z: np.imag( ker_fun( p2, z, y, is_vector ) )
    opts_y = { 'limit': 200 }
    opts_z = { 'limit': 100 }
    opts = [ opts_y, opts_z ]
    quad_r, er = integrate.nquad( func_r, [ [ 0.0, uv_cutoff ], [ 0.0, 1.0 ] ], opts = opts )
    quad_i, ei = integrate.nquad( func_i, [ [ 0.0, uv_cutoff ], [ 0.0, 1.0 ] ], opts = opts )
    print( quad_r, er, quad_i, ei )
    return 2.0 * ( quad_r + 1j * quad_i )

'''
Functions that compute the reference value of the radial integral and the angular
    integral in the self-energy using quad_vec().
p2: momentum square of the quark in GeV^2.
z: angular variable.
y: radial variable.
ker_fun( p2, z, y, is_vector ): kernel function of the integrals.
uv_cutoff: UV cutoff in GeV^2 for radial integrals.
'''
# Because quad_vec is utilized for the integral, both Dirac components are
# evaluated simultaneously.
def vec_ker_map( p2, z, y, ker_fun ):
    fun_v = ker_fun( p2, z, y, True )
    fun_s = ker_fun( p2, z, y, False )
    return np.array( [ np.real( fun_v ), np.imag( fun_v ), np.real( fun_s ), \
                      np.imag( fun_s ) ] )
# Call quad_vec() to compute the radial integrals.
def sde_ker_itg_y( p2, z, ker_fun, uv_cutoff ):
    quad_val, err = integrate.quad_vec( lambda y: vec_ker_map( p2, z, y, \
                                                              ker_fun ), \
                                       0.0, uv_cutoff, \
                                       epsabs = tol_quad_vec, epsrel = tol_quad_vec )
    return quad_val
# Call quad_vec() to compute the angular integrals.
def sde_ker_itg_z( p2, ker_fun, uv_cutoff ):
    quad_val, err = integrate.quad_vec( lambda z: sde_ker_itg_y( p2, z, ker_fun, \
                                                                uv_cutoff ), \
                                       0.0, 1.0, epsabs = tol_quad_vec, \
                                       epsrel = tol_quad_vec )
    quad_val = 2.0 * quad_val
    print( p2, quad_val, 2.0 * err )
    return np.array( [ quad_val[0] + 1j * quad_val[1], \
                      quad_val[2] + 1j * quad_val[3] ] )

'''
Compute the inverse of the quark propagator based on the self-energy.
self_energy_rdd: quark self-energy modulo the constant defined as 1.0 / coe in the function.
is_vector: True returns Dirac vector component. False returns Dirac scalar component.
z_2: wave function renormalization constant.
m_bare: bare mass.
'''
def inv_prop_sde( self_energy_rdd, is_vector, z_2, m_bare ):
    coe = 6.0 * np.pi**3
    if not is_linear_ren_cond:
        self_energy_rdd *= z_2 * z_2
    if is_vector:
        return z_2 + self_energy_rdd / coe
    else:
        return z_2 * m_bare + self_energy_rdd / coe

if __name__ == "__main__":
    ### running parameters ###
    # True if loading previous results for the inverse propagator on the parabola.
    is_load_ab = True
    # True if figures are to be generated.
    is_plot = True
    # True if generating 4 figures, False by default.
    is_4_plot = False
    '''Select if quad_vec() is utilized to compute the reference self-energy,
    otherwise nquad() is used. True by default.'''
    is_quad_vec = True
    # Momentum partition in the BSE, default 0.5.
    eta = 0.5
    # Mass of the bound state that specifies the boundary in the parabola for the
    # accuracy and performance tests.
    mass_bound = 0.13724
    # Number of elements in the geometric series of the loop momentum for the
    # accuracy and performance tests.
    n_kE2_grid = 1152
    ### end running parameters ###
    ### physical paramters ###
    xi = ( eta * mass_bound )**2# Specifies the location and width of the parabola.
    omega = 0.4# IR scale
    d_scale = 0.859# IR strength
    m_rnm_light = 3.6964e-3# renormalized light quark mass
    m_rnm_strange = 0.1# renormalized strange quark mass
    mu2 = 19.0**2# renormalization scale
    uv_cutoff = 1.0e+6# UV cutoff
    n_grid = 2001# number of grid points in the Euclidean-space solution
    c_scale = 10.0# radial variable transforamtion scale
    iter_rel_tol = 1.0e-6# relative tolerance of iterative solver
    iter_abs_tol = 1.0e-6# absolute tolerance of iterative solver
    ### end physical parameters ###
    # parameters for the constructor of SDE_cst_quarks_cpl
    pmts_cpp = ( m_rnm_light, m_rnm_strange, d_scale, omega, \
                c_scale, n_grid, mu2, uv_cutoff, iter_rel_tol, iter_abs_tol )
    sde_quark = SDE_cst_quarks_cpl( *pmts_cpp )
    # Load solution from binary file. False if calling iterative solver.
    sde_quark.iter_solve_sde( True )
    # Get parameters for the constructor of the compiled module for the quark
    # propagator with complex-valued momentum.
    pmts_light = sde_quark.init_prop_pm()
    '''
    pmts_light = ( omega, d_scale, n_grid, c_scale, u_grid, a_grid, b_grid, \
                  uv_cutoff, z_2, z_m * m_rnm, is_linear_ren_cond, sde_n_grid_gc, \
                  is_pv_mt_sde, lambda_pv )
    '''
    # Select parameters needed by the kernel functions.
    args = ( pmts_light[0], pmts_light[1], pmts_light[2], pmts_light[3], \
            pmts_light[4], pmts_light[5], pmts_light[6], pmts_light[7], \
            is_pv_mt_sde, lambda_pv )
    # kernel functions for integrals of quark self-energy with complex momentum
    cp = sde_ker.SDE_SE_cker( *args )
    # momentum grid for the Euclidean-space solution
    u_grid = pmts_light[4]
    p2_SDE = c_scale * u_grid * tfc.vmap( u_grid )
    if is_plot:
        # Plot the Euclidean-space solution.
        fig = plt.figure( figsize = ( 4.0, 3.5 ) )
        plt.semilogx( p2_SDE, pmts_light[5] - 1.0, label = '$A(p^2) - 1.0$' )
        plt.semilogx( p2_SDE, pmts_light[6], '-.', label = '$B(p^2)~\mathrm{GeV}$' )
        plt.legend( loc = 'upper right' )
        plt.xlabel( '$p^2~\mathrm{GeV}^2$', fontsize = font_size )
        plt.ylabel( 'Inverse propagator', fontsize = font_size )
        plt.savefig( 'sf_Euc.eps' , bbox_inches = 'tight', dpi = 300 )
    # Evaluate the parabolic region based on input parameters.
    kE2 = np.logspace( -2, 3, n_kE2_grid )
    p2_grid = np.zeros( np.shape( kE2 ), dtype = complex )
    p2_grid = ( kE2 - xi ) + 1j * 2.0 * np.sqrt( xi * kE2 )
    # Evaluate the inverse propagator using adaptive quadrature.
    p2 = p2_grid
    if is_load_ab:# Load previously saved results.
        temp = np.load('quad_ab.npz')
        p2 = temp['p2']
        a = temp['a']
        b = temp['b']
    else:# Compute the inverse propagator on the parabola.
        t_quad = tfc.Timer('adaptive quadrature')# timer for adptive quadrature
        if mpp_num_processes == 0:# single-thread operation
            a = np.zeros( np.shape( p2 ), dtype = complex )
            b = np.zeros( np.shape( p2 ), dtype = complex )
            if not is_quad_vec:# applies nquad()
                for ind in range( np.size( p2 ) ):
                    a[ind] = sde_ker_itg( cp.self_energy_ker, p2[ind], True, uv_cutoff )
                    b[ind] = sde_ker_itg( cp.self_energy_ker, p2[ind], False, uv_cutoff )
            else:# applies quad_vec()
                for ind in range( np.size( p2 ) ):
                    ker_itg = sde_ker_itg_z( p2[ind], cp.self_energy_ker, uv_cutoff )
                    a[ind] = ker_itg[0]
                    b[ind] = ker_itg[1]
        else:# parallel operations
            if not is_quad_vec:# applies nquad()
                args = itp( p2, [True, False] )
                chuncksize = tfc.mpp_chuncksize( ( np.size( p2 ), 2 ), mpp_num_processes )
                ker_itg = partial( sde_ker_itg, \
                                  ker_fun = cp.self_energy_ker, \
                                  uv_cutoff = uv_cutoff )
                with Pool( mpp_num_processes ) as mpp:
                    results = mpp.starmap( ker_itg, args, chuncksize )
                results = np.array( results )
                results = np.reshape( results, ( np.size( p2 ), 2 ) )
            else:# applies quad_vec()
                chuncksize = tfc.mpp_chuncksize( ( np.size( p2 ), ), mpp_num_processes )
                ker_itg = partial( sde_ker_itg_z,\
                                  ker_fun = cp.self_energy_ker, \
                                  uv_cutoff = uv_cutoff )
                with Pool( mpp_num_processes ) as mpp:
                    results = mpp.map( ker_itg, p2, chuncksize )
                results = np.array( results )
                results = np.reshape( results, ( np.size( p2 ), 2 ) )
            a = results[:,0]
            b = results[:,1]
        t_quad.toc()# Display the time of the computation.
        # Compute the inverse the the propagator.
        a = inv_prop_sde( a, True, pmts_light[8], pmts_light[9] )
        b = inv_prop_sde( b, False, pmts_light[8], pmts_light[9] )
    # Evaluate the inverse propagator using fixed-grid quadrature.
    if is_load_ab:# Load previously saved results.
        temp = np.load('cpp_ab.npz')
        p2_grid = temp['p2_grid']
        a_grid = temp['a_grid']
        b_grid = temp['b_grid']
    else:# Compute the inverse propagator on the parabola.
        t_cpp = tfc.Timer('cpp fixed grid')# timer for fixed-grid quadrature
        if mpp_num_processes == 0:
            a_grid = np.zeros( np.shape( p2 ), dtype = complex )
            b_grid = np.zeros( np.shape( p2 ), dtype = complex )
            for ind in range( np.size( p2 ) ):
                a_grid[ind] = sde_quark.prop_p.inv_prop_ab( p2[ind], True )
                b_grid[ind] = sde_quark.prop_p.inv_prop_ab( p2[ind], False )
        else:
            args = itp( p2, [True, False] )# arguments for starmap()
            chuncksize = tfc.mpp_chuncksize( ( np.size( p2 ), 2 ), mpp_num_processes )
            with Pool( mpp_num_processes ) as mpp:
                results = mpp.starmap( sde_quark.prop_p.inv_prop_ab, args, \
                                      chuncksize )
                results = np.array( results )
                results = np.reshape( results, ( np.size( p2 ), 2 ) )
            a_grid = results[:,0]
            b_grid = results[:,1]
        t_cpp.toc()# Display the time of the computation.
    # Save numberical results for the computed inverse propagator on the parabola.
    if not is_load_ab:
        np.savez( 'quad_ab.npz', p2 = p2, a = a, b = b )
        np.savez( 'cpp_ab.npz', p2_grid = p2_grid, a_grid = a_grid, b_grid = b_grid )
    # Generate figures for comparistion.
    rel_err_a = a_grid / a - 1.0# relative difference in A(p^2)
    rel_err_b = b_grid / b - 1.0# relative difference in B(p^2)
    print( 'Maximum relative error of A function:', 100.0 * np.max( np.abs( rel_err_a ) ), '%' )
    print( 'Maximum relative error of B function:', 100.0 * np.max( np.abs( rel_err_b ) ), '%' )
    if is_plot and is_4_plot:
        fig = plt.figure( figsize = ( 7.75, 6.5 ) )
        ax0 = fig.add_subplot( 2, 2, 1 )
        ax0.semilogx( kE2, np.abs( a_grid ), label = 'Fixed grid' )
        ax0.semilogx( kE2, np.abs( a ), '+', label = 'Adaptive grid', markersize = 10 )
        ax0.set_xlabel( '$k^2~\mathrm{GeV}^2$', fontsize = font_size )
        ax0.set_ylabel( '$|A(p^2)|$', fontsize = font_size )
        ax0.legend( loc = 'upper right' )
        ax1 = fig.add_subplot( 2, 2, 2 )
        ax1.semilogx( kE2, np.abs( b_grid ), label = 'Fixed grid' )
        ax1.semilogx( kE2, np.abs( b ), '+', label = 'Adaptive grid', markersize = 10 )
        ax1.set_xlabel( '$k^2~\mathrm{GeV}^2$', fontsize = font_size )
        ax1.set_ylabel( '$|B(p^2)|~\mathrm{GeV}$', fontsize = font_size )
        ax1.legend( loc = 'upper right' )
        ax2 = fig.add_subplot( 2, 2 ,3 )
        ax2.semilogx( kE2, np.abs( rel_err_a ), '-x' )
        ax2.set_xlabel( '$k^2~\mathrm{GeV}^2$', fontsize = font_size )
        ax2.set_ylabel( 'Relative difference in $A(p^2)$', fontsize = font_size )
        ax2.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
        ax3 = fig.add_subplot( 2, 2 ,4 )
        ax3.semilogx( kE2, np.abs( rel_err_b ), '-x' )
        ax3.set_xlabel( '$k^2~\mathrm{GeV}^2$', fontsize = font_size )
        ax3.set_ylabel( 'Relative difference in $B(p^2)$', fontsize = font_size )
        ax3.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
        plt.savefig( 'test_sf_compiled.eps', dpi = 300, bbox_inches = 'tight' )
    elif is_plot:
        fig = plt.figure( figsize = ( 4.0, 10.5 ) )
        ax0 = fig.add_subplot( 3, 1, 1 )
        ind = np.arange( 0, n_kE2_grid, n_kE2_grid / 48, dtype = int )
        if ind[-1] != n_kE2_grid - 1:
            ind = np.append( ind, np.array( [ n_kE2_grid - 1 ] ) )
        ax0.semilogx( kE2, np.real( a_grid - 1.0 ), label = '$A(p^2)-1.0$ fixed grid' )
        ax0.semilogx( kE2, np.real( b_grid ), '-.', label = '$B(p^2)~\mathrm{GeV}$ fixed grid' )
        ax0.semilogx( kE2[ind], np.real( a[ind] - 1.0 ), '+', \
                     label = '$A(p^2)-1.0$ adaptive', markersize = 10 )
        ax0.semilogx( kE2[ind], np.real( b[ind] ), 'x', \
                     label = '$B(p^2)~\mathrm{GeV}$ adaptive', markersize = 10 )
        ax0.set_xlabel( '$k^2~\mathrm{GeV}^2$', fontsize = font_size )
        ax0.set_ylabel( 'Real part', fontsize = font_size )
        ax0.legend( loc = 'upper right' , frameon = False )
        ax2 = fig.add_subplot( 3, 1, 2 )
        ax2.semilogx( kE2, np.imag( a_grid - 1.0 ), label = '$A(p^2)$ fixed grid' )
        ax2.semilogx( kE2, np.imag( b_grid ), '-.', label = '$B(p^2)~\mathrm{GeV}$ fixed grid' )
        ax2.semilogx( kE2[ind], np.imag( a[ind] - 1.0 ), '+', \
                     label = '$A(p^2)$ adaptive', markersize = 10 )
        ax2.semilogx( kE2[ind], np.imag( b[ind] ), 'x', \
                     label = '$B(p^2)~\mathrm{GeV}$ adaptive', markersize = 10 )
        ax2.set_xlabel( '$k^2~\mathrm{GeV}^2$', fontsize = font_size )
        ax2.set_ylabel( 'Imaginary part', fontsize = font_size )
        ax2.legend( loc = 'lower right' , frameon = False )
        ax1 = fig.add_subplot( 3, 1 ,3 )
        ax1.loglog( kE2[ind], np.abs( rel_err_a[ind] ), '-+', label = '$A(p^2)$' )
        ax1.loglog( kE2[ind], np.abs( rel_err_b[ind] ), '-.x', label = '$B(p^2)$' )
        ax1.set_xlabel( '$k^2~\mathrm{GeV}^2$', fontsize = font_size )
        ax1.set_ylabel( 'Relative difference', fontsize = font_size )
        # ax1.ticklabel_format( axis='y', style='sci', scilimits=(0,0) ) # optional label format
        ax1.legend( loc = 'upper left' )
        plt.savefig( 'test_sf_compiled.eps', dpi = 300, bbox_inches = 'tight' )
