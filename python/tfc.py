# -*- coding: utf-8 -*-
"""
"tfc.py" constains supporting functions for the accuracy and performance tests
of the fixed-grid quadrature algorithm in SDE_quark.

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
import time
import numpy as np

# Evaluate the chunck size for multiprocessing function calls.
def mpp_chuncksize( array_shape, num_processes ):
    chuncksize, residue = divmod( np.prod( array_shape ), num_processes )
    if residue != 0:
        chuncksize += 1
    return chuncksize

# real-time timer class
class Timer:
    def __init__(self, message ):
        self.message = message
        # Record the time of the constructor.
        self.t0 = time.perf_counter()
    def toc(self):
        t1 = time.perf_counter()
        # Display the real time counted from the constructor.
        print( self.message, "took", t1 - self.t0, 's' )

# Radial variable transformation modulo the scale, vmap( u ) = 1.0 / ( 1.0 - u ).
def vmap( u ):
    if np.isscalar( u ):
        try:
            return 1.0 / ( 1.0 - u )
        except ZeroDivisionError:
            return np.inf
    else:
        fun_return = np.zeros( np.shape( u ), dtype = u.dtype )
        fun_return[ u != 1.0 ] = 1.0 / ( 1.0 - u[ u != 1.0 ] )
        fun_return[ u == 1.0 ] = np.inf;
        return fun_return

# Inverse of the radial variable transformation, ivmap( x, c ) = x / ( x + c ).
def ivmap( x, c_scale ):
    if np.isscalar( x ):
        if np.isinf( np.abs( x ) ):
            return 1.0
        else:
            return x / ( x + c_scale )
    else:
        fun_return = np.zeros( np.shape( x ), dtype = x.dtype )
        ind = np.logical_not( np.isinf( np.abs( x ) ) )
        fun_return[ ind ] = x[ ind ] / ( x[ ind ] + c_scale )
        ind = np.logical_not( ind )
        fun_return[ ind ] = 1.0
        return fun_return

# Tolerance tests of f1 and f2 with relative ( rel_tol ) and absolute ( abs_tol )
# tolerances, returns True when either of the differences is lower than the
# corersponding tolerances. If abs_tol is zero, only the condition for the
# absolute difference is counted.
def test_tol_cond( f1, f2, rel_tol, abs_tol ):
    if rel_tol == 0.0:
        abs_cond = np.abs( f1 - f2 ) <= abs_tol
        return np.all( abs_cond )
    else:
        if abs_tol == 0.0:
            return test_rel_diff( f1, f2, rel_tol )
        else:
            abs_cond = np.abs( f1 - f2 ) <= abs_tol
            rel_cond = test_rel_diff( f1, f2, rel_tol )
            return np.all( np.logical_or( rel_cond, abs_cond ) )
# Tolerance test of x with respect to y, returns True when the absolute value
# for the difference is lower than the desired tolerance ( rel_tol ). When the
# reference value of y is zero, returns True only for x = 0.
def test_rel_diff( x, y, rel_tol ):
    if np.isscalar( x ):
        if np.abs( x ) == 0:
            if np.abs( y ) == 0:
                return True
            else:
                return False
        else:
            return 2.0 * abs( ( x - y ) / ( x + y ) ) <= rel_tol
    else:
        zero_x = np.abs( x ) == 0
        zero_y = np.abs( y ) == 0
        fun_return = np.empty( np.shape(x), dtype = bool )
        fun_return[ np.logical_and( zero_x, zero_y ) ] = True
        fun_return[ np.logical_and( zero_x, np.logical_not( zero_y) ) ] = False
        nonzero_x = np.logical_not( zero_x )
        fun_return[ nonzero_x ] = 2.0 * abs( \
                                    ( x[ nonzero_x ] - y[ nonzero_x ] ) / \
                                    ( x[ nonzero_x ] + y[ nonzero_x ] ) \
                                    ) <= rel_tol
        return np.all( fun_return )
