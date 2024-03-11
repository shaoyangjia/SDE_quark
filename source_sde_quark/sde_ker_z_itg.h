/*
Definition of the class Zker that applies adaptive Gauss-Kronrod quadrature
to compute angular integral in quark self-energy. The template parameter T
specifies the shared type of the momentum variable and of the self-energy.

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
// class Zker<typename T>: public SDE_ker<T>, public QUAD_GK<T> \
#include <cmath> #include <complex> \
#include "g_mt.h" #include "quad_gk.h" #include "quad_gc.h"
using namespace std;
template<typename T>
class Zker: public SDE_ker<T>, public QUAD_GK<T>
{
public:
	/* constructor
	omega: IR scale of the Maris-Tandy model in GeV.
	d_scale: IS strength of the Maris-Tandy model in GeV^2.
	is_pv_mt: if Pauli-Villars regularization is applied.
	lambda_pv: mass of the Pauli-Villars regulartor in GeV. */
	Zker( double omega, double d_scale, bool is_pv_mt, double lambda_pv ):
	SDE_ker<T>( omega, d_scale, is_pv_mt, lambda_pv ) {}
	// functions to get model parameters
	double get_omega( void ){ return Gfun::val_omega(); }
  double get_d_scale( void ){ return Gfun::val_d_scale(); }
	bool get_is_pv_mt( void ){ return SDE_ker<T>::Gfun::val_is_pv_mt(); }
	bool get_lambda_pv( void ){ return SDE_ker<T>::Gfun::val_lambda_pv();	}
	/* Function to compute the angular integral for a given external momentum
	x_num and internal loop momentum y_num. is_vector_bool: true returns the
	Dirac vector component, false for the Dirac scalar component.	*/
	T cz_itg( T x_num, double y_num, bool is_vector_bool )
	{
		x = x_num;
		y = y_num;
		bool_is_vector = is_vector_bool;
		QUAD_GK<T>::reset();
	  return 2.0 * QUAD_GK<T>::quad_gk( 0.0, 1.0 );
	}
	// Override the integrand for the adaptive Gauss-Kronrod quadrature.
	T ker( double z ) override
	{ return SDE_ker<T>::sde_ker( x, y, z, bool_is_vector ); }
private:
	T x;// external momentum
	double y;// radial loop momentum
	bool bool_is_vector;// return component selection
};
