/*
Definition of the class Gfun for the dressing functions of the gluon
propagator.

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
using namespace std;
class Gfun// Gluon propagator class applying the Maris-Tandy model
{
public:
	/* full constructor
	omega: IR scale of the Maris-Tandy model in GeV.
	d_num: IR strength of the Maris-Tandy model in GeV^2.
	is_pv_mt: If the Pauli-Villars (PV) regularization is applied.
	lambda_pv: mass of the PV regulartor in GeV.
	n_f: number of active flavors.
	mt: variabl m_t in GeV.
	lQCD: QCD IR scale \Lambda_{\mathrm{QCD}} in GeV. */
	Gfun( double omega, double d_num, bool is_pv_mt, double lambda_pv,
		unsigned n_f = 4, double mt = 0.5, double lQCD = LAMBDA_QCD ):
		d_scale( d_num ), is_pv_mt( is_pv_mt )
	{
		omega2 = pow( omega, 2 );// \omega^2
		coenf = coe_nf( n_f );// coefficient related to n_f
		mt2 = pow( mt, 2 );// m_t^2
		lQCD2 = pow( lQCD, 2 );// \Lambda^2_{\mathrm{QCD}}
		lambda_pv_2 = pow( lambda_pv, 2 );// \Lambda^2_{\mathrm{PV}}
	}
	// Reduced constructor, applies in conjunction with set_gfun_pmts( bool, double ).
	Gfun( double omega, double d_num ): d_scale( d_num )
	{
		unsigned n_f = 4;
		double mt = 0.5;
		double lQCD = LAMBDA_QCD;
		omega2 = pow( omega, 2 );
		coenf = coe_nf( n_f );
		mt2 = pow( mt, 2 );
		lQCD2 = pow( lQCD, 2 );
	}
	// Default constructor, applies in conjunction with set_gfun_pmts( double, double bool, double )
	Gfun( void )
	{
		unsigned n_f = 4;
		double mt = 0.5;
		double lQCD = LAMBDA_QCD;
		coenf = coe_nf( n_f );
		mt2 = pow( mt, 2 );
		lQCD2 = pow( lQCD, 2 );
	}
	void set_gfun_pmts( double omega, double d_num, bool is_pv_mt_in, double lambda_pv )
	{
		d_scale = d_num;
		omega2 = pow( omega, 2 );
		set_gfun_pmts( is_pv_mt_in, lambda_pv );
	}
	void set_gfun_pmts( bool is_pv_mt_in, double lambda_pv )
	{
		is_pv_mt = is_pv_mt_in;
		lambda_pv_2 = pow( lambda_pv, 2 );
	}
	double val_omega( void ){ return sqrt( omega2 ); }// Returns IR scale \omega^2.
  double val_d_scale( void ){ return d_scale; }// Returns IR strength.
	bool val_is_pv_mt( void ){ return is_pv_mt; }// Returns if PV regulator is used.
	double val_lambda_pv( void ){ return sqrt( lambda_pv_2 ); }// mass of the PV regular in GeV
	// g_IR( k2 ) the IR term of the Maris-Tandy model without PV regularization.
	template<typename T>
	T g_mt_ir( T k2 )
	{
		T fun_return = k2 * g_mt_ir_q2( k2 );
		if( is_pv_mt ){ fun_return *= profile_pv( k2 ); }
		return fun_return;
	}
	// default gluon dressing function in the Maris-Tandy model
	template<typename T>
	T gfun_mt( T k2 )
	{
		T fun_return = g_mt_ir_q2( k2 ) * k2 + g_mt_uv( k2 );
		if( is_pv_mt ){ fun_return *= profile_pv( k2 ); }
		return fun_return;
	}
	// gluon dressing function with the option to exclude the IR term
	template<typename T>
	T gfun_mt( T k2, bool is_include_IR )
	{
		T fun_return = g_mt_uv( k2 );
		if( is_include_IR ){ fun_return += g_mt_ir_q2( k2 ) * k2; }
		if( is_pv_mt ){ fun_return *= profile_pv( k2 ); }
		return fun_return;
	}
	// gluon dressing function devided by gluon momentum k^2
	template<typename T>
	T gfun_mt_q2( T k2 )
	{
		T fun_return = g_mt_ir_q2( k2 ) + g_mt_uv( k2 ) / k2;
		if( is_pv_mt ){ fun_return *= profile_pv( k2 ); }
		return fun_return;
	}
	// gluon dressing function devided by gluon momentum k^2 with the option to \
	exclude the IR term
	template<typename T>
	T gfun_mt_q2( T k2, bool is_include_IR )
	{
		T fun_return = g_mt_uv( k2 ) / k2;
		if( is_include_IR ){ fun_return += g_mt_ir_q2( k2 ); }
		if( is_pv_mt ){ fun_return *= profile_pv( k2 ); }
		return fun_return;
	}
private:
	double omega2, d_scale, coenf, mt2, lQCD2;
	bool is_pv_mt;
	double lambda_pv_2;
	// The option to exclude UV term is controlled by IS_PV_MT in consts.h.
	#if IS_UV_MT
		double coe_nf( int nf )
		{	return 48.0 * pow( PI, 2 ) / double( 33 - 2 * nf );	}
	#else
		double coe_nf( int nf ){ return 0.0; }
	#endif
	// PV regulator profile
	template<typename T>
	T profile_pv( T y )
  { return 1.0 / ( y / lambda_pv_2 + 1.0 ); }
	//  ( 1.0 / k2 ) * g_IR( k2 ) as the IR term of the gluon propagator devided \
	by gluon momentum k^2
	template<typename T>
	T g_mt_ir_q2( T k2 )
	{
		T mt_IR = k2 / omega2;
		mt_IR = exp( - mt_IR ) / omega2;
		mt_IR *= 4.0 * pow( PI / omega2, 2 ) * d_scale;
		return mt_IR;
	}
	// g_UV( k2 ) as the UV term of the gluon propagator
	template<typename T>
	T g_mt_uv( T k2 )
	{
		T term_N = 1.0 / ( 4.0 * mt2 );
		if( abs( k2 ) > EPS * mt2 )
		{	term_N *= - k2;
			term_N = 1.0 - exp( term_N );
			term_N /= k2; }
		T term_D = pow( 1.0 + k2 / lQCD2, 2 );
		term_D += pow( EBN, 2 ) - 1.0;
		term_D = 0.5 * log( term_D );
		return coenf * term_N / term_D;
	}
};
