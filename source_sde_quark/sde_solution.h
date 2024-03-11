/* Definitions of struct sde_solution, struct sde_solution_array, and
  class SDE_Solution_IO used in SDE_quark.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3). */
// struct sde_solution; struct sde_solution_array<unsigned num_grid> \
class SDE_Solution_IO ( sde_solution ) \
#include<fstream> #include<string>
using namespace std;
// struct for the Euclidean-space solution of the SDE
struct sde_solution
{
  unsigned n;// number of grid points
  double c;// variable transformation scale in GeV^2
  vector<double> u;// grid for the mapped variable
  vector<double> a;// Dirac vector inverse propagator
  vector<double> b;// Dirac scalar inverse propagator
  double mu2;// renormalization scale in GeV^2
  // m_bare occures in the constructor of Complex_prop<T>.
  double m_rnm;// renormalized mass in GeV. m_bare = m_rnm * z_m
  double uv_cutoff;// UV cutoff of the radial integral in GeV^2
  double z_2;// wave function renormalization constant
  double z_m;// z_4 = z_2 * z_m, mass renormalization constant.
  bool is_linear_ren_cond;// If linear renormalization condition applies.
  bool is_pv_mt;// If Pauli-Vallars regularization applies.
  double lambda_pv;// mass of the Pauli-Villars regulartor
};
// similar to struct sde_solution with u, a, b in array<double,num_grid> \
instead of vector<double>
template<unsigned num_grid>
struct sde_solution_array
{
  unsigned n = num_grid;
  double c;
  array<double, num_grid> u;
  array<double, num_grid> a;
  array<double, num_grid> b;
  double mu2;
  // m_bare occures in the constructor of Complex_prop<T>.
  double m_rnm;// m_bare = m_rnm * z_m
  double uv_cutoff;
  double z_2;
  double z_m;// z_4 = z_2 * z_m
  bool is_linear_ren_cond;
  bool is_pv_mt;
  double lambda_pv;
};

// Class for the input and output of the Euclidean-space solution of the \
SDE for the quark propagator.
class SDE_Solution_IO
{
public:
  // default constructor
  SDE_Solution_IO( void ){}
  /* Full constructor
  n: number of grid points.
  c: variable transformation scale in GeV^2.
  double* u_pointer: address for the grid for the mapped variable.
  double* a_pointer: address for the Dirac vector inverse propagator.
  double* b_pointer: address for Dirac scalar inverse propagator.
  mu2: renormalization scale in GeV^2.
  m_rnm: renormalized mass in GeV
  uv_cutoff: UV cutoff of the radial integral in GeV^2.
  z_2: wave function renormalization constant.
  z_m: mass renormalization constant.
  is_linear_ren_cond: If linear renormalization condition applies.
  is_pv_mt: If Pauli-Vallars regularization applies.
  lambda_pv: Mass of the Pauli-Villars regulartor.
  quark_flavor_in: flavor lable for the quark. */
  SDE_Solution_IO( unsigned n, unsigned c, double* u_pointer, double* a_pointer,
    double* b_pointer, double mu2, double m_rnm, double uv_cutoff, double z_2,
    double z_m, bool is_linear_ren_cond, bool is_pv_mt, double lambda_pv,
    unsigned quark_flavor_in ): quark_flavor( quark_flavor_in ), n( n ), c( c ),
  mu2( mu2 ), m_rnm( m_rnm ), uv_cutoff( uv_cutoff ), z_2( z_2 ), z_m( z_m ),
  is_linear_ren_cond( is_linear_ren_cond ), is_pv_mt( is_pv_mt ),
  lambda_pv( lambda_pv )
  {
    allocate_uab( n );
    for( unsigned ind = 0; ind < n; ind++ )
    {
      u[ind] = u_pointer[ind];
      a[ind] = a_pointer[ind];
      b[ind] = b_pointer[ind];
    }
  }
  // full constructor using struct sde_solution
  SDE_Solution_IO( sde_solution sf_in, unsigned quark_flavor_in ):
  quark_flavor( quark_flavor_in ), n( sf_in.n ), c( sf_in.c ), mu2( sf_in.mu2 ),
  m_rnm( sf_in.m_rnm ), uv_cutoff( sf_in.uv_cutoff ), z_2( sf_in.z_2 ),
  z_m( sf_in.z_m ), is_linear_ren_cond( sf_in.is_linear_ren_cond ),
  is_pv_mt( sf_in.is_pv_mt ), lambda_pv( sf_in.lambda_pv )
  {
    allocate_uab( n );
    for( unsigned ind = 0; ind < n; ind++ )
    {
      u[ind] = sf_in.u[ind];
      a[ind] = sf_in.a[ind];
      b[ind] = sf_in.b[ind];
    }
  }
  // full constructor by loading from binary file
  SDE_Solution_IO( string filename, bool is_print_filename = false )
  { load_sf_bin( filename, true, is_print_filename ); }
  // Save the instance of class to binary file.
  void save_sf_bin( string filename )
  {
    ofstream save_file( filename, ios::binary );
    if( !save_file.is_open() )
    { cout<<"Failed to open save file."<<endl; }
    else
    {
      char* save_data_pointer;
      save_file.seekp(0);
      save_data_pointer = reinterpret_cast<char*>( &quark_flavor );
      save_file.write( save_data_pointer, sizeof quark_flavor );
      save_data_pointer = reinterpret_cast<char*>( &n );
      save_file.write( save_data_pointer, sizeof n );
      save_data_pointer = reinterpret_cast<char*>( &c );
      save_file.write( save_data_pointer, sizeof c );
      //
      save_data_pointer = reinterpret_cast<char*>( u.data() );
      save_file.write( save_data_pointer, streamsize( sizeof( double ) * n ) );
      save_data_pointer = reinterpret_cast<char*>( a.data() );
      save_file.write( save_data_pointer, streamsize( sizeof( double ) * n ) );
      save_data_pointer = reinterpret_cast<char*>( b.data() );
      save_file.write( save_data_pointer, streamsize( sizeof( double ) * n ) );
      //
      save_data_pointer = reinterpret_cast<char*>( &mu2 );
      save_file.write( save_data_pointer, sizeof mu2 );
      save_data_pointer = reinterpret_cast<char*>( &m_rnm );
      save_file.write( save_data_pointer, sizeof m_rnm );
      save_data_pointer = reinterpret_cast<char*>( &uv_cutoff );
      save_file.write( save_data_pointer, sizeof uv_cutoff );
      save_data_pointer = reinterpret_cast<char*>( &z_2 );
      save_file.write( save_data_pointer, sizeof z_2 );
      save_data_pointer = reinterpret_cast<char*>( &z_m );
      save_file.write( save_data_pointer, sizeof z_m );
      save_data_pointer = reinterpret_cast<char*>( &is_linear_ren_cond );
      save_file.write( save_data_pointer, sizeof is_linear_ren_cond );
      save_data_pointer = reinterpret_cast<char*>( &is_pv_mt );
      save_file.write( save_data_pointer, sizeof is_pv_mt );
      save_data_pointer = reinterpret_cast<char*>( &lambda_pv);
      save_file.write( save_data_pointer, sizeof lambda_pv );
      save_file.flush();
      cout<<"Written sde solution to "<<filename<<"."<<endl;
    }
    save_file.close();
  }
  // Load members from a binary file.
  void load_sf_bin( string filename, bool is_constructor = false,
    bool is_print_filename = false )
  {
    char* load_data_pointer;
    ifstream load_file( filename, ios::binary );
    if( !load_file.is_open() )
    { cout<<"Failed to open load file."<<endl; }
    else
    {
      load_file.seekg(0);
      load_data_pointer = reinterpret_cast<char*>( &quark_flavor );
      load_file.read( load_data_pointer, sizeof quark_flavor );
      load_data_pointer = reinterpret_cast<char*>( &n );
      load_file.read( load_data_pointer, sizeof n );
      load_data_pointer = reinterpret_cast<char*>( &c );
      load_file.read( load_data_pointer, sizeof c );
      //
      allocate_uab( n );
      load_data_pointer = reinterpret_cast<char*>( u.data() );
      load_file.read( load_data_pointer, streamsize( sizeof( double ) * n ) );
      load_data_pointer = reinterpret_cast<char*>( a.data() );
      load_file.read( load_data_pointer, streamsize( sizeof( double ) * n ) );
      load_data_pointer = reinterpret_cast<char*>( b.data() );
      load_file.read( load_data_pointer, streamsize( sizeof( double ) * n ) );
      //
      load_data_pointer = reinterpret_cast<char*>( &mu2 );
      load_file.read( load_data_pointer, sizeof mu2 );
      load_data_pointer = reinterpret_cast<char*>( &m_rnm );
      load_file.read( load_data_pointer, sizeof m_rnm );
      load_data_pointer = reinterpret_cast<char*>( &uv_cutoff );
      load_file.read( load_data_pointer, sizeof uv_cutoff );
      load_data_pointer = reinterpret_cast<char*>( &z_2 );
      load_file.read( load_data_pointer, sizeof z_2 );
      load_data_pointer = reinterpret_cast<char*>( &z_m );
      load_file.read( load_data_pointer, sizeof z_m );
      load_data_pointer = reinterpret_cast<char*>( &is_linear_ren_cond );
      load_file.read( load_data_pointer, sizeof is_linear_ren_cond );
      load_data_pointer = reinterpret_cast<char*>( &is_pv_mt );
      load_file.read( load_data_pointer, sizeof is_pv_mt );
      load_data_pointer = reinterpret_cast<char*>( &lambda_pv);
      load_file.read( load_data_pointer, sizeof lambda_pv );
      if( is_print_filename )
      { cout<<"Read sde solution from "<<filename<<"."<<endl; }
    }
    load_file.close();
  }
  // Modify the flavor lable.
  void set_quark_flavor( unsigned quark_flavor_in )
  { quark_flavor = quark_flavor_in; }
  // Display the quark flavor.
  void print_quark_flavor( void )
  {
    switch( quark_flavor )
    {
      case 0:
        cout<<"light"<<endl; break;
      case 1:
        cout<<"up"<<endl; break;
      case 2:
        cout<<"down"<<endl; break;
      case 3:
        cout<<"charm"<<endl; break;
      case 4:
        cout<<"strange"<<endl; break;
      case 5:
        cout<<"top"<<endl; break;
      case 6:
        cout<<"bottom"<<endl; break;
      default:
        cout<<"unspecified"<<endl; break;
    }
  }
  // Display the soluiton.
  void print_solution( void )
  {
    cout<<"SDE solution for the quark of flavor ";
    print_quark_flavor();
    cout<<"number of grid point = "<<n<<endl;
    cout<<"c_scale = "<<c<<endl;
    for( unsigned ind = 0; ind < n; ind ++ )
    { cout<<"("<<u[ind]<<","<<a[ind]<<","<<b[ind]<<") "; }
    cout<<endl;
    cout<<"mu2 = "<<mu2<<endl;
    cout<<"m_rnm = "<<m_rnm<<endl;// m_bare = m_rnm * z_m
    cout<<"uv_cutoff = "<<uv_cutoff<<endl;
    cout<<"z_2 = "<<z_2<<endl;
    cout<<"z_m = "<<z_m<<endl;// z_4 = z_2 * z_m
    cout<<"is_linear_ren_cond = "<<is_linear_ren_cond<<endl;
    cout<<"is_pv_mt = "<<is_pv_mt<<endl;
    cout<<"lambda_pv = "<<lambda_pv<<endl;
  }
  // Generate an instance of sde_solution based on members.
  sde_solution get_sde_solution( void )
  {
    sde_solution sf;
    sf.n = n;
    sf.c = c;
    sf.u.resize( n );
    sf.a.resize( n );
    sf.b.resize( n );
    for( unsigned ind = 0; ind < n; ind++ )
    {
      sf.u[ind] = u[ind];
      sf.a[ind] = a[ind];
      sf.b[ind] = b[ind];
    }
    sf.mu2 = mu2;
    sf.m_rnm = m_rnm;
    sf.uv_cutoff = uv_cutoff;
    sf.z_2 = z_2;
    sf.z_m = z_m;
    sf.is_linear_ren_cond = is_linear_ren_cond;
    sf.is_pv_mt = is_pv_mt;
    sf.lambda_pv = lambda_pv;
    return sf;
  }
  // functions to get model parameters
  bool get_is_pv_mt( void ){ return is_pv_mt; }
  double get_lambda_pv( void ){ return lambda_pv; }
  double get_uv_cutoff( void ){ return uv_cutoff; }
  double get_a_last( void ){ return a[n-1]; }
  double get_b_last( void ){ return b[n-1]; }
  double get_mu2( void ){ return mu2; }
  double get_z_2( void ){ return z_2; }
  double get_z_m( void ){ return z_m; }
  double get_z_4( void ){ return z_2 * z_m; }
  double get_m_bare( void ){ return z_m * m_rnm; }
  bool get_is_linear_ren_cond( void ){ return is_linear_ren_cond; }
private:
  unsigned quark_flavor, n = 0;
  double c, mu2;
  vector<double> u, a, b;
  double m_rnm;// m_bare = m_rnm * z_m
  double uv_cutoff, z_2, z_m;// z_4 = z_2 * z_m
  bool is_linear_ren_cond, is_pv_mt;
  double lambda_pv;
  void allocate_uab( unsigned num_grid )
  { u.resize( num_grid ); a.resize( num_grid ); b.resize( num_grid ); }
};
