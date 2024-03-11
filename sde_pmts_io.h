/* "sde_pmts_io.h" defines the class SDE_PMTS_IO used in SDE_quark.

Copyright (c) 2024 UChicago Argonne, LLC
All Rights Reserved
By: Argonne National Laboratory

This file is part of SDE_quark released under the GNU General Public License
version 3 (GPLv3). */
using namespace std;
class SDE_PMTS_IO
{
public:
  // default constructor
  SDE_PMTS_IO( void ) = default;
  // constructor that loads pmts from binary file with name by string read_filename
  SDE_PMTS_IO( string read_filename ){ read_pmts_bin( read_filename, false ); }
  /* full constructor ( double^3, size_t, unsigned, string, bool, double )
  omega: IR scale of the Maris-Tandy model GeV.
  d_scale: IR strength of the Maris-Tandy model GeV^2.
  c_scale_sde: variable transformation scale of the solution grid in GeV^2.
  n_u_grid: number of grid points for the Euclidean-space solution.
  n_grid_gc: number of grid points for the fixed-grid Gauss-Chebyshev quadrature.
  sf_filename: name of the binary file that saves the solution.
  lambda_pv: mass of the Pauli-Villars regulartor in GeV. */
  SDE_PMTS_IO( double omega, double d_scale, double c_scale_sde, size_t n_u_grid,
    unsigned n_grid_gc, string sf_filename, double lambda_pv ): omega( omega ),
  d_scale( d_scale ), c_scale_sde( c_scale_sde ), n_u_grid( n_u_grid ),
  n_grid_gc( n_grid_gc ), sf_filename( sf_filename ), lambda_pv( lambda_pv ){}
  // public members
  double omega, d_scale, c_scale_sde, lambda_pv;
  size_t n_u_grid;// number of grid points in the iterative solver of the SDE
  unsigned n_grid_gc;// number of grid points in fixed-grid Gauss-Chebyshev \
  quadrature for angular integrals
  string sf_filename;// name for the binary file of the Euclidean-space solution
  // public methods
  array<double,4> get_sde_pmts_double( void )// real-valued parameters
  { return array<double,4> { omega, d_scale, c_scale_sde, lambda_pv }; }
  size_t get_sde_pmts_size_t( void ){ return n_u_grid; }// size_t-valued parameter
  // member function that sets parameters
  void set_sde_pmts( array<double,4> args_double, size_t arg_size_t,
    unsigned n_grid_gc_in, string sf_filename_in )
  {
    set_bse_pmts_double( args_double );
    set_bse_pmts_size_t( arg_size_t );
    n_grid_gc = n_grid_gc_in;
    sf_filename = sf_filename_in;
  }
  // display parameters
  void print_pmts( void )
  {
    cout<<"omega "<<omega<<endl;
    cout<<"d_scale "<<d_scale<<endl;
    cout<<"c_scale_SDE "<<c_scale_sde<<endl;
    cout<<"lambda_pv "<<lambda_pv<<endl;
    cout<<"n_u_grid "<<n_u_grid<<endl;
    cout<<"n_grid_gc "<<n_grid_gc<<endl;
    cout<<"sf_filename "<<sf_filename<<endl;
  }
  // Save parameters in binary file named by filename.
  void write_pmts_bin( string filename )
  {
    ofstream save_file( filename, ios::binary );
    if( !save_file.is_open() )
    { cout<<"Failed to open save file."<<endl; }
    else
    {
      char* data_pointer;
      save_file.seekp(0);
      data_pointer = reinterpret_cast<char*>( &omega );
      save_file.write( data_pointer, sizeof( double ) );
      data_pointer = reinterpret_cast<char*>( &d_scale );
      save_file.write( data_pointer, sizeof( double ) );
      data_pointer = reinterpret_cast<char*>( &c_scale_sde );
      save_file.write( data_pointer, sizeof( double ) );
      data_pointer = reinterpret_cast<char*>( &lambda_pv );
      save_file.write( data_pointer, sizeof( double ) );
      data_pointer = reinterpret_cast<char*>( &n_u_grid );
      save_file.write( data_pointer, sizeof( size_t ) );
      data_pointer = reinterpret_cast<char*>( &n_grid_gc );
      save_file.write( data_pointer, sizeof( unsigned ) );
      size_t sf_filename_size = sf_filename.length();
      data_pointer = reinterpret_cast<char*>( &sf_filename_size );
      save_file.write( data_pointer, sizeof( size_t ) );
      for( size_t ind = 0; ind < sf_filename_size; ind++ )
      { save_file.write( &sf_filename[ind], sizeof( char ) ); }
      save_file.flush();
      cout<<"Written SDE pmts to "<<filename<<"."<<endl;
    }
    save_file.close();
  }
  // Load parameters from binary file named by filename. Display parameters \
  if the second argument is true.
  void read_pmts_bin( string filename, bool is_print_filename = false )
  {
    char* data_pointer;
    ifstream read_file( filename, ios::binary );
    if( !read_file.is_open() )
    { cout<<"Failed to open read file."<<endl; }
    else
    {
      read_file.seekg(0);
      data_pointer = reinterpret_cast<char*>( &omega );
      read_file.read( data_pointer, sizeof( double ) );
      data_pointer = reinterpret_cast<char*>( &d_scale );
      read_file.read( data_pointer, sizeof( double ) );
      data_pointer = reinterpret_cast<char*>( &c_scale_sde );
      read_file.read( data_pointer, sizeof( double ) );
      data_pointer = reinterpret_cast<char*>( &lambda_pv );
      read_file.read( data_pointer, sizeof( double ) );
      data_pointer = reinterpret_cast<char*>( &n_u_grid );
      read_file.read( data_pointer, sizeof( size_t ) );
      data_pointer = reinterpret_cast<char*>( &n_grid_gc );
      read_file.read( data_pointer, sizeof( unsigned ) );
      sf_filename.clear();
      size_t sf_filename_size;
      data_pointer = reinterpret_cast<char*>( &sf_filename_size );
      read_file.read( data_pointer, sizeof( size_t ) );
      for( size_t ind = 0; ind < sf_filename_size; ind++ )
      {
        char sf_char;
        read_file.read( &sf_char, sizeof( char ) );
        sf_filename.push_back( sf_char );
      }
      if( is_print_filename )
      { cout<<"Read BSE pmts from "<<filename<<"."<<endl; }
      read_file.sync();
    }
    read_file.close();
  }
private:
  // set real-valued parameters
  void set_bse_pmts_double( array<double,4> args )
  {
    omega = args[0];
    d_scale = args[1];
    c_scale_sde = args[2];
    lambda_pv = args[3];
  }
  // set the size_t-type parameter
  void set_bse_pmts_size_t( size_t arg ){ n_u_grid = arg; }
};
