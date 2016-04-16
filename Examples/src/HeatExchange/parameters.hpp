#ifndef HH_Parameters_HH
#define HH_Parameters_HH
#include <iosfwd>
#include <string>
struct parameters
{
  //! stationary or non-stationary problem
  int stationary; 
  //! Bar length	
  double L;
  //! First longitudinal dimension
  double a1;
  //! Second longitudinal dimension
  double a2;
  //! Initial temperature (non-stationary problem)
  double Ti;
  //! Dirichlet condition
  double To;
  //! External temperature 
  double Te;
  //! Density (for non-stationary problem)
  double rho;
  //! Heat capacity (for non-stationary problem)
  double Cp; 
  //! Conductivity
  double k; 
  //! Convection coefficient
  double hc;
  //! Final time 
  double T;
  //! Number of elements in space grid
  int M;
  //! Number of elements in time grid
  int N;
  //! System solver type
  int solverType;
  //! max number of iteration for Gauss-Siedel
  int itermax;
  //! Tolerance for stopping criterion
  double toler;
  //! Norm type
  int normType; 
  //! Output mode
  int outputMode;
  //! Results filename
  std::string resultsFilename;
 
  //! Constructor takes default values
  parameters():
    stationary(0),
	L(40.),
    a1(4.),
    a2(50.),
	Ti(20.),
    To(46.),
    Te(20.),
	rho(1.0),
	Cp(1.0),
    k(0.164),
    hc(1.e-6*200.),
	T(1.0),
	M(100),
	N(100),
	solverType(0),
	itermax(1000000),
    toler(1e-8),
    normType(0),
	outputMode(0),
	resultsFilename("result.dat")
  {}
};
//! Prints parameters
std::ostream & operator << (std::ostream &,const parameters &);
#endif
