#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "heatEquationSolver.h"

/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Iterative resolution by Gauss Siedel.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/
//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

//! main program
int main(int argc, char** argv)
{
	using namespace std; // avoid std::
	int status(0); // final program status
	GetPot   cl(argc, argv);
	if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
	// check if we want verbosity
	bool verbose=cl.search(1,"-v");
	// Get file with parameter values
	string inputFilename = cl.follow("parameters.pot","-p");
	cout<<"Reading parameters from "<<inputFilename<<std::endl;
	// read parameters
	const parameters param=readParameters(inputFilename,verbose);
	
	switch (param.stationary)
	{
		case 0:  // Stationary problem
			stationaryHeatEquationSolver(param);
			break;
		case 1:  // Non-stationary problem
			nonStationaryHeatEquationSolver(param);
			break;
	}

	return status;
}
