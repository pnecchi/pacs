#include "GetPot.hpp"
#include "readParameters.hpp"
#include <fstream>

parameters readParameters(std::string const & filename,bool verbose)
{
  // Parameter default constructor fills it with the defaults values
  parameters defaults;
  // checks if file exixts and is readable
  std::ifstream check(filename); 
  if(!check)
    {
      std::cerr<<"ERROR: Parameter file "<<filename<<" does not exist"<<std::endl;
      std::cerr<<"Reverting to default values."<<std::endl;
      if(verbose) std::cout<<defaults;
      check.close();
      return defaults;
    }
  else
    check.close();

  GetPot ifile(filename.c_str());
  parameters values;
  // Read parameters from getpot ddata base
  values.stationary=ifile("stationary",defaults.stationary);
  values.L=ifile("L",defaults.L);
  values.a1=ifile("a1",defaults.a1);
  values.a2=ifile("a2",defaults.a2);
  values.Ti=ifile("Ti",defaults.Ti);
  values.To=ifile("To",defaults.To);
  values.Te=ifile("Te",defaults.Te);
  values.rho=ifile("rho",defaults.rho);
  values.Cp=ifile("Cp",defaults.Cp);
  values.k=ifile("k",defaults.k);
  values.hc=ifile("hc",defaults.hc);
  values.T=ifile("T",defaults.T);
  values.M=ifile("M",defaults.M);
  values.N=ifile("N",defaults.N);
  values.solverType=ifile("solverType",defaults.solverType); 
  values.itermax=ifile("itermax",defaults.itermax);
  values.toler=ifile("toler",defaults.toler);
  values.normType=ifile("normType",defaults.normType);
  values.outputMode=ifile("outputMode", defaults.outputMode);
  values.resultsFilename=ifile("resultsFilename",defaults.resultsFilename.c_str());
   if(verbose)
    {
      std::cout<<"PARAMETER VALUES IN GETPOT FILE"<<"\n";
      ifile.print();
      std::cout<<std::endl;
      std::cout<<"ACTUAL VALUES"<<"\n"<<values;
    }
  return values;
}

