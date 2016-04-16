#include <iostream>
#include "parameters.hpp"
std::ostream & operator << (std::ostream & out,const parameters & p)
{
  out<<"PARAMETER VALUES:"<<"\n";
  out<<"stationary= "<<p.stationary<<"\n";
  out<<"L= "<<p.L<<"\n";
  out<<"a1= "<<p.a1<<"\n";
  out<<"a2= "<<p.a2<<"\n";
  out<<"Ti= "<<p.Ti<<"\n"; 
  out<<"To= "<<p.To<<"\n";
  out<<"Te= "<<p.Te<<"\n";
  out<<"rho= "<<p.rho<<"\n";
  out<<"Cp= "<<p.Cp<<"\n";
  out<<"k= "<<p.k<<"\n";
  out<<"hc= "<<p.hc<<"\n";
  out<<"T= "<<p.T<<"\n";
  out<<"M= "<<p.M<<"\n";
  out<<"N= "<<p.N<<"\n";
  out<<"solverType= "<<p.solverType<<"\n";
  out<<"itermax= "<<p.itermax<<"\n";
  out<<"toler= "<<p.toler<<"\n";
  out<<"normType= "<<p.normType<<"\n";
  out<<"outputMode= "<<p.outputMode<<"\n";
  out<<"resultsFilename= "<<p.resultsFilename<<"\n";
  return out;
}
