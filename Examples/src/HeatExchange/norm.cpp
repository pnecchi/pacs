#include <norm.h>
#include <math.h>	  // sqrt
#include <algorithm>  // std::max
#include <cmath>      // std::abs

double normSup(std::vector<double> const &v)
{
	double normSup = 0.0;
	double temp = 0.0;
	for(size_t i = 0; i < v.size(); ++i)
	{
		temp =  std::abs(v[i]);
		normSup = std::max(normSup, temp);
	}
	return normSup;
}

// pierpaolo - lun 11 apr 2016 15:14:17 CEST
// TODO: The current implementation considers that the vector values 
// correspond to a uniform grid with step h. Write a general version for arbitrary grid. 

double normL2(std::vector<double> const &v, double const h)
{
	double sum = 0.5 * v[0] * v[0];
	for(std::size_t i = 0; i < v.size() - 1; ++i)
	{
		sum += v[i] * v[i];
	}
	sum += 0.5 * v[v.size() - 1] * v[v.size() - 1];
	return sqrt(h * sum);
}

double normH1(std::vector<double> const &v, double const h)
{	
	// Compute first derivative on the grid points
	std::size_t N = v.size();
	std::vector<double> vPrime(N, 0.0);
	vPrime[0] = (v[1] - v[0]) / h;
	vPrime[N-1] = (v[N-1] - v[N-2]) / h;
	for(std::size_t i = 1; i < N-1; ++i)
	{
		vPrime[i] = (v[i+1] - v[i-1]) / (2 * h);
	}
	
	// Compute H_1 norm of v 
	double vNormL2 = normL2(v, h);
	double vPrimeNormL2 = normL2(vPrime, h);
	double vNormH1 = sqrt(vNormL2 * vNormL2 + 
						  vPrimeNormL2 * vPrimeNormL2);
	return vNormH1;
}


