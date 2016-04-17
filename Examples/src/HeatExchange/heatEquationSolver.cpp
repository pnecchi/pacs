#include "heatEquationSolver.h"
#include "parameters.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
#include "norm.h"
#include "TridiagonalMatrix.h"
#include <iostream>

/*! \fn void stationaryHeatEquationSolver(parameters const &param) 
 *  \brief Procedure that solves the steady heat equation in a bar.
 *  \param param The problem parameters 
 */

void stationaryHeatEquationSolver(parameters const &param)
{
	// Read parameters
	const int& itermax= param.itermax;	// Max number of iteration for Gauss-Siedel
	const double& toler=param.toler;	// Tolerance for stopping criterion
	const auto& L= param.L;				// Bar length
	const auto& a1=param.a1;			// First longitudinal dimension
	const auto& a2=param.a2;			// Second longitudinal dimension
	const auto& To=param.To;			// Dirichlet condition
	const auto& Te=param.Te;			// External temperature (Centigrades)
	const auto& k=param.k;				// Thermal conductivity
	const auto& hc=param.hc;			// Convection coefficient
	const auto& M=param.M;				// Number of grid elements
  
	//! Precomputed coefficient for adimensional form of equation
	const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

	// mesh size
	const auto h=1./M;
  
	// Solution vector
	std::vector<double> theta(M+1);
    
	switch (param.solverType)
	{
		case 0:  // Direct solver: Thomas algorithm
			{  // Brackets create scope for local initialization
				// Initialize linear system
				TridiagonalMatrix A(M+1, -1.0, 2.0 + act * h * h, -1.0);
				std::vector<double> rhs(M+1, 0.0);
	
				// Impose boundary conditions 
				A(0, 0) = 1.0;
				A(0, 1) = 0.0;
				A(M, M) = 1.0;
				rhs[0] = (To - Te) / Te;

				// Solve linear system
				solve(A, rhs, theta);
			}
			break;
		case 1:  // Iterative solver: Gauss-Seidel algorithm
			{
				// Gauss Siedel is initialised with a linear variation of T
				for(unsigned int m=0;m <= M;++m)
					theta[m]=(1.-m*h)*(To-Te)/Te;
  
				// Gauss-Seidel
				// epsilon=||x^{k+1}-x^{k}||
				// Stopping criteria epsilon<=toler
				int iter=0;
				double xnew=0.0;
				double epsilon=0.0;
				std::vector<double> thetaDelta(M+1, 0.0); 
				do
				{ 
					// first M-1 row of linear system
					for(int m=1;m < M;m++)
					{   
						xnew  = (theta[m-1]+theta[m+1])/(2.+h*h*act);
						thetaDelta[m] = xnew - theta[m];
						theta[m] = xnew;
					}
				
					//Last row
					xnew = theta[M-1]; 
					thetaDelta[M] = xnew - theta[M];
					theta[M] = xnew;
		
					// Compute norm of the increment
					switch (param.normType)
					{
						case 0:
							epsilon = normSup(thetaDelta);
							break;
						case 1:
							epsilon = normH1(thetaDelta, h);
							break;
						case 2:
							epsilon = normL2(thetaDelta, h);
							break;
					}
		
					// Increment iteration 
					iter=iter+1;     

				}while((epsilon > toler) && (iter < itermax) );
	
				if(iter<itermax)
					std::cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<std::endl;
				else
				{
					std::cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
							"||dx||="<<sqrt(epsilon)<<std::endl;
				}
			}
			break;
	}  /* end switch (param.solverType) */
  
	// Analytic solution
	std::vector<double> thetaa(M+1);
    for(int m=0;m <= M;m++)
		thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));

	// Initialize solution vectors
    std::vector<double> coor(M+1);
    std::vector<double> sol(M+1);
    std::vector<double> exact(M+1);
	
	// Transform results in original dimensions
	for(int m = 0; m <= M; ++m)
	{
		std::tie(coor[m],sol[m],exact[m])=
		std::make_tuple(m*h*L,Te*(1.+theta[m]),thetaa[m]);
	}
	
	// Output	
	if (param.outputMode == 0 || param.outputMode == 2)
	{
		std::cout << "Result file: " << param.resultsFilename << std::endl;
		std::ofstream f(param.resultsFilename);
		for(int m = 0; m<= M; m++)
		{
			f<<coor[m]<<"\t"<<sol[m]<<"\t"<<exact[m]<<std::endl;
		}
		f.close();
	}
	if (param.outputMode == 1 || param.outputMode == 2)
	{
		// Using temporary files (another nice use of tie)
		Gnuplot gp;	
		gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
			"w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
			"w l title 'uex'"<<std::endl;
	}
}


/*! \fn void nonStationaryHeatEquationSolver(parameters const &param) 
 *  \brief Procedure that solves the non-steady heat equation in a bar.
 *  \param param The problem parameters 
 */

void nonStationaryHeatEquationSolver(parameters const &param)
{
	// Read parameters
	const auto& L= param.L;		// Bar length
	const auto& a1=param.a1;	// First longitudinal dimension
	const auto& a2=param.a2;	// Second longitudinal dimension
	const auto& Ti=param.Ti;    // Initial temperature 
	const auto& To=param.To;	// Dirichlet condition
	const auto& Te=param.Te;	// External temperature (Centigrades)
	const auto& rho=param.rho;  // Bar density
	const auto& Cp=param.Cp;    // Bar heat capacity at constant pressure
	const auto& k=param.k;		// Thermal conductivity
	const auto& hc=param.hc;	// Convection coefficient
	const auto& M=param.M;		// Number of grid elements
	const auto& T=param.T;      // Final time. 
	const auto& N=param.N;		// Number of time gridpoints

	// Precompute useful quantities
	const auto hp = 2.0 * hc * (a1 + a2) / (a1 * a2);	// coefficient of heat exchange per unit length
	const auto alpha = k * T / (rho * Cp * L * L);		// adimensional diffusion coefficient 
    const auto beta = hp * T / (rho * Cp);				// adimensional reaction coefficient
	const auto h = 1.0 / M;								// mesh size
	double dt = 1.0 / N;								// time grid size
	double BC0 = (To - Te) / Te;						// Dirichlet boundary condition for x = 0

	// Initialize solution vector and rhs 
	std::vector<double> theta(M+1, Ti / Te - 1.0);	// Solution vector
	std::vector<double> rhs(M+1, 0.0);				// rhs vector of CN linear system

	// Initialize iteration matrices for Crank-Nicholson + FEM scheme
	double FF_diag = 2.0 * h / 3.0;
	double FF_offdiag = h / 6.0;
	double DD_diag = 2.0 / h;
	double DD_offdiag = -1.0 / h;
	double H_lhs_diag = FF_diag + 0.5 * dt * (alpha * DD_diag + beta * FF_diag);
	double H_lhs_offdiag = FF_offdiag + 0.5 * dt * (alpha * DD_offdiag + beta * FF_offdiag);
	double H_rhs_diag = FF_diag - 0.5 * dt * (alpha * DD_diag + beta * FF_diag);
	double H_rhs_offdiag = FF_offdiag - 0.5 * dt * (alpha * DD_offdiag + beta * FF_offdiag);

	TridiagonalMatrix H_lhs(M+1, H_lhs_offdiag, H_lhs_diag, H_lhs_offdiag);
	TridiagonalMatrix H_rhs(M+1, H_rhs_offdiag, H_rhs_diag, H_rhs_offdiag);

/*
 *    // Finite differences 
 *    double G_diag = 2.0 * alpha / (h * h) + beta;
 *    double G_offdiag = - alpha / (h * h);
 *
 *    TridiagonalMatrix H_lhs(M+1, 
 *                            0.5 * dt * G_offdiag, 
 *                            1 + 0.5 * dt * G_diag, 
 *                            0.5 * dt * G_offdiag);
 *
 *    TridiagonalMatrix H_rhs(M+1, 
 *                            - 0.5 * dt * G_offdiag, 
 *                            1 - 0.5 * dt * G_diag, 
 *                            - 0.5 * dt * G_offdiag);
 */

	// Boundary conditions 
	theta[0] = BC0;
	H_lhs(0, 0) = 1.0;
	H_lhs(0, 1) = 0.0;
	H_lhs(M, M-1) = - 1.0;
	H_lhs(M, M) = 1.0;
	
		// Initialize solution vectors for plot
	std::vector<double> tVec((M+1) * (N+1)); 
	std::vector<double> xVec((M+1) * (N+1));	
	std::vector<double> solVec((M+1) * (N+1));

	for(size_t m = 0; m < M+1; ++m)
	{
		std::tie(tVec[m], xVec[m], solVec[m]) =
			std::make_tuple(0.0, m * h * L, Te * (1.0 + theta[m]));
	}

	// Solve ODE with Crank-Nicholson scheme
	for(size_t n = 1; n < N+1; ++n)
	{
		// Compute rhs and fix boundary conditions
		rhs = H_rhs * theta;
		rhs[0] = BC0;
		rhs[M] = 0.0;

		// Solve linear system 
		solve(H_lhs, rhs, theta);
	
		// Store coordinate and solutions 
		for(size_t m = 0; m < M+1; ++m)
		{
			size_t idx = n * (M + 1) + m;
			std::tie(tVec[idx], xVec[idx], solVec[idx]) =
				std::make_tuple(n * dt * T, m * h * L, Te * (1.0 + theta[m]));
		}
	}
		
	// Output	
	if (param.outputMode == 0 || param.outputMode == 2)
	{
		std::cout << "Result file: " << param.resultsFilename << std::endl;
		std::ofstream f(param.resultsFilename);
		for(size_t k = 0; k <= tVec.size(); k++)
		{
			f<<tVec[k]<<"\t"<<xVec[k]<<"\t"<<solVec[k]<<std::endl;
		}
		f.close();
	}
	if (param.outputMode == 1 || param.outputMode == 2)
	{
		// Using temporary files (another nice use of tie)
		Gnuplot gp;	
		gp << "set pm3d" << std::endl;
		gp << "set style data lines" << std::endl;
		gp << "unset key" << std::endl;
	    gp << "set title 'Transient Heat Exchange Problem'" << std::endl;	
		gp << "set xlabel 'Space - x' rotate parallel offset 2,-2" << std::endl;
		gp << "set ylabel 'Time - t' rotate parallel offset -2,-2" << std::endl;
		gp << "set zlabel 'Temperature - T' rotate parallel offset -2,0" << std::endl;
		gp << "set hidden3d" << std::endl;
		gp << "set dgrid3d " << N+1 << "," << M+1 << std::endl; // " qnorm 2" << std::endl;
		gp << "splot"<<gp.file1d(std::tie(xVec,tVec,solVec)) << "pal" << std::endl;
	}
}
