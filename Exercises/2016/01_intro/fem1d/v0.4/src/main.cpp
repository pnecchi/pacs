#include <iostream>
#include <vector>
#include <heatEquation.h>
#include <GetPot>
#include <Matrix.h>

using namespace std;

int main (int argc, char* argv[])
{
	Matrix A(5, 5, 1.0);
	Matrix B(5, 5, 2.0);
	double c = 2.0;

	std::cout << A + B << std::endl;
	std::cout << A / c << std::endl;
	
	//   double xMin = 0.0;
    //double xMax = 1.0;
    //unsigned long NNodes = 10;
    //std::vector<double> solution(NNodes);

    //// Solve heat equation
    //solveHeatEquation(xMin, xMax, NNodes, solution);

    //for (unsigned long i = 0; i < NNodes; i++)
    //{
    //    std::cout << solution[i] << std::endl;
    //}

    return 0;
}

