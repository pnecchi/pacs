#include <iostream>
#include <vector>
#include <heatEquation.h>
#include <GetPot>

int main (int argc, char* argv[])
{
	// Initialize parameters
	double xMin = 0.0;
	double xMax = 1.0;
	unsigned long NNodes = 10;
	std::vector<double> solution(NNodes);

	// Solve heat equation
	solveHeatEquation(xMin, xMax, NNodes, solution);

	for (unsigned long i = 0; i < NNodes; i++)
	{
		std::cout << solution[i] << std::endl;
	}

    return 0;
}

