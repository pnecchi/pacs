#include <Mesh.h>
#include <TridiagonalMatrix.h>
#include <solveTridiagonalSystem.h>

void solveHeatEquation (double xMin,
                        double xMax,
                        unsigned long NNodes,
                        std::vector<double>& solution)
{
    // Initialize finite elements mesh
    Mesh mesh(xMin, xMax, NNodes);
    double h = mesh.getStep();

    // Assemble stiffness matrix
    TridiagonalMatrix A(NNodes);
    std::vector<std::vector<double>> A_loc(2, std::vector<double>(2, 0.0));
    for (unsigned long k = 0; k < NNodes - 1; k++)
    {
        A_loc[0][0] = 1.0 / h;
        A_loc[0][1] = -1.0 / h;
        A_loc[1][0] = -1.0 / h;
        A_loc[1][1] = 1.0 / h;

        std::vector<unsigned long> elementNodes(mesh.getElementNodes(k));
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            {
                A.element(elementNodes[i], elementNodes[j]) += A_loc[j][i];
            }
    }

    // Assemble rhs vector
    std::vector<double> rhs(NNodes, 0.0);
    std::vector<double> rhs_loc(2, 0.0);
    for (unsigned long k = 0; k < NNodes - 1; k++)
    {
        rhs_loc[0] = 0.5 * h;
        rhs_loc[1] = 0.5 * h;

        std::vector<unsigned long> elementNodes(mesh.getElementNodes(k));
        for (int i = 0; i < 2; i++)
        {
            rhs[elementNodes[i]] += rhs_loc[i];
        }
    }

    // Correct for boundary conditions
    A.element(0, 0) = 1.0;
    A.element(0, 1) = 0.0;
    A.element(NNodes-1, NNodes-2) = 0.0;
    A.element(NNodes-1, NNodes-1) = 1.0;
    rhs[0] = 0.0;
    rhs[NNodes-1]= 0.0;

    std::cout << A << std::endl;

    // Solve linear system
    solveTridiagonalSystem(A, rhs, solution);
}
