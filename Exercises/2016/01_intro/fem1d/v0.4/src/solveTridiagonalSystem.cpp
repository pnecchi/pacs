#include <TridiagonalMatrix.h>
#include <vector>

void solveTridiagonalSystem(const TridiagonalMatrix& A,
                            const std::vector<double>& rhs,
                            std::vector<double>& x)
{
    // Initialize auxiliary variables
    unsigned long N = A.size();
    std::vector<double> temp_ud(N, 0.0);
    std::vector<double> temp_rhs(N, 0.0);
    double m;

    // Forward sweep
    m = 1.0 / A.diagonalElement(0);
    temp_ud[0] = m * A.upperDiagonalElement(0);
    temp_rhs[0] = m * rhs[0];
    for (unsigned long i = 1; i < N; i++)
    {
        m = 1.0 / (A.diagonalElement(i) - A.lowerDiagonalElement(i) * temp_ud[i-1]);
        temp_ud[i]  = m * A.upperDiagonalElement(i);
        temp_rhs[i] = m * (rhs[i] - A.lowerDiagonalElement(i) * temp_rhs[i-1]);
    }


    // Backward sweep
    x[N-1] = temp_rhs[N-1];
    for (unsigned long i = N-1; i-- > 0; )
    {
        x[i] = temp_rhs[i] - temp_ud[i] * x[i+1];
    }
}
