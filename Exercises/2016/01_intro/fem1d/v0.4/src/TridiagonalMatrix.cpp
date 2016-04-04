//----------------------------------------------------------------------
// Description: Tridiagonal matrix class - implementation
// Author:      Pierpaolo Necchi
// Email:       pierpaolo.necchi@gmail.com
// Date:        lun 04 apr 2016 17:27:58 CEST
//----------------------------------------------------------------------

#include <TridiagonalMatrix.h>

//----------------------------------------------------------------------
// Frind functions 
//----------------------------------------------------------------------

std::ostream& operator<<(std::ostream &os, const TridiagonalMatrix& TridiagonalMatrix_)
{
    for (size_t i = 0; i < TridiagonalMatrix_.size(); ++i)
    {
        os << TridiagonalMatrix_.lowerDiagonal[i] << " " <<
              TridiagonalMatrix_.diagonal[i] << " " <<
              TridiagonalMatrix_.upperDiagonal[i] << std::endl;
    }
    return os;
}

void swap(TridiagonalMatrix& first, TridiagonalMatrix& second)
{
	using std::swap;
	swap(first.N, second.N);
	swap(first.lowerDiagonal, second.lowerDiagonal);
	swap(first.diagonal, second.diagonal);
	swap(first.upperDiagonal, second.upperDiagonal);
}

void solve(const TridiagonalMatrix& A,
           const std::vector<double>& rhs,
		   std::vector<double>& sol)
{
    // Initialize auxiliary variables
    unsigned long N = A.size();
    std::vector<double> temp_ud(N, 0.0);
    std::vector<double> temp_rhs(N, 0.0);
    double m;

    // Forward sweep
    m = 1.0 / A.diagonal[0];
    temp_ud[0] = m * A.upperDiagonal[0];
    temp_rhs[0] = m * rhs[0];
    for (unsigned long i = 1; i < N; i++)
    {
        m = 1.0 / (A.diagonal[i] - A.lowerDiagonal[i] * temp_ud[i-1]);
        temp_ud[i]  = m * A.upperDiagonal[i];
        temp_rhs[i] = m * (rhs[i] - A.lowerDiagonal[i] * temp_rhs[i-1]);
    }

    // Backward sweep
    sol[N-1] = temp_rhs[N-1];
    for (unsigned long i = N-1; i-- > 0; )
    {
        sol[i] = temp_rhs[i] - temp_ud[i] * sol[i+1];
    }
}

//----------------------------------------------------------------------
// Member functions
//----------------------------------------------------------------------


TridiagonalMatrix::TridiagonalMatrix()
{
    // Standard constructor
}

TridiagonalMatrix::TridiagonalMatrix(unsigned long N_,
				  double lowerDiagonalVal_,
				  double diagonalVal_, 
				  double upperDiagonalVal_)
	: N(N_), 
	  lowerDiagonal(N_, lowerDiagonalVal_),
	  diagonal(N_, diagonalVal_),
	  upperDiagonal(N_, upperDiagonalVal_)
{
	// Nothing to do
}

TridiagonalMatrix::TridiagonalMatrix(const std::vector<double>& lowerDiagonal_,
                                     const std::vector<double>& diagonal_,
                                     const std::vector<double>& upperDiagonal_)
	: N(diagonal_.size()),
	  lowerDiagonal(lowerDiagonal_),
      diagonal(diagonal_),
      upperDiagonal(upperDiagonal_)
{
    // Nothing to do
}

TridiagonalMatrix::TridiagonalMatrix(const TridiagonalMatrix& matrix_)
	: N(matrix_.N), 
	  lowerDiagonal(matrix_.lowerDiagonal),
	  diagonal(matrix_.diagonal),
	  upperDiagonal(matrix_.upperDiagonal)
{
	// Nothing to do
}

TridiagonalMatrix& TridiagonalMatrix::operator=(const TridiagonalMatrix& rhs)
{
	if (this != &rhs)
	{
		TridiagonalMatrix temp(rhs);
		swap(*this, temp);
	}
	return *this;
}

// Access 
unsigned long TridiagonalMatrix::size() const
{
	return N;
}

double& TridiagonalMatrix::operator()(unsigned long rowIdx, unsigned long colIdx)
{
	switch (static_cast<int>(colIdx) - static_cast<int>(rowIdx)) {
		case (-1):
			return lowerDiagonal[rowIdx];
		case 0:
			return diagonal[rowIdx];
		case 1:
			return upperDiagonal[rowIdx];
		default:
			return lowerDiagonal[0];
	}
}

const double& TridiagonalMatrix::operator()(unsigned long rowIdx, unsigned long colIdx) const
{
	switch (static_cast<int>(colIdx) - static_cast<int>(rowIdx)) {
		case (-1):
			return lowerDiagonal[rowIdx];
		case 0:
			return diagonal[rowIdx];
		case 1:
			return upperDiagonal[rowIdx];
		default:
			return lowerDiagonal[0];
	}
}
