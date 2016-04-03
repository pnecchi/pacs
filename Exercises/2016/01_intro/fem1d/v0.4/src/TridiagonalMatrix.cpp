#include <TridiagonalMatrix.h>

std::ostream& operator<<(std::ostream &os, const TridiagonalMatrix& TridiagonalMatrix_)
{
    for (size_t i = 0; i < TridiagonalMatrix_.size(); ++i)
    {
        os << TridiagonalMatrix_.lowerDiagonalElement(i) << " " <<
              TridiagonalMatrix_.diagonalElement(i) << " " <<
              TridiagonalMatrix_.upperDiagonalElement(i) << std::endl;
    }
    return os;
}

TridiagonalMatrix::TridiagonalMatrix()
{
    //ctor
}

TridiagonalMatrix::TridiagonalMatrix(size_t N_)
    : lowerDiagonal(N_), diagonal(N_), upperDiagonal(N_), N(N_)
{
    //
}

TridiagonalMatrix::TridiagonalMatrix(const std::vector<double>& lowerDiagonal_,
                                     const std::vector<double>& diagonal_,
                                     const std::vector<double>& upperDiagonal_)
    : lowerDiagonal(lowerDiagonal_),
      diagonal(diagonal_),
      upperDiagonal(upperDiagonal_),
      N(diagonal_.size())
{
    //
}

TridiagonalMatrix::TridiagonalMatrix(const TridiagonalMatrix& other)
    : lowerDiagonal(other.lowerDiagonal),
      diagonal(other.diagonal),
      upperDiagonal(other.upperDiagonal),
      N(other.size())
{
    //
}

TridiagonalMatrix& TridiagonalMatrix::operator=(const TridiagonalMatrix& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

size_t TridiagonalMatrix::size() const
{
    return N;
}

double TridiagonalMatrix::element(unsigned long rowIdx, unsigned long colIdx) const
{
    double val = 0.0;
	switch (static_cast<int>(colIdx) - static_cast<int>(rowIdx)) {
		case (-1):
			val = lowerDiagonal[rowIdx];
			break;
		case 0:
			val = diagonal[rowIdx];
			break;
		case 1:
			val = upperDiagonal[rowIdx];
			break;
		default:
			break;
	}
	return val; 
}

double TridiagonalMatrix::lowerDiagonalElement(size_t idx) const
{
    return lowerDiagonal[idx];
}

double TridiagonalMatrix::diagonalElement(size_t idx) const
{
    return diagonal[idx];
}

double TridiagonalMatrix::upperDiagonalElement(size_t idx) const
{
    return upperDiagonal[idx];
}

double& TridiagonalMatrix::element(unsigned long rowIdx, unsigned long colIdx)
{
	switch (static_cast<int>(colIdx) - static_cast<int>(rowIdx)) {
		case (-1):
			return lowerDiagonal[rowIdx];
		case 0:
			return diagonal[rowIdx];
		case 1:
			return upperDiagonal[rowIdx];
		default:
			break;
	} 
}

double& TridiagonalMatrix::lowerDiagonalElement(unsigned long idx)
{
    return lowerDiagonal[idx];
}

double& TridiagonalMatrix::diagonalElement(unsigned long idx)
{
    return diagonal[idx];
}

double& TridiagonalMatrix::upperDiagonalElement(unsigned long idx)
{
    return upperDiagonal[idx];
}

