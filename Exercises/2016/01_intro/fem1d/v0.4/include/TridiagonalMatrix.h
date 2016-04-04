//----------------------------------------------------------------------
// Description: Tridiagonal matrix class
// Author:      Pierpaolo Necchi
// Email:       pierpaolo.necchi@gmail.com
// Date:        lun 04 apr 2016 14:21:40 CEST
//----------------------------------------------------------------------

#ifndef TRIDIAGONALMATRIX_H
#define TRIDIAGONALMATRIX_H

#include <iostream>
#include <vector>

class TridiagonalMatrix
{
    friend std::ostream& operator<<(std::ostream &os, const TridiagonalMatrix& TridiagonalMatrix_);
	friend void swap(TridiagonalMatrix& first, TridiagonalMatrix& second);
	friend void solve(const TridiagonalMatrix& A, 
					  const std::vector<double>& rhs, 
					  std::vector<double>& sol);
public:
	// Constructors
    TridiagonalMatrix();
    TridiagonalMatrix(unsigned long N_,
					  double lowerDiagonalVal_=0.0,
					  double diagonalVal_=0.0, 
					  double upperDiagonalVal_=0.0);
	TridiagonalMatrix(const std::vector<double>& lowerDiagonal_,
                      const std::vector<double>& diagonal_,
					  const std::vector<double>& upperDiagonal_);
    TridiagonalMatrix(const TridiagonalMatrix& matrix_);

    // Destructor
	virtual ~TridiagonalMatrix(){};

	// Assignment
	TridiagonalMatrix& operator=(const TridiagonalMatrix& matrix_);

	// Access 
	unsigned long size() const;
	double& operator()(unsigned long rowIdx, unsigned long colIdx);
	const double& operator()(unsigned long rowIdx, unsigned long colIdx) const;

private:
	unsigned long N;
	std::vector<double> lowerDiagonal;
	std::vector<double> diagonal;
    std::vector<double> upperDiagonal;
};

std::ostream& operator<<(std::ostream &os, const TridiagonalMatrix& TridiagonalMatrix_);
void swap(TridiagonalMatrix& first, TridiagonalMatrix& second);
void solve(const TridiagonalMatrix& A, 
      	   const std::vector<double>& rhs, 
		   std::vector<double>& sol);

#endif // TRIDIAGONALMATRIX_H
