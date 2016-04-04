//----------------------------------------------------------------------
// Description: Matrix class for elements of type double
// Author:      Pierpaolo Necchi
// Email:       pierpaolo.necchi@gmail.com
// Date:        lun 04 apr 2016 14:04:17 CEST
//----------------------------------------------------------------------

#ifndef MATRIX_H
#define MATRIX_H

#include <tuple>
#include <vector>
#include <iostream>

//----------------------------------------------------------------------
// Base class
//----------------------------------------------------------------------

class Matrix
{
	friend std::ostream& operator<<(std::ostream& os, const Matrix& A);
	friend void swap(Matrix& first, Matrix& second);
	friend std::vector<double> solve(const Matrix& A, const std::vector<double> rhs);
	
public:
	// Constructors 
	Matrix ();
	Matrix (unsigned long nrows_, 
			unsigned long ncols_, 
			double val_);
	Matrix (const Matrix& matrix_);
	
	// Destuctor
	virtual ~Matrix (){}
	
	// Size
	std::tuple<unsigned long, unsigned long> size() const;
	unsigned long Nrows() const;
	unsigned long Ncols() const;
	
	// Assignment
	virtual Matrix& operator=(const Matrix& matrix_);

	// Access 
	virtual double& operator()(unsigned long rowIdx, unsigned long colIdx);
	virtual const double& operator()(unsigned long rowIdx, unsigned long colIdx) const;

	// Matrix operations
    Matrix transpose();
	Matrix operator+(const Matrix& rhs) const;
	Matrix& operator+=(const Matrix& rhs);
	Matrix operator-(const Matrix& rhs) const;
	Matrix& operator-=(const Matrix& rhs);
	Matrix operator*(const Matrix& rhs) const;
	Matrix& operator*=(const Matrix& rhs);
	
	// Scalar operations
	Matrix operator+(const double& rhs) const;
	Matrix& operator+=(const double& rhs);
	Matrix operator-(const double& rhs) const;
	Matrix& operator-=(const double& rhs);
	Matrix operator*(const double& rhs) const;
	Matrix& operator*=(const double& rhs);
	Matrix operator/(const double& rhs) const;
	Matrix& operator/=(const double& rhs);

	// Matrix-vector multiplication
	virtual std::vector<double> operator*(const std::vector<double>& rhs) const;

protected:
	unsigned long nrows;
	unsigned long ncols;

private:
	std::vector<std::vector<double>> mat;
};

std::ostream& operator<<(std::ostream& os, const Matrix& A);
void swap(Matrix& first, Matrix& second);
std::vector<double> solve(const Matrix& A, const std::vector<double> rhs);

#endif /* end of include guard: MATRIX_H */
