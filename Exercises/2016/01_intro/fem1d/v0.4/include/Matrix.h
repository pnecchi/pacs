#ifndef MATRIX_H
#define MATRIX_H

#include <tuple>
#include <vector>
#include <iostream>

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
	Matrix& operator=(const Matrix& matrix_);

	// Access 
	std::vector<double>& operator[](unsigned long idx);
	const std::vector<double>& operator[](unsigned long idx) const;
	double& operator()(unsigned long rowIdx, unsigned long colIdx);
	const double& operator()(unsigned long rowIdx, unsigned long colIdx) const;

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
	std::vector<double> operator*(const std::vector<double>& rhs) const;
	
private:
	unsigned long nrows;
	unsigned long ncols;
	std::vector<std::vector<double>> mat;
};

std::ostream& operator<<(std::ostream& os, const Matrix& A);
void swap(Matrix& first, Matrix& second);
std::vector<double> solve(const Matrix& A, const std::vector<double> rhs);

#endif /* end of include guard: MATRIX_H */
