#include <Matrix.h>
#include <vector>
#include <iostream>

// Friend functions 
std::ostream& operator<<(std::ostream& os, const Matrix& A)
{
	for (unsigned long i = 0; i < A.Nrows(); ++i) {
		for (unsigned long j = 0; j < A.Ncols(); ++j) {
			std::cout << A(i,j) << ",";
		}
		std::cout << std::endl;
	}
	return os;
}

std::vector<double> solve(const Matrix& A, const std::vector<double> rhs)
{
	std::vector<double> result(A.Ncols(), 0.0);
	// TODO: implement linear system solver
	return result;
}

// Constructors 
Matrix::Matrix ()
{
	// Default constructor
}

Matrix::Matrix (unsigned long nrows_, 
				unsigned long ncols_, 
				double val_)
	: nrows(nrows_), 
	  ncols(ncols_) 
{
	mat.resize(nrows_);
	for (unsigned long i = 0; i < nrows_; i++) {
		mat[i].resize(ncols_, val_);
	}
}

Matrix::Matrix (const Matrix& matrix_)
	: nrows(matrix_.Nrows()), 
	  ncols(matrix_.Ncols()), 
	  mat(matrix_.mat)
{
	// Nothing to do
}

// Size
std::tuple<unsigned long, unsigned long> Matrix::size() const
{
	std::tuple<unsigned long, unsigned long> s = std::make_tuple(this->Nrows(), this->Ncols());
	return s;
}

unsigned long Matrix::Nrows() const
{
	return nrows;
}

unsigned long Matrix::Ncols() const
{
	return ncols;
}

// Assignment
Matrix& Matrix::operator=(const Matrix& matrix_)
{
	if (this != &matrix_)
	{
		Matrix temp(matrix_);
		swap(*this, temp);
	}
	return *this;
}

void swap(Matrix& first, Matrix& second)
{
	using std::swap;
	swap(first.nrows, second.nrows);
	swap(first.ncols, second.ncols);
	swap(first.mat, second.mat);
}

// Access
double& Matrix::operator()(unsigned long rowIdx, unsigned long colIdx)
{
	return (*this).mat[rowIdx][colIdx];
}

const double& Matrix::operator()(unsigned long rowIdx, unsigned long colIdx) const
{
	return (*this).mat[rowIdx][colIdx];
}

// Matrix operations
Matrix Matrix::transpose()
{
	Matrix result(this->Ncols(), this->Nrows(), 0.0);
	for (unsigned long i = 0; i < result.Nrows(); ++i) {
		for (unsigned long j = 0; j < result.Ncols(); ++j) {
			result(i,j) = (*this)(j,i);
		}
	}
	return result;

}

Matrix Matrix::operator+(const Matrix& rhs) const
{
	Matrix result(*this);
	result += rhs;
	return result;
}


Matrix& Matrix::operator+=(const Matrix& rhs)
{
	for (unsigned long i = 0; i < nrows; ++i) {
		for (unsigned long j = 0; j < ncols; ++j) {
			(*this)(i,j) += rhs(i,j);
		}
	}
	return *this;
}

Matrix Matrix::operator-(const Matrix& rhs) const
{
	Matrix result(*this);
	result -= rhs;
	return result;
}

Matrix& Matrix::operator-=(const Matrix& rhs)
{
	for (unsigned long i = 0; i < nrows; ++i) {
		for (unsigned long j = 0; j < ncols; ++j) {
			(*this)(i,j) -= rhs(i,j);
		}
	}
	return *this;
}

Matrix Matrix::operator*(const Matrix& rhs) const
{
	Matrix result(this->Nrows(), rhs.Ncols(), 0.0);
	for (unsigned long i = 0; i < result.Nrows(); ++i) {
		for (unsigned long j = 0; j < result.Ncols(); ++j){
			for (unsigned long k = 0; k < this->Ncols(); k++) {
				result(i,j) += (*this)(i,k) * rhs(k,j);
			}
		}
	}
	return result;
}

Matrix& Matrix::operator*=(const Matrix& rhs)
{
	Matrix result = (*this) * rhs;
	(*this) = result;
	return *this;
}

// Scalar operations
Matrix Matrix::operator+(const double& rhs) const
{
	Matrix result(*this);
	result += rhs;
	return result;
}

Matrix& Matrix::operator+=(const double& rhs)
{
	for (unsigned long i = 0; i < nrows; ++i) {
		for (unsigned long j = 0; j < ncols; ++j) {
			(*this)(i,j) += rhs;
		}
	}
	return *this;
}

Matrix Matrix::operator-(const double& rhs) const
{
	Matrix result(*this);
	result -= rhs;
	return result;
}

Matrix& Matrix::operator-=(const double& rhs)
{
	for (unsigned long i = 0; i < nrows; ++i) {
		for (unsigned long j = 0; j < ncols; ++j) {
			(*this)(i,j) -= rhs;
		}
	}
	return *this;
}

Matrix Matrix::operator*(const double& rhs) const
{
	Matrix result(*this);
	result *= rhs;
	return result;
}

Matrix& Matrix::operator*=(const double& rhs)
{
	for (unsigned long i = 0; i < nrows; ++i) {
		for (unsigned long j = 0; j < ncols; ++j) {
			(*this)(i,j) *= rhs;
		}
	}
	return *this;
}

Matrix Matrix::operator/(const double& rhs) const
{
	Matrix result(*this);
	result /= rhs;
	return result;
}

Matrix& Matrix::operator/=(const double& rhs)
{
	for (unsigned long i = 0; i < nrows; ++i) {
		for (unsigned long j = 0; j < ncols; ++j) {
			(*this)(i,j) /= rhs;
		}
	}
	return *this;
}

// Matrix-vector multiplication
std::vector<double> Matrix::operator*(const std::vector<double>& rhs) const
{
	std::vector<double> result(this->Nrows(), 0.0);
	for (unsigned long i = 0; i < result.size(); ++i) {
		for (unsigned long j = 0; j < this->Ncols(); ++j) {
			result[i] += (*this)(i,j) * rhs[j];
		}
	}
	return result;
}


		
