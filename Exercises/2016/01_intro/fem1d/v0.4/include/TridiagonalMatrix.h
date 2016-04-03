#ifndef TRIDIAGONALMATRIX_H
#define TRIDIAGONALMATRIX_H

#include <iostream>
#include <vector>

class TridiagonalMatrix
{
    friend std::ostream& operator<<(std::ostream &os, const TridiagonalMatrix& TridiagonalMatrix_);

    public:
        // Constructors
        TridiagonalMatrix();
        TridiagonalMatrix(unsigned long N_);
        TridiagonalMatrix(const std::vector<double>& lowerDiagonal_,
                          const std::vector<double>& diagonal_,
                          const std::vector<double>& upperDiagonal_);
        TridiagonalMatrix(const TridiagonalMatrix& other);

        // Destructor
        virtual ~TridiagonalMatrix(){}

        // Assignment
        TridiagonalMatrix& operator=(const TridiagonalMatrix& other);

        // Other methods
        unsigned long size() const;
        double element(unsigned long rowIdx, unsigned long colIdx) const;
        double lowerDiagonalElement(unsigned long idx) const;
        double diagonalElement(unsigned long idx) const;
        double upperDiagonalElement(unsigned long idx) const;
        double& element(unsigned long rowIdx, unsigned long colIdx);
        double& lowerDiagonalElement(unsigned long idx);
        double& diagonalElement(unsigned long idx);
        double& upperDiagonalElement(unsigned long idx);

    private:
        std::vector<double> lowerDiagonal;
        std::vector<double> diagonal;
        std::vector<double> upperDiagonal;
        unsigned long N;
};

std::ostream& operator<<(std::ostream &os, const TridiagonalMatrix& TridiagonalMatrix_);

#endif // TRIDIAGONALMATRIX_H
