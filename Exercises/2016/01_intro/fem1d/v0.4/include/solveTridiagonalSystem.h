#ifndef SOLVETRIDIAGONALSYSTEM_H
#define SOLVETRIDIAGONALSYSTEM_H

void solveTridiagonalSystem(const TridiagonalMatrix& A,
                            const std::vector<double>& rhs,
                            std::vector<double>& x);

#endif

