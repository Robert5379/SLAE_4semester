//
// Created by robert on 07.04.2023.
//

#ifndef SLAE_SOLVER_H
#define SLAE_SOLVER_H

#include "Thridiagonal_matrix.cpp"
#include <vector>

template <typename T=double>
std::vector<T> solver(tridiagonal_matrix &A, std::vector<T> &f);

#endif //SLAE_SOLVER_H
