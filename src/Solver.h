//
// Created by robert on 07.04.2023.
//

#ifndef SLAE_SOLVER_H
#define SLAE_SOLVER_H

#include "Thridiagonal_matrix.h"
#include <vector>

template <typename T=double>
std::vector<T> solver(tridiagonal_matrix &A, std::vector<T> &f);

//std::vector<double> solver(tridiagonal_matrix &A, std::vector<double> &f);

#endif //SLAE_SOLVER_H
