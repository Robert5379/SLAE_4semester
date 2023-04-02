//
// Created by robert on 26.02.23.
//

#ifndef SLAE_HOUSEHOLDER_H
#define SLAE_HOUSEHOLDER_H

#include "Complete_matrix.h"
#include <cmath>

double scalar(const std::shared_ptr<double[]> v, const std::shared_ptr<double[]> x, unsigned int size);
double scalar(const std::shared_ptr<double[]> v, const double* x, unsigned int size);
double scalar(const double* v, const double* x, unsigned int size);
void do_column_R(double*v, std::shared_ptr<double[]> x, unsigned int size, double vv);
std::unique_ptr<Complete_matrix[]> Householder(const Complete_matrix& A);

//#include "Householder.cpp"
#endif //SLAE_HOUSEHOLDER_H
