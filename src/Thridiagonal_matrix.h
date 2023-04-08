//
// Created by robert on 07.04.2023.
//

#ifndef SLAE_THRIDIAGONAL_MATRIX_H
#define SLAE_THRIDIAGONAL_MATRIX_H

#include <iostream>
#include <vector>
#include "Triplex.cpp"

class tridiagonal_matrix{
public:
    tridiagonal_matrix(){};
    void do_matrix(unsigned int N, double a[]);
    unsigned int n;
    std::vector <triplex> data;
};

std::istream& operator >> (std::istream& in, tridiagonal_matrix& p);


#endif //SLAE_THRIDIAGONAL_MATRIX_H
