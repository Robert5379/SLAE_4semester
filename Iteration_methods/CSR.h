//
// Created by robert on 18.02.23.
//

#ifndef SLAE_CSR_H
#define SLAE_CSR_H

#include <memory>
#include <vector>
#include <iostream>
#include <cmath>

class mixed_num_vec{
public:
    mixed_num_vec(std::vector<double> x, unsigned int k);
    std::vector<double> X;
    unsigned int K;
};

class CSR{
public:
    CSR(){};
    CSR(unsigned int N, unsigned int m, const double a[]);
    double get_element(unsigned int line_number, unsigned int column_number)const;
    std::vector<double>  operator*(const std::vector<double>& v)const;
    std::unique_ptr<double[]> multiply(const double x[], unsigned int n) const;
    mixed_num_vec simple_iteration(const std::vector<double> & x0, const std::vector<double> & b, double tau,  double r)const;
    std::vector<double> Jacobi(const std::vector<double> & x0, const std::vector<double> & b, double r) const;
    std::vector<double> Gauss_Seidel(const std::vector<double> & x0, const std::vector<double> & b, double r) const;
    std::istream& operator >> (std::istream& in);
private:
    std::vector<double> data;
    std::vector<unsigned int> column_indexes;
    std::vector<unsigned int> line_indexes;
};
//#include "CSR.cpp"
#endif