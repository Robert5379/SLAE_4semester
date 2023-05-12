//
// Created by robert on 18.02.23.
//

#ifndef SLAE_CSR_H
#define SLAE_CSR_H

#include <vector>
#include <cmath>
#include "Complete_matrix.h"


class mixed_num_vec{
public:
    mixed_num_vec(std::vector<double> x, unsigned int k);
    std::vector<double> X;
    unsigned int K;
};

std::vector<double> operator*(double x, const std::vector<double>& v);
std::vector<double> operator+(const std::vector<double>& x, const std::vector<double>& v);
double operator*(const std::vector<double>& x, const std::vector<double>& v);
std::vector<double> operator/(const std::vector<double>& x, double a);
std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2);
double max_abs(const std::vector<double>& x);
double modul(const std::vector<double>& x);
std::ostream& operator<<(std::ostream & os, const std::vector<double>& v);

class CSR {
    friend std::istream& operator >> (std::istream& in, CSR& A);
public:
    CSR(){};
    CSR(unsigned int lines_number, unsigned int columns_number, const double a[]);
    CSR(const CSR& A);
    double get_element(unsigned int line_number, unsigned int column_number)const;
    std::vector<double>  operator*(const std::vector<double>& v)const;
    void transpose();
    std::unique_ptr<double[]> multiply(const double x[], unsigned int n) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> Simple_iteration(const std::vector<double> & x0, const std::vector<double> & b, double tau,  double r, unsigned int iterations=-1)const;
    std::pair<std::vector<double>, unsigned int> Jacobi(const std::vector<double> & x0, const std::vector<double> & b, double r) const;
    std::pair<std::vector<double>, unsigned int> Gauss_Seidel(const std::vector<double> & x0, const std::vector<double> & b, double r) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> SIM_Chebyshev_acceleration(const std::vector<double> & x0, const std::vector<double> & b,double lambda_max, double lambda_min, double accuracy, double degree =5, unsigned int iterations=-1) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>>SOR(const std::vector<double> & x0, const std::vector<double> & b, double w, double accuracy, unsigned int iterations=-1) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> Symmetrical_Gauss_Seidel(const std::vector<double> & x0, const std::vector<double> & b, double ro, double accuracy, unsigned int iterations=-1) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> SSOR(const std::vector<double> & x0, const std::vector<double> & b, double w, double ro, double accuracy, unsigned int iterations=-1) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> Steepest_descent(const std::vector<double> & x0, const std::vector<double> & b, double accuracy, unsigned int iterations=-1) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> Heavy_ball(const std::vector<double> & x0, const std::vector<double> & b, double accuracy, unsigned int iterations=-1) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> Conjugate_gradient(const std::vector<double> & x0, const std::vector<double> & b, double accuracy, unsigned int iterations=-1) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,unsigned int>> Cholesky_CG(const std::vector<double> & x0, const std::vector<double> & b, double accuracy, unsigned int iterations=-1) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,unsigned int>> GMRES(const std::vector<double> & x0, const std::vector<double> & b, double accuracy) const;
    std::pair<std::vector<double>, std::pair<std::vector<double>,unsigned int>> BiCG(const std::vector<double> & x0, const std::vector<double> & b, double accuracy) const;
private:
    std::vector<double> data;
    std::vector<unsigned int> column_indexes;
    std::vector<unsigned int> line_indexes;
};

std::istream& operator >> (std::istream& in, CSR& A);

//#include "CSR.cpp"
#endif