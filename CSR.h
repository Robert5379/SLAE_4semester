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

std::vector<double> operator*(double x, const std::vector<double>& v);
std::vector<double> operator+(const std::vector<double>& x, const std::vector<double>& v);
double operator*(const std::vector<double>& x, const std::vector<double>& v);
std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2);
double max_abs(const std::vector<double>& x);
double modul(const std::vector<double>& x);
std::ostream& operator<<(std::ostream & os, const std::vector<double>& v);

class CSR {
    friend std::istream& operator >> (std::istream& in, CSR& A);
public:
    CSR(){};
    CSR(unsigned int N, unsigned int m, const double a[]);
    double get_element(unsigned int line_number, unsigned int column_number)const;
    std::vector<double>  operator*(const std::vector<double>& v)const;
    std::unique_ptr<double[]> multiply(const double x[], unsigned int n) const;
    mixed_num_vec Simple_iteration(const std::vector<double> & x0, const std::vector<double> & b, double tau,  double r)const;
    mixed_num_vec Jacobi(const std::vector<double> & x0, const std::vector<double> & b, double r) const;
    mixed_num_vec Gauss_Seidel(const std::vector<double> & x0, const std::vector<double> & b, double r) const;
    std::vector<double> SIM_Chebyshev_acceleration(const std::vector<double> & x0, const std::vector<double> & b,double lambda_max, double lambda_min, double accuracy, double degree =5) const;
    std::vector<double> SOR(const std::vector<double> & x0, const std::vector<double> & b, double w, double accuracy) const;
    std::vector<double> Symmetrical_Gauss_Seidel(const std::vector<double> & x0, const std::vector<double> & b, double ro, double accuracy) const;
    std::vector<double> SSOR(const std::vector<double> & x0, const std::vector<double> & b, double w, double ro, double accuracy) const;
    std::vector<double> Steepest_descent(const std::vector<double> & x0, const std::vector<double> & b, double accuracy) const;
    std::vector<double> Heavy_ball(const std::vector<double> & x0, const std::vector<double> & b, double accuracy) const;
    std::vector<double> Сonjugate_gradient(const std::vector<double> & x0, const std::vector<double> & b, double accuracy) const;
private:
    std::vector<double> data;
    std::vector<unsigned int> column_indexes;
    std::vector<unsigned int> line_indexes;
};

std::istream& operator >> (std::istream& in, CSR& A);

//#include "CSR.cpp"
#endif