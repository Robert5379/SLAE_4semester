//
// Created by robert on 18.02.23.
//

#ifndef SLAE_CSR_H
#define SLAE_CSR_H

#include <memory>
#include <vector>
#include <iostream>


class CSR{
public:
    CSR(){};
    CSR(unsigned int N, unsigned int m, double a[]);
    double get_element(unsigned int i, unsigned int j);
    std::unique_ptr<std::vector<double>>  operator*(std::vector<double>& v);
    std::istream& operator >> (std::istream& in);
private:
    std::vector<double> data;
    std::vector<unsigned int> line_indexes;
    std::vector<unsigned int> column_indexes;
};
#include "CSR.cpp"
#endif