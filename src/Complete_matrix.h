//
// Created by robert on 26.02.23.
//

#ifndef SLAE_COMPLETE_MATRIX_H
#define SLAE_COMPLETE_MATRIX_H

#include <memory>
#include <iostream>

class Complete_matrix {
public:
    Complete_matrix();
    Complete_matrix(const Complete_matrix& A);
    Complete_matrix(unsigned int lines_number, unsigned int columns_number, const double a[]);
    void do_identity(unsigned int n);
    double get_element(unsigned int line_number, unsigned int column_number) const;
    std::unique_ptr<double[]> get_column(unsigned int column_number, unsigned int begin_number=0) const;
    unsigned int get_line_size() const;
    unsigned int get_column_size() const;
    std::unique_ptr<double[]> get_data() const;
    void write_column(const std::shared_ptr<double[]> a, unsigned int column_number, unsigned int begin_number=0);
    void write(double x, unsigned int line_number, unsigned int column_number);
    void transpose();
    Complete_matrix operator-(const Complete_matrix& A) const;
    void operator=(const Complete_matrix& A);
    Complete_matrix operator*(const Complete_matrix& A);
private:
    std::unique_ptr<double[]> data;
    unsigned int n;
    unsigned int m;
};

std::ostream& operator<<(std::ostream & os, const Complete_matrix& A);

#include "Complete_matrix.cpp"
#endif //SLAE_COMPLETE_MATRIX_H
