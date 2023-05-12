//
// Created by robert on 26.02.23.
//

#include "Complete_matrix.h"

Complete_matrix::Complete_matrix(unsigned int lines_number, unsigned int columns_number) {
    this->n=lines_number;
    this->m=columns_number;
    this->data=std::make_unique<double[]>(lines_number*columns_number);
}

Complete_matrix::Complete_matrix(const Complete_matrix &A) {
    unsigned int i;
    this->n=A.get_column_size();
    this->m=A.get_line_size();
    this->data=A.get_data();
}

Complete_matrix::Complete_matrix(unsigned int n0, unsigned int m0, const double *a) {
    unsigned int i;
    this->n=n0;
    this->m=m0;
    this->data=std::make_unique<double[]>(n*m);
    for(i=0;i<n*m;i++){
        data[i]=a[i];
    }
}

void Complete_matrix::do_identity(unsigned int size) {
    n-size;
    m=size;
    data=std::make_unique<double[]>(n*m);
    unsigned int i, j;
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            if(i==j){
                data[i*m+j]=1;
            }else{
                data[i*m+j]=0;
            }
        }
    }
}

double Complete_matrix::get_element(unsigned int line_number, unsigned int column_number) const {
    return data[line_number*m+column_number];
}

std::unique_ptr<double[]> Complete_matrix::get_column(unsigned int column_number, unsigned int begin_number) const {
    std::unique_ptr<double[]> x=std::make_unique<double[]>(n-begin_number);
    unsigned int i;
    for(i=0;i<n-begin_number;i++){
        x[i]=data[(i+begin_number)*m+column_number];
    }
    return x;
}

unsigned int Complete_matrix::get_line_size() const{
    return m;
}

unsigned int Complete_matrix::get_column_size() const{
    return n;
}

std::unique_ptr<double[]> Complete_matrix::get_data() const{
    unsigned int i;
    std::unique_ptr<double[]> x=std::make_unique<double[]>(n*m);
    for(i=0;i<n*m;i++){
        x[i]=data[i];
    }
    return x;
}

void Complete_matrix::write_column(const std::shared_ptr<double[]> a, unsigned int column_number, unsigned int begin_number) {
    unsigned int i;
    for(i=0;i<n-begin_number;i++){
        data[(i+begin_number)*m+column_number]=a[i];
    }
}

void Complete_matrix::write(double x, unsigned int line_number, unsigned int column_number) {
    data[line_number*m+column_number]=x;
}

void Complete_matrix::transpose() {
    std::unique_ptr<double[]> new_data=std::make_unique<double[]>(n*m);
    unsigned int i, j;
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            new_data[j*n+i]=data[i*m+j];
        }
    }
    data=std::move(new_data);
}

void Complete_matrix::operator=(const Complete_matrix &A) {
    unsigned int i;
    this->n=A.get_column_size();
    this->m=A.get_line_size();
    this->data=A.get_data();
}

Complete_matrix Complete_matrix::operator*(const Complete_matrix &A) {
    double *new_data=new double[this->n*A.get_line_size()], t;
    std::unique_ptr<double[]> A_data=A.get_data();
    unsigned int i, j, n1= this->n, m1=A.get_line_size(), i1, j1;
    for(i1=0;i1<n1;i1++){
        for(j1=0;j1<m1;j1++){
            for(i=0;i< this->m;i++){
                t=this->data[i1*n1+i]*A_data[i*m1+j1];
            }
            new_data[i1*m1+j1]=t;
        }
    }
    Complete_matrix B( n1, m1, new_data);
    delete[] new_data;
    return B;
}

std::ostream& operator<<(std::ostream & os, const Complete_matrix& A) {
    unsigned int i, j, n=A.get_column_size(), m=A.get_line_size();
    std::unique_ptr<double[]> data=A.get_data();
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            os << data[i*m+j] << " ";
        }
        os << std::endl;
    }
    return os;
}
