//
// Created by robert on 15.02.23.
//
#include "CSR.h"


CSR::CSR(unsigned int n, unsigned int m, double a[]) {
    unsigned int k=0;
    unsigned int i, j;
    double t1;
    this->line_indexes.push_back(k);
    for(j=0;j<n;j++) {
        for (i = 0; i < m; i++) {
            if (a[i+j*n] != 0) {
                this->data.push_back(a[i+j*m]);
                this->column_indexes.push_back(i);
                k++;
            }
        }
        this->line_indexes.push_back(k);
    }
}

std::istream& CSR::operator>>(std::istream &in)
{
    unsigned int m, k=0, n;
    in >> n >> m;
    std::vector<double> data;
    unsigned int i, j;
    double t1;
    this->line_indexes.push_back(k);
    for(j=0;j<n;j++) {
        for (i = 0; i < m; i++) {
            in >> t1;
            if (t1 != 0) {
                this->data.push_back(t1);
                this->column_indexes.push_back(i);
                k++;
            }
        }
        this->line_indexes.push_back(k);
    }

    return in;
}

std::unique_ptr<std::vector<double>> CSR::operator*(std::vector<double> &v) {
    unsigned int z=0, i;
    std::unique_ptr<std::vector<double>> ans=std::make_unique<std::vector<double>>();
    double t;
    for(i=1;i< this->line_indexes.size();i++) {
        t=0;
        for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
            t+= this->data[z]*v[this->column_indexes[z]];
        }
        ans->push_back(t);
    }
    return ans;
}

double CSR::get_element(unsigned int i, unsigned int j) {
    unsigned int z=0;
    for(z=this->line_indexes[i];z<this->line_indexes[i+1];z++){
        if(this->column_indexes[z]==j){
            return this->data[z];
        }
    }
    return 0;
}

