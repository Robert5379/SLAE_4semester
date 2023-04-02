//
// Created by robert on 15.02.23.
//
#include "CSR.h"


CSR::CSR(unsigned int n, unsigned int m, const double a[]) {
    unsigned int k=0;
    unsigned int i, j;
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

std::vector<double> CSR::operator*(const std::vector<double> &v) const {
    unsigned int z=0, i;
    std::vector<double> ans(v.size());
    double t;
    for(i=1;i< this->line_indexes.size();i++) {
        t=0;
        for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
            t+= this->data[z]*v[this->column_indexes[z]];
        }
        ans[i-1]=t;
    }
    return ans;
}

double CSR::get_element(unsigned int i, unsigned int j) const {
    unsigned int z=0;
    for(z=this->line_indexes[i];z<this->line_indexes[i+1];z++){
        if(this->column_indexes[z]==j){
            return this->data[z];
        }
    }
    return 0;
}

std::unique_ptr<double[]> CSR::multiply(const double *x, unsigned int n) const {
    unsigned int z=0, i;
    std::unique_ptr<double[]> ans=std::make_unique<double[]>(n);
    double t;
    for(i=1;i< this->line_indexes.size();i++) {
        t=0;
        for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
            t+= this->data[z]*x[this->column_indexes[z]];
        }
        ans[i-1]=t;
    }
    return ans;
}

double modul(const std::vector<double>& x){
    double ans=0;
    unsigned int i;
    for(i=0;i<x.size();i++){
        ans+=x[i]*x[i];
    }
    return sqrt(ans);
}

double max_abs(const std::vector<double>& x){
    double ans=0;
    unsigned int i;
    for(i=0;x.size();i++){
        if(ans<abs(x[i])){
            ans=abs(x[i]);
        }
    }
    return ans;
}

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2){
    std::vector<double> ans(v1.size());
    unsigned int i;
    for(i=0;i<v1.size();i++){
        ans[i]=v1[i]-v2[i];
    }
    return ans;
}

std::vector<double> operator*(double x, const std::vector<double>& v){
    std::vector<double> ans(v.size());
    unsigned int i;
    for(i=0;i<v.size();i++){
        ans[i]=v[i]*x;
    }
    return ans;
}

mixed_num_vec::mixed_num_vec(std::vector<double> x, unsigned int k):X(x), K(k) {}


mixed_num_vec CSR::simple_iteration(const std::vector<double> & x0, const std::vector<double> &b, double tau, double r) const {
    std::vector<double>x=x0;
    std::vector<double>r_n=(*this)*x-b;
    unsigned int k=0;
    while(modul(r_n)>r){
        x=x-tau*(r_n);
        r_n=(*this)*x-b;
        k++;
    }
    return mixed_num_vec(x, k);
}

std::vector<double> CSR::Jacobi(const std::vector<double> &x0, const std::vector<double> &b, double r) const {
    std::vector<double>x=x0;
    std::vector<double>x_next(x0.size());
    std::vector<double>D_(x0.size());
    double t;
    unsigned int i, z, n=x0.size();
    for(i=0;i<n;i++){
        D_[i]= 1/this->get_element(i, i);
    }
    while(modul(b-(*this)*x)>r){
        for(i=1;i< this->line_indexes.size();i++) {
            t = 0;
            for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
                if(i-1!=column_indexes[z]) {
                    t += this->data[z] * x[this->column_indexes[z]];
                }
            }
            x_next[i-1]=D_[i-1]*(b[i-1]-t);
        }
        x=x_next;
    }
    return x;
}

std::vector<double> CSR::Gauss_Seidel(const std::vector<double> &x0, const std::vector<double> &b, double r) const {
    std::vector<double>x=x0;
    std::vector<double>D_(x0.size());
    double t;
    unsigned int i, z, n=x0.size();
    for(i=0;i<n;i++){
        D_[i]= 1/this->get_element(i, i);
    }
    while(modul(b-(*this)*x)>r){
        for(i=1;i< this->line_indexes.size();i++) {
            t = 0;
            for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
                if(i-1!=column_indexes[z]) {
                    t += this->data[z] * x[this->column_indexes[z]];
                }
            }
            x[i-1]=D_[i-1]*(b[i-1]-t);
        }
    }
    return x;
}

