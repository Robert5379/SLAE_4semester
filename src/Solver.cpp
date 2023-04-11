
#include "Solver.h"

template <typename T>
std::vector<T> solver(tridiagonal_matrix &A, std::vector<T> &f){
    unsigned int i, n=A.n;
    double p[n];
    double q[n];
    std::vector<T> x(A.n);
    p[1]=-A.data[0].c/A.data[0].b;
    q[1]=f[0]/A.data[0].b;

    for(i=2;i<n;i++){
        p[i]=-A.data[i-1].c/(A.data[i-1].b+A.data[i-1].a*p[i-1]);
        q[i]=(f[i-1]-A.data[i-1].a*q[i-1])/(A.data[i-1].b+A.data[i-1].a*p[i-1]);
    }
    x[A.n-1]=(f[A.n-1]-A.data[A.n-1].a*q[n-1])/(A.data[n-1].a*p[n-1]+A.data[n-1].b);

    for(i=n-1;i>0;i--){
        x[i-1]=p[i]*x[i]+q[i];
    }
    return x;
}

template std::vector<double> solver(tridiagonal_matrix &A, std::vector<double> &f);
template std::vector<float> solver(tridiagonal_matrix &A, std::vector<float> &f);