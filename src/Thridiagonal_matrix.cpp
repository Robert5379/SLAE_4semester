//#include "Triplex.cpp"
//#include <iostream>
//#include <vector>

#include "Thridiagonal_matrix.h"

//template <typename T=double>

void tridiagonal_matrix::do_matrix(unsigned int N, double *a) {
        n = N;
        unsigned int i = 0;
        triplex tr;
        tr.b = a[0];
        tr.c = a[1];
        data.push_back(tr);
        for (i = 2; i < 3 * n - 6; i = i + 3) {
            tr.a = a[i];
            tr.b = a[i + 1];
            tr.c = a[i + 2];
            data.push_back(tr);
        }
        tr.a = a[3 * n - 4];
        tr.b = a[3 * n - 3];
        data.push_back(tr);
    }



std::istream& operator >> (std::istream& in, tridiagonal_matrix& p)
{
    in >> p.n;
    unsigned int i, j;
    triplex tr;
    double t1, t2;
    in >> t1>>t2;
    tr.b=t1;
    tr.c=t2;
    p.data.push_back(tr);
    for(i=1;i<p.n-1;i++){
        in >> tr.a >>tr.b >> tr.c;
        p.data.push_back(tr);
    }
    in >> tr.a >>tr.b;

    p.data.push_back(tr);
    return in;
}

