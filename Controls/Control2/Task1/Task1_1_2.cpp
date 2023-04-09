//
// Created by robert on 09.04.2023.
//
#include "src/CSR.h"
#include <fstream>
using namespace std;

int main(){
    double *a=new double[83521];
    unsigned int i,j, n=289, l=17;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                a[n*i+j]=8;
            }else{
                if(i<n-1&&j==i+1){
                    a[n*i+j]=1;
                }else{
                    if(i>0&&j==i-1){
                        a[n*i+j]=1;
                    }else{
                        if(i<n-l&&j==i+l){
                            a[n*i+j]=1;
                        }else{
                            if(i>l-1&&j==i-l){
                                a[n*i+j]=1;
                            }else{
                                a[n*i+j]=0;
                            };
                        }
                    }
                }
            }
        }
    }

    CSR A(n, n, a);
    delete[] a;
    std::vector<double> x0(289), b(289), x;
    std::vector<unsigned int> k;
    double lambda_max=8+4*cos(M_PI/18), lambda_min=8+4*cos(17*M_PI/18);
    for(i=0;i<289;i++){
        x0[i]=0;
        b[i]=1;
    }
    std::ofstream out;
    out.open("/home/robert/CLionProjects/untitled/Control1/Task1_1_2.txt");
    double accuracy=1e-13;

    auto x_MPI= (get<1>(A.Simple_iteration(x0, b, 2/(lambda_max+lambda_min), accuracy)));
    x=get<0>(x_MPI);
    k=get<1>(x_MPI);
    for(i=0;i<x.size();i++){
        out << x[i] << " ";
    }
    out << endl;
    for(i=0;i<k.size();i++){
        out << k[i] << " ";
    }

    out.close();
    return 0;
}