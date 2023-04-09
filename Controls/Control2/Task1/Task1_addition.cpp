//
// Created by robert on 08.04.2023.
//
#include "src/CSR.h"
#include <fstream>
using namespace std;

int main(){
    double a[83521];
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


    std::vector<double> x0(289), b(289), v2;

    for(i=0;i<289;i++){
        x0[i]=0;
        b[i]=1;
    }
    std::ofstream out;
    out.open("/home/robert/CLionProjects/untitled/Control1/Task1.txt", std::ios::app);
    double accuracy=1e-13;

    v2=get<0>(A.Symmetrical_Gauss_Seidel(x0, b, 0.5, accuracy));
    out << "Symmetrical_Gauss_Seidel" << endl;
    for(i=0;i<289;i++){
        out << v2[i] << " ";
    }

    out.close();
    return 0;
}