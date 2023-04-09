//
// Created by robert on 08.04.2023.
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
    std::vector<double> x_MPI, x0(289), b(289), x_MPI_best, x_Sym_GS;
    std::vector<double> Chebyshev;
    double lambda_max=8+4*cos(M_PI/18), lambda_min=8+4*cos(17*M_PI/18);
    for(i=0;i<289;i++){
        x0[i]=0;
        b[i]=1;
    }
    std::ofstream out;
    out.open("/home/robert/CLionProjects/untitled/Control1/Task1.txt");
    double accuracy=1e-13;

    x_MPI= (get<0>(A.Simple_iteration(x0, b, 1/lambda_max, accuracy)));
    out << "Simple_iteration1" << endl;
    for(i=0;i<289;i++){
        out << x_MPI[i] << " ";
    }
    out << endl;

    x_MPI_best=get<0>(A.Simple_iteration(x0, b, 2/(lambda_max+lambda_min), accuracy));
    out << "Gradient Descent with Optimal Parameter" << endl;
    for(i=0;i<289;i++){
        out << x_MPI_best[i] << " ";
    }
    out << endl;

    Chebyshev=get<0>(A.SIM_Chebyshev_acceleration(x0, b, lambda_max, lambda_min, accuracy, 3));
    out << "SIM_Chebyshev_acceleration" << endl;
    for(i=0;i<289;i++){
        out << Chebyshev[i] << " ";
    }
    out << endl;

    out.close();
    return 0;
}