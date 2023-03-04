//
// Created by robert on 04.03.23.
//
#include "/home/robert/CLionProjects/SLAE/CSR.h"
#include <fstream>
using namespace std;

int main(){
    double a[]={10, 1, 0, 1, 7, 0, 0, 0.1, 1};
    CSR A(3, 3, a);
    std::vector<double>x0={0, 0, 0}, b={20, 30, 1};
    double tau=0.1;
    unsigned int i;
    double tau_data[21];
    unsigned int k_data[21];

    for(i=0;i<21;i++){
        tau_data[i]=tau;
        mixed_num_vec ans=A.simple_iteration(x0, b, tau, 0.000000000001);
        k_data[i]=ans.K;
        tau/=1.3;
    }

    std::ofstream out;
    out.open("/home/robert/CLionProjects/untitled/tau_data.txt");
    for(i=0;i<21;i++){
        out << tau_data[i] << " ";
    }
    out.close();
    out.open("/home/robert/CLionProjects/untitled/k_data.txt");
    for(i=0;i<21;i++){
        out << k_data[i] << " ";
    }
    out.close();
    return 0;
}