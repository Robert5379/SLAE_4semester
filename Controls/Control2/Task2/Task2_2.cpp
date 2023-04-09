//
// Created by robert on 09.04.2023.
//
#include "src/CSR.h"
#include <fstream>
using namespace std;

int main(){
    double a[]={1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1.5};
    CSR A(4, 4, a);
    std::vector<double> x0={0, 0, 0, 0}, b={3, 3, 3, 3}, b1={1.5, 1.5, 1.5, 1.5}, x1;
    double accuracy=1e-13, lambda_max=1.5, lambda_min=1, c=4, nev_prev;
    unsigned int i;
    std::ofstream out;
    out.open("/home/robert/CLionProjects/untitled/Control1/Task2_2.txt");
    auto x=(A.Simple_iteration(x0, b1, 0.9*2/lambda_max, accuracy));
    i=1;
    while(i<100) {
        nev_prev=get<0>(get<1>(x))[i-1];
        x = (A.Simple_iteration(x0, b1, 0.9 * 2 / lambda_max, accuracy, i));
        out << get<0>(x)[0] << " " << get<0>(x)[3] << " ";
        i++;
    }
    out << endl;
    i=1;
    while(i<100) {
        nev_prev=get<0>(get<1>(x))[i-1];
        x = (A.Simple_iteration(x0, b1, 2 / (lambda_max+lambda_min), accuracy, i));
        out << get<0>(x)[0] << " " << get<0>(x)[3] << " ";
        i++;
    }
    out << endl;
    i=1;
    while(i<100) {
        nev_prev=get<0>(get<1>(x))[i-1];
        x = (A.Steepest_descent(x0, b1,  accuracy, i));
        out << get<0>(x)[0] << " " << get<0>(x)[3] << " ";
        i++;
    }
    out << endl;
    i=1;
    while(i<100) {
        nev_prev=get<0>(get<1>(x))[i-1];
        x = (A.SIM_Chebyshev_acceleration(x0, b1, lambda_max, lambda_min, accuracy, 3, i));
        out << get<0>(x)[0] << " " << get<0>(x)[3] << " ";
        i++;
    }
    out << endl;
    i=1;
    while(i<100) {
        nev_prev=get<0>(get<1>(x))[i-1];
        x = (A.Conjugate_gradient(x0, b1, accuracy, i));
        out << get<0>(x)[0] << " " << get<0>(x)[3] << " ";
        i++;
    }
    out.close();
    return 0;
}