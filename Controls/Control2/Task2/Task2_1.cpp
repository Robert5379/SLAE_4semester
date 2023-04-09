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
    double accuracy=1e-13, lambda_max=1.5, lambda_min=1, c=4;
    unsigned int i;
    std::ofstream out;
    out.open("/home/robert/CLionProjects/untitled/Control1/Task2_1.txt");
    out << "SIM" << endl;
    out << "x_min" << endl;
    x1=get<0>(A.Simple_iteration(x0, b1, 0.9*2/lambda_max, accuracy));
    out << x1 << endl;
    out << "f_min" << endl;
    out << x1*(A*x1)-b*x1+c << endl;
    out << "SIM_optimal" << endl;
    out << "x_min" << endl;
    x1=get<0>(A.Simple_iteration(x0, b1, 2/(lambda_max+lambda_min), accuracy));
    out << x1 << endl;
    out << "f_min" << endl;
    out << x1*(A*x1)-b*x1+c << endl;
    out << "Steepest_descent" << endl;
    out << "x_min" << endl;
    x1=get<0>(A.Steepest_descent(x0, b1,  accuracy));
    out << x1 << endl;
    out << "f_min" << endl;
    out << x1*(A*x1)-b*x1+c << endl;
    out << "SIM_Chebyshev" << endl;
    out << "x_min" << endl;
    x1=get<0>(A.SIM_Chebyshev_acceleration(x0, b1, lambda_max, lambda_min, accuracy, 3));
    out << x1 << endl;
    out << "f_min" << endl;
    out << x1*(A*x1)-b*x1+c << endl;
    out << "Conjugate_gradient" << endl;
    out << "x_min" << endl;
    x1=get<0>(A.Conjugate_gradient(x0, b1, accuracy));
    out << x1 << endl;
    out << "f_min" << endl;
    out << x1*(A*x1)-b*x1+c << endl;
    return 0;
}