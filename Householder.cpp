//
// Created by robert on 26.02.23.
//

#include "Householder.h"


double scalar(const std::shared_ptr<double[]> v, const std::shared_ptr<double[]> x, unsigned int size) {
    unsigned int i;

    double ans=0;
    for(i=0;i<size;i++){
        ans+=v[i]*x[i];
    }
    return ans;
}

double scalar(const std::shared_ptr<double[]> v, const double* x, unsigned int size) {
    unsigned int i;

    double ans=0;
    for(i=0;i<size;i++){
        ans+=v[i]*x[i];
    }
    return ans;
}

double scalar(const double* v, const double* x, unsigned int size) {
    unsigned int i;

    double ans=0;
    for(i=0;i<size;i++){
        ans+=v[i]*x[i];
    }
    return ans;
}

std::unique_ptr<double[]> count_v();

void do_column_R(double*v, std::shared_ptr<double[]> x, unsigned int size, double vv){
    unsigned int i;
    double alpha=2* scalar(x, v, size)/ vv;
    for(i=0;i<size;i++){
        x[i]-=alpha*v[i];
    }

}

std::unique_ptr<Complete_matrix[]> Householder(const Complete_matrix& A){
    Complete_matrix R(A);
    unsigned int i=0, j, n=A.get_column_size();
    double* v1=new double[n];
    std::shared_ptr<double[]> x=A.get_column(0);
    double modul=0, alpha, vv;
    for(i=0;i<A.get_column_size();i++){
        modul+=x[i]*x[i];
        v1[i]=x[i];
    }

    if(v1[0]>=0) {
        v1[0] += sqrt(modul);
    }else{
        v1[0] -= sqrt(modul);
    }
    modul= scalar(v1, v1, n);
    double* q=new double[n*n];
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                q[i*n+j]=1;
            }else{
                q[i*n+j]=0;
            }
            q[i*n+j]-=2*v1[i]*v1[j]/modul;
        }
    }
    Complete_matrix Q(n, n, q);
    Q.transpose();
    delete[] q;
    //do matrix R
    for(i=0;i<n;i++){
        v1=new double[n-i];
        x=R.get_column(i, i);
        for(j=0;j<n-i;j++){
            v1[j]=x[j];
        }
        if(v1[0]>=0) {
            v1[0] += sqrt(scalar(x, x, n - i));
        }else{
            v1[0] -= sqrt(scalar(x, x, n - i));
        }
        vv=scalar(v1, v1, n-i);
        alpha=2* scalar(x, v1, n-i)/ vv;
        R.write(x[0]-alpha*v1[0], i, i);
        for(j=i+1;j<n;j++){
            R.write(0, j, i);
        }
        for(j=i+1;j<n;j++){
            x=R.get_column(j, i);
            do_column_R(v1, x, n-i, vv);
            R.write_column(x, j, i);
        }
        if(i>0){
            for(j=0;j<n;j++){
                x=Q.get_column(j, i);
                do_column_R(v1, x, n-i, vv);
                Q.write_column(x, j, i);
            }
        }

    }
    delete[] v1;
    Q.transpose();
    std::unique_ptr<Complete_matrix[]> ans=std::make_unique<Complete_matrix[]>(2);
    ans[0]=Q;
    ans[1]=R;
    return ans;
};
