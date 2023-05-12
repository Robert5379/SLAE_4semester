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
    this->line_indexes.resize(n+1);
    this->data.resize(k);
    this->column_indexes.resize(k);
}

CSR::CSR(const CSR &A): data(A.data), column_indexes(A.column_indexes), line_indexes(A.line_indexes) {}

std::istream& operator>>(std::istream &in, CSR& A)
{
    unsigned int m, k=0, n;
    in >> n >> m;
    std::vector<double> data;
    unsigned int i, j;
    double t1;
    A.line_indexes.push_back(k);
    for(j=0;j<n;j++) {
        for (i = 0; i < m; i++) {
            in >> t1;
            if (t1 != 0) {
                A.data.push_back(t1);
                A.column_indexes.push_back(i);
                k++;
            }
        }
        A.line_indexes.push_back(k);
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

double operator*(const std::vector<double>& x, const std::vector<double>& v){
    double ans=0;
    unsigned int i;
    for(i=0;i<x.size();i++){
        ans+=x[i]*v[i];
    }
    return ans;

}

std::vector<double> operator/(const std::vector<double>& x, double a){
    std::vector<double> ans=x;
    unsigned int i;
    for(i=0;i<x.size();i++){
        ans[i]/=a;
    }
    return ans;
}

std::vector<double> operator+(const std::vector<double>& x, const std::vector<double>& v){
    std::vector<double> ans=x;
    for(unsigned int i=0;i<x.size();i++){
        ans[i]+=v[i];
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

std::ostream& operator<<(std::ostream & os, const std::vector<double>& v){
    unsigned int i;
    for(i=0;i<v.size();i++){
        os << v[i] << " ";
    }
    return os;
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


void CSR::transpose() {
    unsigned int N= this->data.size(), W= this->line_indexes.size(), i, S=0, t, j1, j2, RIndex, IIndex, j;
    std::vector<double> tVals(N);
    std::vector<unsigned int> tCols(N);
    std::vector<unsigned int> tRows(W);
    for(i=0;i<N;i++){
        tRows[this->column_indexes[i]+1]++;
    }

    for(i=1;i<W;i++){
        t=tRows[i];
        tRows[i]=S;
        S+=t;
    }

    for(i=0;i<W-1;i++){
        j1= this->line_indexes[i];
        j2= this->line_indexes[i+1];
        for(j=j1;j<j2;j++){
            RIndex= this->column_indexes[j];
            IIndex= tRows[RIndex+1];
            tVals[IIndex]= this->data[j];
            tCols[IIndex]=i;
            tRows[RIndex+1]++;
        }
    }

    this->data=tVals;
    this->column_indexes=tCols;
    this->line_indexes=tRows;
}


std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> CSR::Simple_iteration(const std::vector<double> & x0, const std::vector<double> &b, double tau, double r, unsigned int iterations) const {
    std::vector<double>x=x0;
    std::vector<double>r_n=(*this)*x-b;
    unsigned int k=0;
    std::vector<unsigned int>k_vector;
    std::vector<double> nev;
    double tolerance=modul(r_n);
    while(tolerance>r&&k<iterations){
        x=x-tau*(r_n);
        r_n=(*this)*x-b;
        k++;
        k_vector.push_back(k);
        tolerance= modul(r_n);
        nev.push_back(tolerance);
    }
    return std::make_pair(x, std::make_pair(nev, k_vector));
}

std::pair<std::vector<double>, unsigned int> CSR::Jacobi(const std::vector<double> &x0, const std::vector<double> &b, double r) const {
    std::vector<double>x=x0;
    std::vector<double>x_next(x0.size());
    double diag;
    double t;
    unsigned int i, z, n=x0.size(), k=0;
    while(modul(b-(*this)*x)>r){
        for(i=1;i< this->line_indexes.size();i++) {
            t = 0;
            for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
                if(i-1!=column_indexes[z]) {
                    t += this->data[z] * x[this->column_indexes[z]];
                } else{
                    diag=this->data[z];
                }
            }
            x_next[i-1]=(b[i-1]-t)/diag;
        }
        x=x_next;
        k++;
    }
    return std::make_pair(x, k);
}

std::pair<std::vector<double>, unsigned int> CSR::Gauss_Seidel(const std::vector<double> &x0, const std::vector<double> &b, double r) const {
    std::vector<double>x=x0;
    double diag;
    double t;
    unsigned int i, z, n=x0.size(), k=0;
    while(modul(b-(*this)*x)>r){
        for(i=1;i< this->line_indexes.size();i++) {
            t = 0;
            for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
                if(i-1!=column_indexes[z]) {
                    t += this->data[z] * x[this->column_indexes[z]];
                } else{
                    diag=this->data[z];
                }
            }
            x[i-1]=(b[i-1]-t)/diag;
        }
        k++;
    }
    return std::make_pair(x, k);
}

std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> CSR::SIM_Chebyshev_acceleration(const std::vector<double> &x0, const std::vector<double> &b,
                                                    double lambda_max, double lambda_min,  double accuracy, double degree, unsigned int iterations) const {
    unsigned int n=1, i, k=0;
    int j;
    for(i=0;i<degree;i++){
        n*=2;
    }
    double *tau=new double[n];
    tau[0]= cos(M_PI/(2*n));
    double sin_b=sin(M_PI/(2*n)), cos_a=cos(M_PI/n), sin_a=sin(M_PI/n), a1=(lambda_max-lambda_min)/2, b1=(lambda_max+lambda_min)/2;
    for(i=1;i<n;i++){
        tau[i]=tau[i-1]*cos_a-sin_b*sin_a;
        sin_b=tau[i]*sin_a+cos_a*sin_b;
    }
    for(i=1;i<n;i++){
        tau[i]=1/(a1*tau[i]+b1);
    }

    std::vector<double>x=x0;
    std::vector<double>accuracy_n=(*this)*x-b;
    unsigned int *tau_num=new unsigned int[n], n_i=1;
    tau_num[0]=0;
    tau_num[1]=1;
    for(i=0;i<degree;i++){
        j=n_i;
        n_i*=2;
        for(;j>=0;j--){
            tau_num[2*j+1]=n_i-1-tau_num[j];
            tau_num[2*j]=tau_num[j];
        }
    }
    std::vector<unsigned int>k_vector;
    std::vector<double> nev;
    double tolerance=modul(accuracy_n);
    while(tolerance>accuracy&&k<iterations){
        for(i=0;i<n;i++) {
            x = x - tau[tau_num[i]] * (accuracy_n);
            accuracy_n = (*this) * x - b;
            tolerance= modul(accuracy_n);
            k++;
            k_vector.push_back(k);
            nev.push_back(tolerance);
        }
    }
    delete [] tau;
    delete [] tau_num;
    return std::make_pair(x, std::make_pair(nev, k_vector));
}

std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> CSR::SOR(const std::vector<double> &x0, const std::vector<double> &b, double w,
                             double accuracy, unsigned int iterations) const {
    std::vector<double>x=x0;
    double diag;
    double t;
    unsigned int i, z, k=0;
    std::vector<unsigned int>k_vector;
    std::vector<double> nev;
    double tolerance=modul(b-(*this)*x);
    while(tolerance>accuracy&&k<iterations){
        for(i=1;i< this->line_indexes.size();i++) {
            t = 0;
            for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
                if(i-1!=column_indexes[z]) {
                    t += this->data[z] * x[this->column_indexes[z]];
                } else{
                    diag=this->data[z];
                }
            }
            x[i-1]=(w*b[i-1]-w*t-(w-1)*diag*x[i-1])/diag;
        }
        tolerance=modul(b-(*this)*x);
        k++;
        k_vector.push_back(k);
        nev.push_back(tolerance);
    }
    return std::make_pair(x, std::make_pair(nev, k_vector));
}

std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> CSR::Symmetrical_Gauss_Seidel(const std::vector<double> &x0, const std::vector<double> &b, double ro,
                                                  double accuracy, unsigned int iterations) const {
    std::vector<double>x_prev=x0;
    std::vector<double>x=x0;
    std::vector<double>vector_for_swap=x0;
    double diag, mu_prev=1, mu=1/ro, mu_next;
    unsigned int i, z, k=0;

    mu_next=2*mu/ro-mu_prev;
    // обычная итерация симметричного ГЗ (над вектором x)
    for(i=1;i< this->line_indexes.size();i++) {
        x[i-1]=b[i-1];
        for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
            if(i-1!=column_indexes[z]) {
                x[i-1] -= this->data[z] * x[this->column_indexes[z]];
            } else{
                diag=this->data[z];
            }
        }
        x[i-1]/=diag;
    }
    for(i=this->line_indexes.size()-1;i>0;i--) {
        x[i-1]=b[i-1];
        for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
            if(i-1!=column_indexes[z]) {
                x[i-1] -= this->data[z] * x[this->column_indexes[z]];
            } else{
                diag=this->data[z];
            }
        }
        x[i-1]/=diag;
    }
    //конец обычной итерации симметричного ГЗ

    std::vector<unsigned int>k_vector;
    std::vector<double> nev;
    double tolerance=modul(b-(*this)*x);

    while(tolerance>accuracy&&k<iterations){
        vector_for_swap=x;
        for(i=1;i< this->line_indexes.size();i++) {
            x[i-1]=b[i-1];
            for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
                if(i-1!=column_indexes[z]) {
                    x[i-1] -= this->data[z] * x[this->column_indexes[z]];
                } else{
                    diag=this->data[z];
                }
            }
            x[i-1]/=diag;
        }
        for(i=this->line_indexes.size()-1;i>0;i--) {
            x[i-1]=b[i-1];
            for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
                if(i-1!=column_indexes[z]) {
                    x[i-1] -= this->data[z] * x[this->column_indexes[z]];
                } else{
                    diag=this->data[z];
                }
            }
            x[i-1]/=diag;
        }
        for(i=0;i< this->data.size();i++){
            x[i]=2*mu*x[i]/(mu_next*ro)-mu_prev*x_prev[i]/mu_next;
        }
        x_prev=vector_for_swap;
        mu_prev=mu;
        mu=mu_next;
        mu_next=2*mu/ro-mu_prev;
        tolerance=modul(b-(*this)*x);
        k++;
        nev.push_back(tolerance);
        k_vector.push_back(k);
    }
    return std::make_pair(x, std::make_pair(nev, k_vector));
}

std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> CSR::SSOR(const std::vector<double> &x0, const std::vector<double> &b, double w, double ro,
                              double accuracy, unsigned int iterations) const {
    std::vector<double>x_prev=x0;
    std::vector<double>x=x0;
    std::vector<double>vector_for_swap=x0;
    double diag, mu_prev=1, mu=1/ro, mu_next, temp;
    unsigned int i, z, k=0;

    mu_next=2*mu/ro-mu_prev;
    // обычная итерация SSOR (над вектором x)
    for(i=1;i< this->line_indexes.size();i++) {
        temp=x[i-1];
        x[i-1]=b[i-1];
        for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
            if(i-1!=column_indexes[z]) {
                x[i-1] -= this->data[z] * x[this->column_indexes[z]];
            } else{
                diag=this->data[z];
            }
        }
        x[i-1]*=w;
        x[i-1]/=diag;
        x[i-1]+=(1-w)*temp;
    }
    for(i=this->line_indexes.size()-1;i>0;i--) {
        temp=x[i-1];
        x[i-1]=b[i-1];
        for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
            if(i-1!=column_indexes[z]) {
                x[i-1] -= this->data[z] * x[this->column_indexes[z]];
            } else{
                diag=this->data[z];
            }
        }
        x[i-1]*=w;
        x[i-1]/=diag;
        x[i-1]+=(1-w)*temp;
    }
    //конец обычной итерации симметричного ГЗ

    std::vector<unsigned int>k_vector;
    std::vector<double> nev;
    double tolerance=modul(b-(*this)*x);

    while(tolerance>accuracy&&k<iterations){
        vector_for_swap=x;
        for(i=1;i< this->line_indexes.size();i++) {
            temp=x[i-1];
            x[i-1]=b[i-1];
            for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
                if(i-1!=column_indexes[z]) {
                    x[i-1] -= this->data[z] * x[this->column_indexes[z]];
                } else{
                    diag=this->data[z];
                }
            }
            x[i-1]*=w;
            x[i-1]/=diag;
            x[i-1]+=(1-w)*temp;
        }
        for(i=this->line_indexes.size()-1;i>0;i--) {
            temp=x[i-1];
            x[i-1]=b[i-1];
            for (z = this->line_indexes[i-1]; z < this->line_indexes[i]; z++) {
                if(i-1!=column_indexes[z]) {
                    x[i-1] -= this->data[z] * x[this->column_indexes[z]];
                } else{
                    diag=this->data[z];
                }
            }
            x[i-1]*=w;
            x[i-1]/=diag;
            x[i-1]+=(1-w)*temp;
        }
        for(i=0;i< this->data.size();i++){
            x[i]=2*mu*x[i]/(mu_next*ro)-mu_prev*x_prev[i]/mu_next;
        }
        x_prev=vector_for_swap;
        mu_prev=mu;
        mu=mu_next;
        mu_next=2*mu/ro-mu_prev;
        tolerance=modul(b-(*this)*x);
        k++;
        nev.push_back(tolerance);
        k_vector.push_back(k);
    }
    return std::make_pair(x, std::make_pair(nev, k_vector));
}

std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> CSR::Steepest_descent(const std::vector<double> &x0, const std::vector<double> &b,
                                          double accuracy, unsigned int iterations) const {
    std::vector<double>x=x0;
    std::vector<double>r_i=(*this)*x-b;
    std::vector<double>Ar=(*this)*r_i;
    double alpha=r_i*r_i/(r_i*Ar);
    unsigned int k=0;
    std::vector<unsigned int>k_vector;
    std::vector<double> nev;
    double tolerance=modul(r_i);
    while(tolerance>accuracy&&k<iterations){
        x=x-alpha*(r_i);
        r_i=r_i-alpha*Ar;
        Ar=(*this)*r_i;
        alpha=r_i*r_i/(r_i*Ar);
        tolerance=modul(r_i);
        k++;
        k_vector.push_back(k);
        nev.push_back(tolerance);
    }
    return std::make_pair(x, std::make_pair(nev, k_vector));
}

std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> CSR::Heavy_ball(const std::vector<double> &x0, const std::vector<double> &b,
                                     double accuracy, unsigned int iterations) const {
    std::vector<double>x_prev=x0;
    std::vector<double>x=x0;
    std::vector<double>Ad(x0.size());
    std::vector<double>r_i=(*this)*x-b;
    std::vector<double>delta(x0.size());

    double alpha=(r_i * r_i) / (r_i * ((*this) * r_i)), beta, rAr, rAd, rr;
    x = x - alpha * r_i;
    unsigned int k=0;
    std::vector<unsigned int>k_vector;
    std::vector<double> nev;
    double tolerance=modul(r_i);
    while(tolerance>accuracy&&k<iterations){
        delta = x - x_prev;
        x_prev=x;
        Ad=(*this) * delta;
        rAd = r_i * Ad;
        rAr = r_i * ((*this) * r_i);
        rr = r_i * r_i;
        beta = (rr * rAd - r_i * delta * rAr) / (delta * Ad* rAr - rAd * rAd);
        alpha = (rr + beta * rAd) / rAr;
        x = x - alpha * r_i + beta * delta;
        r_i = (*this) * x - b;
        tolerance=modul(r_i);
        k++;
        k_vector.push_back(k);
        nev.push_back(tolerance);
    }
    return std::make_pair(x, std::make_pair(nev, k_vector));
}

std::pair<std::vector<double>, std::pair<std::vector<double>,std::vector<unsigned int>>> CSR::Conjugate_gradient(const std::vector<double> &x0, const std::vector<double> &b,
                                            double accuracy, unsigned int iterations) const {
    std::vector<double>x=x0;
    std::vector<double>r_i=(*this)*x-b;
    std::vector<double>r_prev=r_i;
    std::vector<double>d=r_i;

    double alpha, dr=d*r_i;
    alpha=dr/(d*((*this)*d));
    unsigned int k=0;
    std::vector<unsigned int>k_vector;
    std::vector<double> nev;
    double tolerance=modul(r_i);
    while(tolerance>accuracy&&k<iterations){
        x=x-alpha*d;
        r_prev=r_i;
        r_i=(*this)*x-b;
        d=r_i+(r_i*r_i/(d*r_prev))*d;
        dr=d*r_i;
        alpha=dr/(d*((*this)*d));
        tolerance=modul(r_i);
        k++;
        k_vector.push_back(k);
        nev.push_back(tolerance);
    }
    return std::make_pair(x, std::make_pair(nev, k_vector));
}
/*
std::pair<std::vector<double>, std::pair<std::vector<double>,unsigned int>> CSR::Cholesky_CG(
        const std::vector<double> &x0, const std::vector<double> &b, double accuracy, unsigned int iterations) const {
    std::vector<double>x=x0;
    std::vector<double>r_i=(*this)*x-b;
    std::vector<double>r_prev=r_i;
    std::vector<double>d=r_i;

    double alpha, dr=d*r_i;
    alpha=dr/(d*((*this)*d));
    unsigned int k=0;
    std::vector<double> nev;
    double tolerance=modul(r_i);
    while(tolerance>accuracy&&k<iterations){
        x=x-alpha*d;
        r_prev=r_i;
        r_i=(*this)*x-b;
        d=r_i+(r_i*r_i/(d*r_prev))*d;
        dr=d*r_i;
        alpha=dr/(d*((*this)*d));
        tolerance=modul(r_i);
        k++;
        nev.push_back(tolerance);
    }
}
*/
std::pair<std::vector<double>, std::pair<std::vector<double>,unsigned int>> CSR::GMRES(const std::vector<double> &x0,
                                                                                       const std::vector<double> &b,
                                                                                       double accuracy) const {
    std::vector<double> x = x0;
    unsigned int n = this->line_indexes.size()-1;
    std::vector<std::vector<double>> v(n+1, std::vector<double>(n));
    Complete_matrix H(n+1, n);
    std::vector<std::pair<double, double>> SinCos(n);
    std::vector<double> r_0 = (*this) * x0 - b;
    std::vector<double> nev{std::sqrt(r_0 * r_0)};
    r_0.clear();
    unsigned int iterations_number=0, i, j, k;
    double t, t1, t2;
    long int i1, j1;
    std::vector<double> e(n+1);
    std::vector<double> y(n);
    while (nev[iterations_number] > accuracy) {
        std::vector<double> r = (*this) * x - b;
        e[0] = std::sqrt(r * r);
        v[0] = r / e[0];
        for (i = 0; i < n; i++) {
            v[i+1]=(*this)*v[i];
            for(j=0;j<i+1;j++){
                H.write(v[j]*v[i+1], j, i);
                v[i+1]=v[i+1]-H.get_element(j, i)*v[j];
            }
            H.write(std::sqrt(v[i+1]*v[i+1]), i+1, i);

            for(j=0;j<i;j++){
                t=H.get_element(j, i);
                H.write(SinCos[j].second * H.get_element(j, i) + SinCos[j].first * H.get_element(j+1, i), j, i);
                H.write(-SinCos[j].first * t + SinCos[j].second * H.get_element(j+1, i), j, i);
            }
            t1=H.get_element(i, i);
            t2=H.get_element(i + 1, i);
            double t = std::sqrt(t2 * t2 + t1 * t1);
            SinCos[i].first = t2/ t;
            SinCos[i].second = t1 / t;
            H.write(SinCos[i].second * t1 + SinCos[i].first * t2, i, i);
            H.write(0, i+1, i);

            t=e[i];
            e[i] = SinCos[i].second * e[i];
            e[i + 1] = -SinCos[i].first * t;
            nev.push_back(std::abs(e[i + 1]));
            iterations_number++;
        }

        for ( i1 = n-1; i1 >= 0; i1--) {
            y[i1] = e[i1];
            for ( j1= n-1; j1 > i1; j1--) {
                y[i1] -= H.get_element(i1, j1) * y[j1];
            }
            y[i1] /= H.get_element(i1, i1);
            x = x - y[i1] * v[i1];
        }
    }

    return std::make_pair(x, std::make_pair(nev, iterations_number));
}

std::pair<std::vector<double>, std::pair<std::vector<double>,unsigned int>> CSR::BiCG(const std::vector<double> &x0,
                                                                                      const std::vector<double> &b,
                                                                                      double accuracy) const {
    CSR A_T(*this);
    A_T.transpose();
    std::vector<double> x = x0;
    std::vector<double> r = (*this) * x - b;
    std::vector<double> nev{std::sqrt(r * r)};
    std::vector<double> r_wave, p, p_wave;
    std::vector<double> Ap;
    unsigned int i, iterations_number=0, j, k;
    double q, t, rr, rr_prev;
    while(nev[iterations_number] > accuracy){
        r = (*this) * x - b;
        r_wave=r;
        p=r;
        p_wave=r;
        t=0;
        rr_prev=1;
        for(i=0;i<this->line_indexes.size();i++){
            rr = r * r_wave;
            if (rr == 0) break;
            t = rr / rr_prev;
            p = r + t * p;
            p_wave = r_wave + t * p_wave;
            rr_prev = rr;
            Ap = (*this) * p;
            q = rr / (r_wave * Ap);
            x = x - q * p;
            r = r - q * Ap;
            r_wave = r_wave - q * (A_T * p_wave);
            nev.push_back(std::sqrt(r * r));
            iterations_number++;
        }
    }
    return std::make_pair(x, std::make_pair(nev, iterations_number));
}
