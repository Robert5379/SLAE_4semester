#include <gtest/gtest.h>
#include "../src/Solver.cpp"

TEST(SolverTest, SolverCorrectness1) {
    std::vector<double> X(3), f(3), test(3);
    f={3, 6, 2};
    X[0]=1.49;
    X[1]=-0.02;
    X[2]=-0.68;
    double a[7]={2, -1, 5, 4, 2, 1, -3};

    tridiagonal_matrix A;
    A.do_matrix(3, a);
    test=solver(A, f);

    for (int i = 0; i <3;++i) {
        EXPECT_TRUE(abs(X[i]-test[i])<0.01) << "Vectors X and test differ at index " << i;
    }

}

TEST(SolverTest, SolverCorrectness2) {
    unsigned int const n=5;
    std::vector<double> X(n), f(n), test(n);
    f={-25, 72, -69, -156, 20};
    X={-10, 5, -2, -10, -3};

    double a[3*n-2]={2, -1, -3, 8, -1, -5, 12, 2, -6, 18, -4, -5, 10};

    tridiagonal_matrix A;
    A.do_matrix(n, a);
    test=solver(A, f);

    for (int i = 0; i <3;++i) {
        EXPECT_TRUE(abs(X[i]-test[i])<0.01) << "Vectors X and test differ at index " << i;
    }

}

TEST(SolverTest, SolverCorrectness3) {
    unsigned int const n=4;
    std::vector<double> X(n), f(n), test(n);
    f={9, 11, 14, 13};
    X={3.65, -0.49, 2.78, 0.32};

    double a[3*n-2]={3, 4, 1, 2, 3, 23, 9, 1, 4, 6};

    tridiagonal_matrix A;
    A.do_matrix(n, a);
    test=solver(A, f);

    for (int i = 0; i <3;++i) {
        EXPECT_TRUE(abs(X[i]-test[i])<0.01) << "Vectors X and test differ at index " << i;
    }

}