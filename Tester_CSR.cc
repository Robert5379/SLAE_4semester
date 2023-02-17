#include <gtest/gtest.h>
#include "CSR.cpp"

TEST(CSRTest, CSRCorrectness1) {
    double a[]={1, 2, 3, 4};
    CSR A(2, 2, a);
    EXPECT_EQ(3, A.get_element(1, 0)) << "Wrong element";
}

TEST(CSRTest, CSRCorrectness2) {
    double a[]={1, 2, 3, 4};
    CSR A(2, 2, a);
    std::vector<double> v={1, 2}, *new_v=A*v, test={5, 11};

    EXPECT_EQ(test, *new_v) << "Wrong * №2";
}

TEST(CSRTest, CSRCorrectness3) {
    double a[]={1, 2, 3, 4, 5, 6, 7, 8, 5};
    CSR A(3, 3, a);
    std::vector<double> v={1, 2, 1}, *new_v=A*v, test={8, 20, 28};

    EXPECT_EQ(test, *new_v) << "Wrong * №3";
}

TEST(CSRTest, CSRCorrectness4) {
    double a[]={1, 2, 3, 4, 5, 6, 7, 8, 5};
    CSR A(3, 3, a);
    std::vector<double> v={0, 0, 0}, *new_v=A*v, test={0, 0, 0};
    EXPECT_EQ(test, *new_v) << "Wrong * №4";
}

TEST(CSRTest, CSRCorrectness5) {
    double a[]={0, 0, 0, 0, 0, 0, 0, 0, 0};
    CSR A(3, 3, a);
    std::vector<double> v={9, 7, 56}, *new_v=A*v, test={0, 0, 0};
    EXPECT_EQ(test, *new_v) << "Wrong * №5";
}