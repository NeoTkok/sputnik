#ifndef TEST_CPP
#define TEST_CPP

#include <gtest/gtest.h>

#include "Class_sputnik.hpp"
#include "Orbit.hpp"

#include <cmath>


TEST(Orbit, Rp) {
    orbit A1(2,1/sqrt(2),0.,0.,0.);
    EXPECT_NEAR(2-sqrt(2), A1.getRp(), 1e-15);
}

TEST(Orbit, Ra) {
    orbit A1(2,1/sqrt(2),0.,0.,0.);
    EXPECT_NEAR(2+sqrt(2), A1.getRa(), 1e-15);
}

TEST(Orbit, p) {
    orbit A1(2,1/sqrt(2),0.,0.,0.);
    EXPECT_NEAR(A1.a * (1-A1.e * A1.e), A1.getP() ,1e-15);
}

TEST(Orbit, a_xyz1) {
    orbit A1(2, 1/sqrt(2), M_PI/3 , 0., M_PI/4);
    std::array<double, 3> r = {sqrt(6)/4, -sqrt(2)/4, sqrt(2)/2};
    for(int i = 0; i < 3; ++i)
        EXPECT_NEAR(r[i], A1.getPlane()[i] ,1e-15);
}

TEST(Orbit, a_xyz2) {
    orbit A1(2, 1/sqrt(2), M_PI , 0., M_PI);
    std::array<double, 3> r = {0, 0, -1};
    for(int i = 0; i < 3; ++i)
        EXPECT_NEAR(r[i], A1.getPlane()[i] ,1e-15);
}






int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



#endif