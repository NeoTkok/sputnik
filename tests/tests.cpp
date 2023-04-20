#ifndef TEST_CPP
#define TEST_CPP

#include <gtest/gtest.h>

#include "Sputnik.hpp"
#include "Orbit.hpp"

#include <cmath>


const double GM = 3.98600441888888 * 1e14;  

TEST(Orbit, Rp) {
    orbit A(2,1/sqrt(2),0.,0.,0.);
    EXPECT_NEAR(2-sqrt(2), A.getRp(), 1e-15);
}

TEST(Orbit, Ra) {
    orbit A(2,1/sqrt(2),0.,0.,0.);
    EXPECT_NEAR(2+sqrt(2), A.getRa(), 1e-15);
}

TEST(Orbit, P) {
    orbit A(2,1/sqrt(2),0.,0.,0.);
    EXPECT_NEAR(A.getRp()*(1 + A.e), A.getP() ,1e-15);
}

TEST(Orbit, a_xyz1) {
    orbit A(2, 1/sqrt(2), M_PI/3 , 0., M_PI/4);
    std::array<double, 3> r = {sqrt(6)/4, -sqrt(2)/4, sqrt(2)/2};
    for(int i = 0; i < 3; ++i)
        EXPECT_NEAR(r[i], A.getPlane()[i] ,1e-15);
}

TEST(Orbit, a_xyz2) {
    orbit A(2, 1/sqrt(2), M_PI , 0., M_PI);
    std::array<double, 3> r = {0, 0, -1};
    for(int i = 0; i < 3; ++i)
        EXPECT_NEAR(r[i], A.getPlane()[i] ,1e-15);
}


TEST(Orbit, peresech_1){
    orbit A(2, 1/sqrt(2), M_PI/4 , 0., M_PI / 4);
    orbit B(2, 1/sqrt(2), M_PI/4 , 0., 3 * M_PI / 4);
    std::array<double, 3> r = {1/sqrt(2), 1/sqrt(2), 0};
    for(int i = 0; i < 3; ++i)
        EXPECT_NEAR(r[i], peresech(A, B)[i] ,1e-15);
}

TEST(Orbit, peresech_2){
    orbit A(2, 1/sqrt(2), 0 , 0., M_PI / 4);
    orbit B(2, 1/sqrt(2), M_PI/2 , 0., 3 * M_PI / 4);
    std::array<double, 3> r = {1/sqrt(3), 1/sqrt(3), 1/sqrt(3)};
    for(int i = 0; i < 3; ++i)
        EXPECT_NEAR(r[i], peresech(A, B)[i] ,1e-15);
}


TEST(Sputnik, h) {
    orbit A(2, 1/sqrt(2), M_PI , 0., M_PI);
    sputnik S(A, GM);   
    EXPECT_NEAR(S.get_h(), - GM * (1 - A.e * A.e) / A.getP() ,1e-15);
}


TEST(Sputnik, vr) {
    orbit A(2, 1/sqrt(2), M_PI , 0., M_PI);
    sputnik S(A,GM);
    EXPECT_NEAR(0., S.v_r(0) ,1e-15);
}

TEST(Sputnik, vtau) {
    orbit A(2, 1/sqrt(2), M_PI , 0., M_PI);
    sputnik S(A,GM);
    EXPECT_NEAR(S.v_tau(-M_PI/4), S.v_tau(M_PI/4) ,1e-15);
}


TEST(Orbit, f) {
    orbit A(10., 1./2., 0., 0., 0.);
    orbit B(10., 1./2., 0., -M_PI / 2., 0.);
    EXPECT_NEAR(-B.getP()+A.getRp(), delta_f(A,B,0) ,1e-15);
    EXPECT_NEAR(A.getRa()-B.getP(), delta_f(A,B,M_PI) ,1e-15);
}


TEST(Orbit, null_f) {
    orbit A(10., 1./2., 0., 0., 0.);
    orbit B(10., 1./2., 0., -M_PI / 2., 0.);
    std::vector<double> Y =  null_f(A,B);
    EXPECT_NEAR(M_PI*3/4, exact_null_f(A, B, Y[1], Y[2], 1e-13) ,1e-12);
    EXPECT_NEAR(M_PI*7/4, exact_null_f(A, B, Y[5], Y[6], 1e-13) ,1e-12);
}


// аналитичекски протестировал df; null_df; exact_null_df
/*
TEST(Sputnik, sign) {
    orbit A(10., 1./2., 0., 0., 0.);
    orbit B(10., 1./2., 0., -M_PI / 2., 0.);

    std::vector<double> Z = null_df(A, B);

    for(int i = 0; i < Z.size(); ++i)
        std::cout << i << " - " << Z[i] << "\n";
    
//    std::cout << std::setprecision(15) << exact_null_df(A, B, Z[1], Z[2], 1e-13) -  3*M_PI/4   << std::endl;
}
*/




int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



#endif