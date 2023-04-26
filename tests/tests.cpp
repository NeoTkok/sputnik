#ifndef TEST_CPP
#define TEST_CPP

#include <gtest/gtest.h>

#include "Sputnik.hpp"
#include "Orbit.hpp"

#include <cmath>


const double GM = 3.98600441888888 * 1e14;  

//************* проверяю правильность вывода расстояния до перицентра *********************
//
TEST(Orbit, Rp) {
    orbit A(1.,1./sqrt(2),0.,0.,0.);
    ASSERT_NEAR(2-sqrt(2), A.getRp(), 1e-15);
}

//************* проверю правильность вывода расстояния до апоцентра ***********************
TEST(Orbit, Ra) {
    orbit A(1.,1./sqrt(2),0.,0.,0.);
    ASSERT_NEAR(2+sqrt(2), A.getRa(), 1e-15);
}


//***************определение уравнения плоскости (Ах+Ву+Сz+D = 0, D = 0)*******************
// 
TEST(Orbit, a_xyz1) {
    orbit A(1., 1/sqrt(2), M_PI/3 , 0., M_PI/4);
    std::array<double, 3> r = {sqrt(6)/4, -sqrt(2)/4, sqrt(2)/2};
    for(int i = 0; i < 3; ++i)
        ASSERT_NEAR(r[i], A.getPlane()[i] ,1e-15);
}

TEST(Orbit, a_xyz2) {
    orbit A(1., 1/sqrt(2), M_PI , 0., M_PI);
    std::array<double, 3> r = {0, 0, -1};
    for(int i = 0; i < 3; ++i)
        ASSERT_NEAR(r[i], A.getPlane()[i] ,1e-15);
}

//********* пересечение плоскостей. получение нормированного направляющего вектора ************
//************** при совпадении плоскостей получаем (0 0 0) ***********************************
TEST(Orbit, peresech_1){
    orbit A(1., 1/sqrt(2), M_PI/4 , 0., M_PI / 4);
    orbit B(1., 1/sqrt(2), M_PI/4 , 0., 3 * M_PI / 4);
    std::array<double, 3> r = {1/sqrt(2), 1/sqrt(2), 0};
    for(int i = 0; i < 3; ++i)
        ASSERT_NEAR(r[i], peresech(A, B)[i] ,1e-15);
}

TEST(Orbit, peresech_2){
    orbit A(1., 1/sqrt(2), 0 , 0., M_PI / 4);
    orbit B(1., 1/sqrt(2), M_PI/2 , 0., 3 * M_PI / 4);
    std::array<double, 3> r = {1/sqrt(3), 1/sqrt(3), 1/sqrt(3)};
    for(int i = 0; i < 3; ++i)
        ASSERT_NEAR(r[i], peresech(A, B)[i] ,1e-15);
}

// наложение двух плоскостей))
TEST(Orbit, peresech_3){
    orbit A(1., 1/sqrt(2), 0 , 0., 0);
    orbit B(1., 1/sqrt(2), M_PI/2 , 0., M_PI);
    std::array<double, 3> r = {0., 0., 0.};
    for(int i = 0; i < 3; ++i)
        ASSERT_NEAR(r[i], peresech(A, B)[i] ,1e-13);
}

// *********** тестирование функции разности значений при заданном угле в Rp и Ra ********** 
TEST(Orbit, f) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., -M_PI / 2., 0.);
    ASSERT_NEAR(-B.p+A.getRp(), delta_f(A,B,0) ,1e-15);
    ASSERT_NEAR(A.getRa()-B.p, delta_f(A,B,M_PI) ,1e-15);
}


// ************* определение точных точек на отрезке **********
TEST(Orbit, exact_null_f_on_intersection) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., -M_PI / 2., 0.);
    std::vector<double> Y = intersection(A, B, 1e-12);

    ASSERT_NEAR(M_PI*3/4, Y[0] ,1e-13);
    ASSERT_NEAR(M_PI*7/4, Y[1] ,1e-13);
}

//  -//- наложение  
TEST(Orbit, intersection_2) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., 0., 0.);

    std::vector<double> Z = intersection(A, B, 1e-12); // получение extr
    
    ASSERT_LT(2, Z.size());
}   

// *********** тестирование разности значений df при заданном угле в Rp и Ra ********** 
TEST(Orbit, df) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., -M_PI / 2., 0.);
    ASSERT_LT(0, delta_df(A,B,M_PI));
    ASSERT_GT(0, delta_df(A,B,0));
}

// **** определение отрезков на которых df меняет знак, т. е отрезки где существует нуль ******
// null_df - вектор {-1, phi_1, phi_2, -1, -1, phi_3, phi_4, -1}
TEST(Orbit, null_df) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., -M_PI / 2., 0.);
    ASSERT_LT(delta_df(A, B, null_df(A, B)[1] - 1e-10), delta_df(A,B, null_df(A, B)[1]));
    ASSERT_GT(delta_df(A, B, null_df(A, B)[2] + 1e-10), delta_df(A,B, null_df(A, B)[2]));
    ASSERT_GT(delta_df(A, B, null_df(A, B)[5] - 1e-10), delta_df(A,B, null_df(A, B)[5]));
    ASSERT_LT(delta_df(A, B, null_df(A, B)[6] + 1e-10), delta_df(A,B, null_df(A, B)[6]));
}

// ************* определение точных точек на отрезке где f меняет знак (точный нуль) **********
TEST(Orbit, exact_null_df) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., -M_PI / 2., 0.);
    std::vector<double> Y =  null_df(A,B);
    ASSERT_GT(delta_df(A, B, exact_null_df(A, B, Y[1], Y[2], 1e-13) + 1e-10), delta_df(A, B, exact_null_df(A, B, Y[1], Y[2], 1e-13)));
    ASSERT_LT(delta_df(A, B, exact_null_df(A, B, Y[1], Y[2], 1e-13) - 1e-10), delta_df(A, B, exact_null_df(A, B, Y[1], Y[2], 1e-13)));
}




// + протестировал df; null_df; exact_null_df:
// аналитически проверил данные функции в geogebra 
// всё сошлось!!! скрины при запросе)
// также тестировались пересечения орбит (функция f)
// пример получения значений заскринины 
/*
TEST(Orbita, plus) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., -M_PI / 2., 0.);
    std::vector<double> Z = null_df(A, B);
    for(int i = 0; i < Z.size(); ++i)
        std::cout << i << " - " << Z[i] << "\n";    
    std::cout << std::setprecision(15) << exact_null_df(A, B, Z[1], Z[2], 1e-13)   << std::endl;
}
*/


// также аналитически протестирован угол между пересечениями кривых
// получил хороший параметр, delta(phi) от которого я собираюсь отсупить 1е-4
// при сравнении пересечения орбит под острым углом, среднем и близким к прямому
// далее -> сверлся с geogebra (фото имеется) 
// (слишком маленькие phi - плохо, машинные е, слишкоб большие - тоже плохо
// в итоге добиваюсь ошибки косинуса в 4(при очень малых углах) - 6 порядке - не идеально, но с производными было бы хуже

TEST(Orbit, cos_beta_1) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., -M_PI / 2., 0.);

    std::vector<double> Y = intersection(A, B, 1e-12);

    ASSERT_NEAR(0.539504, cos_beta(A, B, Y[0]) ,1e-6);
}
// большой угол пересечения
TEST(Orbit, cos_beta_2) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    //x1 = 7.08; y1 = -2.08;
    double a = 2.08 * sqrt(2) + 5 * sqrt(2) / 2;
    double e = 5 * sqrt(2. / 4.) / a;
    orbit B(a * (1 - e * e), e, 0., -M_PI / 4., 0.);

    std::vector<double> Y = intersection(A, B, 1e-8);

    ASSERT_NEAR(0.99923897, cos_beta(A,B, Y[1]) ,1e-5);
}
// и острый угол
TEST(Orbit, cos_beta_3) {
    orbit A(10. * 3. / 4., 1./2., 0., 0., 0.);
        
    double a = 100;
    double e = sqrt(1. - 25. / 10000.); 
    orbit B(a * (1 - e*e), e, 0., -M_PI, 0.);

    std::vector<double> Y = intersection(A, B, 1e-10);
    ASSERT_NEAR(0.056178, cos_beta(A,B,Y[1]), 1e-5);
}


// тестирование экстримальных точек функции f
TEST(Orbit, extr_1) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., -M_PI / 2., 0.);

    std::vector<double> Z = extr(A, B, 1e-12); // получение extr

    ASSERT_NEAR(1.3599277029, Z[0] ,1e-11);
    ASSERT_NEAR(3.352461277485, Z[1] ,1e-11);

}
// тест где кривые накладывеются пересекаются(значит экетремумов больше 4)))
TEST(Orbit, extr_2) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., 0., 0.);

    std::vector<double> Z = extr(A, B, 1e-12); // получение extr
    
    ASSERT_LT(2, Z.size());
}



TEST(Orbit, touch_1) {
    orbit A(10., 1./2., 0., 0., 0.);
    orbit B(10., 1./2., 0., - M_PI / 2., 0.);

    std::vector<double> Z = extr(A, B, 1e-12); // получение extr
    
    ASSERT_NEAR(touch(A,B, Z[0], 1e-9), 0 ,1e-11);

} 

TEST(Orbit, touch_2) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    double a = 15./2. + 1e-8;
    double e =  1./3.;
    orbit B(a * (1. - e*e), e, 0., 0., 0.);
    
    std::vector<double> Z = extr(A, B, 1e-14); // получение extr

    ASSERT_NEAR(touch(A,B, Z[0], 1e-9), 0 ,1e-4);
    ASSERT_NEAR(touch(A,B, Z[0], 1e-3), 1 ,1e-4);

}


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//************* тестирование полной энергии спутника ***************
TEST(Sputnik, h) {
    orbit A(1., 1/sqrt(2), M_PI , 0., M_PI);
    sputnik S(A, GM);   
    ASSERT_NEAR(S.get_h(), - GM * (1 - A.e * A.e) / A.p ,1e-15);
    ASSERT_GT(0,S.get_h());
}

//******* тестирование радиальной скорости спутника ********************
TEST(Sputnik, v_r) {
    orbit A(1, 1/sqrt(2), M_PI , 0., M_PI);
    sputnik S(A,GM);
    ASSERT_NEAR(0., S.v_r(0) ,1e-15);
}

// ************ тестирование тангенциальной скорости спутника *******
TEST(Sputnik, v_tau) {
    orbit A(1, 1/sqrt(2), M_PI , 0., M_PI);
    sputnik S(A,GM);
    ASSERT_NEAR(S.v_tau(-M_PI/4), S.v_tau(M_PI/4) ,1e-15);
}

// ****** тестирование полной скорости спутника через проекции ************
TEST(Sputnik, v_1) {
    orbit A(10'000'000, 1/sqrt(2), M_PI , 0., M_PI);
    sputnik S(A, GM);
    ASSERT_NEAR(S.v_tau(M_PI) * S.v_tau(M_PI) + S.v_r(M_PI) * S.v_r(M_PI), S.get_v(M_PI) * S.get_v(M_PI), 1e-8);
}

// получение разумной скорости в апоцентре
TEST(Sputnik, v_2) {
    orbit A(42'164'000, 0, M_PI , 0., M_PI);
    sputnik S(A,GM);
    ASSERT_NEAR(3074.66628, S.get_v(M_PI) ,1e-5);
}


// тестирование маневра в некой точке пересечения

TEST(Sputnik, DELTA_V) {
    orbit A(10'000'000. * 3. / 4., 1./2., 0., 0., 0.);
    orbit B(10'000'000. * 3. / 4., 1./2., 0., -M_PI / 2., 0.);

    sputnik S(A, GM);

    std::vector<double> Z = intersection(A,B,1e-8);// вектор пересечений
 
    ASSERT_NEAR(5154.9359445, S.delta_V(B, Z[0]), 1e-6);
}

TEST(Sputnik, DELTA_VVVV) {
    orbit A(20'000'000. * 3./4., 1./2., 0., 0., 0.);
    double c = 5'000'000.;
    double a = c + 3'819'660.;
    orbit B(a * (1 - c * c / a / a), c/a, 0., -M_PI / 2., 0.);

    sputnik S(A, GM);
    double x = S.manevr(B,1e-8);

}






int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



#endif