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
        ASSERT_NEAR(r[i], peresech(A, B, 1e-8)[i] ,1e-15);
}

TEST(Orbit, peresech_2){
    orbit A(1., 1/sqrt(2), 0 , 0., M_PI / 4);
    orbit B(1., 1/sqrt(2), M_PI/2 , 0., 3 * M_PI / 4);
    std::array<double, 3> r = {1/sqrt(3), 1/sqrt(3), 1/sqrt(3)};
    for(int i = 0; i < 3; ++i)
        ASSERT_NEAR(r[i], peresech(A, B, 1e-8)[i] ,1e-15);
}

// наложение двух плоскостей))
TEST(Orbit, peresech_3){
    orbit A(1., 1/sqrt(2), 0 , 0., 0);
    orbit B(1., 1/sqrt(2), M_PI/2 , 0., M_PI);
    std::array<double, 3> r = {0., 0., 0.};
    for(int i = 0; i < 3; ++i)
        ASSERT_NEAR(r[i], peresech(A, B, 1e-8)[i] ,1e-13);
}

TEST(Orbit, ugol_1) {
    orbit A(20'000'000. * 3./4., 1./2., M_PI/5., -M_PI/3., M_PI/4.);
    
    std::array<double,3> a = A.a_r(0);
    ASSERT_NEAR(a[0]*a[0]+a[1]*a[1]+a[2]*a[2],A.getR(0)*A.getR(0),1e-10);
}

TEST(Orbit, ugol_2) {
    orbit A(20'000'000. * 3./4., 1./2., M_PI/2., 0, M_PI/3.);

    std::array<double,3> a = A.a_r(M_PI/2.);
    ASSERT_NEAR(-0.5 * A.getR(M_PI/2.),A.a_r(M_PI/2.)[0],1e-8);
    ASSERT_NEAR(A.getR(0),A.a_r(0)[1],1e-8);
    ASSERT_NEAR(-sqrt(3.)/2. * A.getR(-M_PI/2.),A.a_r(-M_PI/2.)[2],1e-8);
}
// *********** тестирование функции разности значений при заданном угле в Rp и Ra ********** 
TEST(Orbit, f) {
    orbit A(10. * (3./4.), 1./2., 0., 0., 0.);
    orbit B(10. * (3./4.), 1./2., 0., -M_PI / 2., 0.);
    ASSERT_NEAR(-B.p+A.getRp(), delta_f(A,B,0) ,1e-1);
    ASSERT_NEAR(A.getRa()-B.p, delta_f(A,B,M_PI) ,1e-1);
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

// а теперь добавляем пространство
TEST(Orbit, cos_gama_1) {
    orbit A(30'000'000., 1./sqrt(2), M_PI/4. , M_PI / 3., M_PI / 2.);
    orbit B(30'000'000., 1./sqrt(2), M_PI/4. , 0., M_PI/3.);

//    ASSERT_NEAR(sqrt(3) / 2, cos_gama(A, B, 1e-5) ,1e-4);
}

TEST(Orbit, cos_gama_2) {
    orbit A(10'000'000. * (3./4.), 1./2., M_PI/4., 0., M_PI/2.);
    orbit B(10'000'000. * (3./4.), 1./2., M_PI/4., 0., M_PI / 4.);

    ASSERT_NEAR(sqrt(1./2.), cos_gama(A, B, 1e-5) ,1e-6);
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


// проверка точного касания
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
    ASSERT_NEAR(touch(A,B, Z[0], 1e-8), 1 ,1e-4);
}

// угол одной oрбиты при котором происходит пересечение плоскостей
// т е векторы совподают(радиальный и плоскость) 
TEST(Orbit, ugol_vect){
    orbit A(1., 1./sqrt(2), M_PI/4. , M_PI / 4., M_PI / 4.);
    orbit B(1., 1./sqrt(2), M_PI/4. , 0., 3. * M_PI / 4.);
    std::array<double, 3> r = {1./sqrt(2), 1./sqrt(2), 0.};

    ASSERT_NEAR(PHI(r, A, 1e-5), 5.38257531 - M_PI / 4. , 1e-6);
    ASSERT_NEAR(PHI(r, B, 1e-5), 0., 1e-6);
}

// угол между орбитами в пространстве
TEST(Orbit, ugol_peresech_1){
    orbit A(30'000'000., 1./sqrt(2), M_PI/4., 0., M_PI / 4.);
    orbit B(30'000'000., 1./sqrt(2), M_PI/4., 0., M_PI/ 4.1);

    ASSERT_NEAR(cos_gama(A, B, 1e-7), 0.999817, 1e-6);
}

TEST(Orbit, ugol_peresech_2){
    orbit A(30'000'000., 1./sqrt(2), M_PI/4., 0., M_PI / 5.);
    orbit B(30'000'000., 1./sqrt(2), M_PI/4., 0, 2. * M_PI / 5.);

    ASSERT_NEAR(cos_gama(A, B, 1e-7), cos(M_PI/5), 1e-3);

}

TEST(Orbit, ugol_peresech_3){
    orbit A(30'000'000., 1./sqrt(2), M_PI/4., 0., M_PI / 4.);
    orbit B(30'000'000., 1./sqrt(2), M_PI/4., 0., 3 * M_PI/ 4.);

    ASSERT_NEAR(cos_gama(A, B, 1e-7), 0., 1e-6);
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



// далее идет тестирование функции v_min_p для различных конфигурация

// целевая внутри(но снаружи окружности, см описание функции)
TEST(Sputnik, kasanie_1) {
    orbit A(20'000'000. * 3./4., 1./2., 0., 0., 0.);

    orbit B(30'000'000. * 3./4., 1./2., 0., 0, 0.);

    sputnik S(B, GM);
    ASSERT_NEAR(S.v_min_p(A,1e-7).SPEED, -361.07397, 1e-4);

}

// точное касание с окружностью
TEST(Sputnik, kasanie_2_tochn) {
    orbit A(10'000'000. * 3./4., 1./2., 0., 0., 0.);

    orbit B(30'000'000. * 3./4., 1./2., 0., 0, 0.);

    sputnik S(B, GM);

    ASSERT_NEAR((S.v_min_p(A,1e-7)).SPEED, -1158.5454, 1e-4);
}

// близ этой оркружности но только снаружи
TEST(Sputnik, kasanie_3) {
    orbit A(10'010'000. * 3./4., 1./2., 0., 0., 0.);
    orbit B(30'000'000. * 3./4., 1./2., 0., 0, 0.);

    sputnik S(B, GM);

    ASSERT_NEAR((S.v_min_p(A,1e-7)).SPEED, -1157.25739, 1e-4);
}

// близ этой окружности, но только внутри неё
// важно заметить, что у новой орбиты S.v_min_p(A,1e-7)).ORB
// перецинтр развернулся на 180 градусов
TEST(Sputnik, kasanie_4) {
    orbit A(9'900'000. * 3./4., 1./2., 0., 0., 0.);
    orbit B(30'000'000. * 3./4., 1./2., 0., 0, 0.);
    
    sputnik S(B, GM);

    ASSERT_NEAR((S.v_min_p(A,1e-7)).SPEED, -1171.51367, 1e-4);
}

// добавление скорости
TEST(Sputnik, kasanie_5) {
    orbit A(26'000'000. * 3./4., 1./2., 0., 0., 0.);
    orbit B(110'000'000. * 3./4., 1./2., 0., 0., 0.);

    sputnik S(A, GM);
    
    ASSERT_NEAR((S.v_min_p(B,1e-6)).SPEED, 757.7616, 1e-4);
}

// добавление совсем маленькой скорости
TEST(Sputnik, kasanie_6) {
    orbit A(20'000'000. * 3./4., 1./2., 0., 0., 0.);
    orbit B(20'100'000. * 3./4., 1./2., 0., 0., 0.);
    sputnik S(A, GM);

    ASSERT_NEAR((S.v_min_p(B,1e-6)).SPEED, 4.81267, 1e-4);
}


// далее тесты чекают фунцию поворота, которая 
// выдаёт скорости для маневра и новые орбиты
TEST(Sputnik, povorot_1){
    orbit A(25'000'000., 1./sqrt(2), M_PI/4. , M_PI / 40., M_PI / 4.);
    orbit B(25'000'000., 1./sqrt(2), M_PI/4. , 0.,  M_PI / 4.);
    sputnik S(A, GM);

    ASSERT_NEAR(S.povorot(B,1e-5).SPEED , 221.698 , 1e-3);
}

TEST(Sputnik, povorot_2){
    orbit A(20'000'000., 1./sqrt(2), M_PI/4. , M_PI / 4., M_PI / 4.);
    orbit B(25'000'000., 1./sqrt(2), M_PI/4. , 0., 2 * M_PI / 4.);
    sputnik S(A, GM);

    ASSERT_NEAR(S.povorot(B,1e-5).SPEED , 5'633.44 , 1e-2);

}
TEST(Sputnik, povorot_3){
    orbit A(30'000'000., 1./sqrt(2), M_PI/3. , M_PI / 4., M_PI / 4.);
    orbit B(25'000'000., 1./sqrt(2), M_PI/4. , M_PI /5., 2. * M_PI / 4.);
    sputnik S(A, GM);

    ASSERT_NEAR(S.povorot(B,1e-5).SPEED , 4'203.72 , 1e-2);

}
TEST(Sputnik, povorot_4){
    orbit A(20'000'000., 1./sqrt(2), M_PI/4. , M_PI / 4., M_PI / 4.);
    orbit B(25'000'000., 1./sqrt(2), M_PI/4. , 0., 3. * M_PI / 4.);
    sputnik S(A, GM);

    ASSERT_NEAR(S.povorot(B,1e-5).ORB.OMEGA, M_PI/4 , 1e-6);
    ASSERT_NEAR(S.povorot(B,1e-5).ORB.omega, M_PI/4 , 1e-6);
    ASSERT_NEAR(S.povorot(B,1e-5).ORB.i, 3 * M_PI/4 , 1e-6);
}


TEST(Sputnik, manevr_1){
    orbit A(26'000'000., 1./sqrt(2), M_PI/4. , 0., M_PI / 4.);
    orbit B(25'000'000., 1./sqrt(2), M_PI/4. , 0., M_PI / 4.1);
    sputnik S(A, GM);

    ASSERT_NEAR(S.manevr(B, 1e-7), 309.318 , 1e-2);
}

TEST(Sputnik, manevr_2){
    orbit A(26'000'000., 1./sqrt(2), M_PI/4. , 0., M_PI / 4.);
    orbit B(25'000'000., 1./sqrt(2), M_PI/4. , 0., M_PI / 4.);
    sputnik S(A, GM);

    ASSERT_NEAR(S.manevr(B, 1e-7), 124.275 , 1e-2);
}

TEST(Sputnik, manevr_3){
    orbit A(26'000'000., 1./sqrt(2), M_PI/4. , 0., M_PI / 4.);
    orbit B(35'000'000., 1./sqrt(2), -M_PI/4. , M_PI/6., M_PI / 3.);
    sputnik S(A, GM);

    ASSERT_NEAR(S.manevr(B, 1e-7), 5208.662 , 1e-2);
}

TEST(Sputnik, manevr_4){
    orbit A(26'000'000., 1./sqrt(6), M_PI/2. , 0., M_PI / 4.);
    orbit B(15'000'000., 1./sqrt(2), M_PI/4. , 0., M_PI / 4.5);
    sputnik S(A, GM);

    ASSERT_NEAR(S.manevr(B, 1e-7), 6237.26 , 1e-2);
}

TEST(Sputnik, manevr_5){
    orbit A(15'000'000., 1./sqrt(2), 0 , M_PI/2., 0.);
    orbit B(25'000'000., 1./sqrt(2), 0 , 0., 0.);
    sputnik S(A, GM);

    ASSERT_NEAR(S.manevr(B, 1e-7), 3299.67 , 1e-2);
}
TEST(Sputnik, manevr_6){
    orbit A(15'000'000., 1./sqrt(2), 0. , M_PI/2., 0.);
    orbit B(25'000'000., 1./sqrt(2), 0. , 0., 0.);
    sputnik S(B, GM);

    ASSERT_NEAR(S.manevr(A, 1e-7), 5111.45 , 1e-2);
}

TEST(Sputnik, manevr_7){
    orbit A(26'000'000., 1./sqrt(2), M_PI/2. , 5., M_PI / 4.);
    orbit B(125'000'000., 1./sqrt(2), -M_PI/4. , 0., M_PI / 2.5);
    sputnik S(A, GM);

    ASSERT_NEAR(S.manevr(B, 1e-6), 5057.88 , 1e-2);
}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



#endif